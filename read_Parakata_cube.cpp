// *********************************************************
// Purpose: Read Parakata CUBE files
//
// 2003-2021, Roberto Flores-Moreno
// *********************************************************

#include <stdlib.h>

#include <iostream>
#include <fstream>

#include <global.h>
#include <fcns.h>

using std::ifstream;

static bool readbox(const char* filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return false;
  };
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  while(f.getline(line,MAX_STR_SIZE))
  {
    if (strcasestr(line,"@box")!=NULL)
    {
      if (strstr(line,"{")==NULL)
      {
        f >> str;
        if (str[0]!='{')
        {
          errmsg(__FILE__,__LINE__,"@box must be followed by {",0);
          goto CLOSE_RETURN;
        }
      }
      // Load box
      int i;
      for (i=0;i<12;i++) f >> LATTICE[i];
      //rfm for (i=0;i<12;i++) LATTICE[i] = AngstromToBohr(LATTICE[i]); 
      for (i=1;i<4;i++) f >> NLATTICE[i]; 
      NLATTICE[0] = NLATTICE[1]*NLATTICE[2]*NLATTICE[3];
      f.close();
      return true;
    }
  }
CLOSE_RETURN:
  f.close();
  return false;
}

static void readscalar(const char* filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  // free memory previously used by scalar and allocate 
  // enough space for new field
  if (SCALAR_FIELD) free(SCALAR_FIELD);
  while(f.getline(line,MAX_STR_SIZE))
  {
    if (strcasestr(line,"@scalar")!=NULL)
    {
      if (strstr(line,"{")==NULL)
      {
        f >> str;
        if (str[0]!='{')
        {
          errmsg(__FILE__,__LINE__,"@scalar must be followed by {",0);
          goto CLOSE_RETURN;
        }
      }
      SCALAR_FIELD = (double*) malloc(NLATTICE[0]*sizeof(double));
      double v;
      int i,j,k,l;
      l = 0;
      for (k=0;k<NLATTICE[3];k++)
      {
        for (j=0;j<NLATTICE[2];j++)
        {
          for (i=0;i<NLATTICE[1];i++)
          {
            if (f.eof()) 
              errmsg(__FILE__,__LINE__,"Not enough data, check file",0);
            f >> v;
            SCALAR_FIELD[l] = v;
            l++;
          }
        }
      }
    }
  }
CLOSE_RETURN:
  f.close();
}

static void readgeo(const char* filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  while(f.getline(line,MAX_STR_SIZE))
  {
    if (strcasestr(line,"@molecule")!=NULL)
    {
      if (strstr(line,"{")==NULL)
      {
        f >> str;
        if (str[0]!='{')
        {
          errmsg(__FILE__,__LINE__,"@molecule must be followed by {",0);
          goto CLOSE_RETURN;
        }
      }
      NATOM = 0;
      f >> str;
      while (str[0] != '}')
      {
        if ( str[0] >= '0' && str[0] <= '9' )
          ATOMNO[NATOM] = atoi(str);
        else
          ATOMNO[NATOM] = symtoan(str);
        f >> COORD[1][NATOM][0];
        f >> COORD[1][NATOM][1];
        f >> COORD[1][NATOM][2];
        NATOM++;
        f >> str;
      } 
    }
  }
CLOSE_RETURN:
  f.close();
}

void read_Parakata_cube(const char*filename)
{
  // Load box
  if (!readbox(filename)) 
  {
    errmsg(__FILE__,__LINE__,"Box not found, loading aborted",0);
    return;
  }

  // Load scalar field
  readscalar(filename);

  // Load structure
  readgeo(filename);
};
