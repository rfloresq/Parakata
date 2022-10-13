// *********************************************************
// Purpose: Read XYZ files
//
// 2003-2012, Roberto Flores-Moreno
// *********************************************************

#include <stdlib.h>

#include <fstream>

#include <fcns.h>
#include <global.h>

using std::ifstream;

void readxyz(const char* filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };

  int i,iatom;
  char str[MAX_STR_SIZE];

  f >> NATOM;
  if (NATOM>MAXATOM) 
  {
    NATOM = MAXATOM;
    sprintf(str,"Loading only %d atoms",NATOM);
    errmsg(__FILE__,__LINE__,str,0);
  };
  f.getline(str,MAX_STR_SIZE);
  f.getline(str,MAX_STR_SIZE);
  for ( iatom = 0 ; iatom < NATOM ; iatom++ )
  {
    f >> str;
    if ( str[0] >= '0' && str[0] <= '9' )
    {
      ATOMNO[iatom] = atoi(str);
    }
    else
    {
      ATOMNO[iatom] = symtoan( str );
    };
    for ( i = 0 ; i < 3 ; i++ )
    {
      f >> COORD[1][iatom][i];
      COORD[1][iatom][i] = AngstromToBohr(COORD[1][iatom][i]);
    }
  }


  f.close();

  center();
};
