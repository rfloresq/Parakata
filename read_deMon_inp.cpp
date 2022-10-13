// *********************************************************
// Purpose: Read deMon2k INP files
//
// Roberto Flores-Moreno
// *********************************************************

#include <fstream>

#include <fcns.h>
#include <global.h>

using std::ifstream;

void read_deMon_inp(const char*filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };

  int i;
  bool ToBohr = true;
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];

  while(f.getline(line,MAX_STR_SIZE))
  {
    if (strncasecmp(line,"GEOME",5)==0)
    {
      if (strcasestr(line,"MAT")!=NULL)
      {
        errmsg(__FILE__,__LINE__,"Cannot read deMon2k INP Z-Matrix yet",0);
        return;
      }
      if (strcasestr(line,"BOHR")!=NULL) ToBohr = false;
      NATOM = 0;
      f >> str;
      while ((i=symtoan(str))>0)
      {
        ATOMNO[NATOM] = i;
        f >> COORD[1][NATOM][0];
        f >> COORD[1][NATOM][1];
        f >> COORD[1][NATOM][2];
        f.getline(line,MAX_STR_SIZE);
        if (ToBohr)
        {
          COORD[1][NATOM][0] = AngstromToBohr(COORD[1][NATOM][0]);
          COORD[1][NATOM][1] = AngstromToBohr(COORD[1][NATOM][1]);
          COORD[1][NATOM][2] = AngstromToBohr(COORD[1][NATOM][2]);
        }
        NATOM++;
        if (f.eof()) goto CLOSE;
        f >> str;
      }
      break;
    }
  }
CLOSE:
  f.close();
};
