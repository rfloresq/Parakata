// *********************************************************
// Purpose: Point Group operations applier driver
//
// Nov 2004, Roberto Flores-Moreno
// Aug 2005, Roberto Flores-Moreno
// Mar 2007, Roberto Flores-Moreno
// *********************************************************

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fcns.h>
#include <global.h>

#include <Vector.h>

int pgbdrv(char* pg,int* ref) 
{
  int n;

  if ( strcmp(pg,"Cs") == 0 )
  {
    strcpy(pg,"C1h");
  }
  else if ( strcmp(pg,"Ci") == 0 )
  {
    strcpy(pg,"S2");
  }
  else if ( strcmp(pg,"Vh") == 0 )
  {
    strcpy(pg,"D2h");
  } 
  else if ( strcmp(pg,"Vd") == 0 )
  {
    strcpy(pg,"D2d");
  }
  else if ( strcmp(pg,"S6v") == 0 )
  {
    strcpy(pg,"D3d");
  }
  else if ( strcmp(pg,"S8v") == 0 )
  {
    strcpy(pg,"D4d");
  };

  n = 0;
  if ((pg[0] != 'T')&&(pg[0] != 'O')&&(pg[0] != 'I')&&
      (pg[1] != '*')&&(pg[1] != 's'))
  {
    char s[10];
    sprintf(s,"%c",pg[1]);
    n = atoi(s);
  };

  Vector v(&COORD[1][ref[0]][0]);
  for ( int iatom = 0 ; iatom < NATOM ; iatom++ )
  {
    for (int i = 0 ; i < 3 ; i++ )
    {
      COORD[1][iatom][i] -= v[i];
    }
  }

  if ( pg[0] == 'C' )
  {
    appcgroup(pg,n,&ref[1]);
  }
  else if ( pg[0] == 'D' )
  {
    appdgroup(pg,n,&ref[1]);
  }
  else if ( pg[0] == 'S' )
  {
    appsgroup(n,&ref[1]);
  }
  else if ( pg[0] == 'T' )
  {
    apptgroup(&ref[1]);
  }
  else if ( pg[0] == 'O' )
  {
    appogroup(&ref[1]);
  }
  else if ( pg[0] == 'I' )
  {
    appigroup(&ref[1]);
  };

  center();

  return 0;
};
