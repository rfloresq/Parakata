// *********************************************************
//
// Purpose: Symmetry Point/Spatial Group operations caller 
//
// Aug 2007, Roberto Flores-Moreno, Zeferino Gomez-Sandoval
//
// *********************************************************

#include <fcns.h>

int sgdrv(char* gname,int* ref) 
{
  if ( gname[0] == 's' )
  {
    return sgbdrv(gname,ref);
  }
  else 
  {
    return pgbdrv(gname,ref);
  };
};
