// *********************************************************
// Purpose: Make a copy of a non-redundant atom if required
//
// Nov 2004, Roberto Flores-Moreno
// *********************************************************

#include <global.h>

#include <Vector.h>

void clonatom(int from,Vector newpos)
{
  int i,iatom;  
  Vector v;  
  
  for ( iatom = 0 ; iatom < NATOM ; iatom++ )
  {
    if ( ATOMNO[iatom] == 0 ) continue;
    v = newpos;
    v -= Vector(&COORD[1][iatom][0]);
    if ( v.Norm() < 0.5 )
    {
      return;
    };
  };

  for ( i = 0 ; i < 3 ; i++ )
  {
    COORD[1][NATOM][i] = newpos[i];
  }
  ATOMNO[NATOM] = ATOMNO[from];
  NATOM++;
};
