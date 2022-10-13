// *********************************************************
// Purpose: Apply translation
//
// Aug 2007, Roberto Flores Moreno: Spatial groups for polymers
// *********************************************************

#include <fcns.h>
#include <global.h>

#include <Vector.h>

void apptrans(Vector axis,double shift, int nr)
{
  int iatom;
  Vector normal,v;

  normal = axis;
  normal ^= 1.0;
  normal *= shift;

  for ( iatom = 0 ; iatom < nr ; iatom++ )
  {
    if ( ATOMNO[iatom] == 0 ) continue;

    v = Vector(&COORD[1][iatom][0]);
    v += normal;
    clonatom(iatom,v);
  };
};
