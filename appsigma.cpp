// *********************************************************
// Purpose: Apply sigma plane reflexion 
//
// Nov 2004, Roberto Flores Moreno 
// Aug 2005, Roberto Flores Moreno 
// Mar 2007, Roberto Flores Moreno 
// *********************************************************

#include <fcns.h>
#include <global.h>

#include <Vector.h>

void appsigma(Vector normal,int nr)
{
  int iatom;
  double sp;
  Vector va,vb;

  normal ^= 1.0;

  for ( iatom = 0 ; iatom < nr ; iatom++ )
  {
    if ( ATOMNO[iatom] == 0 ) continue;

    sp = Vector(&COORD[1][iatom][0]).Dot(normal);
    va = Vector(&COORD[1][iatom][0]);
    vb = normal;
    vb *= 2.0*sp;
    va  -= vb;
    clonatom(iatom,va);
  };
};
