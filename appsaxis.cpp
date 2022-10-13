// *********************************************************
// Purpose: Apply Sn rotation operation 
//
// Nov 2004, Roberto Flores-Moreno 
// Aug 2005, Roberto Flores-Moreno 
// Mar 2007, Roberto Flores-Moreno 
// *********************************************************

#include <fcns.h>
#include <global.h>

#include <Vector.h>

void appsaxis(int n,int i,Vector axis,int nr)
{
  int iatom;
  double sp,phi;
  Vector v,va,vb;

  if ( i%2 == 0 )
  {
    appcaxis(n,i,axis,nr);
  }
  else
  {
    phi = 360.0*((double)i)/((double)n);
    axis ^= 1.0;

    for ( iatom = 0 ; iatom < nr ; iatom++ )
    {
      if ( ATOMNO[iatom] == 0 ) continue;

      v = Vector(&COORD[1][iatom][0]);
      v.Rotate(axis,phi); 
      sp = Vector( v ).Dot( axis );
      va = Vector( v );
      vb = axis;
      vb *= 2.0*sp;
      va  -= vb;
      clonatom(iatom,va);
    };
  };
}; 
