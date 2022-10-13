// *********************************************************
// Purpose: Apply T group operations 
//
// Nov 2004, Roberto Flores Moreno
// Aug 2005, Roberto Flores Moreno
// Mar 2007, Roberto Flores Moreno
// *********************************************************

#include <math.h>

#include <fcns.h>
#include <global.h>

#include <Vector.h>

int apptgroup(int *ref)
{
  int       i,iaxis,nr;
  double    tangle; 
  Vector axisb;  
  Vector vector;
  Vector axis[4]; 

  nr = NATOM;

  axis[0] = Vector(&COORD[1][ref[0]][0]);
  axis[0] /= axis[0].Norm();
  axisb = Vector(&COORD[1][ref[1]][0]);
  axis[2] = Vector( axis[0] );
  axis[2] *= axisb.Dot( axis[0] )/axis[0].Dot( axis[0] );
  axisb -= axis[2];
  axisb /= axisb.Norm();

  tangle = acos(-1.0/3.0);
  axis[2] = axis[0] > axisb;
  axis[1] = axis[0];
  axis[1].Rotate( axis[2] , 180.0*tangle/M_PI );
  axis[2] = axis[1];
  axis[2].Rotate( axis[0] , 120.0 );
  axis[3] = axis[2];
  axis[3].Rotate( axis[0] , 120.0 );

  for ( iaxis = 0 ; iaxis  < 4 ; iaxis++ )
  {
    appcaxis(3,1, axis[iaxis], nr );
    appcaxis(3,2, axis[iaxis], nr );
    if ( iaxis == 0 )
    {
      vector = axis[0] > axisb;
    }
    else 
    {
      vector = axis[iaxis-1] > axis[iaxis];
    };
    for ( i = 0 ; i < 3 ; i++ )
    {
      appsigma( vector , nr );
      vector.Rotate( axis[iaxis] , 120.0 );
    };
  };

  for ( iaxis = 0 ; iaxis  < 3 ; iaxis++ )
  {
    if ( iaxis == 2 )
    {
      vector = axis[iaxis] + axis[1];
    }
    else
    {
      vector = axis[iaxis] + axis[iaxis+1];
    };
    for ( i = 0 ; i < 3 ; i++ )
    {
      appsaxis(4, i + 1, vector , nr );
    };
  };

  return 0;
};
