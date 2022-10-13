// *********************************************************
// Purpose: Apply O group operations 
//
// Nov 2004, Roberto Flores Moreno
// Aug 2005, Roberto Flores Moreno
// Feb 2007, Roberto Flores Moreno
// *********************************************************

#include <math.h>

#include <fcns.h>
#include <global.h>

int appogroup(int *ref)
{
  int i,j,iaxis;    
  int nr;            
  Vector axisa;    
  Vector c2axis;  
  Vector axis[5];

  nr = NATOM;

  axis[0] = Vector(&COORD[1][ref[0]][0]);
  axis[0] /= axis[0].Norm();
  axisa = axis[0];
  axis[1] = Vector(&COORD[1][ref[1]][0]);
  axis[2] = Vector( axis[0] );
  axis[2] *= axis[1].Dot( axis[0] )/axis[0].Dot( axis[0] );
  axis[1] -= axis[2];
  axis[1] /= axis[1].Norm();
  axis[2] = axis[0] > axis[1];
  axis[3] = axis[0];
  axis[4] = axis[1];

  for ( iaxis = 0 ; iaxis  < 3 ; iaxis++ )
  {
    c2axis = axis[iaxis+1] + axis[iaxis+2];
    c2axis /= c2axis.Norm(); 
    for ( i = 0 ; i < 3 ; i++ )
    {
      appsaxis( 4 , i+1, axis[iaxis] , nr );
      appcaxis( 4 , i+1, axis[iaxis] , nr );
      appcaxis( 2 , 1, c2axis , nr );
      c2axis.Rotate( axis[iaxis] , 90.0 );
      appcaxis( 2 , 1, c2axis , nr );
      for ( j = 0 ;  j < 3 ; j++ )
      {
        appsigma( axis[iaxis+j] , nr );
      };
    };
  };

  appsaxis( 2 , 1 , axisa , nr );

  axis[0] += axis[1];
  axis[0] += axis[2];
  axis[0] /= axis[0].Norm();
  for ( iaxis = 1 ; iaxis < 4 ; iaxis++ )
  {
    axis[iaxis] = axis[iaxis-1];
    axis[iaxis].Rotate( axisa , 90.0 );
  };
 
  for ( iaxis = 0 ; iaxis < 4 ; iaxis++ )
  {
    for ( i = 0 ; i < 5 ; i++ )
    {
      appsaxis( 6 , i+1, axis[iaxis] , nr );
    };
  };

  return 0;
};
