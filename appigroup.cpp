// *********************************************************
// Purpose: Apply I group operations 
//
// 2004-2012: Roberto Flores Moreno 
// *********************************************************

#include <math.h>

#include <fcns.h>
#include <global.h>

int appigroup(int *ref)
{
  int i,iaxis,nr;

  nr = NATOM;

  Vector vector,axisa,axisb,axis[12];
  Vector axis2[20],axis3[15],plane[15];

  axisa = Vector(&COORD[1][ref[0]][0]);
  axisa /= axisa.Norm();
  axisb = Vector(&COORD[1][ref[1]][0]);
  vector = Vector( axisa );
  vector *= axisb.Dot( axisa )/axisa.Dot( axisa );
  axisb -= vector;
  axisb /= axisb.Norm();
  vector = axisa > axisb;

  axis[0] = axisa;
  axis[1] = axisa;
  axis[1].Rotate( vector , atan(2.0)*180.0/M_PI );
  for ( iaxis = 2 ; iaxis < 6 ; iaxis++ )
  {
    axis[iaxis] = axis[iaxis-1];
    axis[iaxis].Rotate( axisa , 72.0 );
  };
  axis[6] = axis[4];
  axis[6] *= -1.0; 
  axis[7] =  axis[5];
  axis[7] *= -1.0; 
  axis[8] = axis[1];
  axis[8] *= -1.0; 
  axis[9] = axis[2];
  axis[9] *= -1.0; 
  axis[10] = axis[3];
  axis[10] *= -1.0; 
  axis[11] = axis[0];
  axis[11] *= -1.0; 

  for ( iaxis = 0 ; iaxis < 12 ; iaxis++ )
  {
    for ( i = 0 ; i < 4 ; i++ )
    {
      appsaxis( 10 , i + 1 , axis[iaxis] , nr );
    };
  };

  axis2[0] = axis[0] + axis[1] + axis[2];
  axis2[0] /= axis2[0].Norm();
  axis2[5] = axis[1] + axis[2] + axis[6];
  axis2[5] /= axis2[5].Norm();
  for ( iaxis = 1 ; iaxis < 5 ; iaxis++ )
  {
    axis2[iaxis] = axis2[iaxis-1];
    axis2[iaxis].Rotate( axisa , 72.0 );
    axis2[iaxis+5] = axis2[iaxis+4];
    axis2[iaxis+5].Rotate( axisa , 72.0 );
  };
  for ( iaxis = 10 ; iaxis < 20 ; iaxis++ )
  {
    axis2[iaxis] = axis2[iaxis-10];
    axis2[iaxis] *= -1.0; 
  };

  for ( iaxis = 0 ; iaxis < 20 ; iaxis++ )
  {
    for ( i = 0 ; i < 2 ; i++ ) 
    {
      appsaxis( 6 , i + 1 , axis2[iaxis] , nr );
    };
  };

  for ( iaxis = 0 ; iaxis < 4 ; iaxis++ )
  {
    axis3[iaxis] = axis[0] + axis[iaxis+1];
    axis3[iaxis+5] = axis[iaxis+1] + axis[iaxis+2];
  };
  axis3[4] = axis[0] + axis[5];
  axis3[9] = axis[5] + axis[1];
  axis3[10] = axis[1] + axis[6];
  for ( iaxis = 11 ; iaxis < 15 ; iaxis++ )
  {
    axis3[iaxis] = axis3[iaxis-1];
    axis3[iaxis].Rotate( axisa , 72.0 ); 
  };

  for ( iaxis = 0 ; iaxis < 15 ; iaxis++ )
  {
    appcaxis( 2 , 1 , axis3[iaxis] , nr );
  };

  for (iaxis = 0 ; iaxis < 5 ; iaxis++ )
  {
    plane[iaxis] = axisa > axis[iaxis+1];
    plane[iaxis+10] = axis[iaxis+1] > axis[iaxis+6];
  };
  for ( iaxis = 0 ; iaxis < 4 ; iaxis++ )
  {
    plane[iaxis+5] = axis[iaxis+1] > axis[iaxis+2];
  };
  plane[9] = axis[5] > axis[1];

  for ( iaxis = 0 ; iaxis < 15 ; iaxis++ )
  {
    appsigma( plane[iaxis] , nr );
  };

  appsaxis( 2 , 1 , axisa , nr );

  return 0;
};


