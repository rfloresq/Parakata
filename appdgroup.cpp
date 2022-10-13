// *********************************************************
// Purpose: Apply D group operations 
//
// 2004-2012: Roberto Flores Moreno 
// *********************************************************

#include <fcns.h>
#include <global.h>

int appdgroup(char* pg,int n,int *ref)
{
  int i,top,step,nr; 

  nr = NATOM;

  Vector axisa(COORD[1][ref[0]][0]);
  axisa /= axisa.Norm();

  if ( pg == "D*h" )
  {
    appsaxis(2,1,axisa,nr);
    return 0;
  };
  
  for ( i = 1 ; i < n ; i++ )
  {
    appcaxis( n , i , axisa , nr );
  };

  Vector axisb(&COORD[1][ref[1]][0]);
  Vector vector( axisa );
  vector *= axisb.Dot( axisa )/axisa.Dot( axisa );
  axisb -= vector;
  axisb /= axisb.Norm();
  vector = axisb;
  for ( i = 1 ; i <= n ; i++ )
  {
    appcaxis( 2 , 1 , vector , nr );
    vector.Rotate( axisa , 180.0/((double)n) );
  };

  if ( pg[2] == 'h' )
  {
    if ( n%2 == 0 )
    {
      top = n - 1;
      step = 1;
    }
    else
    {
      top = 2*n-1;
      step = 2;
    };
    for ( i = 1 ; i <= top ; i += step )
    {
      appsaxis( n , i , axisa , nr );
    };
  }
  else if ( pg[2] == 'd' )
  {
    for ( i = 1 ; i <= 2*n-1 ; i += 2 )
    {
      appsaxis( n , i , axisa , nr );
    };
  };

  vector = axisa > axisb;
  if ( pg[2] == 'd' )
  {
    vector.Rotate( axisa , 180.0/((double)n) );
  };
  for ( i = 1 ; i <= n ; i++ )
  {
    appsigma( vector , nr );
    vector.Rotate( axisa , 360.0/((double)n) );
  };

  if ( pg[2] == 'h' )
    appsigma( axisa , nr );

  return 0;
};
