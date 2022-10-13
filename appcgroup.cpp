// *********************************************************
// Purpose: Apply C group operations 
//
// 2004-2007: Roberto Flores Moreno 
// *********************************************************

#include <string.h>

#include <fcns.h>
#include <global.h>
#include <Vector.h>

int appcgroup(char* pg,int n,int *ref)
{
  int i,nr;  

  if ( n < 2 ) return 0;

  Vector axisa(&COORD[1][ref[0]][0]);

  nr = NATOM;
  for ( i = 1 ; i < n ; i++ )
  {
    appcaxis(n,i,axisa,nr);
  };
  
  if ( strlen(pg) == 2 ) return 0;
  
  if ( pg[2] == 'h' )
  {
    for ( i = 1 ; i < n ; i++ )
    {
      appsaxis( n , i , axisa , nr );
    };
    appsigma( axisa , nr );
  }
  else 
  {
    Vector axisb(&COORD[1][ref[1]][0]);
    Vector normal( axisa );
    normal *= axisb.Dot( axisa )/axisa.Dot( axisa );
    axisb -= normal;
    normal = axisa > axisb;
    for ( i = 1 ; i <= n ; i++ )
    {
      appsigma(normal,nr);
      normal.Rotate(axisa,360.0/((double)n));
    };
  };
  return 0;
};

