// *********************************************************
// Purpose: Apply S group operations 
//
// Nov 2004, Roberto Flores-Moreno
// Aug 2005, Roberto Flores-Moreno
// Mar 2007, Roberto Flores-Moreno
// *********************************************************

#include <stdio.h>

#include <fcns.h>
#include <global.h>

int appsgroup(int n,int *ref)
{
  int  i;    
  int nr;    
  char s[10];
  
  if ( n%2 == 1 )
  {
    sprintf(s,"C%dh",n);
    return appcgroup(s,n,ref);
  }
  else
  {
    nr = NATOM;

    Vector axisa(&COORD[1][ref[0]][0]);
    axisa /= axisa.Norm();

    if ( n >= 4 )
    {
      for ( i = 1 ; i <= n/2-1 ; i++ )
      {
        appcaxis( n , i , axisa , nr );
      };
    };

    for ( i = 1 ; i < n ; i++ )
    {
      appsaxis( n , i , axisa , nr );
    };

  };

  return 0;
};
