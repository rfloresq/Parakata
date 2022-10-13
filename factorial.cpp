// ***********************************************************
// Purpose: Factorial and related functions
//
// Roberto Flores Moreno
// 2003-2005: Cinvestav-IPN
// 2007-2008: Universidad de Guanajuato
// 2008-2009: Cinvestav-IPN
// 2009-2010: Universidad de Guadalajara
// ***********************************************************

#include <macros.h>

double factorial(int n)
{
  double r;

  if ( n < 0 )
  {
    r = 0.0;
  }
  else
  {
    r = 1.0;
    for ( int i = 1; i <= n; i++ )
    {
      r *= (double)i;
    }
  }
  return r;
}

double double_factorial(int n)
{
  int i;
  double r;

  if ( n <= 1 )
  {
    r = 1.0;
  }
  else if ( MOD(n,2) == 1 )
  {
    r = 1.0;
    for ( i = 1 ; i <= n ; i += 2 )
    {
      r *= (double)i;
    }
  }
  else
  {
    r = 1.0;
    for ( i = 2 ; i <= n ; i += 2 )
    {
      r *= (double)i;
    }
  }

  return r;
}

double noverk(int n,int k)
{
  double r;

  if (k>=0&&k<=n)
  {
    r = factorial(n)/(factorial(n-k)*factorial(k));
  }
  else
  {
    r = 0.0;
  }
  return r;
}
