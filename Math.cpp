#include <Math.h>

double Factorial(int n)
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


