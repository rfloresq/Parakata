// ******************************************************************
//
// Purpose: Vector (3D) management
//
// Roberto Flores-Moreno
// 2003-2005: Cinvestav-IPN
// 2007-2008: Universidad de Guanajuato
// 2008-2009: Cinvestav-IPN
// 2009-2010: Universidad de Guadalajara
//
// ******************************************************************

#include <math.h> 
#include <stdio.h>

#include <Vector.h>

// Constructors
Vector::Vector( double *v )
{
  x = v[0];
  y = v[1];
  z = v[2];
};

Vector::Vector( double xx , double yy , double zz ) 
{
  x = xx;
  y = yy;
  z = zz;
};

Vector::~Vector()
{
}

// Access
double Vector::operator [] ( short n )
{
  if ( n == 0 )
  {
    return x;
  }
  else if ( n == 1 )
  {
    return y;
  }
  else if ( n == 2 )
  {
    return z;
  } 
  else
  {
    return 0.0;
  };
};

// this = result 

void Vector::operator += ( const Vector &cv )
{
  Vector v = cv;
  x += v[0];
  y += v[1];
  z += v[2];
};

void Vector::operator -= ( const Vector &cv )
{
  Vector v = cv; 
  x -= v[0];
  y -= v[1];
  z -= v[2];
};

void Vector::operator *= ( double s )
{
  x *= s;
  y *= s;
  z *= s;
};
void Vector::operator *= ( int n )
{
  x *= double(n);
  y *= double(n);
  z *= double(n);
};

void Vector::operator /= ( double s )
{
  x /= s;
  y /= s;
  z /= s;
};

void Vector::operator <= ( const Vector &cv )
{
  double m[3][3];
  Vector v = cv;

  m[0][0] = 0.0;   m[0][1] = -v[2]; m[0][2] = v[1];
  m[1][0] = v[2];  m[1][1] = 0.0;   m[1][2] = -v[0];
  m[2][0] = -v[1]; m[2][1] = v[0];  m[2][2] = 0.0;

  operator %= ( m );
};

void Vector::operator >= ( const Vector &v )
{
  operator <= ( v );
  operator *= ( -1.0 );
};

void Vector::operator %= ( double m[3][3] )
{
  double nx,ny,nz;

  nx = Dot( Vector( m[0][0] , m[0][1] , m[0][2] ) ); 
  ny = Dot( Vector( m[1][0] , m[1][1] , m[1][2] ) );
  nz = Dot( Vector( m[2][0] , m[2][1] , m[2][2] ) );

  x = nx;
  y = ny;
  z = nz;
};

void Vector::operator ^= ( double s ) 
{
  if ( Norm() > 0.0 )
  {
    operator *= ( s/Norm() );
  };
};


// this != result

Vector Vector::operator + ( const Vector &v )
{
  Vector u( x , y , z );
  u += v;
  return u;
};

Vector Vector::operator - ( const Vector &v )
{
  Vector u( x , y , z );
  u -= v;
  return u;
};

Vector Vector::operator < ( const Vector &v )
{
  Vector u( x , y , z );
  u <= v;
  return u;
};

Vector Vector::operator > ( const Vector &v )
{
  Vector u( x , y , z );
  u >= v;
  return u;
};

Vector Vector::operator % ( double m[3][3] )
{
  Vector u( x , y , z );
  u %= m;
  return u;
};

Vector Vector::operator ^ ( double s ) 
{
  Vector u( x , y , z );
  u ^= s;
  return u;
};

Vector Vector::operator * ( double s ) 
{
  Vector u( x , y , z );
  u *= s;
  return u;
};

Vector Vector::operator * ( int n ) 
{
  Vector u( x , y , z );
  u *= n;
  return u;
};

Vector Vector::operator / ( double s ) 
{
  Vector u( x , y , z );
  u /= s;
  return u;
};

// functions

double Vector::Dot( const Vector &cv )
{
  Vector v = cv;
  return( x*v[0] + y*v[1] + z*v[2] );
};

double Vector::Norm() 
{
  return  sqrt( Dot( Vector( x , y , z ) ) );
};

void Vector::Rotate( const Vector &v , double angle )
{
  double cosp,sinp,cost;
  double m[3][3];

  Vector u = v;

  if ( u.Norm() == 0.0 ) 
  {
    return;
  }
  u ^= 1.0;

  angle = angle*M_PI/180.0;

  cosp = cos(angle);
  sinp = sin(angle); 
  cost = 1.0 - cosp;

  m[0][0] = u[0]*u[0]*cost + cosp;
  m[0][1] = u[0]*u[1]*cost - u[2]*sinp;
  m[0][2] = u[0]*u[2]*cost + u[1]*sinp;
  m[1][0] = u[0]*u[1]*cost + u[2]*sinp;
  m[1][1] = u[1]*u[1]*cost + cosp;
  m[1][2] = u[1]*u[2]*cost - u[0]*sinp;
  m[2][0] = u[0]*u[2]*cost - u[1]*sinp;
  m[2][1] = u[1]*u[2]*cost + u[0]*sinp;
  m[2][2] = u[2]*u[2]*cost + cosp;

  operator %= ( m );
};

