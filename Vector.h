// ******************************************************************
// Purpose: Vector (3D) management
//
// 2003-2012: Roberto Flores-Moreno
// ******************************************************************

#ifndef VECTOR_H
#define VECTOR_H

class Vector
{
public:
  Vector( double xx = 0.0 ,
          double yy = 0.0 ,
          double zz = 0.0 );
  Vector( double* );
  virtual ~Vector();

  double x;
  double y;
  double z; 

  Vector operator + ( const Vector & );
  Vector operator - ( const Vector & );
  Vector operator * ( double );
  Vector operator * ( int );
  Vector operator / ( double );
  Vector operator % ( double m[3][3] );
  Vector operator ^ ( double );
  Vector operator < ( const Vector & );
  Vector operator > ( const Vector & );

  void operator += ( const Vector & );
  void operator -= ( const Vector & );
  void operator >= ( const Vector & );
  void operator <= ( const Vector & );
  void operator *= ( double );
  void operator *= ( int );
  void operator /= ( double );
  void operator %= ( double m[3][3] );
  void operator ^= ( double );

  virtual double operator [] ( short n );

  double Dot( const Vector & );
  double Norm(void);
  void Rotate( const Vector & , double );


};

#endif

