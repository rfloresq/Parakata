// *********************************************************
// Purpose: Compute distance between two atoms a-b
//
// Roberto Flores Moreno
// 2003-2005: Cinvestav-IPN
// 2007-2008: Universidad de Guanajuato
// 2008-2009: Cinvestav-IPN
// 2009-2010: Universidad de Guadalajara
// *********************************************************

#include <global.h>

#include <Vector.h>

double distance(int a, int b)
{
  double val;
  Vector va,vb,vc;

  va = Vector(&COORD[1][a][0]);
  vb = Vector(&COORD[1][b][0]);
  vc = va - vb;
  val = vc.Norm(); 
  return val;
}
