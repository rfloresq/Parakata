// *********************************************************
// Purpose: Center coordinates
//
// Roberto Flores Moreno
// 2003-2005: Cinvestav-IPN
// 2007-2008: Universidad de Guanajuato
// 2008-2009: Cinvestav-IPN
// 2009-2021: Universidad de Guadalajara
// *********************************************************

#include <global.h>
#include <Vector.h>

void center(void)
{
    
  Vector ref(0.0,0.0,0.0);
  for (int iatom = 0 ; iatom < NATOM ; iatom++ )
  {
    ref += Vector(&COORD[1][iatom][0]);
  }
  ref /= (double)NATOM;
  for (int iatom = 0 ; iatom < NATOM ; iatom++ )
  {
    COORD[1][iatom][0] -= ref[0];
    COORD[1][iatom][1] -= ref[1];
    COORD[1][iatom][2] -= ref[2];
  }
  
};
