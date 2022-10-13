// *********************************************************
//
// Purpose: Space Group operations applier driver
//
// Aug 2007, Roberto Flores-Moreno:  Space groups for polymers
//
// *********************************************************

#include <iostream>

#include <QtWidgets>

#include <Vector.h>

#include <fcns.h>
#include <global.h>

int sgbdrv(char* sg,int* ref) 
{
  int nr;

  nr = NATOM;

  if ( sg[1] == '1' )
  {
    Vector axis(&COORD[1][ref[1]][0]);
    axis -= Vector(&COORD[1][ref[0]][0]);
    double sd = QInputDialog::getDouble(0,
                "Translation Unit Length","R = ",1.0,0.0,10000.0,3);
    double nt = QInputDialog::getInt(0,
                "Number on Monomers","NM = ",2,2,10000);
    nt = nt - 1;
    double s = 0.0;
    for ( int it = 0 ; it < nt ; it++ )
    {
      s = s + sd;
      apptrans(axis,s,nr);
    }
  }
  else 
  {
    cout << "sgbdrv: Unsupported spatial group"<<endl;
    exit(EXIT_FAILURE);
  };

  center();

  return 0;
};
