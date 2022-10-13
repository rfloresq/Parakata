// *******************************************************************
// Purpose: Read external files
//
// 2003-2021, Roberto Flores Moreno
// *******************************************************************

#include <iostream>
#include <ctype.h>
#include <string.h>

#include <fcns.h>
#include <global.h>

#include <Parakata.h>
#include <Viewer.h>
#include <PlotManager.h>
#include <SCFWin.h>
#include <OPTWin.h>
#include <FREQWin.h>

void reader(const char*filename,const char*fmt)
{
  int iatom,l;

  NSCFCYC = 0;
  NVIB = 0;
  NORB = 0;
  MOLDIPOLE[0] = 0.0;
  MOLDIPOLE[1] = 0.0;
  MOLDIPOLE[2] = 0.0;
  for (iatom=0;iatom<NATOM;iatom++)
  {  
     mb.ab[iatom].nshl = 0;
     ATOMCHARGE[iatom] = 0.0;
     ATOMFUKUI_H[iatom] = 0.0;
     ATOMFUKUI_L[iatom] = 0.0;
     ATOMFORCE[iatom][0] = 0.0;
     ATOMFORCE[iatom][1] = 0.0;
     ATOMFORCE[iatom][2] = 0.0;
     ATOMNUCFUKUI_H[iatom][0] = 0.0;
     ATOMNUCFUKUI_H[iatom][1] = 0.0;
     ATOMNUCFUKUI_H[iatom][2] = 0.0;
     ATOMNUCFUKUI_L[iatom][0] = 0.0;
     ATOMNUCFUKUI_L[iatom][1] = 0.0;
     ATOMNUCFUKUI_L[iatom][2] = 0.0;
  }
  if (DMR1) free(DMR1);
  if (DMR2) free(DMR2);
  if (DMR3) free(DMR3);
  if (FUKUIL_MATRIX) free(FUKUIL_MATRIX);
  if (FUKUIR_MATRIX) free(FUKUIR_MATRIX);
  
  if ( strcasecmp( fmt , "XYZ" ) == 0 )
  {
    readxyz(filename);
  }
  else if ( strcasecmp( fmt , "deMon2k" ) == 0 )
  {
    l = strlen(filename);
    if ( filename[l-1]=='p' || filename[l-1] == 'w') read_deMon_inp(filename);
    else read_deMon_out(filename);
  }
  else if ( strcasecmp( fmt , "Nagual" ) == 0 )
  {
    read_Nagual_out(filename);
  }
  else if ( strcmp( fmt , "Parakata" ) == 0 )
  {
    read_Parakata_cube(filename);
  }
  // Prepare basis and basis/orbital settings
  setup_basis();
  if (prk) 
  {  
    prk->viewer->Redraw();
    if ( strcasecmp( fmt , "xyz" ) != 0 )
      prk->viewer->pmgr->Launch();
  }
};
