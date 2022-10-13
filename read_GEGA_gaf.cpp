// *********************************************************
// Purpose: Read GAF files
//
// Roberto Flores-Moreno
// 2003-2005: Cinvestav-IPN
// 2007-2008: Universidad de Guanajuato
// 2008-2009: Cinvestav-IPN
// 2009-2010: Universidad de Guadalajara
//
// *********************************************************

#include <stdlib.h>

#include <string>

#include <fstream>

#include <fcns.h>
#include <global.h>

using namespace std;
using std::ifstream;

void sort_gaf_by_energy(void)
{
  int i,iatom,j,k;
  double hold;
  int atomno[MAXATOM];
  double coord[MAXATOM][3];

  for (i=0;i<NOPTCYC;i++)
  {
    k = i;
    for (j=i+1;j<NOPTCYC;j++)
    {
      if (OPT_ENERGY[j]<OPT_ENERGY[k]) k = j;
    }
    hold = OPT_ENERGY[i];
    get_structure("OPT",i);
    for (iatom=0;iatom<NATOM;iatom++)
    {
      atomno[iatom] = ATOMNO[iatom];
      coord[iatom][0] = COORD[1][iatom][0];
      coord[iatom][1] = COORD[1][iatom][1];
      coord[iatom][2] = COORD[1][iatom][2];
    }
    get_structure("OPT",k);
    save_structure("OPT",i);
    OPT_ENERGY[i] = OPT_ENERGY[k];
    for (iatom=0;iatom<NATOM;iatom++)
    {
      ATOMNO[iatom] = atomno[iatom];
      COORD[1][iatom][0] = coord[iatom][0];
      COORD[1][iatom][1] = coord[iatom][1];
      COORD[1][iatom][2] = coord[iatom][2];
    }
    save_structure("OPT",k);
    OPT_ENERGY[k] = hold;
  }
}

void read_GEGA_gaf(const char* filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  NATOM = 0;
  NOPTCYC = 0;
  while(f.getline(line,MAX_STR_SIZE))
  {
    string stringLine = line;

    int sindex=stringLine.find("NUMBER OF ATOMS");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      sindex=stringLine.find("=")+1;
      NATOM = atoi(&line[sindex]);
    }
    sindex=stringLine.find("THE FOLLOWING IS A MINIMUM");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      f >> str;
      f >> OPT_ENERGY[NOPTCYC];
      f.getline(line,MAX_STR_SIZE);
      for (int iatom=0;iatom<NATOM;iatom++)
      {
        f >> str;
        ATOMNO[iatom] = symtoan(str);
        f >> COORD[1][iatom][0];
        f >> COORD[1][iatom][1];
        f >> COORD[1][iatom][2];
        COORD[1][iatom][0] = AngstromToBohr(COORD[1][iatom][0]);
        COORD[1][iatom][1] = AngstromToBohr(COORD[1][iatom][1]);
        COORD[1][iatom][2] = AngstromToBohr(COORD[1][iatom][2]);
      }
      save_structure("OPT",NOPTCYC);
      NOPTCYC++;
      if (NOPTCYC==MAXOPTCYC) 
      {
        errmsg(__FILE__,__LINE__,"NOt all structures loaded",0);
        goto CLOSE;
      }
    }
  }

CLOSE:
  f.close();

  save_structure("SCF",NOPTCYC-1);
  sort_gaf_by_energy();
};
