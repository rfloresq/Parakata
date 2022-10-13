// **************************************************************
// Purpose: Read OUT (Nagyual 1.3) files
//
// Roberto Flores-Moreno
//
// **************************************************************

#include <stdlib.h>
#include <malloc.h>

#include <iostream>
#include <fstream>
#include <string>

#include <param.h>
#include <macros.h>
#include <fcns.h>
#include <global.h>

using namespace std;
using std::ifstream;

bool read_Nagual_out_geo(const char*filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return false;
  };
  
  NOPTCYC = 0;
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  bool ret = false;
  while (f.getline(line,MAX_STR_SIZE))
  {
    int sindex;

    string stringLine=line;
    sindex=stringLine.find(" Number of atoms:");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      NATOM = atoi(&line[17]);
      NELECTRON = 0; // FIXME
      if (NATOM>MAXATOM) 
      {
        NATOM = MAXATOM;
        sprintf(str,"Loading only %d atoms",NATOM);
        errmsg(__FILE__,__LINE__,str,0);
      }
      continue;
    }
    sindex=stringLine.find("Geometry in angstroms:");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      for( int iatom = 0; iatom < NATOM ; iatom++ )
      {
        f >> str;
        f >> COORD[1][iatom][0];
        f >> COORD[1][iatom][1];
        f >> COORD[1][iatom][2];
        COORD[1][iatom][0] = AngstromToBohr(COORD[1][iatom][0]);
        COORD[1][iatom][1] = AngstromToBohr(COORD[1][iatom][1]);
        COORD[1][iatom][2] = AngstromToBohr(COORD[1][iatom][2]);
        ATOMNO[iatom] = symtoan(str);
        f.getline(line, MAX_STR_SIZE);
      }
      ret = true;
    }
  }
  f.close();
  center();
  save_structure("SCF",0);
  return ret;
}

bool read_Nagual_out_bas(const char*filename)
{
  char line[MAX_STR_SIZE];
  char basis[MAXATOM][MAX_STR_SIZE];
  int iatom;

  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return false;
  };

  while(f.getline(line,MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex=stringLine.find("Basis sets:");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      f.getline(line,MAX_STR_SIZE);
      for (iatom = 0 ; iatom < NATOM ; iatom++)
      {
        f.getline(line,MAX_STR_SIZE);
        f.getline(line,MAX_STR_SIZE);
        f.getline(line,MAX_STR_SIZE);

        AtomBasis ab;
        ShellBasis sb;
        int n,iblk,nblk;
        ab.nshl = 0;
        ab.ncao = 0;
        ab.nsao = 0;
        f >> nblk;
        for (iblk=0;iblk<nblk;iblk++)
        {
          f >> n;
          f >> sb.l;
          f >> sb.k;
          for (n=0;n<sb.k;n++)
          {
            f >> sb.z[n];
            f >> sb.d[n];
          }
          ab.sb[ab.nshl] = sb;
          ab.ncao += (sb.l+1)*(sb.l+2)/2;
          ab.nsao += 2*sb.l+1;
          ab.nshl++;
        }
        mb.ab[iatom] = ab;
      }
    }
  } 
  return true;
}

bool read_Nagual_out_mos(const char* filename)
{
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];

  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return false;
  };

  NORB = 0;
  NBASIS = 0;
  int nblk = 5;
  while(f.getline(line,MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex;
    sindex=stringLine.find(" Number of basis functions:");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      NBASIS = atoi(&line[31]);
    sindex=stringLine.find(" Coefficients for set   1");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      cout << NBASIS << endl;
      int i,imo,ibas,ll,ul,n;
      bool betas;
      MO mo[5];
      int TOP = NBASIS;
      strcpy(mo[0].spin,"ALPHA");
      mo[0].nao = TOP;
      for (i=1;i<nblk;i++) 
      {  
        strcpy(mo[i].spin,mo[0].spin);
        mo[i].nao = NBASIS;
      }
      ll = 0;
      ul = MIN(ll + nblk-1,TOP-1);
      while(ul>=ll)
      {
        for (imo = ll; imo <= ul; imo++)
        {
          f>>n;
          if (n-imo!=1) break;
        }
        ul = imo-1;
        for (ibas=0;ibas<TOP;ibas++)
        {
          for (imo = ll; imo <= ul; imo++) 
          {
            f >> mo[imo-ll].c[ibas];
          }
        }
        for (imo = ll; imo <= ul; imo++)
        {
          MOS[NORB] = mo[imo-ll];
          NORB++;
        }
        ll = ul+1;
        ul = MIN(ll + nblk-1,TOP-1);
      }
    }
  }

  f.close();

  if (NORB>0) return true;
  else return false;

}


void read_Nagual_out(const char*filename)
{
  if (!read_Nagual_out_geo(filename)) 
    errmsg(__FILE__,__LINE__,"GEOMETRY loading failed",0);


  if (!read_Nagual_out_bas(filename))
    errmsg(__FILE__,__LINE__,"BASIS loading failed",0);
  else 
  {
    if (!read_Nagual_out_mos(filename))
      errmsg(__FILE__,__LINE__,"Basis found but MOs missed",0);
  }
}
