// **************************************************************
// Purpose: Read OUT (deMon2k 6.1.2) files
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

bool read_deMon_out_geo(const char*filename)
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
  bool grad,orig;
  bool xw = true;
  bool ret = false;
  orig = false;
  while (f.getline(line,MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex=stringLine.find("CONDENSED FUKUI"); 
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    { int iatom;
      f.getline(line,MAX_STR_SIZE);
      f.getline(line,MAX_STR_SIZE);
      f.getline(line,MAX_STR_SIZE);
      f.getline(line,MAX_STR_SIZE);
      for(iatom=0;iatom<NATOM;iatom++)
      {
        f >> line;
        f >> ATOMFUKUI_H[iatom];
        f >> ATOMFUKUI_L[iatom];
        f.getline(line,MAX_STR_SIZE);
      }
      continue;
    }  
    sindex=stringLine.find("MULLIKEN POPULATION ANALYSIS");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) xw = false;
    sindex=stringLine.find(" ATOMIC NET CHARGES"); 
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    { int iatom;
      for(iatom=0;iatom<NATOM;iatom++)
      {
        f >> line;
        if (xw) f >> line;
        f >> ATOMCHARGE[iatom];
        cout << ATOMCHARGE[iatom] << endl;
      }
      continue;
    } 
    sindex=stringLine.find(" NUMBER OF ATOMS:");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      NATOM = atoi(&line[17]);
      f >> str;
      f >> str;
      f >> str;
      f >> str;
	NELECTRON = atoi(str);
      if (NATOM>MAXATOM) 
      {
        NATOM = MAXATOM;
        sprintf(str,"Loading only %d atoms",NATOM);
        errmsg(__FILE__,__LINE__,str,0);
      }
      continue;
    }
    sindex=stringLine.find("TOTAL ENERGY");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      sindex=stringLine.find("=")+1;
      OPT_ENERGY[NOPTCYC] = strtod(&line[sindex],0);
    }
    sindex=stringLine.find("INPUT ORIENTATION IN");
    if(!((sindex>=0)&&(sindex<(signed)stringLine.size())))
      sindex=stringLine.find("COORDINATES OF OPTIMIZATION STEP");
    if(!((sindex>=0)&&(sindex<(signed)stringLine.size())))
      sindex=stringLine.find("CARTESIAN GRADIENTS IN A.U.");
    if(!((sindex>=0)&&(sindex<(signed)stringLine.size())))
      sindex=stringLine.find("MOLECULE ORIENTATION");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      sindex=stringLine.find("CARTESIAN GRADIENTS IN A.U.");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) { grad = true; }
      else { grad = false; }

      sindex=stringLine.find("ORIGINAL INPUT ORIENTATION");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) { orig = true; }

      bool AA = false;
      if ( strstr( line, "ANGSTROM" ) != NULL ) AA = true;
      f.getline(line, MAX_STR_SIZE);
      f.getline(line, MAX_STR_SIZE);
      f.getline(line, MAX_STR_SIZE); 
      for( int iatom = 0; iatom < NATOM ; iatom++ )
      {
        f >> str;
        f >> str;
        f >> COORD[1][iatom][0];
        f >> COORD[1][iatom][1];
        f >> COORD[1][iatom][2];
        if (AA) 
        {
          COORD[1][iatom][0] = AngstromToBohr(COORD[1][iatom][0]);
          COORD[1][iatom][1] = AngstromToBohr(COORD[1][iatom][1]);
          COORD[1][iatom][2] = AngstromToBohr(COORD[1][iatom][2]);
        }
        ATOMNO[iatom] = symtoan(str);
        f.getline(line, MAX_STR_SIZE);
      }
      ret = true;
      if (grad) 
      {
        save_structure("GRAD",NOPTCYC+1);
      }
      else
      {
        center();
        save_structure("OPT",NOPTCYC+1);
        NOPTCYC++;
      }
    }
  }
  f.close();
  if (orig) 
  {
    center();
    save_structure("SCF",0);
  }
  return ret;
}

void read_out_cps(const char*filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };
  
  NCPs = 0;
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  while (f.getline(line,MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex=stringLine.find("(RANK,SIGNATURE):"); 
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    { 
      f >> str;
      f >> str;
      f >> str;
      f >> CPCOORD[NCPs][0];
      f >> CPCOORD[NCPs][1];
      f >> CPCOORD[NCPs][2];
      CPCOORD[NCPs][0] = AngstromToBohr(CPCOORD[NCPs][0]);
      CPCOORD[NCPs][1] = AngstromToBohr(CPCOORD[NCPs][1]);
      CPCOORD[NCPs][2] = AngstromToBohr(CPCOORD[NCPs][2]);
      NCPs++;
    }
  }
  if (NCPs>0)
  {
    Vector ref(0.0,0.0,0.0);
    for (int i = 0 ; i < NCPs ; i++ )
    {
      ref += Vector(&CPCOORD[i][0]);
    }
    ref /= (double)NCPs;
    for (int i = 0 ; i < NCPs ; i++ )
    {
      CPCOORD[i][0] -= ref[0];
      CPCOORD[i][1] -= ref[1];
      CPCOORD[i][2] -= ref[2];
    }
  }
  f.close();
}
  

bool read_out_scf(const char*filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return false;
  };
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  bool ret = false;
  while (f.getline(line,MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex=stringLine.find(" CYCLE:");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
       NSCFCYC = 0;
       while (strstr(line,"CYCLE:")!=NULL)
       {
         f.getline(line,MAX_STR_SIZE);
         f.getline(line,MAX_STR_SIZE);
         f >> str; f >> str; f >> str; f >> str;
         f >> SCF_ENERGY[NSCFCYC];
         NSCFCYC++;
         f.getline(line,MAX_STR_SIZE);
         f.getline(line,MAX_STR_SIZE);
         f.getline(line,MAX_STR_SIZE);
       }
       ret = true;
       break;
    }
  }
  f.close();
  return ret;
}

bool read_out_freq(const char*filename)
{
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  string stringLine;

  ifstream f(filename);
  if (!f)
  {
    errmsg(__FILE__,__LINE__,"Could not open file ",0);
    return false;
  }

  NVIB = 0;
  while(f.getline(line, MAX_STR_SIZE))
  {
    stringLine=line;
    int startingIndex=stringLine.find("MODE:");
    if(startingIndex>=0)
    {
      f >> str;
      f >> VIB[NVIB][0];
      f >> str;
      f >> str;
      f >> VIB[NVIB][1];
      f.getline(line, MAX_STR_SIZE);
      f.getline(line, MAX_STR_SIZE);
      f.getline(line, MAX_STR_SIZE);
      f.getline(line, MAX_STR_SIZE);
      for (int iatom=0;iatom<NATOM;iatom++)
      { 
        f >> str;
        ATOMNO[iatom] = symtoan(str);
        f >> COORD[1][iatom][0];
        f >> COORD[1][iatom][1];
        f >> COORD[1][iatom][2];
      }
      save_structure("VIB",NVIB);
      NVIB++;
    }
  }
  f.close();
  get_structure("SCF",0);
  if (NVIB>0) return true;
  else return false;
}

bool read_out_bas(const char*filename)
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
    int sindex=stringLine.find(" GENERAL BASIS SET:");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      for (iatom = 0 ; iatom < NATOM ; iatom++)
      {
        strcpy(basis[iatom],&line[20]);
      }
      f.getline(line,MAX_STR_SIZE);
      stringLine=line;
      sindex=stringLine.find("BASIS");
      while((sindex>=0)&&(sindex<(signed)stringLine.size()))
      {
        int n;
        char str[3];
        n = 0;
        str[n++] = line[1];
        if (line[2]!=' ') str[n++] = line[2];
        str[n++] = '\0';
        n = symtoan(str);
        for (iatom = 0 ; iatom < NATOM ; iatom++)
        {
          if (ATOMNO[iatom]==n)
          {
            if (line[12]==':') 
            {
              strcpy(basis[iatom],&line[14]);
            }
            else
            {
              strcpy(basis[iatom],&line[15]);
            }
          }
        }
        f.getline(line,MAX_STR_SIZE);
        stringLine=line;
        sindex=stringLine.find("BASIS");
      }
      break;
    }
  }
  for (iatom = 0 ; iatom < NATOM ; iatom++)
  {
    if (!loadbasis(iatom,basis[iatom])) 
    {
      f.close();
      return false;
    }
  } 
  return true;
}

bool read_out_mos(const char* filename)
{
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];

  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return false;
  };

  SPHORB = false;
  NORB = 0;
  NBASIS = 0;
  bool OpenShell = false;
  bool gotcycle = false;
  int nblk = 5;
  while(f.getline(line,MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex;
    if (!gotcycle)
    {
      sindex=stringLine.find(" CYCLE:");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) gotcycle = true;
    }
    sindex=stringLine.find(" SPHERICAL ORBITALS WILL BE USED");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    {
      SPHORB = true;
      NSPH = 0;
      for (int iatom=0;iatom<NATOM;iatom++)
        NSPH += mb.ab[iatom].nsao;
    }
    if (NBASIS == 0)
    {
      sindex=stringLine.find(" NUMBER OF ORBITALS:");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) NBASIS = atoi(&line[20]);
    }
    sindex=stringLine.find(" PRINT BASIS SET CONTRACTION");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      nblk = 3;
    cout << line << endl;
    sindex=stringLine.find(" MO COEFFICIENTS");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())&&gotcycle)
    {
      cout << "Encontramos coeficientes" << endl;
      int i,imo,ibas,ll,ul,n;
      bool betas;
      MO mo[5];
      int TOP = NBASIS;
      if (SPHORB) TOP = NSPH;
      strcpy(mo[0].spin,"ALPHA");
      sindex=stringLine.find(" BETA MO COEFFICIENTS");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      { strcpy(mo[0].spin,"BETA");OpenShell = true; betas = true;}
      else betas =false;
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
        for (imo = ll; imo <= ul; imo++) 
	{
	  f>>mo[imo-ll].energy;
	  if (betas) EBEB[0][imo] = mo[imo-ll].energy;
	  else EBEA[0][imo] = mo[imo-ll].energy;
	}
        for (imo = ll; imo <= ul; imo++) f>>mo[imo-ll].occ;
        for (ibas=0;ibas<TOP;ibas++)
        {
          f >> str;
          f >> str;
          f >> str;
          f >> str;
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

void read_out_kmosr(const char* filename)
{
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];

  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };

  SPHORB = false;
  KNORB = 0;
  NBASIS = 0;
  int k;
  bool OpenShell = false;
  bool gotcycle = false;
  bool gotk = false;
  int nblk = 5;
  while(f.getline(line,MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex;
    if (!gotcycle)
    {
      sindex=stringLine.find(" CYCLE:");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) gotcycle = true;
    }
    sindex=stringLine.find(" K INDEX");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    {
      k = atoi(&line[9]);
      if (k==1) gotk = true;
    }
    sindex=stringLine.find(" SPHERICAL ORBITALS WILL BE USED");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    {
      SPHORB = true;
      NSPH = 0;
      for (int iatom=0;iatom<NATOM;iatom++)
        NSPH += mb.ab[iatom].nsao;
    }
    if (NBASIS == 0)
    {
      sindex=stringLine.find(" NUMBER OF ORBITALS:");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) NBASIS = atoi(&line[20]);
    }
    sindex=stringLine.find(" PRINT BASIS SET CONTRACTION");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      nblk = 3;
    sindex=stringLine.find(" MO COEFFICIENTS REAL");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())&&gotcycle&&gotk)
    {
      int i,imo,ibas,ll,ul,n;
      bool betas;
      MO mo[5];
      int TOP = NBASIS;
      if (SPHORB) TOP = NSPH;
      strcpy(mo[0].spin,"ALPHA");
      sindex=stringLine.find(" BETA MO COEFFICIENTS REAL");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      { strcpy(mo[0].spin,"BETA");OpenShell = true; betas = true;}
      else betas =false;
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
        for (imo = ll; imo <= ul; imo++) 
	{
	  f>>mo[imo-ll].energy;
	  if (betas) KEBEB[0][imo] = mo[imo-ll].energy;
	  else KEBEA[0][imo] = mo[imo-ll].energy;
	}
        for (imo = ll; imo <= ul; imo++) f>>mo[imo-ll].occ;
        for (ibas=0;ibas<TOP;ibas++)
        {
          f >> str;
          f >> str;
          f >> str;
          f >> str;
          for (imo = ll; imo <= ul; imo++) 
          {
            f >> mo[imo-ll].c[ibas];
          }
        }
        for (imo = ll; imo <= ul; imo++)
        {
          KMOSR[KNORB] = mo[imo-ll];
          KNORB++;
        }
        ll = ul+1;
        ul = MIN(ll + nblk-1,TOP-1);
      }
      break;
    }
  }

  f.close();


  return;
}

void read_out_kmosi(const char* filename)
{
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];

  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };

  SPHORB = false;
  KNORB = 0;
  NBASIS = 0;
  int k;
  bool OpenShell = false;
  bool gotcycle = false;
  bool gotk = false;
  int nblk = 5;
  while(f.getline(line,MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex;
    if (!gotcycle)
    {
      sindex=stringLine.find(" CYCLE:");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) gotcycle = true;
    }
    sindex=stringLine.find(" K INDEX");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    {
      k = atoi(&line[9]);
      if (k==1) gotk = true;
    }
    sindex=stringLine.find(" SPHERICAL ORBITALS WILL BE USED");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    {
      SPHORB = true;
      NSPH = 0;
      for (int iatom=0;iatom<NATOM;iatom++)
        NSPH += mb.ab[iatom].nsao;
    }
    if (NBASIS == 0)
    {
      sindex=stringLine.find(" NUMBER OF ORBITALS:");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) NBASIS = atoi(&line[20]);
    }
    sindex=stringLine.find(" PRINT BASIS SET CONTRACTION");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      nblk = 3;
    sindex=stringLine.find(" MO COEFFICIENTS IMAGINARY");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())&&gotcycle&&gotk)
    {
      int i,imo,ibas,ll,ul,n;
      bool betas;
      MO mo[5];
      int TOP = NBASIS;
      if (SPHORB) TOP = NSPH;
      strcpy(mo[0].spin,"ALPHA");
      sindex=stringLine.find(" BETA MO COEFFICIENTS IMAGINARY");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      { strcpy(mo[0].spin,"BETA");OpenShell = true; betas = true;}
      else betas =false;
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
        for (imo = ll; imo <= ul; imo++) 
	{
	  f>>mo[imo-ll].energy;
	  if (betas) KEBEB[0][imo] = mo[imo-ll].energy;
	  else KEBEA[0][imo] = mo[imo-ll].energy;
	}
        for (imo = ll; imo <= ul; imo++) f>>mo[imo-ll].occ;
        for (ibas=0;ibas<TOP;ibas++)
        {
          f >> str;
          f >> str;
          f >> str;
          f >> str;
          for (imo = ll; imo <= ul; imo++) 
          {
            f >> mo[imo-ll].c[ibas];
          }
        }
        for (imo = ll; imo <= ul; imo++)
        {
          KMOSI[KNORB] = mo[imo-ll];
          KNORB++;
        }
        ll = ul+1;
        ul = MIN(ll + nblk-1,TOP-1);
      }
      break;
    }
  }

  f.close();


  return;
}

void read_deMon_out_fukui(const char* filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  int ibas,jbas,k,ll,n,nk,ul;
  while (f.getline(line,MAX_STR_SIZE))
  {
    int sindex = -1;
    int rindex = -1;
    string stringLine=line;
    sindex=stringLine.find(" NUMBER OF ORBITALS:");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    {
      NBASIS = atoi(&line[20]);
      continue;
    }
    sindex=stringLine.find(" FUKUI^- MATRIX");
    if (sindex<0) sindex=stringLine.find(" FUKUI^+ MATRIX");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
#define NBLK 4
      int nblk = NBLK;
      double *F;
      int SQBASIS1 = NBASIS*NBASIS;
      int SQBASIS2 = (NBASIS*(NBASIS+1))/2;
      F = (double*)malloc(SQBASIS1*sizeof(double));
      ll = 0;
      ul = MIN(ll + nblk-1,NBASIS-1);
      while(ul>=ll)
      {
        for (jbas = ll; jbas <= ul; jbas++)
        {
          f>>n;
          if (n-jbas!=1) break;
        }
        ul = jbas-1;
        for (ibas=0;ibas<NBASIS;ibas++)
        {
          f>>n;
   	  f>>str;
	  f>>str;
	  f>>str;
          if (n-ibas!=1) {
            cout << "UPS: something went wrong on reading"<<endl;
            cout << str << endl;
            cout << "n= "<<n<<" ibas = "<<ibas<<endl;
            exit(0);
          }
          for (jbas = ll; jbas <= ul; jbas++) 
          {
            f >> F[jbas*NBASIS+ibas];
          }
        }
        ll = ul+1;
        ul = MIN(ll + nblk-1,NBASIS-1);
      }

      rindex=stringLine.find(" FUKUI^- MATRIX");
      if((rindex>=0)&&(rindex<(signed)stringLine.size()))
      {
        FUKUIL_MATRIX = (double*) malloc(SQBASIS2*sizeof(double));
        k = 0;
        for (ibas=0;ibas<NBASIS;ibas++)
          for (jbas=ibas;jbas<NBASIS;jbas++)
          {
            FUKUIL_MATRIX[k] = F[jbas*NBASIS+ibas];
            if (jbas!=ibas)
              FUKUIL_MATRIX[k] = 2.0*FUKUIL_MATRIX[k];
            k++;
          }
      }

      rindex=stringLine.find(" FUKUI^+ MATRIX");
      if((rindex>=0)&&(rindex<(signed)stringLine.size()))
      {
        FUKUIR_MATRIX = (double*) malloc(SQBASIS2*sizeof(double));
        k = 0;
        for (ibas=0;ibas<NBASIS;ibas++)
          for (jbas=ibas;jbas<NBASIS;jbas++)
          {
            FUKUIR_MATRIX[k] = F[jbas*NBASIS+ibas];
            if (jbas!=ibas)
              FUKUIR_MATRIX[k] = 2.0*FUKUIR_MATRIX[k];
            k++;
          }
      }

      free(F);
    }
  }
  f.close();
}

void read_deMon_out_response(const char* filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  int ibas,jbas,k,ll,n,nk,ul;
  while (f.getline(line,MAX_STR_SIZE))
  {
    int sindex = -1;
    int rindex = -1;
    string stringLine=line;
    sindex=stringLine.find(" NUMBER OF ORBITALS:");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    {
      NBASIS = atoi(&line[20]);
      continue;
    }
    sindex=stringLine.find(" DENSITY MATRIX RESPONSE");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
#define NBLK 4
      int nblk = NBLK;
      double *F;
      int SQBASIS1 = NBASIS*NBASIS;
      int SQBASIS2 = (NBASIS*(NBASIS+1))/2;
      F = (double*)malloc(SQBASIS1*sizeof(double));
      ll = 0;
      ul = MIN(ll + nblk-1,NBASIS-1);
      while(ul>=ll)
      {
        for (jbas = ll; jbas <= ul; jbas++)
        {
          f>>n;
          if (n-jbas!=1) break;
        }
        ul = jbas-1;
        for (ibas=0;ibas<NBASIS;ibas++)
        {
          f>>n;
   	  f>>str;
	  f>>str;
	  f>>str;
          if (n-ibas!=1) {
            cout << "UPS: something went wrong on reading"<<endl;
            cout << str << endl;
            cout << "n= "<<n<<" ibas = "<<ibas<<endl;
            exit(0);
          }
          for (jbas = ll; jbas <= ul; jbas++) 
          {
            f >> F[jbas*NBASIS+ibas];
          }
        }
        ll = ul+1;
        ul = MIN(ll + nblk-1,NBASIS-1);
      }

      rindex=stringLine.find(" DENSITY MATRIX RESPONSE X");
      if((rindex>=0)&&(rindex<(signed)stringLine.size()))
      {
        DMR1 = (double*) malloc(SQBASIS2*sizeof(double));
        k = 0;
        for (ibas=0;ibas<NBASIS;ibas++)
          for (jbas=ibas;jbas<NBASIS;jbas++)
          {
            DMR1[k] = F[jbas*NBASIS+ibas];
            if (jbas!=ibas)
              DMR1[k] = 2.0*DMR1[k];
            k++;
          }
      }
      rindex=stringLine.find(" DENSITY MATRIX RESPONSE Y");
      if((rindex>=0)&&(rindex<(signed)stringLine.size()))
      {
        DMR2 = (double*) malloc(SQBASIS2*sizeof(double));
        k = 0;
        for (ibas=0;ibas<NBASIS;ibas++)
          for (jbas=ibas;jbas<NBASIS;jbas++)
          {
            DMR2[k] = F[jbas*NBASIS+ibas];
            if (jbas!=ibas)
              DMR2[k] = 2.0*DMR2[k];
            k++;
          }
      }
      rindex=stringLine.find(" DENSITY MATRIX RESPONSE Z");
      if((rindex>=0)&&(rindex<(signed)stringLine.size()))
      {
        DMR3 = (double*) malloc(SQBASIS2*sizeof(double));
        k = 0;
        for (ibas=0;ibas<NBASIS;ibas++)
          for (jbas=ibas;jbas<NBASIS;jbas++)
          {
            DMR3[k] = F[jbas*NBASIS+ibas];
            if (jbas!=ibas)
              DMR3[k] = 2.0*DMR3[k];
            k++;
          }
      }

      free(F);
    }
  }
  f.close();
}

void read_AZD(const char*filename)
{
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  string stringLine;

  ifstream f(filename);
  if (!f)
  {
    errmsg(__FILE__,__LINE__,"Could not open file ",0);
    return;
  }

  while(f.getline(line, MAX_STR_SIZE))
  {
    stringLine=line;
    int startingIndex=stringLine.find("AZD FIELD");
    if(startingIndex>=0)
    {
      f >> str;
      cout << str << endl;
    }
  }
  f.close();
}

void read_deMon_ebe(const char* filename)
{
  for (int i=0;i<NBASIS;i++)
  {
    for (int j=0;j<2;j++)
    {
      EBEA[j][i] = -9999.99/TOEV(1.0);     // Undefined value
      EBEB[j][i] = -9999.99/TOEV(1.0);     // Undefined value
      PSA[j][i] = -1.0;
      PSB[j][i] = -1.0;
    }
  }

  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };
  char line[MAX_STR_SIZE];
  char str[MAX_STR_SIZE];
  int ibas,jbas,k,ll,n,nk,ul;
  bool ebealpha;
  while (f.getline(line,MAX_STR_SIZE))
  {
    int sindex = -1;
    int rindex = -1;
    string stringLine=line;
    sindex=stringLine.find(" DIAGONAL GENERALIZED ONE-PARTICLE PROPAGATOR");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    {
      f.getline(line,MAX_STR_SIZE);
      stringLine=line;
      sindex=stringLine.find(" ALPHA");
      if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
	ebealpha = true;
      else
	ebealpha = false;
      f.getline(line,MAX_STR_SIZE);
      f.getline(line,MAX_STR_SIZE);
      f.getline(line,MAX_STR_SIZE);
      f.getline(line,MAX_STR_SIZE);
      stringLine=line;
      sindex=stringLine.find(" ORBITAL NUMBER:");
      while ((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      {
        int n = atoi(&line[16]);
	double ebe;
        n--;
	// KT
        f >> str;
	f >> str;
        f >> ebe;
	if (ebealpha) 
	{
	  EBEA[0][n] = ebe/TOEV(1.0);
	  f >> PSA[0][n];
	}
	else
	{
	  EBEB[0][n] = ebe/TOEV(1.0);
	  f >> PSB[0][n];
	}
        f.getline(line,MAX_STR_SIZE);
	// EP2
        f >> str;
	f >> str;
        f >> ebe;
	if (ebealpha) 
	{
	  EBEA[1][n] = ebe/TOEV(1.0);
	  f >> PSA[1][n];
	}
	else
	{
	  EBEB[1][n] = ebe/TOEV(1.0);
	  f >> PSB[1][n];
	}
        f.getline(line,MAX_STR_SIZE);
        f.getline(line,MAX_STR_SIZE);
        stringLine=line;
        sindex=stringLine.find(" ORBITAL NUMBER:");
      }
      continue;
    }
  }
  f.close();
}

void read_deMon_out_nffez(const char*filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };
  
  char line[MAX_STR_SIZE];
  int iatom,iindex,sindex;

  iindex = 0;
  while (f.getline(line,MAX_STR_SIZE))
  {
    string stringLine=line;

    sindex=stringLine.find("RESPONSE IN ATOMIC FORCES"); 
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
    {
      iindex++;
      for(iatom=0;iatom<NATOM;iatom++)
      {
        f >> line;
	if (iindex==1)
        {
          f >> ATOMNFFEX[iatom][0];
          f >> ATOMNFFEX[iatom][1];
          f >> ATOMNFFEX[iatom][2];
	}
	else if (iindex==2)
        {
          f >> ATOMNFFEY[iatom][0];
          f >> ATOMNFFEY[iatom][1];
          f >> ATOMNFFEY[iatom][2];
	}
	else if (iindex==3)
        {
          f >> ATOMNFFEZ[iatom][0];
          f >> ATOMNFFEZ[iatom][1];
          f >> ATOMNFFEZ[iatom][2];
	}
      }
      if (iindex==3||iindex==NATOM) break;
    }  
  }
  f.close();
}

void read_deMon_out(const char*filename)
{
  if (!read_deMon_out_geo(filename)) 
    errmsg(__FILE__,__LINE__,"GEOMETRY loading failed",0);

  read_out_cps(filename);

  (void) read_out_scf(filename);


 // (void) read_out_freq(filename);


  read_AZD(filename);


  if (!read_out_bas(filename))
    errmsg(__FILE__,__LINE__,"BASIS loading failed",0);
  else 
  {
    if (!read_out_mos(filename))
      errmsg(__FILE__,__LINE__,"Basis found but MOs missed",0);

    read_out_kmosr(filename);
    read_out_kmosi(filename);

    // Load Fukui functions
    read_deMon_out_fukui(filename);

    // Load electrical response 
    read_deMon_out_response(filename);
  }

  read_deMon_ebe(filename);
  read_deMon_out_nffez(filename);

}
