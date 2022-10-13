// *****************************************************************
// Purpose: Basis set related routines
//
// Roberto Flores-Moreno
// *****************************************************************

#include <math.h>

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <param.h>
#include <types.h>
#include <macros.h>
#include <global.h>
#include <fcns.h>

void printnormbas(MoleculeBasis* mb)
{
  AtomBasis ab;
  ShellBasis sb;
  int iatom,ico,ishl,igto,jgto,lx,ly,lz;
  double factor,norm;

  for ( iatom = 0 ; iatom < NATOM ; iatom++ )
  {
    ab = mb->ab[iatom];
    for ( ishl = 0 ; ishl < ab.nshl ; ishl++ )
    {
      sb = ab.sb[ishl];
      norm = pow(2.0,double(sb.l))*pow(2.0/M_PI,0.75);
      // Make sure to apply this only once
      for ( igto = 0 ; igto < sb.k ; igto++ )
      {
        sb.d[igto] = norm*pow(sb.z[igto],(2.0*sb.l+3.0)*0.25)*sb.d[igto];
      }
      factor = pow(M_PI,1.5)/pow(2.0,(double)sb.l);
      norm = 0.0;
      for ( igto = 0 ; igto < sb.k ; igto++ )
      {
        for ( jgto = 0 ; jgto < sb.k ; jgto++ )
        {
          norm += sb.d[igto]*sb.d[jgto]*factor/
                  pow(sb.z[igto]+sb.z[jgto],sb.l+1.5);
        }
      }
      ico = 0;
      for ( lx = sb.l ; lx >= 0 ; lx-- )
      {
        for ( ly = sb.l - lx ; ly >= 0 ; ly-- )
        {
          lz = sb.l - lx - ly;
          sb.ncsto[ico] = 1.0/sqrt(norm*
                                   double_factorial(2*lx-1)*
                                   double_factorial(2*ly-1)*
                                   double_factorial(2*lz-1));
          ico++;
        }
      }
      ab.sb[ishl] = sb;
    }
    mb->ab[iatom] = ab;
  }
}


// Normalization of basis sets
void normbas(MoleculeBasis* mb)
{
  AtomBasis ab;
  ShellBasis sb;
  int iatom,ico,ishl,igto,jgto,lx,ly,lz;
  double factor,norm;

  for ( iatom = 0 ; iatom < NATOM ; iatom++ )
  {
    ab = mb->ab[iatom];
    for ( ishl = 0 ; ishl < ab.nshl ; ishl++ )
    {
      sb = ab.sb[ishl];
      norm = pow(2.0,double(sb.l))*pow(2.0/M_PI,0.75);
      // Make sure to apply this only once
      for ( igto = 0 ; igto < sb.k ; igto++ )
      {
        sb.d[igto] = norm*pow(sb.z[igto],(2.0*sb.l+3.0)*0.25)*sb.d[igto];
      }
      factor = pow(M_PI,1.5)/pow(2.0,(double)sb.l);
      norm = 0.0;
      for ( igto = 0 ; igto < sb.k ; igto++ )
      {
        for ( jgto = 0 ; jgto < sb.k ; jgto++ )
        {
          norm += sb.d[igto]*sb.d[jgto]*factor/
                  pow(sb.z[igto]+sb.z[jgto],sb.l+1.5);
        }
      }
      ico = 0;
      for ( lx = sb.l ; lx >= 0 ; lx-- )
      {
        for ( ly = sb.l - lx ; ly >= 0 ; ly-- )
        {
          lz = sb.l - lx - ly;
          sb.ncsto[ico] = 1.0/sqrt(norm*
                                   double_factorial(2*lx-1)*
                                   double_factorial(2*ly-1)*
                                   double_factorial(2*lz-1));
          ico++;
        }
      }
      ab.sb[ishl] = sb;
    }
    mb->ab[iatom] = ab;
  }
}
//
// Purpose: Generate the orbital transformation matrices for the
//          transformation of orbitals from the Cartesian representation 
//          to the spherical representation.
//
// Lit.: H.B. Schlegel, M.J. Frisch, IJQC. 54, 83 (1995)
//
void gbtc( void )
{
  int expo,i,ic,is,j,k,l,lx,ly,lz,m,ma,ncar,nsph;
  double s,s1,s2;

  for (l=0;l<=MAXLBAS;l++)
  {
    nsph = 2*l+1;
    ncar = (l+1)*(l+2)/2;
    for (is=0;is<nsph;is++)
    {
      for (ic=0;ic<ncar;ic++)
      {
        CTOSTM[l][is][ic] = 0.0;
      }
    }
  }

  for (l=0;l<=MAXLBAS;l++)
  {
    for (lx=l;lx>=0;lx--)
    {
      for (ly=l-lx;ly>=0;ly--)
      {
        lz = l - lx - ly;
        ic = lz + (ly+lz)*(ly+lz+1)/2;
        for (m=-l;m<=l;m++)
        {
          is = l + m;
          ma = ABS(m);
          j = lx + ly - ma;
          if (j>=0&&(MOD(j,2)==0))
          {
            j /= 2;
            s1 = 0.0;
            for (i=0;i<=(l-ma)/2;i++)
            {
              s2 = 0.0;
              for (k=0;k<=j;k++)
              {
                if (((m<0)&&(MOD(ABS(ma-lx),2)==1))||
                    ((m>0)&&(MOD(ABS(ma-lx),2)==0))) 
                {
                  expo = (ma - lx + 2*k)/2;
                  s = pow(-1.0,double(expo))*sqrt(2.0);
                }
                else if ((m==0)&&(MOD(lx,2)==0))
                {
                  expo = k - lx/2;
                  s = pow(-1.0,double(expo));
                }
                else s = 0.0;
                s2 += noverk(j,k)*noverk(ma,(lx-2*k))*s;
              }
              s1 += noverk(l,i)*noverk(i,j)*pow(-1.0,double(i))*
                    factorial(2*l-2*i)/
                    factorial(l-ma-2*i)*s2;
            }
            CTOSTM[l][is][ic] = sqrt((factorial(2*lx)*factorial(2*ly)*
                                      factorial(2*lz)*factorial(l)*
                                      factorial(l-ma))/(factorial(lx)*
                                      factorial(ly)*factorial(lz)*
                                      factorial(2*l)*factorial(l+ma)))*
                                s1/(pow(2.0,double(l))*factorial(l));
          }
          else CTOSTM[l][is][ic] = 0.0;
        }
      }
    }
  }
}

static double gtorad(double z,int l,double tol)
{
  // A couple of parameters here
  int maxit = 100;
  double eps = 1.0e-12;

  int i;
  double delta,dg,g,rl,zr;

  double ldt,rad;

  ldt = log(1.0/tol);
  rad = ldt/z;
  if (l!=0)
  {
    rad = sqrt(rad);
    rl = (double)l;
    g = sqrt(rl/2.0);
    if (g>rad) rad = 0.5*(rad+g);
    for(i=0;i<maxit;i++)
    {
      //j = i;
      zr = z*rad;
      g = ldt + rl*log(rad) - zr*rad;
      if (ABS(rad)<eps) 
      {
        rad = 100.0;
        break;
      }
      dg = rl/rad - 2.0*zr;
      if (ABS(dg)<=eps) delta = 0.0; 
      else delta = g/dg; 
      rad = rad-delta;
      if (g<=0.0) break;
    } 
    rad = rad*rad;
  }
  return rad;
}

// Evaluate effective shell radius
void evaluate_atomrad(void)
{
  double tol = 1.0e-10;

  int i,iatom,ishl;
  double atomrad,rad,z;

  for (iatom=0;iatom<NATOM;iatom++)
  {
    AtomBasis ab = mb.ab[iatom];
    atomrad = 0.0;
    for (ishl=0;ishl<ab.nshl;ishl++)
    {
      ShellBasis sb = ab.sb[ishl];
      z = sb.z[0];
      for (i=1;i<sb.k;i++)
        if (z>sb.z[i]) 
          z = sb.z[i];
      rad = gtorad(z,sb.l,tol);
      sb.maxrad = rad;
      ab.sb[ishl] = sb;
      if (atomrad<rad) atomrad = rad;
    }
    ab.maxrad = atomrad;
    mb.ab[iatom] = ab;
  }
}

// Trasform spherical MO coefficients to cartesian MO coefficients.
void sc2cc(void)
{
  int floor_car,floor_sph,iatom,ic,imo,is,ishl,l,ncar,nsph;
  AtomBasis ab;
  ShellBasis sb;

  for (imo=0;imo<NORB;imo++)
  {
    // Load spherical
    MO mosph = MOS[imo];

    // Evaluate cartesian
    MO mocar = mosph;
    floor_car = 0;
    floor_sph = 0;
    for (iatom=0;iatom<NATOM;iatom++)
    {
      ab = mb.ab[iatom];
      for (ishl=0;ishl<ab.nshl;ishl++)
      {
        sb = ab.sb[ishl];
        l = sb.l;
        ncar = (l+1)*(l+2)/2;
        nsph = 2*l+1;
        for (ic=0;ic<ncar;ic++)
        {
          mocar.c[floor_car+ic] = 0.0;
          for (is=0;is<nsph;is++)
          {
            mocar.c[floor_car+ic] += CTOSTM[l][is][ic]*mosph.c[floor_sph+is];
          }
        }
        floor_car += ncar;
        floor_sph += nsph;
      }
    }
    // Save new coefficients
    MOS[imo] = mocar;
  }
}

// Interfaz routine
void setup_basis(void)
{
  // Evaluate NBASIS and NSPH
  NBASIS = 0;
  NSPH = 0;
  for (int iatom = 0; iatom < NATOM ; iatom++)
  {
    NBASIS += mb.ab[iatom].ncao;
    NSPH += mb.ab[iatom].nsao;
  }

  //Normalize basis sets
//  normbas(&mb);

  // Cartesian to spherical transformation
  if (SPHORB) 
  {
    gbtc();
    sc2cc();
    SPHORB = false;
  }

  // Screening radius (squared)
  evaluate_atomrad();
printnormbas(&mb);
}

bool loadbasis(int atom, char* basisname)
{
  char line[MAX_STR_SIZE];

  if (strncasecmp(basisname,"READ IN",7)==0)  
  {
    errmsg(__FILE__,__LINE__,"Please, do not use READ IN basis option",0);
    return false;
  }

  const char* local_basis = "BASIS";
  ifstream f(local_basis);
  if ( f.fail() )
  {
    const char* master_basis = "/usr/local/deMon/BASIS";
    char msg[MAX_STR_SIZE];
    sprintf(msg,"Error opening %s, Trying %s",local_basis,master_basis);
    //RFM errmsg(__FILE__,__LINE__,msg,0);
    f.open(master_basis);
    if (f.fail() )
    {
      sprintf(msg,"Error opening %s",master_basis);
      errmsg(__FILE__,__LINE__,msg,0);
      return false;
    }
  };
  string flag = "O-";
  flag += ELEMENT_NAME[ATOMNO[atom]];
  while(f.getline(line,MAX_STR_SIZE))
  {
    string stringLine=line;
    int sindex=stringLine.find(flag);
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      if (strstr(line,basisname)!=NULL)
      {
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
        mb.ab[atom] = ab;
        return true;
      }
    }
  }
  f.close();
  return false;
}

