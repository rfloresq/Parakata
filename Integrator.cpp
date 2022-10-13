// Evaluation of gaussian AO integrals 
// Roberto Flores-Moreno, 2008
//
// Lit: S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
//      M. Head-Gordon and J. A. Pople, J. Chem. Phys. 89, 5777 (1988)


#include <iostream>

#include <types.h>
#include <Integrator.h>
#include <Math.h>

using namespace std;

int gaop[MAXLBAS+3][MAXLBAS+3][MAXLBAS+3];
int raop[MAXLBAS+3][MAXLBAS+3][MAXLBAS+3];
double ftab[2*MAXLBAS+6+1][MAX_TAB_GAM+1];

// Initialize integration routines
Integrator::Integrator()
{
  int gcount,l,lx,ly,lz,rcount;

  gcount = 0;
  for (l=0;l<=MAXLBAS+2;l++)
  {
    rcount = 0;
    for (lx=l;lx>=0;lx--)
    {
      for (ly=l-lx;ly>=0;ly--)
      {
        lz = l - lx - ly;
        raop[lx][ly][lz] = rcount;
        gaop[lx][ly][lz] = gcount;
        rcount++;
        gcount++;
      }
    }
  }
  // Initialize Gamma Function
  int i,j,k;
  int nitermax = 30;
  int nmax = 2*MAXLBAS+6;
  double eps = 1.0e-15;
  double bessel,expterm,prefak,preterm,produkt,serie,sumterm,term,ttab;
  double r[nitermax+11];

  for (i=0;i<=nmax;i++)
    ftab[i][0] = 1.0/(2*i+1);
  for (i=1;i<=MAX_TAB_GAM;i++)
  {
    ttab = double(i)/10.0;
    r[nitermax+10] = 0.0;
    for (j=1;j<=nitermax+9;j++)
      r[nitermax+10-j] = -ttab/(4*(nitermax+10-j) + 2.0- ttab*r[nitermax+11-j]);
    bessel = (2*sinh(ttab/2))/ttab;
    prefak = exp(-ttab/2)*bessel;
    term = 1.0;
    serie = prefak*(1.0/(2.0*nmax + 1.0));
    for (k=1;k<=nitermax;k++)
    {
      preterm = (2.0*k + 1.0)/(2.0*nmax + 1.0);
      term = term*(2.0*nmax - 2.0*k + 1.0)/(2.0*nmax + 2.0*k + 1.0);
      produkt = 1.0;
      for (l=1;l<=k;l++)
        produkt = produkt*r[l];
      sumterm = prefak*preterm*term*produkt;
      if (ABS(sumterm)<=eps) goto TABOK;
      else serie = serie + sumterm;
    }
    cout << "Fatal integrals_initialize, contact support" << endl;
    exit(EXIT_FAILURE);
TABOK:
    ftab[nmax][i] = serie;
    expterm = exp(-ttab);
    for (j=1;j<=nmax;j++)
      ftab[nmax-j][i] = 1.0/(2*(nmax-j)+1)*(2*ttab*ftab[nmax+1-j][i] + expterm);
  }
}

/*

// Evaluate overlap integrals
void Integrator::Overlap(double *ra,double *rb,GaussianShell *sa,
 GaussianShell *sb, double **ib,int nsd, bool apply_normalization)
{
  int ax,ay,az,daop,i,j,k,lab;
  double zp,zp2,xi,r2;
  double r[3],rr[3];

  // Initialize
  lab = sa->l + sb->l;
  r[0] = rb[0]-ra[0];
  r[1] = rb[1]-ra[1];
  r[2] = rb[2]-ra[2];
  r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];

  daop = PSDAOPointer(lab);
  double v[daop];
  double **h;

  h = new double*[daop];
  for (i=0;i<daop;i++)
    h[i] = new double[daop];
  for (i=0;i<daop;i++)
    for (j=0;j<daop;j++)
      h[i][j] = 0.0;

  // Main loop on primitives
  for (i=0;i<(signed)sa->z.size();i++)
  {
    for (j=0;j<(signed)sb->z.size();j++)
    {
      zp = sa->z[i] + sb->z[j];
      zp2 = 2.0*zp;
      xi = sa->z[i]*sb->z[j]/zp;
      rr[0] = sb->z[j]/zp*r[0];
      rr[1] = sb->z[j]/zp*r[1];
      rr[2] = sb->z[j]/zp*r[2];
      // s-s type
      v[gaop[0][0][0]] = pow(GREEN_PI/zp,1.5)*exp(-xi*r2);
      // vertical recurrence relation x
      if (lab>0) v[gaop[1][0][0]] = rr[0]*v[gaop[0][0][0]];
      for (ax=2;ax<=lab;ax++)
      {
        v[gaop[ax][0][0]] = rr[0]*v[gaop[ax-1][0][0]] +             
                    (double(ax-1)/zp2)*v[gaop[ax-2][0][0]];
      }
      // vertical recurrence relation y
      for (ax=0;ax<=lab;ax++)
      {
        if (lab-ax>0) v[gaop[ax][1][0]] = rr[1]*v[gaop[ax][0][0]];
        for (ay=2;ay<=lab-ax;ay++)
        {
          v[gaop[ax][ay][0]] = rr[1]*v[gaop[ax][ay-1][0]] +               
                            (double(ay-1)/zp2)*v[gaop[ax][ay-2][0]];
        }
      }
      // vertical recurrence relation z
      for (ax=0;ax<=lab;ax++)
      {
        for (ay=0;ay<=lab-ax;ay++)
        {
          if (lab-ax-ay>0) v[gaop[ax][ay][1]] = rr[2]*v[gaop[ax][ay][0]];
          for (az=2;az<=lab-ax-ay;az++)
          {
            v[gaop[ax][ay][az]] = rr[2]*v[gaop[ax][ay][az-1]] +  
                               (double(az-1)/zp2)*v[gaop[ax][ay][az-2]];
          }
        }
      }
      // contraction
      double factor = (pow(sa->z[i],nsd))*sa->d[i]*sb->d[j];
      for (k=0;k<daop;k++)
        h[k][0] += factor*v[k];
    }
  }
  // horizontal recurrence relation
  HRR(sa->l,sb->l,r,h);
  NormalizeTwoIndex(sa,sb,ib,h,apply_normalization);
  

  for (i=0;i<daop;i++)
    delete[] h[i];
  delete[] h;
}
*/

// save and normalize
void Integrator::NormalizeTwoIndex(ShellBasis *sa, ShellBasis *sb, 
double **ib,double **h,bool apply_normalization)
{
  int ax,ay,az,bx,by,bz,i,ii,j;

  for (ax=0;ax<=sa->l;ax++)
  {
    for (ay=0;ay<=sa->l-ax;ay++)
    {
      az = sa->l - ax - ay;
      i = raop[ax][ay][az];
      ii = gaop[ax][ay][az];
      for (bx=0;bx<=sb->l;bx++)
      {
        for (by=0;by<=sb->l-bx;by++)
        {
          bz = sb->l-bx-by;
          j = raop[bx][by][bz];
          if (apply_normalization)
            ib[i][j] = sa->ncsto[i]*sb->ncsto[j]*h[ii][gaop[bx][by][bz]];
          else ib[i][j] = h[ii][gaop[bx][by][bz]];
        }
      }
    }
  }
}

// Calculation of integrals over horizontal recurrence relation.
void Integrator::HRR(int la,int lb,double* r,double **h)
{
  int ax,ay,az,bx,by,bz,lab,l1,l2;

  lab = la + lb;
  for (l1=1;l1<=lb;l1++)
  {
    for (l2=la;l2<=lab-l1;l2++)
    {
      for (ax=0;ax<=l2;ax++)
      {
        for (ay=0;ay<=l2-ax;ay++)
        {
          az = l2 - ax - ay;
          // bz up
          h[gaop[ax][ay][az]][gaop[0][0][l1]] =  
          h[gaop[ax][ay][az+1]][gaop[0][0][l1-1]] - 
          r[2]*h[gaop[ax][ay][az]][gaop[0][0][l1-1]];
          // by up
          for (by=1;by<=l1;by++)
          {
            bz = l1 - by;
            h[gaop[ax][ay][az]][gaop[0][by][bz]] =
            h[gaop[ax][ay+1][az]][gaop[0][by-1][bz]] - 
            r[1]*h[gaop[ax][ay][az]][gaop[0][by-1][bz]];
          }
          // bx up
          for (bx=1;bx<=l1;bx++)
          {
            for (by=0;by<=l1-bx;by++)
            {
              bz = l1 - bx - by;
              h[gaop[ax][ay][az]][gaop[bx][by][bz]] = 
              h[gaop[ax+1][ay][az]][gaop[bx-1][by][bz]] -
              r[0]*h[gaop[ax][ay][az]][gaop[bx-1][by][bz]];
            }
          }
        } 
      }
    }
  }
}

/*
// Evaluate kinetic energy integrals
void Integrator::Kinetic(double *ra,double *rb,GaussianShell *sa,
 GaussianShell *sb, double **ib,int nsd, bool apply_normalization)
{
  int ap,ax,ay,az,i,j,na,nb;
  double **ob;
  GaussianShell *wsa;

  int l = GREEN_MAX(sa->l,sb->l) + nsd + 2;
  int dover = ((l+1)*(l+2))/2;
  na = ((sa->l+1)*(sa->l+2))/2;
  nb = ((sb->l+1)*(sb->l+2))/2;

  ob = new double*[dover];
  for(i=0;i<dover;i++)
    ob[i] = new double[dover];
  
  wsa = new GaussianShell();
  for (i=0;i<(signed)sa->z.size();i++)
  {
    wsa->z.push_back( sa->z[i] );
    wsa->d.push_back( sa->d[i] );
  }

  wsa->l = sa->l;
  Overlap(ra,rb,wsa,sb,ob,1+nsd,false);

  for (i=0;i<na;i++)
    for (j=0;j<nb;j++)
      ib[i][j] = (2*sa->l+3)*ob[i][j];

  wsa->l = sa->l + 2;
  Overlap(ra,rb,wsa,sb,ob,2+nsd,false);

  for (ax=0;ax<=sa->l;ax++)
  {
    for (ay=0;ay<=sa->l-ax;ay++)
    {
      az = sa->l - ax - ay;
      ap = raop[ax][ay][az];
      for (j=0;j<nb;j++)
      {
        ib[ap][j] += - 2.0*ob[raop[ax+2][ay][az]][j]       
                     - 2.0*ob[raop[ax][ay+2][az]][j]      
                     - 2.0*ob[raop[ax][ay][az+2]][j];              
      }
    }
  }
  if (sa->l>1)
  {
    wsa->l = sa->l - 2;
    Overlap(ra,rb,wsa,sb,ob,0+nsd,false);

    for (ax=0;ax<=sa->l;ax++)
    {
      for (ay=0;ay<=sa->l-ax;ay++)
      {
        az = sa->l - ax - ay;
        ap = raop[ax][ay][az];
        for (j=0;j<nb;j++)
        {
          if (ax>1) ib[ap][j] += -(ax*(ax-1)/2.0)*ob[raop[ax-2][ay][az]][j];
          if (ay>1) ib[ap][j] += -(ay*(ay-1)/2.0)*ob[raop[ax][ay-2][az]][j]; 
          if (az>1) ib[ap][j] += -(az*(az-1)/2.0)*ob[raop[ax][ay][az-2]][j]; 
        }
      }
    }
  }
  for(i=0;i<dover;i++)
    delete[] ob[i];
  delete[] ob;

  delete wsa;

  if (apply_normalization)
  {
    for (i=0;i<na;i++)
      for (j=0;j<nb;j++)
        ib[i][j] *= sa->ncsto[i]*sb->ncsto[j];
  }
}
*/

// Evaluate nuclear atraction integrals
void Integrator::Core(double *ra,double *rb,double *rc,ShellBasis *sa,
ShellBasis *sb, double **ib)
{
  int ax,ay,az,daop,i,j,k,lab,ma,mm;
  double factor,r2,sf,t,xi,zp;
  double p[3],q[3],r[3];
  double f[MAXLBAS+1];
  double **h,**v;

  lab = sa->l + sb->l;
  daop = PSDAOPointer(lab);

  v = new double*[daop];
  for (i=0;i<daop;i++)
    v[i] = new double[lab+1];

  h = new double*[daop];
  for (i=0;i<daop;i++)
    h[i] = new double[daop];
  for (i=0;i<daop;i++)
    for (j=0;j<daop;j++)
      h[i][j] = 0.0;

  // Initialize
  r[0] = rb[0]-ra[0];
  r[1] = rb[1]-ra[1];
  r[2] = rb[2]-ra[2];
  r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];

  for (i=0;i<sa->k;i++)
  {
    for (j=0;j<sb->k;j++)
    {
      zp = sa->z[i] + sb->z[j];
      xi = sa->z[i]*sb->z[j]/zp;
      sf = sb->z[j]/zp;
      p[0] = (ra[0]*(sa->z[i]) + rb[0]*(sb->z[j]))/zp;
      p[1] = (ra[1]*(sa->z[i]) + rb[1]*(sb->z[j]))/zp;
      p[2] = (ra[2]*(sa->z[i]) + rb[2]*(sb->z[j]))/zp;
      factor = 2.0*M_PI/zp*exp(-xi*r2)*sa->d[i]*sb->d[j];
      q[0] = p[0] - rc[0];
      q[1] = p[1] - rc[1];
      q[2] = p[2] - rc[2];
      t = zp*(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
      // s-s type
      Gamma(lab,t,f);
      for (k=0;k<=lab;k++)
        v[gaop[0][0][0]][k] = factor*f[k];
      // ax up
      mm = lab;
      for (k=0;k<mm;k++)
        v[gaop[1][0][0]][k] = sf*r[0]*v[gaop[0][0][0]][k] 
                              - q[0]*v[gaop[0][0][0]][k+1];
      for (ax=2;ax<=mm;ax++)
      {
        ma = mm - ax;
        for (k=0;k<=ma;k++)
          v[gaop[ax][0][0]][k] = sf*r[0]*v[gaop[ax-1][0][0]][k] 
          -q[0]*v[gaop[ax-1][0][0]][k+1] + (ax-1)/(2.0*zp)*   
          (v[gaop[ax-2][0][0]][k]-v[gaop[ax-2][0][0]][k+1]);
      }
      // ay up
      for (ax=0;ax<=lab;ax++)
      {
        mm = lab - ax;
        for (k=0;k<mm;k++)
          v[gaop[ax][1][0]][k] = sf*r[1]*v[gaop[ax][0][0]][k]      
                               - q[1]*v[gaop[ax][0][0]][k+1];
        for (ay=2;ay<=mm;ay++)
        {
          ma = mm - ay;
          for (k=0;k<=ma;k++)
            v[gaop[ax][ay][0]][k] = sf*r[1]*v[gaop[ax][ay-1][0]][k]      
            -q[1]*v[gaop[ax][ay-1][0]][k+1] + (ay-1)/(2.0*zp)*           
            (v[gaop[ax][ay-2][0]][k]-v[gaop[ax][ay-2][0]][k+1]);
        }
      }
      // az up
      for (ax=0;ax<=lab;ax++)
      {
        for (ay=0;ay<=lab-ax;ay++)
        {
          mm = lab - ax - ay;
          for (k=0;k<mm;k++)
            v[gaop[ax][ay][1]][k] = sf*r[2]*v[gaop[ax][ay][0]][k]  
                                  - q[2]*v[gaop[ax][ay][0]][k+1];
          for (az=2;az<=mm;az++)
          {
            ma = mm - az;
            for (k=0;k<=ma;k++)
              v[gaop[ax][ay][az]][k] = sf*r[2]*v[gaop[ax][ay][az-1]][k]  
              -q[2]*v[gaop[ax][ay][az-1]][k+1] + (az-1)/(2.0*zp)*        
              (v[gaop[ax][ay][az-2]][k]-v[gaop[ax][ay][az-2]][k+1]);
          }
        }
      }
      // contraction
      factor = 1.0;
      for (k=0;k<daop;k++)
        h[k][0] += factor*v[k][0];
    }
  }
  for (i=0;i<daop;i++)
    delete[] v[i];
  delete[] v;
  HRR(sa->l,sb->l,r,h);
  NormalizeTwoIndex(sa,sb,ib,h,true);
  for (i=0;i<daop;i++)
    delete[] h[i];
  delete[] h;
}

// Calculation of the incomplete gamma function F(t)
// for multicenter integrals over Gaussian functions.
//
// L. E. McMurchie and E.R. Davidson, 
// J. Comp. Phys. 26, 218 (1978).
//
// Original version from Andreas M. Koester, 1996
#define GAMMA_EPS 1.0e-13
void Integrator::Gamma(int m,double t,double* f)
{
  int i,k,ttab;
  double a,b,c,d,expterm;

  if (t<0.0) t = GAMMA_EPS;
  if (t<=GAMMA_EPS) 
  {
    f[m] = 1.0/(2.0*m + 1.0);
    for (i=1;i<=m;i++)
      f[m-i] = 1.0/(2.0*(m-i) + 1.0);
    return;
  }
  else if (t<=12.0)
  {
    ttab = int(10*t+0.5);
    f[m] = ftab[m][ttab];
    for (k=1;k<=6;k++)
      f[m] += ftab[m+k][ttab]*(pow(double(ttab)/10.0 - t,k))/Factorial(k);
    if (m>0) expterm = exp(-t);
    for (i=1;i<=m;i++)
      f[m-i] = 1.0/(2*(m-i)+1) * (2*t*f[m+1-i] + expterm);
    return;
  }
  else if (t<=15.0)
  {
    a = 0.4999489092;
    b = 0.2473631686;
    c = 0.3211809090;
    d = 0.3811559346;
    f[0] = 0.5*sqrt(M_PI/t) - (exp(-t)/t)*(a - b/t + c/(t*t) - d/(t*t*t));
  }
  else if (t<=18.0)
  {
    a = 0.4998436875;
    b = 0.2424943800;
    c = 0.2464284500;
    f[0] = 0.5*sqrt(M_PI/t) - (exp(-t)/t)*(a - b/t + c/(t*t));
  }
  else if (t<=24.0)
  {
    a = 0.4990931620;
    b = 0.2152832000;
    f[0] = 0.5*sqrt(M_PI/t) - (exp(-t)/t)*(a - b/t);
  }
  else if (t<=30.0)
  {
    a = 0.49000000;
    f[0] = 0.5*sqrt(M_PI/t) - (exp(-t)/t)*a;
  }
  else
  {
    f[0] = 0.5*sqrt(M_PI/t);
  }
  if (t>(2.0*m+36)) 
  {
    for (i=1;i<=m;i++)
      f[i] = (2*i-1)/(2*t)*f[i-1];
  }
  else
  {
    expterm = exp(-t);
    for (i=1;i<=m;i++)
      f[i] = 1/(2*t) * ((2*i-1)*f[i-1] - expterm);
  }
}

/*
// Evaluate coulomb integrals
void Integrator::Interaction2(double *ra,double *rc,GaussianShell *sa,GaussianShell *sc,double **blk, int nsda)
{
  int i,k;
  double oblk[GREEN_MAX_NCO][GREEN_MAX_NCO][GREEN_MAX_NCO][GREEN_MAX_NCO];
  bool skipd[2];
  GaussianShell *sb,*sd;

  sb = new GaussianShell();
  sb->l = 0;
  sb->nco = 1;
  sb->z.push_back( 0.0 );
  sb->d.push_back( 1.0 );
  sb->ncsto.push_back( 1.0 );

  sd = new GaussianShell( *sb );

  Interaction(ra,ra,rc,rc,sa,sb,sc,sd,oblk,GREEN_TOL_NUM,skipd,nsda,0,0);
  delete sb;
  delete sd;

  for (i=0;i<sa->nco;i++)
    for (k=0;k<sc->nco;k++)
      blk[i][k] = oblk[i][0][k][0];

}

void Integrator::Interaction3(double *ra,double *rb,double *rc,GaussianShell *sa, GaussianShell *sb, GaussianShell *sc, double blk[GREEN_MAX_NCO][GREEN_MAX_NCO][GREEN_MAX_NCO],bool skip, int nsda,int nsdb)
{
  int i,j,k;
  double oblk[GREEN_MAX_NCO][GREEN_MAX_NCO][GREEN_MAX_NCO][GREEN_MAX_NCO];
  bool skipd[2];
  GaussianShell *sd;

  sd = new GaussianShell();
  sd->l = 0;
  sd->nco = 1;
  sd->z.push_back( 0.0 );
  sd->d.push_back( 1.0 );
  sd->ncsto.push_back( 1.0 );
  Interaction(ra,rb,rc,rc,sa,sb,sc,sd,oblk,GREEN_TOL_NUM,skipd,nsda,nsdb,0);
  delete sd;
  skip = skipd[0];
  for (i=0;i<sa->nco;i++)
    for (j=0;j<sb->nco;j++)
      for (k=0;k<sc->nco;k++)
        blk[i][j][k] = oblk[i][j][k][0];
}

void Integrator::Interaction(double *ra,double *rb,double *rc,double *rd,
GaussianShell *sa, GaussianShell *sb, GaussianShell *sc, GaussianShell *sd,
double blk[GREEN_MAX_NCO][GREEN_MAX_NCO][GREEN_MAX_NCO][GREEN_MAX_NCO],double eps, 
bool *skip, int nsda,int nsdb,int nsdc)
{
  int ax,ay,az,bm1,bx,by,bz,cx,cy,cz,dc,dd,dleft,dright,dx,dy,dz;
  int i,id,ii,im,im1,im2,j,jd,jj,jm1,k,l,lab,lcd,ll,lm1,lm2,lt,lx,ly,lz,m;
  double ab2,cd2,f1,f2,f3,zp,zp2,zq,zq2,zw,zw2;
  double xip,xiq,xiw,ecd,ep,eq,ew,pa,qc,wp,wq;
  double p[3],q[3],rab[3],rcd[3],w[3];
  double f[2*GREEN_MAX_L_I+1];
  double **h;
  double ***v;

  lab = sa->l + sb->l;
  lcd = sc->l + sd->l;
  lt = lab + lcd;

  rab[0] = rb[0] - ra[0];
  rab[1] = rb[1] - ra[1];
  rab[2] = rb[2] - ra[2];
  ab2 = rab[0]*rab[0]+rab[1]*rab[1]+rab[2]*rab[2];

  rcd[0] = rd[0] - rc[0];
  rcd[1] = rd[1] - rc[1];
  rcd[2] = rd[2] - rc[2];
  cd2 = rcd[0]*rcd[0]+rcd[1]*rcd[1]+rcd[2]*rcd[2];

  dc = ((sc->l+1)*(sc->l+2))/2;
  dd = ((sd->l+1)*(sd->l+2))/2;

  dleft = PSDAOPointer(lab);
  dright = PSDAOPointer(lcd);

  h = new double*[dleft];
  for (i=0;i<dleft;i++)
    h[i] = new double[dright];
  for (i=0;i<dleft;i++)
    for (j=0;j<dright;j++)
       h[i][j] = 0.0;

  v = new double**[dleft];
  for (i=0;i<dleft;i++) 
  {
    v[i] = new double*[dright];
    for (j=0;j<dright;j++)
      v[i][j] = new double[lt+1];
  }

  skip[0] = true;
  skip[1] = true;
  for (i=0;i<(signed)sa->z.size();i++)
  {
    for (j=0;j<(signed)sb->z.size();j++)
    {
      zp = sa->z[i] + sb->z[j];
      zp2 = 2.0*zp;
      xip = sa->z[i]*sb->z[j]/zp;
      ep = xip*ab2;
      p[0] = (ra[0]*(sa->z[i]) + rb[0]*(sb->z[j]))/zp;
      p[1] = (ra[1]*(sa->z[i]) + rb[1]*(sb->z[j]))/zp;
      p[2] = (ra[2]*(sa->z[i]) + rb[2]*(sb->z[j]))/zp;
      f1 = 2.0*exp(-ep)*sa->d[i]*sb->d[j];
      if (GREEN_ABS(2.0*GREEN_PI*f1/zp)<eps) goto DOSMIL;
      f1 = f1*(pow(sa->z[i],nsda))*(pow(sb->z[j],nsdb));
      skip[0] = false;
      for (k=0;k<(signed)sc->z.size();k++)
      {
        for (l=0;l<(signed)sd->z.size();l++)
        {
          zq = sc->z[k] + sd->z[l];
          zq2 = 2.0*zq;
          xiq = sc->z[k]*sd->z[l]/zq;
          eq = xiq*cd2;
          q[0] = (rc[0]*(sc->z[k]) + rd[0]*(sd->z[l]))/zq;
          q[1] = (rc[1]*(sc->z[k]) + rd[1]*(sd->z[l]))/zq;
          q[2] = (rc[2]*(sc->z[k]) + rd[2]*(sd->z[l]))/zq;
          zw = zp + zq;
          zw2 = 2*zw;
          xiw = zp*zq/zw;
          w[0] = p[0] - q[0];
          w[1] = p[1] - q[1];
          w[2] = p[2] - q[2];
          ew = xiw*(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
          w[0] = (p[0]*zp + q[0]*zq)/zw;
          w[1] = (p[1]*zp + q[1]*zq)/zw;
          w[2] = (p[2]*zp + q[2]*zq)/zw;
          ecd = exp(-eq)*sc->d[k]*sd->d[l];
          if (GREEN_ABS(2.0*GREEN_PI*ecd/zq)<eps) goto MIL;
          f2 = f1*ecd*sqrt(xiw/GREEN_PI)*(pow(GREEN_PI*GREEN_PI/(zp*zq),1.5))*
               (pow(sc->z[k],nsdc));
          skip[1] = false;
          Gamma(lt,ew,xiw,f);  // ==================
          for (im=0;im<=lt;im++)
            f[im] = f2*f[im];
          // s-s-s-s type
          for (im=0;im<=lt;im++)
            v[0][0][im] = f[im];
          // vrr right side
          for (ll=1;ll<=lcd;ll++)
          {
            lm1 = ll - 1;
            lm2 = GREEN_MAX(0,ll-2);
            m = lt - ll;
            qc = q[2] - rc[2];
            wq = w[2] - q[2];
            ii = gaop[0][0][ll];
            im1 = gaop[0][0][lm1];
            im2 = gaop[0][0][lm2];
            for (im=0;im<=m;im++)
              v[0][ii][im] = qc*v[0][im1][im] + wq*v[0][im1][im+1] +        
              (double(lm1)/zq2)*(v[0][im2][im]-zp/zw*v[0][im2][im+1]);
            qc = q[1] - rc[1];
            wq = w[1] - q[1];
            for (ly=1;ly<=ll;ly++)
            {
              lm1 = ly - 1;
              lm2 = GREEN_MAX(0,ly-2);
              lz = ll - ly;
              ii = gaop[0][ly][lz];
              im1 = gaop[0][lm1][lz];
              im2 = gaop[0][lm2][lz];
              for (im=0;im<=m;im++)
                v[0][ii][im] = qc*v[0][im1][im] + wq*v[0][im1][im+1] +      
                (double(lm1)/zq2)*(v[0][im2][im]-zp/zw*v[0][im2][im+1]);
            }
            qc = q[0] - rc[0];
            wq = w[0] - q[0];
            for (lx=1;lx<=ll;lx++)
            {
              lm1 = lx - 1;
              lm2 = GREEN_MAX(0,lx-2);
              for (ly=0;ly<=ll-lx;ly++)
              {
                lz = ll - lx - ly;
                ii = gaop[lx][ly][lz];
                im1 = gaop[lm1][ly][lz];
                im2 = gaop[lm2][ly][lz];
                for (im=0;im<=m;im++)
                  v[0][ii][im] = qc*v[0][im1][im] + wq*v[0][im1][im+1] +    
                  (double(lm1)/zq2)*(v[0][im2][im]-zp/zw*v[0][im2][im+1]);
              }
            }
          }
          // vrr left side
          for (ll=1;ll<=lab;ll++)
          {
            lm1 = ll - 1;
            lm2 = GREEN_MAX(0,ll-2);
            m = lab - ll;
            pa = p[2]-ra[2];
            wp = w[2]-p[2];
            ii = gaop[0][0][ll];
            im1 = gaop[0][0][lm1];
            im2 = gaop[0][0][lm2];
            for (im=0;im<=m;im++)
              for (id=0;id<dright;id++)
                v[ii][id][im] = pa*v[im1][id][im] + wp*v[im1][id][im+1] + 
                (double(lm1)/zp2)*(v[im2][id][im]-zq/zw*v[im2][id][im+1]);
            for (bz=1;bz<=lcd;bz++)
            {
              f3 = double(bz)/zw2;
              bm1 = bz - 1;
              for (bx=0;bx<=lcd-bz;bx++)
              {
                for (by=0;by<=lcd-bx-bz;by++)
                {
                  jj = gaop[bx][by][bz];
                  jm1 = gaop[bx][by][bm1];
                  for (im=0;im<=m;im++)
                    v[ii][jj][im] = v[ii][jj][im]+ f3*v[im1][jm1][im+1];
                }
              }
            }
            pa = p[1]-ra[1];
            wp = w[1]-p[1];
            for (ly=1;ly<=ll;ly++)
            {
              lm1 = ly - 1;
              lm2 = GREEN_MAX(0,ly-2);
              lz = ll - ly;
              ii = gaop[0][ly][lz];
              im1 = gaop[0][lm1][lz];
              im2 = gaop[0][lm2][lz];
              for (im=0;im<=m;im++)
                for (id=0;id<dright;id++)
                  v[ii][id][im] = pa*v[im1][id][im] + wp*v[im1][id][im+1] + 
                  (double(lm1)/zp2)*(v[im2][id][im]-zq/zw*v[im2][id][im+1]);
              for (by=1;by<=lcd;by++)
              {
                f3 = double(by)/zw2;
                bm1 = by - 1;
                for (bx=0;bx<=lcd-by;bx++)
                {
                  for (bz=0;bz<=lcd-bx-by;bz++)
                  {
                    jj = gaop[bx][by][bz];
                    jm1 = gaop[bx][bm1][bz];
                    for (im=0;im<=m;im++)
                      v[ii][jj][im] += f3*v[im1][jm1][im+1];
                  }
                }
              }
            }
            pa = p[0]-ra[0];
            wp = w[0]-p[0];
            for (lx=1;lx<=ll;lx++)
            {
              lm1 = lx - 1;
              lm2 = GREEN_MAX(0,lx-2);
              for (ly=0;ly<=ll-lx;ly++)
              {
                lz = ll - lx - ly;
                ii = gaop[lx][ly][lz];
                im1 = gaop[lm1][ly][lz];
                im2 = gaop[lm2][ly][lz];
                for (im=0;im<=m;im++)
                  for (id=0;id<dright;id++)
                    v[ii][id][im] = pa*v[im1][id][im] + wp*v[im1][id][im+1] + 
                    (double(lm1)/zp2)*(v[im2][id][im]-zq/zw*v[im2][id][im+1]);
                for (bx=1;bx<=lcd;bx++)
                {
                  f3 = double(bx)/zw2;
                  bm1 = bx - 1;
                  for (by=0;by<=lcd-bx;by++)
                  {
                    for (bz=0;bz<=lcd-bx-by;bz++)
                    {
                      jj = gaop[bx][by][bz];
                      jm1 = gaop[bm1][by][bz];
                      for (im=0;im<=m;im++)
                        v[ii][jj][im] += f3*v[im1][jm1][im+1];
                    }
                  }
                }
              }
            }
          }
          // Contraction
          for (id=0;id<dleft;id++)
            for (jd=0;jd<dright;jd++)
              h[id][jd] += v[id][jd][0];
MIL:
          continue;
        }
      }
DOSMIL:
      continue;
    }
  }
  for (i=0;i<dleft;i++) 
  {
    for (j=0;j<dright;j++)
      delete[] v[i][j];
    delete[] v[i];
  }
  delete[] v;

  if (skip[1]) 
  {
    for (i=0;i<dleft;i++)
      delete[] h[i];
    delete[] h;

    for (k=0;k<sc->nco;k++)
      for (l=0;l<sd->nco;l++)
        for (i=0;i<sa->nco;i++)
          for (j=0;j<sb->nco;j++)
            blk[i][j][k][l] = 0.0;
    return;
  }
  double **hr;
  hr = new double*[dright];
  for (i=0;i<dright;i++)
    hr[i] = new double[dright];

  double hh[dleft][dc][dd];
  for (i=0;i<dleft;i++)
  {
    for (id=0;id<dright;id++)
      hr[id][0] = h[i][id];
    HRR(sc->l,sd->l,rcd,hr);
    for (cx=0;cx<=sc->l;cx++)
    {
      for (cy=0;cy<=sc->l-cx;cy++)
      {
        cz = sc->l - cx - cy;
        k = raop[cx][cy][cz];
        for (dx=0;dx<=sd->l;dx++)
        {
          for (dy=0;dy<=sd->l-dx;dy++)
          {
            dz = sd->l-dx-dy;
            l = raop[dx][dy][dz];
            hh[i][k][l] = hr[gaop[cx][cy][cz]][gaop[dx][dy][dz]];
          }
        }
      }
    }
  }
  for (i=0;i<dleft;i++)
    delete[] h[i];
  delete[] h;

  for (i=0;i<dright;i++)
    delete[] hr[i];
  delete[] hr;

  double **hl;
  hl = new double*[dleft];
  for (i=0;i<dleft;i++)
    hl[i] = new double[dleft];
  for (k=0;k<dc;k++)
  {
    for (l=0;l<dd;l++)
    {
      for (id=0;id<dleft;id++)
        hl[id][0] = hh[id][k][l];
      HRR(sa->l,sb->l,rab,hl);
      for (ax=sa->l;ax>=0;ax--)
      {
        for (ay=sa->l-ax;ay>=0;ay--)
        {
          az = sa->l - ax - ay;
          i = raop[ax][ay][az];
          for (bx=sb->l;bx>=0;bx--)
          {
            for (by=sb->l-bx;by>=0;by--)
            {
              bz = sb->l-bx-by;
              j = raop[bx][by][bz];
              blk[i][j][k][l] = hl[gaop[ax][ay][az]][gaop[bx][by][bz]]*     
              sa->ncsto[i]*sb->ncsto[j]*sc->ncsto[k]*sd->ncsto[l];
            }
          }
        }
      }
    }
  }
  for (i=0;i<dleft;i++)
    delete[] hl[i];
  delete[] hl;
}


// Dipole integrals: <a|r|b>
void Integrator::Dipole(double *ra,double *rb,GaussianShell *sa,
 GaussianShell *sb, double **ib,int indx)
{
  int ap,ap1,ax,ay,az,i,j,na,nb;
  double **ob;
  GaussianShell *wsa;

  int l = GREEN_MAX(sa->l,sb->l) + 1;
  int dover = ((l+1)*(l+2))/2;
  na = ((sa->l+1)*(sa->l+2))/2;
  nb = ((sb->l+1)*(sb->l+2))/2;

  ob = new double*[dover];
  for(i=0;i<dover;i++)
    ob[i] = new double[dover];
  
  wsa = new GaussianShell();
  wsa->z.clear();
  wsa->d.clear();
  for (i=0;i<(signed)sa->z.size();i++)
  {
    wsa->z.push_back( sa->z[i] );
    wsa->d.push_back( sa->d[i] );
  }
  wsa->l = sa->l + 1;

  Overlap(ra,rb,wsa,sb,ob,0,false);

  for (ax=0;ax<=sa->l;ax++)
  {
    for (ay=0;ay<=sa->l-ax;ay++)
    {
      az = sa->l - ax - ay;
      ap = raop[ax][ay][az];
      if (indx==0) ap1 = raop[ax+1][ay][az];
      else if (indx==1) ap1 = raop[ax][ay+1][az];
      else ap1 = raop[ax][ay][az+1];
      for (j=0;j<nb;j++)
        ib[ap][j] = ob[ap1][j];       
    }
  }
  delete wsa;

  Overlap(ra,rb,sa,sb,ob,0,false);

  for (i=0;i<na;i++)
    for (j=0;j<nb;j++)
      ib[i][j] += ra[indx]*ob[i][j];

  for(i=0;i<dover;i++)
    delete[] ob[i];
  delete[] ob;

  for (i=0;i<na;i++)
    for (j=0;j<nb;j++)
      ib[i][j] *= sa->ncsto[i]*sb->ncsto[j];
}

void Integrator::SetInteraction( void(*f)(int,double,double,double*) )
{
  Gamma = f;
}

void Integrator::SetInteractionPotential(GaussianShell *sa)
{
  pot = new GaussianShell( *sa );
}

void GammaSolvent(int m,double t,double rho,double* J)
{
  int i,k;
  double expf,fk;

  if (!pot) 
  {
    cout << "There is no potential for GammaSolvent"<<endl;
    return;
  }

  // Since we are now fitting to Morse
  for (i=0;i<=m;i++)
    J[i] = 0.0;

  for (k=0;k<(signed)pot->z.size();k++)
  {
    fk = pot->z[k]/(pot->z[k]+rho);
    expf = pot->d[k]*exp(-fk*t)*0.5*sqrt(GREEN_PI/rho)*pow(1.0-fk,1.5);
    J[0] += expf;
    for (i=1;i<=m;i++)
    {
      expf *= fk;
      J[i] += expf;
    }
  }
}

void Integrator::CoreAlt(double *ra,double *rb,double *rc,GaussianShell *sa,
 GaussianShell *sb, double **ib,int nsda,int nsdb,bool apply_normalization)
{
  int ax,ay,az,daop,i,j,k,lab,mm,ma;
  double factor,zp,xi,sf,r2,t;
  double r[3],p[3],q[3],f[GREEN_MAX_L_I+1];
  double **v,**h;

  lab = sa->l + sb->l;
  daop = PSDAOPointer(lab);

  v = new double*[daop];
  for (i=0;i<daop;i++)
    v[i] = new double[lab+1];
  for (i=0;i<daop;i++)
    for (j=0;j<=lab;j++)
      v[i][j] = 0.0;

  h = new double*[daop];
  for (i=0;i<daop;i++)
    h[i] = new double[daop];
  for (i=0;i<daop;i++)
    for (j=0;j<daop;j++)
      h[i][j] = 0.0;

  // Initialize
  r[0] = rb[0]-ra[0];
  r[1] = rb[1]-ra[1];
  r[2] = rb[2]-ra[2];
  r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];

  for (i=0;i<(signed)sa->z.size();i++)
  {
    for (j=0;j<(signed)sb->z.size();j++)
    {
      zp = sa->z[i] + sb->z[j];
      xi = sa->z[i]*sb->z[j]/zp;
      sf = sb->z[j]/zp;
      p[0] = (ra[0]*(sa->z[i]) + rb[0]*(sb->z[j]))/zp;
      p[1] = (ra[1]*(sa->z[i]) + rb[1]*(sb->z[j]))/zp;
      p[2] = (ra[2]*(sa->z[i]) + rb[2]*(sb->z[j]))/zp;
      factor = 2.0*GREEN_PI/zp*exp(-xi*r2)*sa->d[i]*sb->d[j];
      q[0] = p[0] - rc[0];
      q[1] = p[1] - rc[1];
      q[2] = p[2] - rc[2];
      t = zp*(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
      // s-s type
      GammaCoulomb(lab,t,t,f);
      for (k=0;k<=lab;k++)
        v[gaop[0][0][0]][k] = factor*f[k];
      // ax up
      mm = lab;
      if (mm>0)
        for (k=0;k<mm;k++)
          v[gaop[1][0][0]][k] = 0.1;//sf*r[0]*v[gaop[0][0][0]][k] 
                            //  - q[0]*v[gaop[0][0][0]][k+1];
      // contraction
      factor = 1.0;
      if (nsda>0) factor *= pow(sa->z[i],nsda);
      if (nsdb>0) factor *= pow(sb->z[j],nsdb);
      for (k=0;k<daop;k++)
        h[k][0] += factor*v[k][0];
    }
  }
  cout << "LAB = " << lab<<" DAOP = "<<daop<<endl;
  for (k=0;k<daop;k++)
  {
    h[k][0] = 0.1;
    h[k][1] = 0.0;
    h[k][2] = 0.0;
    h[k][3] = 0.0;
 }
  for (i=0;i<daop;i++)
    delete[] v[i];
  delete[] v;
  //rfm HRR(sa->l,sb->l,r,h);
  // NormalizeTwoIndex(sa,sb,ib,h,apply_normalization);
  int bx,by,bz,ii;

  for (ax=0;ax<=sa->l;ax++)
  {
    for (ay=0;ay<=sa->l-ax;ay++)
    {
      az = sa->l - ax - ay;
      i = raop[ax][ay][az];
      ii = gaop[ax][ay][az];
      for (bx=0;bx<=sb->l;bx++)
      {
        for (by=0;by<=sb->l-bx;by++)
        {
          bz = sb->l-bx-by;
          j = raop[bx][by][bz];
          ib[i][j] = h[ii][gaop[bx][by][bz]];
        }
      }
    }
  }
  for (i=0;i<daop;i++)
    delete[] h[i];
  delete[] h;
}

*/
