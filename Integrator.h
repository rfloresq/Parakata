
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#define MAX_TAB_GAM 120

#define PSDAOPointer(l) ((l)*(l)*(l)-(l))/6 + (l)*(l) + 2*(l) + 1;

// +2 for kinetic energy integrals and avoid shifting full program
extern int gaop[MAXLBAS+3][MAXLBAS+3][MAXLBAS+3];
extern int raop[MAXLBAS+3][MAXLBAS+3][MAXLBAS+3];

extern double ftab[2*MAXLBAS+6+1][MAX_TAB_GAM+1];

class Integrator
{
  public:
    Integrator(void);

    void Core(double*,double*,double*,ShellBasis*,ShellBasis*,double**);

  protected:

    void HRR(int,int,double*,double**);
    void NormalizeTwoIndex(ShellBasis*,ShellBasis*,double**,double**,bool);
    void Gamma(int,double,double*);
};


#endif

