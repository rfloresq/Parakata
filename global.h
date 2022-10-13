// *******************************************************************
// Purpose: Definition of global variables.
//
// 2003-2012: Roberto Flores-Moreno
// *******************************************************************

#ifndef GLOBAL_H
#define GLOBAL_H

#define ZMatrix 0
#define Cartesians 1

#include "types.h"
#include "Parakata.h"
#include <vector>

using namespace std;

extern Parakata* prk;

extern MoleculeBasis mb;

extern int NORB;
extern MO MOS[MAXNORB];

extern int KNORB;
extern int PBCN1;
extern int PBCN2;
extern int PBCN3;
extern double PBCK[3];
extern double PBCCD1[3];
extern double PBCCD2[3];
extern double PBCCD3[3];
extern MO KMOSR[MAXNORB];
extern MO KMOSI[MAXNORB];

extern int NATOM;
extern int NSCFCYC;
extern int NOPTCYC;
extern int NVIB;
extern int NBASIS;
extern int NSPH;
extern int NCPs;
extern int selectedMethod;
extern int selectedeBasis;
extern int selectednucBasis;
extern int NELECTRON;
extern int HOMON[2];
extern int LUMON[2];
extern int NDEGHOMO[2];
extern int NDEGLUMO[2];

extern int CCSET;
extern int ATOMNO[MAXATOM];
extern int NLATTICE[4];
extern int NA[MAXATOM];
extern int NB[MAXATOM];
extern int NC[MAXATOM];

extern bool drawmol;
extern bool Gaussbool;
extern bool Demonbool;
extern bool nuclearCheckbool;

extern vector<double> rhog;
extern vector<double> rhod;

extern double MOLDIPOLE[3];

extern double ATOMCHARGE[MAXATOM];
extern double ATOMFUKUI_H[MAXATOM];
extern double ATOMFUKUI_L[MAXATOM];
extern double ATOMFORCE[MAXATOM][3];
extern double ATOMNUCFUKUI_H[MAXATOM][3];
extern double ATOMNUCFUKUI_L[MAXATOM][3];
extern double ATOMNFFEX[MAXATOM][3];
extern double ATOMNFFEY[MAXATOM][3];
extern double ATOMNFFEZ[MAXATOM][3];
extern double LATTICE[12];
extern double COORD[2][MAXATOM][3];
extern double CPCOORD[MAXATOM][3];
// Order KT, EP2
extern double PSA[2][MAXNBAS];
extern double PSB[2][MAXNBAS];
extern double EBEA[2][MAXNBAS];
extern double EBEB[2][MAXNBAS];
extern double KEBEA[2][MAXNBAS];
extern double KEBEB[2][MAXNBAS];
extern double SCF_ENERGY[MAXSCFCYC];
extern double OPT_ENERGY[MAXOPTCYC];
extern double VIB[3*MAXATOM][2];
extern double CTOSTM[MAXLBAS][2*MAXLBAS+1][(MAXLBAS+1)*(MAXLBAS+2)/2];
extern double geosrs;
extern double geocrs;
extern double* SCALAR_FIELD;

// Quadratic fields: matrices
extern double* DENSITY_MATRIX;          // Fukui f^-
extern double* DMR1;                    
extern double* DMR2;                    
extern double* DMR3;                    
extern double* FUKUIL_MATRIX;          // Fukui f^-
extern double* FUKUIR_MATRIX;          // Fukui f^+

extern double ELEMENT_COV_R[MAXEL+1];
extern double ELEMENT_VDW_R[MAXEL+1];
extern double ELEMENT_COLOR[MAXEL+1][3];
extern double ELEMENT_WEIGHT[MAXEL+1];
extern char* ELEMENT_SYMBOL[MAXEL+1];
extern char* ELEMENT_NAME[MAXEL+1];
extern bool SPHORB;


#endif // GLOBAL_H
