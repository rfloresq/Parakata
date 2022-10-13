// ******************************************************************
// Purpose: Definition of some variable types
//
// 2003-2012: Roberto Flores Moreno
// ******************************************************************

#ifndef TYPES_H
#define TYPES_H

#include "param.h"

typedef struct ShellBasis
{
  int l;
  int k;
  double maxrad;
  double z[MAXCON];
  double d[MAXCON];
  double ncsto[MAXNCO];
} ShellBasis;

typedef struct AtomBasis
{
  int nshl;
  int ncao;
  int nsao;
  double maxrad;
  ShellBasis sb[MAXASHL];
} AtomBasis;

typedef struct MoleculeBasis 
{
  AtomBasis ab[MAXATOM];
} MoleculeBasis;

typedef struct MO
{
  char spin[6];
  char sym[15];
  double energy;
  double occ;
  int nao;
  double c[MAXNBAS];
} MO;

typedef struct PlotEntry
{
  char name[256];
  int type;
  int style;
  int shiny;
  int point_size;
  int line_width;
  double iso;
  double color[4];
  unsigned long list;
} PlotEntry;

typedef struct ImageLabel 
{
  char text[1024];
  int x;
  int y;
} ImageLabel;

#endif // TYPES_H
