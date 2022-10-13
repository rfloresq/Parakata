// *******************************************************************
// Purpose: Definition of parameters
//
// 2003-2012: Roberto Flores-Moreno
// *******************************************************************

#ifndef PARAM_H
#define PARAM_H

#include "macros.h"

#define PROGRAM_NAME "Parakata 0.1"

#define MAX_STR_SIZE 256

#define MAXANIM    10
#define MAXASHL    20
#define MAXATOM   360
#define MAXCON     25
#define MAXEL     103
#define MAXGEO   1000
#define MAXLBAS     4
#define MAXNBAS  1000 
#define MAXSCFCYC 100
#define MAXOPTCYC 500
#define MAXTYP    100

#define MAXNCO  (MAXLBAS+1)*(MAXLBAS+2)/2
#define MAXNORB  4*MAXNBAS
#define MAXNSTRU  MAX(100,MAXOPTCYC)

#endif // PARAM_H
