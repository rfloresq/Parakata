// *******************************************************************
// Purpose: Management of structures
//
// Roberto Flores Moreno
// *******************************************************************
#include <string.h>

#include<iostream>

#include <global.h>

#define SCF_OFFSET  0
#define OPT_OFFSET  SCF_OFFSET+1
#define VIB_OFFSET OPT_OFFSET+MAXOPTCYC
#define ANIM_OFFSET VIB_OFFSET+3*MAXATOM

typedef struct Structure
{
  int natom;
  int atomno[MAXATOM];
  double x[MAXATOM];
  double y[MAXATOM];
  double z[MAXATOM];
} Structure;

Structure STRUCT_BUF[ANIM_OFFSET+MAXANIM+1];
Structure STRUCT_AUX[ANIM_OFFSET+MAXANIM+1];

int modid(char* type, int id)
{
  int mid;

  if (strcasecmp(type,"SCF")==0) mid = id + SCF_OFFSET;
  else if (strcasecmp(type,"OPT")==0) mid = id + OPT_OFFSET;
  else if (strcasecmp(type,"VIB")==0) mid = id + VIB_OFFSET;
  else if (strcasecmp(type,"ANIM")==0) mid = id + ANIM_OFFSET;

  return mid;
}

void save_structure(char*type, int id)
{
  int mid = modid(type,id);

  Structure s;

  s.natom = NATOM;
  for (int iatom=0;iatom<NATOM;iatom++)
  { 
    s.atomno[iatom] = ATOMNO[iatom];
    s.x[iatom] = COORD[1][iatom][0];
    s.y[iatom] = COORD[1][iatom][1];
    s.z[iatom] = COORD[1][iatom][2];
  }

  if ( strcmp( type , "GRAD") == 0 ) STRUCT_AUX[mid] = s;
  else STRUCT_BUF[mid] = s;
}

void get_structure(char*type,int id)
{
  int mid = modid(type,id);

  Structure s;

  s = STRUCT_BUF[mid];
  NATOM = s.natom;
  for (int iatom=0;iatom<NATOM;iatom++)
  { 
    ATOMNO[iatom] = s.atomno[iatom];
    COORD[1][iatom][0] = s.x[iatom];
    COORD[1][iatom][1] = s.y[iatom];
    COORD[1][iatom][2] = s.z[iatom];
  }

  s = STRUCT_AUX[mid];
  NATOM = s.natom;
  for (int iatom=0;iatom<NATOM;iatom++)
  { 
    ATOMFORCE[iatom][0] = s.x[iatom];
    ATOMFORCE[iatom][1] = s.y[iatom];
    ATOMFORCE[iatom][2] = s.z[iatom];
  }
}
