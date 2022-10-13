// *********************************************************
// Purpose: Function prototypes
//
// 2003-2012: Roberto Flores Moreno
// *********************************************************

#ifndef FCNS_H
#define FCNS_H

#include <QtCore>
#include "types.h"

#include "Vector.h"

void center(void);
void errmsg(const char*,int,const char*,int);
void gbtc(void);
void get_structure(char*,int);
bool loadbasis(int,char*);

void reader(const char*,const char*);
void readcif(const char*);
void readxyz(const char*);
void readgcube(const char*);
void readDrhocube(const char*);
void read_MOLDEN_mol(const char*);
void read_MOLEKEL_mkl(const char*);
void read_deMon_out(const char*);
void read_deMon_inp(const char*);
void read_Nagual_out(const char*);
void read_GEGA_gaf(const char*);
void read_Parakata_cube(const char*);

void save_structure(char*,int);
void setup_basis(void);
void setup_fields(void);

void writer(const char*,const char*);
void write_deMon_inp(const char*,QString(quantum_type));
void writexyz(const char*);
void write_MOLDEN_mol(const char*);

void zctocc(int);
int sgdrv(char*,int*); 
int sgbdrv(char*,int*); 
int pgbdrv(char*,int*); 
int appcgroup(char*,int,int*);
int appdgroup(char*,int,int*);
int appsgroup(int,int*);
int apptgroup(int*);
int appogroup(int*);
int appigroup(int*);
void build_deMon_input(char*,const QString& );
void clonatom(int,Vector);
void apptrans(Vector,double,int);
void appcaxis(int,int,Vector,int); 
void appsaxis(int,int,Vector,int);
void appsigma(Vector,int);


double distance(int,int);
double double_factorial(int);
double factorial(int);
double noverk(int,int);

int symtoan(char*);

#endif  // FCNS_H
