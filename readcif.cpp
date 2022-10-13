// **************************************************************
// Purpose: Read CIF files
//
// Roberto Flores-Moreno
// 2003-2005: Cinvestav-IPN
// 2007-2008: Universidad de Guanajuato
// 2008-2009: Cinvestav-IPN
// 2009-2010: Universidad de Guadalajara
//
// **************************************************************

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include <fcns.h>
#include <global.h>

#define MAX_STR_SIZE_CIF 80

using namespace std;
using std::ifstream;

double getreal(const char* istr)
{
  char str[MAX_STR_SIZE];

  strcpy((char*)istr,str);
  for (int i=0;i<(signed)strlen(str);i++) 
  {
    if(str[i]=='(') 
    {
      str[i] = '\0';
      break;
    }
  }
  return strtod(str,0);
}

void readcif(const char*filename)
{
  ifstream f(filename);
  if ( f.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
    return;
  };

  int sindex;
  string stringLine;
  char line[MAX_STR_SIZE_CIF];
  char str[MAX_STR_SIZE_CIF];
  double a,b,c;
  while(f.getline(line,MAX_STR_SIZE_CIF))
  {
    stringLine=line;
    sindex=stringLine.find("_cell_length_a");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      a = getreal(&line[16]);
    sindex=stringLine.find("_cell_length_b");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      b = getreal(&line[16]);
    sindex=stringLine.find("_cell_length_c");
    if((sindex>=0)&&(sindex<(signed)stringLine.size())) 
      c = getreal(&line[16]);
    sindex=stringLine.find("_atom_site_disorder_group");
    if((sindex>=0)&&(sindex<(signed)stringLine.size()))
    {
      NATOM = 0;
      f >> str;
      f >> str;
      while (strstr(str,"_atom")==NULL)
      {
        ATOMNO[NATOM] = symtoan(str);
        f >> str; COORD[1][NATOM][0] = AngstromToBohr(a*getreal(str));
        f >> str; COORD[1][NATOM][1] = AngstromToBohr(b*getreal(str));
        f >> str; COORD[1][NATOM][2] = AngstromToBohr(c*getreal(str));
        f.getline(line,MAX_STR_SIZE_CIF);
        f >> str;
        f >> str;
        NATOM++;
        if (NATOM==MAXATOM) 
        {
          sprintf(str,"Loading only %d atoms",NATOM);
          errmsg(__FILE__,__LINE__,str,0);
          break;
        };
      }
      break;
    }
  }

  f.close();

  center();
};
