// *********************************************************
//
// Purpose: symbol-to-atomic number interpreter
//
// Roberto Flores-Moreno 
// 2003-2005: Cinvestav-IPN
// 2007-2008: Universidad de Guanajuato
// 2008-2009: Cinvestav-IPN
// 2009-2010: Universidad de Guadalajara
//
// *********************************************************

#include <ctype.h>
#include <string.h>

#include <global.h>

int symtoan( char *sym )
{
  int an = 0;
  char tst[3];

  if ( strlen(sym) == 0 )
  {
    return 0;
  }
  else
  {
    tst[0] = toupper(sym[0]);
    if ( strlen(sym) == 1 )
    {
      tst[1] = '\0';
    }
    else
    {
      if ( sym[1] >= '0' && sym[1] <= '9' )
      {
        tst[1] = '\0';
      }
      else
      {
        tst[1] = tolower(sym[1]);
        tst[2] = '\0';
      };
    };
  };

  for ( int i = 0 ; i < 104 ; i++ )
  {
    if ( strcmp( tst , ELEMENT_SYMBOL[i] ) == 0 )
    {
      an = i;
      break;
    };
  };
  return an;

}
