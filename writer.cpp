// **************************************************************
// Purpose: Write external files 
//
// 2003-2012: Roberto Flores-Moreno
// **************************************************************

#include <string.h>

#include <fcns.h>

void writer(const char*filename,const char*fmt)
{
  if ( strncasecmp( fmt , "xyz", 3 ) == 0 )
  {
    writexyz(filename);
  }
};
