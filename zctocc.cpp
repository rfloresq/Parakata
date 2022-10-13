// *********************************************************
//
// Purpose: Obtain Cartesian coordinates fom Z-Matrix
//
// Dec 2004, Roberto Flores Moreno 
// Feb 2007, Roberto Flores Moreno
//
// *********************************************************

#include <fcns.h>
#include <global.h>
#include <Vector.h>

void zctocc( int llatom ) 
{
  int i,iatom,bondref,angleref;

  for ( iatom = llatom ; iatom < NATOM ; iatom++ )
  {
    bondref = NA[iatom];

    for ( i = 0 ; i < 3 ; i++ )
    {
      COORD[1][iatom][i] = COORD[1][bondref][i];
    }

    if ( iatom > 0 )
    {
      angleref = NB[iatom];
      Vector shift(1.0,0.0,0.0);
      if ( iatom != 1 )
      {
        shift = Vector(&COORD[1][angleref][0]) - 
                Vector(&COORD[1][bondref][0]);
        shift ^= (1.0);    
      };
      shift *= COORD[0][iatom][0];
      for ( i = 0 ; i < 3 ; i++ )
      {
        COORD[1][iatom][i] += shift[i];
      }

      if ( iatom > 1 )
      {
        Vector axis;
        shift = Vector(&COORD[1][angleref][0]) - 
                Vector(&COORD[1][bondref][0]);
        if ( iatom == 2 )
        {
          axis = Vector(0.0,1.0,0.0) < shift;
        }
        else
        {
          axis = Vector(&COORD[1][NC[iatom]][0]) - 
                 Vector(&COORD[1][bondref][0]);
          axis <= shift;
        };
        shift = Vector(&COORD[1][iatom][0]) - 
                Vector(&COORD[1][bondref][0]);
        shift.Rotate(axis,COORD[0][iatom][1]);
        for ( i = 0 ; i < 3 ; i++ )
        {
          COORD[1][iatom][i] = COORD[1][bondref][i] + shift[i];
        }

        if ( iatom > 2 )
        {
          axis = Vector(&COORD[1][angleref][0]) - 
                 Vector(&COORD[1][bondref][0]);
          shift = Vector(&COORD[1][iatom][0]) - 
                  Vector(&COORD[1][bondref][0]);
          shift.Rotate(axis,COORD[0][iatom][2]);
          for ( i = 0 ; i < 3 ; i++ )
          {
            COORD[1][iatom][i] = COORD[1][bondref][i] + shift[i];
          }
        };
      };
    };
  };
  center();
};

