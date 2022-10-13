// ******************************************************************
// Purpose: EBEWin class functions.
//
// 2003-2021, Roberto Flores Moreno
// ******************************************************************

#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <QtGui>

#include <global.h>
#include <fcns.h>

#include <EBEWin.h>
#include <Viewer.h>
#include <Parakata.h>
#include <Vector.h>
#include <macros.h>


using namespace std;

EBEWin::EBEWin( Parakata *is)
           : QWidget( 0 )
{
  prk = is;
  tabWidget = new QTabWidget;

  QGroupBox *settBox = new QGroupBox(tr("Settings"));
  QHBoxLayout *settLayout = new QHBoxLayout;
  settLayout->addStretch(1);

  QPushButton *doneB = new QPushButton(tr("Done"));
  connect( doneB, SIGNAL(clicked()), this, SLOT(close()));

  settLayout->addWidget( doneB );

  settBox->setLayout( settLayout );

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget( tabWidget );
  mainLayout->addWidget( settBox );
  setLayout( mainLayout );
  setWindowTitle(tr( "Parakata: EBEs" ));

  tabs.clear();
  listados.clear();
}



EBEWin::~EBEWin()
{
  tabs.clear();
  listados.clear();
}


void EBEWin::NewList( char** entriesa , char** entriesb, int nentries , 
		      char* name ) 
{
  tabs.push_back( new QWidget( 0 ) );
  QListWidget *slListBox = new QListWidget;
  QStringList sllist;

  //connect( slListBox, SIGNAL(currentRowChanged(int)), 
  //              this, SLOT(ChangeStructure()));

  for ( int i = 0; i < nentries; i++ )
  {
    sllist.append(tr(entriesa[i]));
    sllist.append(tr(entriesb[i]));
  }
  slListBox->insertItems( 0, sllist);

  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget( slListBox );
  tabs[tabs.size()-1]->setLayout( layout );
  listados.push_back( slListBox );
  tabWidget->addTab( tabs[tabs.size()-1] , QString(name) );
}


void EBEWin::Launch( void )
{
  if ( isVisible() )
  {
    hide();
    return;
  }

  tabWidget->clear();
  tabs.clear();
  listados.clear();

  char *texta[MAXNBAS];
  char *textb[MAXNBAS];

  for ( int i = 0 ; i  < MAXNBAS ; i++ )
  { 
    texta[i] = (char*) malloc(256*sizeof(char));
    if ( ! texta[i] )
    {
      errmsg(__FILE__,__LINE__, 
              "Memory allocation failed, aborting" , 2);
      exit( EXIT_FAILURE );
    }
    textb[i] = (char*) malloc(256*sizeof(char));
    if ( ! textb[i] )
    {
      errmsg(__FILE__,__LINE__, 
              "Memory allocation failed, aborting" , 2);
      exit( EXIT_FAILURE );
    }
  }

  for ( int i = 0; i < NBASIS; i++ )
  {
    sprintf( texta[i],"Alpha %d: %20.3f\n", i+1,TOEV(EBEA[0][i]) );
    sprintf( textb[i],"Beta %d: %20.3f\n", i+1,TOEV(EBEB[0][i]) );
  };
  NewList( texta , textb, NBASIS, "KT" );
  for ( int i = 0; i < NBASIS; i++ )
  {
    sprintf( texta[i],"Alpha %d: %20.3f\n", i+1,TOEV(EBEA[1][i]) );
    sprintf( textb[i],"Beta %d: %20.3f\n", i+1,TOEV(EBEB[1][i]) );
  };
  NewList( texta , textb, NBASIS, "G1PP2D" );

  show();
}

