// ******************************************************************
// Purpose: OPTWin class functions.
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

#include <OPTWin.h>
#include <Viewer.h>
#include <Parakata.h>
#include <Vector.h>
#include <macros.h>


using namespace std;

OPTWin::OPTWin( Parakata *is)
           : QWidget( 0 )
{
  prk = is;
  tabWidget = new QTabWidget;

  QGroupBox *settBox = new QGroupBox(tr("Settings"));
  QHBoxLayout *settLayout = new QHBoxLayout;
  settLayout->addStretch(1);

  anim_delay = 100;
  QSpinBox *adS = new QSpinBox( this );
  adS->setRange( 1 , 1000000 );
  adS->setSingleStep( 5 );
  adS->setValue( anim_delay );
  QLabel* lab = new QLabel(tr("Delay:"));
  lab->setBuddy(adS);
  connect( adS , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangeAnimSpeed( int ) ) );

  animButton =  new QToolButton( this );
  animButton->setText( "ANIMATE" ); 
  connect(animButton,SIGNAL(clicked()), this,SLOT(SetAnimation())); 

  QPushButton *doneB = new QPushButton(tr("Done"));
  connect( doneB, SIGNAL(clicked()), this, SLOT(close()));

  settLayout->addWidget( adS );
  settLayout->addWidget( animButton );
  settLayout->addWidget( doneB );

  settBox->setLayout( settLayout );

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget( tabWidget );
  mainLayout->addWidget( settBox );
  setLayout( mainLayout );
  setWindowTitle(tr( "Parakata: Geometries" ));

  tabs.clear();
  listados.clear();

  timer = new QTimer(this);
}



OPTWin::~OPTWin()
{
  tabs.clear();
  listados.clear();
}


void OPTWin::NewList( char** entries , int nentries , char* name ) 
{
  tabs.push_back( new QWidget( 0 ) );
  QListWidget *slListBox = new QListWidget;
  QStringList sllist;

  connect( slListBox, SIGNAL(currentRowChanged(int)), 
                this, SLOT(ChangeStructure()));

  for ( int i = 0; i < nentries; i++ )
  {
    sllist.append(tr(entries[i]));
  };
  slListBox->insertItems( 0, sllist);

  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget( slListBox );
  tabs[tabs.size()-1]->setLayout( layout );
  listados.push_back( slListBox );
  tabWidget->addTab( tabs[tabs.size()-1] , QString(name) );
}


void OPTWin::Launch( void )
{
  if ( isVisible() )
  {
    hide();
    return;
  }

  tabWidget->clear();
  tabs.clear();
  listados.clear();

  char *text[MAXOPTCYC];

  for ( int i = 0 ; i  < MAXOPTCYC ; i++ )
  { 
    text[i] = (char*) malloc(256*sizeof(char));
    if ( ! text[i] )
    {
      errmsg(__FILE__,__LINE__, 
              "Memory allocation failed, aborting" , 2);
      exit( EXIT_FAILURE );
    }
  }

  QStringList optlist;
  for ( int i = 0; i < NOPTCYC-1; i++ )
  {
    sprintf( text[i],"%20.10f\n", OPT_ENERGY[i+1]);
    optlist.append(tr(text[i]) );
  };

  double id[MAXOPTCYC];
  for (int i=0;i<NOPTCYC-1;i++)
  {
    id[i] = (double)i;
  }
  NewList( text , NOPTCYC-1 , "Optimization" );

  show();
}

void OPTWin::ChangeStructure()
{
  int table = tabWidget->currentIndex();
  if ( listados[table]->currentItem() == NULL ) 
  {
    errmsg(__FILE__,__LINE__, 
           "You must choose a field first" , 0 );
    return;
  };
  int n = listados[table]->currentRow();
  get_structure("OPT",n+1);
  prk->viewer->Redraw();
}

void OPTWin::Advance(void)
{
  int table = tabWidget->currentIndex();
  int n = listados[table]->currentRow();
  if (n+1<listados[table]->count()) n++;
  else n = 0;
  listados[table]->setCurrentRow(n);
  timer->start(anim_delay);
}

void OPTWin::SetAnimation(void)
{
  if ( animButton->text() == "ANIMATE" )
  {
    connect(timer,SIGNAL(timeout()),this,SLOT(Advance()));
    animButton->setText("STOP");
    Advance();
  }
  else
  {
    disconnect(timer,SIGNAL(timeout()),this,SLOT(Advance()));
    animButton->setText("ANIMATE");
  }
}

void OPTWin::ChangeAnimSpeed(int ns)
{
  anim_delay = ns;
}
