// ******************************************************************
//
// Purpose: SCF Window object
//
// Roberto Flores-Moreno
// 2003-2005: Cinvestav-IPN
// 2007-2008: Universidad de Guanajuato
// 2008-2009: Cinvestav-IPN
// 2009-2010: Universidad de Guadalajara
//
// ******************************************************************
//
#include <global.h>

#include <QtGui>

#include <SCFWin.h>   
#include <Parakata.h>
#include <CurvePlotter.h>

SCFWin::SCFWin( Parakata *is )
      : QWidget( 0 )
{
  prk = is;
  setWindowTitle( "SCF" );
  tabWidget = new QTabWidget;

  QGridLayout *mainLayout = new QGridLayout;

  QToolButton *startButton =  new QToolButton( this );
  startButton->setText( "Start" ); 
//  connect(startButton,SIGNAL(clicked()), this,SLOT(Start())); 

  QToolButton *stopButton =  new QToolButton( this );
  stopButton->setText( "Stop" ); 
//  connect(stopButton,SIGNAL(clicked()), this,SLOT(Stop())); 

  QToolButton *doneButton =  new QToolButton( this );
  doneButton->setText( "Done" ); 
  connect(doneButton,SIGNAL(clicked()), this,SLOT(Launch())); 

  frame = new QWidget( 0 );
  tabWidget->addTab( frame , tr("SCF Energy (a.u.)") );

  SCFPlotter = new CurvePlotter( this );
  SCFPlotter->resize( 400 , 200 ); 

  mainLayout->addWidget( tabWidget       , 1 , 0 , 6, 4 );
  mainLayout->addWidget( startButton     , 7 , 0 , 1, 1 );
  mainLayout->addWidget( stopButton      , 7 , 1 , 1, 1 );
  mainLayout->addWidget( doneButton      , 7 , 2 , 1, 1 );
  mainLayout->addWidget( SCFPlotter       , 8 , 0 , 6, 4 );
  setLayout( mainLayout );

  hide();
}

SCFWin::~SCFWin()
{
}

void SCFWin::Launch()
{
  if ( isVisible() )
  {
    hide();
  }
  else
  {
    Load();
    show();
  };
}

void SCFWin::Load()
{
  int i;
  char name[MAX_STR_SIZE];
  QListWidget *SCFListBox = new QListWidget;
  QStringList SCFlist;

  for ( i = 1; i < NSCFCYC ; i++ )
  {
    sprintf( name,"Cycle %d: %20.6f\n", i +1 ,SCF_ENERGY[i]);
    SCFlist.append(tr( name) );
  };

  SCFListBox->insertItems( 0, SCFlist);

  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget( SCFListBox );
  frame->setLayout( layout );

  double emin,emax;
  double cval[MAXSCFCYC];
  double eval[MAXSCFCYC];
  emin = SCF_ENERGY[0];
  emax = SCF_ENERGY[0];
  for ( i = 1 ; i < NSCFCYC ; i++ )
  {
    cval[i-1] = (double) i;
    eval[i-1] = (double) SCF_ENERGY[i];
    if ( SCF_ENERGY[i] < emin ) emin = SCF_ENERGY[i];
    if ( SCF_ENERGY[i] > emax ) emax = SCF_ENERGY[i];
  }
  SCFPlotter->SetData( cval , eval ,  NSCFCYC-1 );
  //SCFPlotter->SetBounds( (double)0 , (double)(NSCFCYC+1) ,  
  //                       emin , emax );
}

