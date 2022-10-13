// ******************************************************************
// Purpose: Listed Curves
//
// Roberto Flores-Moreno 
// ******************************************************************

#include <iostream>
#include <ctype.h>
#include <math.h>

#include <QtGui>

#include <macros.h>
#include <global.h>
#include <ListedCurve.h>   
#include <Parakata.h>
#include <CurvePlotter.h>

ListedCurve::ListedCurve( Parakata *is, const char* title )
           : QWidget( 0 )
{
  prk = is;
  setWindowTitle( title );

  QGridLayout *mainLayout = new QGridLayout;

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

  QToolButton *doneButton =  new QToolButton( this );
  doneButton->setText( "DONE" ); 
  connect(doneButton,SIGNAL(clicked()), this,SLOT(Launch())); 

  frame = new QWidget( 0 );
  QVBoxLayout *layout = new QVBoxLayout;
  frame->setLayout( layout );

  plotter = new CurvePlotter( this );
  plotter->resize( 400 , 200 ); 
  connect(this,SIGNAL(selectionChanged(int)),
          plotter,SLOT(HighlightPoint(int)));

  mainLayout->addWidget( frame           , 1 , 0 , 6, 4 );
  mainLayout->addWidget( plotter         , 7 , 0 , 6, 4 );
  mainLayout->addWidget( lab             , 15 , 0 , 1, 1 );
  mainLayout->addWidget( adS             , 15 , 1 , 1, 1 );
  mainLayout->addWidget( animButton      , 15 , 2 , 1, 1 );
  mainLayout->addWidget( doneButton      , 15 , 3 , 1, 1 );
  setLayout( mainLayout );

  theListBox = new QListWidget;
  connect(theListBox,SIGNAL(currentRowChanged(int)),
          this,SLOT(ChangeSelection(int)));

  timer = new QTimer(this);
}

ListedCurve::~ListedCurve()
{
}

void ListedCurve::Launch()
{
  if ( isVisible() )
  {
    hide();
  }
  else
  {
    show();
  };
}

void ListedCurve::Load(QStringList thelist, double* xdata, double* ydata, int n)
{
  theListBox->clear();
  theListBox->insertItems( 0, thelist);
  (frame->layout())->addWidget( theListBox );
  plotter->SetData( xdata , ydata ,  n );
  update();
}

void ListedCurve::ChangeSelection(int n)
{
  emit selectionChanged(n);
  update();
}

void ListedCurve::Advance(void)
{
  int n = theListBox->currentRow();
  if (n+1<theListBox->count()) n++;
  else n = 0;
  theListBox->setCurrentRow(n);
  timer->start(anim_delay);
}

void ListedCurve::SetAnimation(void)
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

void ListedCurve::ChangeAnimSpeed(int ns)
{
  anim_delay = ns;
}
