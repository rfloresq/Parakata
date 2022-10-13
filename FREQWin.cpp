// ******************************************************************
// Purpose: Frequency window
//
// 2003-2012, Roberto Flores-Moreno
// ******************************************************************

#include <iostream>
#include <ctype.h>
#include <math.h>

#include <QtGui>

#include <macros.h>
#include <global.h>
#include <fcns.h>

#include <Parakata.h>
#include <Viewer.h>
#include <FREQWin.h>
#include <CurvePlotter.h>

FREQWin::FREQWin( Parakata *is )
       : QWidget( 0 )
{
  prk = is;
  setWindowTitle( "Normal Modes Window" );

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
  num_anim_steps = MAXANIM;
  forward = true;
}

FREQWin::~FREQWin()
{
}

void FREQWin::Launch()
{
  if ( isVisible() )
  {
    hide();
  }
  else
  {
    Update();
    show();
  };
}

void FREQWin::ChangeSelection(int n)
{
  int iatom,istep;
  double DISP[MAXATOM][3];
  char filename[MAX_STR_SIZE];

  // Get displacements
  get_structure("VIB",n);
  for (iatom=0;iatom<NATOM;iatom++)
  {
    DISP[iatom][0] = 2.0*COORD[1][iatom][0];
    DISP[iatom][1] = 2.0*COORD[1][iatom][1];
    DISP[iatom][2] = 2.0*COORD[1][iatom][2];
  }
  
  get_structure("SCF",0);
  save_structure("ANIM",0);
  for (istep=1;istep<=num_anim_steps;istep++)
  {
    for (iatom=0;iatom<NATOM;iatom++)
    {
      COORD[1][iatom][0] += DISP[iatom][0]/num_anim_steps;
      COORD[1][iatom][1] += DISP[iatom][1]/num_anim_steps;
      COORD[1][iatom][2] += DISP[iatom][2]/num_anim_steps;
    }
    save_structure("ANIM",istep);
    // Save structure for re-start of optimizations 
    sprintf(filename,"freq%dfrom%d.xyz",istep,num_anim_steps);
    writexyz(filename);
  }
  get_structure("SCF",0);
  prk->viewer->Redraw();
  anim_pos = 0;
  //update();
}

void FREQWin::Advance(void)
{
  if (forward)
  {
    anim_pos++;
    if (anim_pos>=num_anim_steps) 
    {
      forward = false;
      anim_pos--;
    }
  }
  else
  {
    anim_pos--;
    if (anim_pos<=0) 
    {
      forward = true;
      anim_pos++;
    }
  }
  get_structure("ANIM",anim_pos);
  prk->viewer->Redraw();
  timer->start(anim_delay);
}

void FREQWin::SetAnimation(void)
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
    get_structure("SCF",0);
    prk->viewer->Redraw();
  }
}

void FREQWin::ChangeAnimSpeed(int ns)
{
  anim_delay = ns;
}

void FREQWin::Update( void )
{
  char name[MAX_STR_SIZE];

  QStringList freqlist;
  for ( int imode = 0; imode < NVIB; imode++ )
  {
    sprintf( name,"%10.1f\n", VIB[imode][0]);
    cout << ":  "<<imode << " "<<VIB[imode][0]<< " "<<VIB[imode][1]<<endl;
    freqlist.append(tr( name) );
  };
  theListBox->clear();
  theListBox->insertItems( 0, freqlist);
  (frame->layout())->addWidget( theListBox );
  double fval[4000];
  double fint[4000];
  double maxinten,a;
  int i,j;
  double f,inten,swidth,sshift;
  swidth = 5*5;
  f = 0.0;
  j = 0;
  do 
  {
    f += 10.0;
    inten = 0.0;
    for ( i = 0 ; i < NVIB ; i++ )
    {
      sshift = f - VIB[i][0];
      sshift = sshift*sshift;
      a = VIB[i][1]; 
      inten += a*(swidth/(swidth+sshift));
    }
    fval[j] = f;
    fint[j] = inten;
    j += 1;
  } while ( f < 4000.0 );
/*
  maxinten = 0.0;
  for ( i = 0 ; i < j ; i++ )
  {
    maxinten = MAX(maxinten,fint[i]);
  }
  for ( i = 0 ; i < j ; i++ )
  {
    fint[i] = maxinten - fint[i];
  }
*/
  
  plotter->SetData( fval , fint ,  j );
  //update();
}
