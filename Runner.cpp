// ******************************************************************
// Purpose: class Runner functions. 
//
// 2012: Roberto Flores-Moreno
// ******************************************************************

#include <fstream>
#include <iomanip>

#include <Runner.h>   

#include <macros.h>
#include <fcns.h>
#include <global.h>

using namespace std;

Runner::Runner(void)
      : QWidget( 0 )
{

  programBox = new QComboBox;
  programBox->addItem(QString("deMon2k"));

  theoryBox = new QComboBox;
  theoryBox->addItem(QString("HF"));
  theoryBox->addItem(QString("DFT"));

  basisBox = new QComboBox;
  basisBox->addItem(QString("(DZVP)"));
  basisBox->addItem(QString("(TZVP)"));
  basisBox->addItem(QString("(RECP|SD)"));

  optCheck = new QCheckBox(this);
  optCheck->setText("Optimize");

  freqCheck = new QCheckBox(this);
  freqCheck->setText("Frequencies");

  fukuiCheck = new QCheckBox(this);
  fukuiCheck->setText("Fukui");

  chargeSpin = new QDoubleSpinBox( this );
  chargeSpin->setRange( -10.0 , 10.0 );
  chargeSpin->setSingleStep( 1.0 );
  chargeSpin->setDecimals( 1 );
  chargeSpin->setValue( 0.0 );
  QLabel* chargelab = new QLabel(tr("Charge:"));
  chargelab->setBuddy(chargeSpin);

  QGridLayout *mainLayout = new QGridLayout;

  runButton =  new QToolButton( this );  
  runButton->setText( "Run" );
  connect( runButton,SIGNAL(clicked()),this,SLOT(Run()));

  QToolButton* doneButton =  new QToolButton( this );  
  doneButton->setText( "Done" );
  connect( doneButton,SIGNAL(clicked()),this,SLOT(hide()));

  mainLayout->addWidget( programBox   ,      0 , 0 , 1 , 1);
  mainLayout->addWidget( theoryBox    ,      1 , 0 , 1 , 1);
  mainLayout->addWidget( basisBox     ,      2 , 0 , 1 , 1);
  mainLayout->addWidget( optCheck     ,      3 , 0 , 1 , 1);
  mainLayout->addWidget( freqCheck    ,      4 , 0 , 1 , 1);
  mainLayout->addWidget( fukuiCheck   ,      5 , 0 , 1 , 1);
  mainLayout->addWidget( runButton    ,      6 , 0 , 1 , 1);
  mainLayout->addWidget( chargelab    ,      7 , 0 , 1 , 1);
  mainLayout->addWidget( chargeSpin   ,      7 , 1 , 1 , 1);
  mainLayout->addWidget( doneButton   ,      8 , 0 , 1 , 1);

  setLayout( mainLayout );

  setWindowTitle( "Run calculation" );
  hide();
}

Runner::~Runner()
{
}


void Runner::Run()
{
  QString quantum_program;

  quantum_program = programBox->currentText();

  if (quantum_program == "deMon2k" )
  {
    Prepare_deMon();
    system("cp $SNP_JABALI_BASIS ./BASIS "); 
    system("$SNP_JABALI_BINARY ");
    reader("deMon.out","deMon2k");
  }
}

void Runner::Prepare_deMon(void)
{
  ofstream fout("deMon.inp");
  if ( fout.fail() )
  {
    errmsg(__FILE__,__LINE__,"Error opening file",0);
  };

  int i,iatom;

  fout <<"#\n";
  fout <<"# Generated by "<<PROGRAM_NAME<<"\n";
  fout <<"#\n";
  fout <<"GUESS CORE"<<endl; 
  fout <<"ERIS DIRECT"<<endl; 
  fout <<"CHARGE "<<chargeSpin->value()<<endl; 
  fout <<"AUXIS (GEN-A2)"<<endl; 
  fout <<"BASIS "<< (basisBox->currentText()).toStdString()<<endl;
  
  if ( optCheck->isChecked() ) fout <<"OPTIMIZE"<<endl;
  if ( freqCheck->isChecked() ) fout <<"FREQU"<<endl;
  if ( fukuiCheck->isChecked() )
  {
    fout << "PRINT MOS FUKUI" <<endl;
    fout << "FUKUI" <<endl;
  }
  else 
  {
    fout <<"PRINT MOS"<<endl; 
  }

  fout <<"GEOMETRY CARTESIAN ANGSTROMS"<<endl;
 
  for ( iatom = 0 ; iatom < NATOM ; iatom++ )
  {
    fout << setw(2) << left << ELEMENT_SYMBOL[ATOMNO[iatom]] << " ";
    for ( i = 0 ; i < 3 ; i++ )
    {
      fout << fixed << setw(10) << setprecision(5) << right 
           << (double)BohrToAngstrom(COORD[1][iatom][i]) << " ";
    }
    fout << "\n";
  }
  fout << "END\n";

  fout.close();
};

