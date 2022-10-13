// *********************************************************
// Purpose: Geometry Editor object
//
// 2005-2007, 2012: Roberto Flores Moreno
// *********************************************************

#include <math.h>

#include <QtOpenGL/QGLWidget>
#include <macros.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <QtWidgets>
#include <QStringList>
#include <fcns.h>
#include <global.h>
#include <GeometryEditor.h>   
#include <Parakata.h>
#include <Viewer.h>   
#include <ElementSelector.h>
#include <Vector.h>

using namespace std;
using std::ifstream;


GeometryEditor::GeometryEditor( Parakata *is )
              : QWidget( 0 )
{
  prk = is;
  setWindowTitle( "Parakata: Geometry Editor" );
  
  QList<QString> elebas;
  QList<QString> nucbas;

  QStringList Headers;
  Headers << "Atom" 
          << "Bond" << " " 
          << "Angle" << " " 
          << "Dihedral" << " " 
          << "e-Basis" << "Nuc-Basis";
  QGridLayout *mainLayout = new QGridLayout;

  tw = new QTableWidget( MAXATOM , 9 );
  connect( tw , SIGNAL( itemChanged(QTableWidgetItem*) ), 
           this , SLOT( Edited( QTableWidgetItem*) ));
  tw->setHorizontalHeaderLabels(Headers);

  QToolButton *addButton =  new QToolButton( this );
  addButton->setText( "Add Atom" ); 
  connect(addButton,SIGNAL(clicked()),this,SLOT(AddAtomStart())); 

  QToolButton *delButton =  new QToolButton( this );
  delButton->setText( "Delete Atom" ); 
  connect(delButton,SIGNAL(clicked()),this,SLOT(DeleteAtom())); 

  QToolButton *doneButton =  new QToolButton( this );
  doneButton->setText( "Done" ); 
  connect(doneButton,SIGNAL(clicked()),this,SLOT(hide())); 

  QToolButton *fdButton =  new QToolButton( this );
  fdButton->setText( "Fit Distance" );
  connect(fdButton,SIGNAL(clicked()),this,SLOT(FitDistance()));

  QCheckBox* cb = new QCheckBox( this );
  cb->setText( "Def. Atom" );
  connect( cb , SIGNAL( clicked() ) , this , SLOT( UseDefault() ) );


  QCheckBox* cz = new QCheckBox( this );
  cz->setText( "Cartesians" ); 
  connect( cz , SIGNAL( clicked() ) , this , SLOT( ChangeCoordinates() ) );

  QComboBox *monitorComboBox = new QComboBox;
  monitorComboBox->addItem( QString("Monitor"));
  monitorComboBox->addItem( QString("Distance"));
  monitorComboBox->addItem( QString("Angle"));
  monitorComboBox->addItem( QString("Dihedral"));
  monitorComboBox->addItem( QString("Clean"));
  connect(monitorComboBox, SIGNAL(activated(const QString &)),
          this, SLOT( SetupMonitor(const QString &)));

  QComboBox *pgComboBox = new QComboBox;
  pgComboBox->addItem( QString("Symmetry Group"));
  pgComboBox->addItem( QString("Td"));
  pgComboBox->addItem( QString("Oh"));
  pgComboBox->addItem( QString("Ih"));
  pgComboBox->addItem( QString("Cs"));
  pgComboBox->addItem( QString("Ci"));
  pgComboBox->addItem( QString("C1"));
  pgComboBox->addItem( QString("C2h"));
  pgComboBox->addItem( QString("C3h"));
  pgComboBox->addItem( QString("C4h"));
  pgComboBox->addItem( QString("C5h"));
  pgComboBox->addItem( QString("C2v"));
  pgComboBox->addItem( QString("C3v"));
  pgComboBox->addItem( QString("C4v"));
  pgComboBox->addItem( QString("C5v"));
  pgComboBox->addItem( QString("C6v"));
  pgComboBox->addItem( QString("C7v"));
  pgComboBox->addItem( QString("C8v"));
  pgComboBox->addItem( QString("C9v"));
  pgComboBox->addItem( QString("C*v"));
  pgComboBox->addItem( QString("S2"));
  pgComboBox->addItem( QString("S4"));
  pgComboBox->addItem( QString("S6"));
  pgComboBox->addItem( QString("S8"));
  pgComboBox->addItem( QString("D2h"));
  pgComboBox->addItem( QString("D3h"));
  pgComboBox->addItem( QString("D4h"));
  pgComboBox->addItem( QString("D5h"));
  pgComboBox->addItem( QString("D6h"));
  pgComboBox->addItem( QString("D7h"));
  pgComboBox->addItem( QString("D8h"));
  pgComboBox->addItem( QString("D9h"));
  pgComboBox->addItem( QString("D*h"));
  pgComboBox->addItem( QString("D3d"));
  pgComboBox->addItem( QString("D5d"));
  pgComboBox->addItem( QString("D7d"));
  pgComboBox->addItem( QString("D9d"));
  pgComboBox->addItem( QString("s1c1"));
  connect(pgComboBox, SIGNAL(activated(const QString &)),
          this, SLOT( BuildGroup(const QString &)));

  QSpinBox *risS = new QSpinBox( this );
  risS->setRange( 1 , 10000 );
  risS->setSingleStep( 1 );
  risS->setValue( 4 );
  QLabel* rislab = new QLabel(tr("RIS (alkane):"));
  rislab->setBuddy(risS);
  connect( risS , SIGNAL(valueChanged( int )) , this , SLOT( GoRIS( int ) ) );

  QSpinBox *tubeS = new QSpinBox( this );
  tubeS->setRange( 1 , 1000 );
  tubeS->setSingleStep( 1 );
  tubeS->setValue( 1 );
  QLabel* tubelab = new QLabel(tr("Nanotube:"));
  tubelab->setBuddy(tubeS);
  connect( tubeS , SIGNAL(valueChanged( int )) , this , SLOT( BuildTube(int)));

  mainLayout->addWidget( tw             , 0 , 0 , 3, -1 );
  mainLayout->addWidget( addButton      , 3 , 0 , 1, 1 );
  mainLayout->addWidget( delButton      , 3 , 1 , 1, 1 );
  mainLayout->addWidget( cb             , 3 , 2 , 1, 1 );
  mainLayout->addWidget( rislab         , 3 , 3 , 1, 1 );
  mainLayout->addWidget( risS           , 3 , 4 , 1, 1 );
  mainLayout->addWidget( monitorComboBox, 4 , 0 , 1, 1 );
  mainLayout->addWidget( pgComboBox     , 4 , 1 , 1, 1 );
  mainLayout->addWidget( cz             , 4 , 2 , 1, 1 );
  mainLayout->addWidget( tubelab        , 4 , 3 , 1, 1 );
  mainLayout->addWidget( tubeS          , 4 , 4 , 1, 1 );
  mainLayout->addWidget( fdButton       , 5 , 0 , 1, 1 );
  mainLayout->addWidget( doneButton     , 5 , 1 , 1, 1 );
  setLayout( mainLayout );

  default_element = 6;
  use_default = false;
}

GeometryEditor::~GeometryEditor()
{
}

void GeometryEditor::Launch()
{
  if ( isVisible() )
  {
    hide();
  }
  else
  {
    show();
    Update();
  };
}

void GeometryEditor::AddAtomStart( )
{

  if ( !use_default )
  {
    prk->es->exec(); 
    default_element = prk->es->GetSelectedElement();
  };
    Enabled( false );
  if ( NATOM < 2 )
  {
    AddAtom( NATOM, 0, 0 );
  }
  else 
  {
    if ( NATOM < 3 )
    {
      prk->Report( "Select one atom" ); 
     }
    else if ( NATOM < 4 )
    {
      prk->Report( "Select two atoms" );
    }
    else
    {
      prk->Report( "Select three atoms" );  
};
  Enabled( false );
    prk->status |= SFB;
    prk->SetCursor( "Selection" );

  };
}

void GeometryEditor::Update()
{
  int iatom;

  disconnect( tw , SIGNAL( itemChanged(QTableWidgetItem*) ), 
             this , SLOT( Edited( QTableWidgetItem*) ));

  for ( iatom = 0 ; iatom < NATOM ; iatom++ )
  {
    if ( ! tw->item(iatom,0) )
    {
      tw->setItem(iatom,0, new QTableWidgetItem( " " ));
      tw->setItem(iatom,1, new QTableWidgetItem( " " ));
      tw->setItem(iatom,2, new QTableWidgetItem( " " ));
      tw->setItem(iatom,3, new QTableWidgetItem( " " ));
      tw->setItem(iatom,4, new QTableWidgetItem( " " ));
      tw->setItem(iatom,5, new QTableWidgetItem( " " ));
      tw->setItem(iatom,6, new QTableWidgetItem( " " ));
      tw->setItem(iatom,7, new QTableWidgetItem( " cc-pvdz " ));
      tw->setItem(iatom,8, new QTableWidgetItem( " dirac " ));
    }
    (tw->item(iatom,0))->setText(QString(ELEMENT_SYMBOL[ATOMNO[iatom]]));
    if (CCSET==ZMatrix) 
    {
      if ( iatom > 0 ) 
      {
        (tw->item(iatom,1))->setText( QString::number( NA[iatom] ) );
        (tw->item(iatom,2))->setText( QString::number( COORD[0][iatom][0] ) );
        if ( iatom > 1 )
        {
          (tw->item(iatom,3))->setText( QString::number( NB[iatom] ) );
          (tw->item(iatom,4))->setText( QString::number( COORD[0][iatom][1] ) );
          if ( iatom > 2 )
          {
            (tw->item(iatom,5))->setText(QString::number(NC[iatom] ) );
            (tw->item(iatom,6))->setText(QString::number(COORD[0][iatom][2] ) );
          }
          else
          {
            (tw->item( iatom , 5 ))->setText( " " ); 
            (tw->item( iatom , 6 ))->setText( " " );
          }
        }
        else
        {
          (tw->item( iatom , 3 ))->setText( " " ); 
          (tw->item( iatom , 4 ))->setText( " " );
          (tw->item( iatom , 5 ))->setText( " " );
          (tw->item( iatom , 6 ))->setText( " " );
        }
      }
      else
      {
        (tw->item( iatom , 1 ))->setText( " " ); 
        (tw->item( iatom , 2 ))->setText( " " );
        (tw->item( iatom , 3 ))->setText( " " );
        (tw->item( iatom , 4 ))->setText( " " );
        (tw->item( iatom , 5 ))->setText( " " );
        (tw->item( iatom , 6 ))->setText( " " );
      }
    }
    else
    { 
      (tw->item( iatom , 1 ))->setText( QString( "x" ) );
      (tw->item( iatom , 2 ))->setText( QString::number( COORD[1][iatom][0] ) );
      (tw->item( iatom , 3 ))->setText( QString( "y" ) );
      (tw->item( iatom , 4 ))->setText( QString::number( COORD[1][iatom][1] ) );
      (tw->item( iatom , 5 ))->setText( QString( "z" ) );
      (tw->item( iatom , 6 ))->setText( QString::number( COORD[1][iatom][2] ) );
    }
  }

  connect( tw , SIGNAL( itemChanged(QTableWidgetItem*) ), 
           this , SLOT( Edited( QTableWidgetItem*) ));

  Enabled( true );
}

void GeometryEditor::DeleteAtom()
{
  if ( NATOM > 0 )
  {
    NATOM--;
    (tw->item( NATOM , 0 ))->setText( QString( " " ) );
    (tw->item( NATOM , 1 ))->setText( QString( " " ) );
    (tw->item( NATOM , 2 ))->setText( QString( " " ) );
    (tw->item( NATOM , 3 ))->setText( QString( " " ) );
    (tw->item( NATOM , 4 ))->setText( QString( " " ) );
    (tw->item( NATOM , 5 ))->setText( QString( " " ) );
    (tw->item( NATOM , 6 ))->setText( QString( " " ) );
    (tw->item( NATOM , 7 ))->setText( QString( " " ) );
    (tw->item( NATOM , 8 ))->setText( QString( " " ) );   
    elebas.removeLast(); 
    nucbas.removeLast();
   prk->viewer->Redraw();
    Enabled( true );
  };
}

void GeometryEditor::Enabled( bool enabled )
{
  setEnabled( enabled );
}

void GeometryEditor::Edited( QTableWidgetItem* item ) 
{
  int iatom = tw->currentRow();
  int col = tw->currentColumn();

  if ( col == 0 )
  {
    int an = symtoan( (char*)(item->text().toStdString()).c_str() );
    if ( an != ATOMNO[iatom] )
    {
      ATOMNO[iatom] = an;
      item->setText( QString(ELEMENT_SYMBOL[ATOMNO[iatom]])); 
      prk->viewer->Redraw();
    };
  }
  else if ( col == 2 )
  {
    COORD[CCSET][iatom][0] = item->text().toDouble();
    if (CCSET==ZMatrix) zctocc(iatom);
    prk->viewer->Redraw();
  }
  else if ( col == 4 )
  {
    COORD[CCSET][iatom][1] = item->text().toDouble();
    if (CCSET==ZMatrix) zctocc(iatom);
    prk->viewer->Redraw();
  }
  else if ( col == 6 )
  {
    COORD[CCSET][iatom][2] = item->text().toDouble();
    if (CCSET==ZMatrix) zctocc(iatom);
    prk->viewer->Redraw();
  } 
  else if ( col == 7 )
  {
   elebas[iatom] = item->text();
  }
  else if ( col == 8 )
  {
  nucbas[iatom] = item->text();
  }
}

void GeometryEditor::BuildGroup( const QString &s )
{
  //rfm if ( s == "Symmetry Group" ) return;
  sprintf( gname , "%s", (s.toStdString()).c_str() );
  prk->status |= SFPG;
  Enabled ( false );
}

void GeometryEditor::SetupMonitor( const QString &s )
{
  if ( s == "Monitor" ) return;
  if ( s == "Distance" ) 
  {
    prk->SetCursor( "Selection" );
    prk->status |= SFRV;
    prk->Report("Select two atoms");
  }
  else if ( s == "Angle" )
  {
    prk->SetCursor( "Selection" );
    prk->status |= SFAV;
    prk->Report("Select three atoms");
  }
  else if ( s == "Dihedral" )
  {
    prk->SetCursor( "Selection" );
    prk->status |= SFDV;
    prk->Report("Select four atoms");
  }
  else if ( s == "Clean" ) 
  {
    prk->viewer->monitored_bonds.clear();
    prk->viewer->monitored_angles.clear();
    prk->viewer->monitored_dihedrals.clear();
    prk->viewer->Redraw();
  };
}


void GeometryEditor::ChangeCoordinates()
{
  if (CCSET==Cartesians)
  {
    CCSET = ZMatrix;
  }
  else
  {
    CCSET = Cartesians;
  }
  Update();
}

void GeometryEditor::UseDefault()
{
  use_default =  (use_default ? false : true );
}

void GeometryEditor::FitDistance( void )
{
  prk->SetCursor( "Selection" );
  prk->status |= SFFV;
  prk->Report("Select two atoms");
}


void GeometryEditor::AddAtom( int na, int nb, int nc )
{
  if (NATOM==MAXATOM) 
  {
    errmsg("GeometryEditor::AddAtom",0,"Maximum number of atoms reached",0);
    return;
  };
  elebas.append("cc-pvdz");
  nucbas.append("dirac");
  ATOMNO[NATOM] = default_element;
  NA[NATOM] = na-1;
  NB[NATOM] = nb-1;
  NC[NATOM] = nc-1;
  COORD[0][NATOM][0] = ELEMENT_COV_R[ATOMNO[NATOM]] + 
                       ELEMENT_COV_R[ATOMNO[NA[NATOM]]];
  COORD[0][NATOM][1] = 109.467;
  COORD[0][NATOM][2] = 180.0;
  NATOM++;
  zctocc(NATOM-1);
  Update();
  prk->viewer->Redraw();
}

void GeometryEditor::GoRIS(int n)
{
  NATOM = -n;
  writexyz("RIS.xyz");
}
                                                                                                                                                                                     
void GeometryEditor::BuildTube(int n)
{
  int i,k,nc;
  double angle;
  double bond_length,h,l,r;

  nc = 9;
  bond_length = AngstromToBohr(1.3);

  h = bond_length*sin(M_PI/6.0);
  l = sqrt(4.0*(bond_length*bond_length-h*h));
  r = l/(2.0*sin(M_PI/nc));
  angle = 360.0/nc;
  
  Vector z(0.0,0.0,1.0);
  Vector pos(r,0.0,0.0);
  NATOM = 0;
  for (k=0;k<n;k++)
  {
    for (i=0;i<nc;i++)
    {
      ATOMNO[NATOM] = default_element;
      COORD[1][NATOM][0] = pos[0];
      COORD[1][NATOM][1] = pos[1];
      COORD[1][NATOM][2] = pos[2];
      pos.Rotate(z,angle);
      NATOM++;
    }
    angle /= 2.0;
    pos.Rotate(z,angle);
    angle *= 2.0;
    pos += Vector(0.0,0.0,h);
    for (i=0;i<nc;i++)
    {
      ATOMNO[NATOM] = default_element;
      COORD[1][NATOM][0] = pos[0];
      COORD[1][NATOM][1] = pos[1];
      COORD[1][NATOM][2] = pos[2];
      pos.Rotate(z,angle);
      NATOM++;
    }
    pos += Vector(0.0,0.0,bond_length);
  }
  center();
}
