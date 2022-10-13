// ******************************************************************
// Purpose: class Parakata functions. 
//
// 2003-2021: Roberto Flores-Moreno
// ******************************************************************

#include <stdarg.h>
#include <pthread.h>

#include <iostream>

#include <QtCore>
#include <QtGui>

#include <QMainWindow>
#include <global.h>
#include <fcns.h>
#include <Parakata.h>   
#include <Viewer.h>
#include <EBEWin.h>
#include <SCFWin.h>
#include <OPTWin.h>
#include <FREQWin.h>
#include <ElementSelector.h>
#include <PlotManager.h>
#include <ImageRender.h>
#include <GeometryEditor.h>
#include <Runner.h>

static const char  *READ_FORMATS[] = {"deMon2k","Nagual","XYZ",NULL};

static const char *WRITE_FORMATS[] = {"XYZ", NULL};

Parakata::Parakata( QApplication *iapp, int argc, char** argv )
     : QMainWindow( 0 )
{
  app = iapp;

  status = 0;

  tabWidget = new QTabWidget;

  freqwin = new FREQWin(this);
  scfwin = new SCFWin(this);

  optwin = new OPTWin(this);
  optwin->hide();

  ebewin = new EBEWin(this);
  ebewin->hide();

  es = new ElementSelector(this,0);

  geoeditor = new GeometryEditor( this );
  geoeditor->hide();

  viewer = new Viewer( this );
  runner = new Runner();

  statusbar = new QStatusBar( this );

 int i;


  i = 0;
  while ( READ_FORMATS[i] != NULL )
  {
    RFMTA.push_back( new QAction( tr( READ_FORMATS[i] ), this ) );
    RFMTA[i]->setStatusTip(tr("Open file in %1 format")
                           .arg(READ_FORMATS[i]));
    RFMTA[i]->setCheckable( true );
    connect( RFMTA[i], SIGNAL(triggered()), this, SLOT( ReadFile() ));
    i++;
  };

  i = 0;
  while ( WRITE_FORMATS[i] != NULL )
  {
    WFMTA.push_back( new QAction( tr( WRITE_FORMATS[i] ), this ) );
    WFMTA[i]->setStatusTip(tr("Write file in %1 format")
                           .arg(WRITE_FORMATS[i]));
    WFMTA[i]->setCheckable( true );
    connect( WFMTA[i], SIGNAL(triggered()), this, SLOT( WriteFile() ));
    i++;
  };

  buildAct = new QAction(tr("&Build Molecule"), this);
  buildAct->setShortcut(tr("Ctrl+B"));
  buildAct->setStatusTip(tr("Edit molecular geometry"));
  connect(buildAct, SIGNAL(triggered()),
          geoeditor, SLOT(show()));

  imgAct = new QAction(tr("&Image"), this);
  imgAct->setStatusTip(tr("Write image file"));
  connect(imgAct, SIGNAL(triggered()), this, SLOT(WriteImage()));


  runAct = new QAction(tr("Run"), this);
  runAct->setShortcut(tr("Ctrl+R"));
  runAct->setStatusTip(tr("Run calculation"));
  connect(runAct, SIGNAL(triggered()), runner, SLOT(show()));

  plotAct = new QAction(tr("Plot"), this);
  plotAct->setShortcut(tr("Ctrl+P"));
  plotAct->setStatusTip(tr("Plot a molecular entity"));
  connect(plotAct, SIGNAL(triggered()), viewer->pmgr, SLOT(Launch()));

  optAct = new QAction(tr("&Geometry"), this);
  optAct->setShortcut(tr("Ctrl+G"));
  optAct->setStatusTip(tr("Structural changes"));
  connect(optAct,SIGNAL(triggered()), optwin, SLOT(Launch()));

  ebeAct = new QAction(tr("&EBEs"), this);
  ebeAct->setShortcut(tr("Ctrl+E"));
  ebeAct->setStatusTip(tr("Electron binding energies"));
  connect(ebeAct,SIGNAL(triggered()), ebewin, SLOT(Launch()));

  coloratomAct = new QAction(tr("Atom type"), this);
  coloratomAct->setStatusTip(tr("Change atom type color"));
  connect(coloratomAct, SIGNAL(triggered()),this,SLOT(ChangeAtomColor()));

  colorbgAct = new QAction(tr("Background"), this);
  colorbgAct->setStatusTip(tr("Change background color"));
  connect(colorbgAct, SIGNAL(triggered()), viewer,SLOT(ChangeBackgroundColor()));

  colortagAct = new QAction(tr("Tags"), this);
  colortagAct->setStatusTip(tr("Change tags color"));
  connect(colortagAct, SIGNAL(triggered()),viewer,SLOT(ChangeTagColor()));

  colorarrowAct = new QAction(tr("Arrows"), this);
  colorarrowAct->setStatusTip(tr("Change arrows color"));
  connect(colorarrowAct, SIGNAL(triggered()),viewer,SLOT(ChangeArrowColor()));

  colorsurfAct = new QAction(tr("Surface"), this);
  colorsurfAct->setStatusTip(tr("Change surface color"));
  connect(colorsurfAct, SIGNAL(triggered()),
          this,SLOT(ChangeSurfaceColor()));

  quitAct = new QAction(tr("&Quit"), this);
  quitAct->setShortcut(tr("Ctrl+Q"));
  quitAct->setStatusTip(tr("Finish the application"));
  connect(quitAct, SIGNAL(triggered()),
          app, SLOT(closeAllWindows()));

  aboutAct = new QAction(tr("&About"), this);
  aboutAct->setStatusTip(tr("Show the application's About box"));
  connect(aboutAct, SIGNAL(triggered()), this, SLOT(About()));

  fileMenu = new QMenu(tr("&File"),this);

  readMenu = new QMenu(tr("&Read"),this);


  for( size_t i = 0 ; i < RFMTA.size() ; i++ )
  {
    readMenu->addAction( RFMTA[i] );
  };

  fileMenu->addMenu(readMenu);

  writeMenu = new QMenu(tr("&Write"),this);
  for( size_t i = 0 ; i < WFMTA.size() ; i++ )
  {
    writeMenu->addAction( WFMTA[i] );
  };
  fileMenu->addMenu(writeMenu);

  fileMenu->addAction( imgAct );
  fileMenu->addSeparator();

  fileMenu->addAction( quitAct );

  editMenu = new QMenu(tr("&Edit"),this);
  editMenu->addAction( buildAct );

  computeMenu = new QMenu(tr("&Analize"),this);
  computeMenu->addAction(runAct);
  computeMenu->addAction(plotAct);
  computeMenu->addAction(optAct);
  computeMenu->addAction(ebeAct);

  colorMenu = new QMenu(tr("Color"),this);
  colorMenu->addAction(coloratomAct);
  colorMenu->addAction(colorbgAct);
  colorMenu->addAction(colortagAct) ;
  colorMenu->addAction(colorarrowAct) ;
  colorMenu->addAction(colorsurfAct) ;

  helpMenu = new QMenu(tr("&Help"),this);
  helpMenu->addAction(aboutAct);

  menuBar()->addMenu(fileMenu);
  menuBar()->addMenu(editMenu);
  menuBar()->addMenu(computeMenu);
  menuBar()->addMenu(colorMenu);
  menuBar()->addMenu(helpMenu);

  //QGridLayout *mainLayout = new QGridLayout;

  scfButton =  new QToolButton( this );
  scfButton->setText( "SCF" ); 
  connect(scfButton,SIGNAL(clicked()),
          this,SLOT(LaunchSCFWin()));

  freqButton = new QToolButton( this );
  freqButton->setText( "Frequencies" );
  connect(freqButton,SIGNAL(clicked()), 
          this,SLOT(LaunchFREQWin()));

  inButton = new QToolButton( this );
  inButton->setText( "Zoom In" );
  connect( inButton,SIGNAL(clicked()),viewer,SLOT(ZoomIn()));
  
  outButton = new QToolButton( this );
  outButton->setText( "Zoom Out" );
  connect( outButton,SIGNAL(clicked()),viewer,SLOT(ZoomOut()));

  giraButton = new QToolButton( this );
  giraButton->setText( "Start Rotation" );
  connect( giraButton,SIGNAL(clicked()),this,SLOT(SetSpin()));

  modelComboBox = new QComboBox;
  modelComboBox->addItem( "Model" );
  modelComboBox->addItem( QString("None"));
  modelComboBox->addItem( QString("Balls & Sticks"));
  modelComboBox->addItem( QString("Balls"));
  modelComboBox->addItem( QString("Sticks"));
  modelComboBox->addItem( QString("Wireframe"));
  modelComboBox->addItem( QString("Van der Waals"));
  modelComboBox->addItem( QString("Custom"));
  connect(modelComboBox, SIGNAL(activated(const QString &)),
          this, SLOT(ChangeModel(const QString &)));

  tagComboBox = new QComboBox;
  tagComboBox->addItem( "Labels" );
  tagComboBox->addItem( "None");
  tagComboBox->addItem( "Symbols");
  tagComboBox->addItem( "Numbers");
  tagComboBox->addItem( "Sym+Num");
  tagComboBox->addItem( "Charges");
  tagComboBox->addItem( "Fukui HOMO");
  tagComboBox->addItem( "Fukui LUMO");
  tagComboBox->addItem( "Fukui Average");
  tagComboBox->addItem( "Fukui Difference");
  connect(tagComboBox, SIGNAL(activated(const QString &)),
          viewer, SLOT(ChangeTags(const QString &)));

  arrowComboBox = new QComboBox;
  arrowComboBox->addItem( "Arrows" );
  arrowComboBox->addItem( "None");
  arrowComboBox->addItem( "Dipole");
  arrowComboBox->addItem( "Forces");
  arrowComboBox->addItem( "NFF HOMO");
  arrowComboBox->addItem( "NFF LUMO");
  arrowComboBox->addItem( "NFFEX");
  arrowComboBox->addItem( "NFFEY");
  arrowComboBox->addItem( "NFFEZ");
  connect(arrowComboBox, SIGNAL(activated(const QString &)),
          viewer, SLOT(ChangeArrows(const QString &)));

  QComboBox *mComboBox = new QComboBox;
  mComboBox->addItem( QString("Measurement"));
  mComboBox->addItem( QString("Distance"));
  mComboBox->addItem( QString("Angle"));
  mComboBox->addItem( QString("Dihedral"));
  mComboBox->addItem( QString("Clean"));
  connect(mComboBox, SIGNAL(activated(const QString &)),
          this, SLOT( SetupMeasurement(const QString &)));

  styleComboBox = new QComboBox;
  styleComboBox->addItems(QStyleFactory::keys());
  connect(styleComboBox, SIGNAL(activated(const QString &)),
          this, SLOT(MakeStyle(const QString &)));


  QWidget *visualizer = new QWidget(0);
  QVBoxLayout *visuallayout = new QVBoxLayout;

  visuallayout->addWidget( modelComboBox,1,Qt::AlignCenter);
  visuallayout->addWidget( scfButton,1,Qt::AlignCenter);
  visuallayout->addWidget( freqButton,1,Qt::AlignCenter);
  visuallayout->addWidget( inButton,1,Qt::AlignCenter);
  visuallayout->addWidget( outButton,1,Qt::AlignCenter);
  visuallayout->addWidget( giraButton,1,Qt::AlignCenter);
  visuallayout->addWidget( styleComboBox,1,Qt::AlignCenter);
  visuallayout->addWidget( tagComboBox,1,Qt::AlignCenter);
  visuallayout->addWidget( arrowComboBox,1,Qt::AlignCenter);
  visuallayout->addWidget( statusbar,1,Qt::AlignCenter);
  visuallayout->addWidget( mComboBox,1,Qt::AlignCenter);

  // styleComboBox->hide();

  visualizer->setLayout(visuallayout);

  setCentralWidget( tabWidget );

  //tabWidget->clear();
  tabWidget->addTab( viewer , QString("Display") );
  tabWidget->addTab( visualizer , QString("Control") );
  //tabWidget->addTab( winplot , QString("Plotter") );
  //tabWidget->addTab( geoeditor , QString("Geometry") );

  //setLayout( mainLayout );

  timer = new QTimer(this);
  drawmol = true;

  setWindowTitle( PROGRAM_NAME );
  move( 0 , 0  );
  resize(640 , 100);
}

Parakata::~Parakata()
{
}


void Parakata::Report( const char *fmt, ... )
{
  va_list args;   
  char    s[256];

  va_start( args , fmt );
  vsprintf(s , fmt , args );
  va_end( args );
 
  if (statusbar) statusbar->showMessage( QString( s ) );
}

void Parakata::About()
{
  char str[MAX_STR_SIZE];

  sprintf(str,"%s\n\nCopyright 2021 Roberto Flores-Moreno\nr.flores@academicos.udg.mx",PROGRAM_NAME);
  QMessageBox::about( this, "About Parakata",str);
}

void Parakata::ChangeModel( const QString &newmodel )
{
  if ( newmodel == "Model" ) return;

  if ( newmodel == "None" ) drawmol = false;
  else drawmol = true;

  if ( newmodel == "Wireframe" )
  {
    geosrs = 0.0;
    geocrs = 0.0;
  }
  else if ( newmodel == "Balls" )
  {
    geosrs = 0.3;
    geocrs = 0.0;
  }
  else if ( newmodel == "Balls & Sticks" )
  {
    geosrs = 0.3;
    geocrs = 0.12;
  }
  else if ( newmodel == "Sticks" )
  {
    geosrs = -0.1;
    geocrs = -0.1;
  }
  else if ( newmodel == "Van der Waals" )
  {
    geosrs = -0.0001;
    geocrs = 0.0;
  }
  else if ( newmodel == "Custom" )
  {
    bool ok;
    geosrs = QInputDialog::getDouble(this, "B&S Model editor", 
                                     "Atom SF: ",0.3, 0.0, 1.0,2,&ok);
    geocrs = QInputDialog::getDouble(this, "B&S Model editor", 
                                     "Bond SF: ",0.12, 0.0, 1.0,2,&ok);
  };
  viewer->Redraw();
}

void Parakata::MakeStyle( const QString &style )
{
  QApplication::setStyle(QStyleFactory::create(style));
}

void Parakata::ReadFile()
{
  char fileformat[MAX_STR_SIZE];
  size_t i;

  for ( i = 0 ; i < RFMTA.size() ; i++ )
  {
    if ( RFMTA[i]->isChecked() ) 
    {
      RFMTA[i]->setChecked( false );
      sprintf(fileformat,"%s",((RFMTA[i]->text()).toStdString()).c_str());
      break;
    };
  };

  QString filter("All Files (*.*)");
  filter +=      "\nOnly    (*.out)";
  filter +=      "\nOnly    (*.inp)";
  filter +=      "\nOnly    (*.mol)";
  filter +=      "\nOnly    (*.mkl)";
  filter +=      "\nOnly    (*.log)";
  filter +=      "\nOnly    (*.new)";
  filter +=      "\nOnly    (*.com)";
  filter +=      "\nOnly    (*.xyz)";
  filter +=      "\nOnly    (*.cif)";
  filter +=      "\nOnly    (*.gaf)";
  filter +=      "\nOnly    (*.nmr)";
  QString fn = QFileDialog::getOpenFileName(this,
                 QString("Read file"),
                 QDir::currentPath(),
                 filter);
  if ( ! ( fn.isEmpty() || fn.isNull() ) )
  {
    char filename[MAX_STR_SIZE];

    sprintf(filename,"%s",(fn.toStdString()).c_str());
    Report( "Reading file %s ...", filename );
    string tit = PROGRAM_NAME;
    tit += ":";
    tit += filename;
    setWindowTitle( tit.c_str() );
    DoReadFile(filename,fileformat);
    Report( "Reading file %s ... DONE", filename );
  };
}

void Parakata::DoReadFile(char* filename, char* format)
{
  reader(filename,format);
  ChangeModel("Balls & Sticks");
  viewer->Redraw();
}

void Parakata::WriteFile()
{
  char sffmt[MAX_STR_SIZE];
  char fileformat[MAX_STR_SIZE];
  QString ffmt;
  size_t i;

  for ( i = 0 ; i < WFMTA.size() ; i++ )
  {
    if ( WFMTA[i]->isChecked() ) 
    {
      WFMTA[i]->setChecked( false );
      sprintf(sffmt,"%s",((WFMTA[i]->text()).toStdString()).c_str());
      strcpy(fileformat,sffmt);
      // File extension comes only from first 3 characters
      sffmt[3] = '\0';
      ffmt = sffmt;
      break;
    };
  };

  QString initialPath = QDir::currentPath() + "/untitled.";
  initialPath = initialPath + ffmt.toLower();
  QString fn = QFileDialog::getSaveFileName( this,
               tr( "Write %1 file").arg(ffmt.toUpper()),
               initialPath, tr("All Files (*.%1)").arg(ffmt.toLower()) );
  if ( ! ( fn.isEmpty() || fn.isNull() ) )
  {
    char filename[MAX_STR_SIZE];

    sprintf(filename,"%s",(fn.toStdString()).c_str());
    Report( "Writing file %s ...", filename );
    writer(filename , fileformat );
    Report( "Writing file %s ... DONE", filename );
  };
}

void Parakata::WriteImage()
{ viewer->image_render->show(); }

void Parakata::NewMolecule()
{
  NATOM = 0;
  viewer->Redraw();
  geoeditor->Update();
}

void Parakata::ChangeAtomColor()
{
  prk->es->exec(); 
  int sel = 54; //prk->es->GetSelectedElement();
  QColor color; 
  color = QColorDialog::getColor( color , this );
  ELEMENT_COLOR[sel][0] = color.red()/255.0;
  ELEMENT_COLOR[sel][1] = color.green()/255.0;
  ELEMENT_COLOR[sel][2] = color.blue()/255.0;
  viewer->Redraw();
}

void Parakata::ChangeSurfaceColor()
{
  prk->SetCursor( "Selection" );
  prk->status |= SFSURF;
  Report("Select a surface");
}

void Parakata::LaunchSCFWin( void ) 
{ 
  if (prk)
  {
    if (prk->scfwin) 
    {
      prk->scfwin->Launch(); 
    }
  }
}
  

void Parakata::LaunchFREQWin( void ) 
{ 
  if (prk)
  {
    if (prk->freqwin) 
    {
      prk->freqwin->Launch(); 
    }
  }
}

void Parakata::SetupMeasurement( const QString &s )
{
  if ( s == "Measurement" ) return;
  if ( s == "Distance" && NATOM > 1 ) 
  {
    prk->SetCursor( "Selection" );
    prk->status |= SFRV;
    Report("Select two atoms");
  }
  else if ( s == "Angle" && NATOM > 2)
  {
    prk->SetCursor( "Selection" );
    prk->status |= SFAV;
    Report("Select three atoms");
  }
  else if ( s == "Dihedral" && NATOM > 3 )
  {
    prk->SetCursor( "Selection" );
    prk->status |= SFDV;
    Report("Select four atoms");
  }
  else if ( s == "Clean" && viewer ) 
  {
    viewer->monitored_bonds.clear();
    viewer->monitored_angles.clear();
    viewer->monitored_dihedrals.clear();
    viewer->Redraw();
  };
  viewer->image_render->show(); 
}

void Parakata::SetSpin(void)
{
  if (giraButton->text()=="Start Rotation")
  {
    giraButton->setText("Stop Rotation");
    connect(timer,SIGNAL(timeout()),
            this,SLOT(AdvanceSpin()));
    AdvanceSpin();
  }
  else if (giraButton->text()=="Stop Rotation")
  {
    giraButton->setText("Start Rotation");
    disconnect(timer,SIGNAL(timeout()),
               this,SLOT(AdvanceSpin()));
  }
}

void Parakata::AdvanceSpin(void)
{
  viewer->Rotate( 5.0, 0.0, 1.0, 0.0); 
  viewer->Redraw();
  timer->start(40);
}

void Parakata::LaunchGeoEditor() { geoeditor->Launch(); }


void Parakata::SetCursor( const char *type )
{
  QString typestr = type;
  QCursor cursor;
  if ( typestr == "Normal" )
  {
    cursor = QCursor( Qt::ArrowCursor );
  }
  else if ( typestr == "Wait" )
  {
    cursor = QCursor( Qt::WaitCursor );
  }
  else if ( typestr == "Selection" )
  {
    cursor = QCursor( Qt::PointingHandCursor );
  };
  setCursor( cursor );
} 
