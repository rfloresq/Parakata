// ******************************************************************
// Purpose: PlotManager class functions.
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

#include <PlotManager.h>
#include <Viewer.h>
#include <Parakata.h>
#include <Integrator.h>
#include <Vector.h>
#include <macros.h>


#define TOEV(a) ((a)*27.21138)

#define ISO             0
#define CONTOUR         1
#define COLOURED_ISO    2
#define COLOURED_PLANE  3
#define CURVE           4
#define NPLOT_TYPES     5

#define SOLID_STYLE         0
#define POINTS_STYLE        1
#define TRANSLUCENT_STYLE   2
#define MESH_STYLE          3
#define NPLOT_STYLES        4

#define MAXGPL 128

using namespace std;

static int PLOT_TYPE;
static int EUCLID = 3;
static char* typenames[NPLOT_TYPES] = {"Isosurface","Contours",
"Coloured Iso-Density","Coloured Plane","Curve"}; 
static char* stylenames[NPLOT_STYLES] = {"Solid","Points", "Mesh"};
static char *flagbox_names[PLOT_TAB_NROW][PLOT_TAB_NCOL] = {
{"Shiny","Reference","Draw Box"},
{"Hide all","Hold Previous","Nothing"},
{"Approximate","Average","Nothing"},
{"Nothing","Nothing","Nothing"}};

double *DENSITY_MATRIX = 0;
double *AUX_MATRIX = 0;
double *FUKUIL_MATRIX = 0;
double *FUKUIR_MATRIX = 0;
double *DMR1 = 0;
double *DMR2 = 0;
double *DMR3 = 0;

void add_to_density_matrix(double *P, int from, int to)
{
  int imo,ibas,jbas,k;

  for (imo=from;imo<=to;imo++)
  {
    MO mo = MOS[imo];
    k = 0;
    for (ibas=0;ibas<NBASIS;ibas++)
    {
      for (jbas=ibas;jbas<NBASIS;jbas++)
      {
        P[k] += mo.c[ibas]*mo.c[jbas]*mo.occ;
        // Off diagonal elements multiplied by 2
        if (jbas!=ibas)
          P[k] += mo.c[ibas]*mo.c[jbas]*mo.occ;
        k++;
      }
    }
  }
}

void build_density_matrix(double *P)
{
  // Refresh memory
  int homo,ibas,imo,jbas,k,lumo;

  k = 0;
  for (ibas=0;ibas<NBASIS;ibas++)
    for (jbas=ibas;jbas<NBASIS;jbas++)
    {
      P[k] = 0.0;
      k++;
    }

  // alpha
  lumo = 0;
  while (MOS[lumo].occ>0.0) lumo++;
  homo = lumo - 1;
  add_to_density_matrix(P, 0, homo);

  // beta
  for (imo=0;imo<NORB;imo++)
  {
    if (strcasestr(MOS[imo].spin,"beta")!=NULL) 
    {
      lumo = imo;
      while (MOS[lumo].occ>0.0) lumo++;
      homo = lumo - 1;
      add_to_density_matrix(P, imo, homo);
      break;
    }
  }
}



void build_rod_matrix(double *matrix, int id)
{
  int ll,ul,ibas,imo,jbas,k;

  k = 0;
  for (ibas=0;ibas<NBASIS;ibas++)
    for (jbas=ibas;jbas<NBASIS;jbas++)
    {
      matrix[k] = 0.0;
      k++;
    }

  ll = id;
  ul = id;
  while (ABS(MOS[ll-1].energy-MOS[id].energy)<0.001) ll--;
  while (ABS(MOS[ul+1].energy-MOS[id].energy)<0.001) ul++;
  for (imo=ll;imo<=ul;imo++)
  {
    MO mo = MOS[imo];
    k = 0;
    for (ibas=0;ibas<NBASIS;ibas++)
    {
      for (jbas=ibas;jbas<NBASIS;jbas++)
      {
        matrix[k] += mo.c[ibas]*mo.c[jbas]/((double)(ul-ll+1));
        if (jbas!=ibas)
          matrix[k] += mo.c[ibas]*mo.c[jbas]/((double)(ul-ll+1));
        k++;
      }
    }
  }
}
void build_spin_matrix(double *matrix)
{
  int ibas,imo,jbas,k;
  double w;

  k = 0;
  for (ibas=0;ibas<NBASIS;ibas++)
    for (jbas=ibas;jbas<NBASIS;jbas++)
    {
      matrix[k] = 0.0;
      k++;
    }

  for (imo=0;imo<NORB;imo++)
  {
    MO mo = MOS[imo];
    if ( MOS[imo].occ > 0.00000001 )
    {
      if (strcasestr(MOS[imo].spin,"alpha")!=NULL) 
      {
        w = 1.0;
      }
      else if (strcasestr(MOS[imo].spin,"beta")!=NULL) 
      {
        w = -1.0;
      }
      else
      {
        w = 0.0;
      }
      w = w*MOS[imo].occ;
      k = 0;
      for (ibas=0;ibas<NBASIS;ibas++)
      {
        for (jbas=ibas;jbas<NBASIS;jbas++)
        {
          matrix[k] += w*mo.c[ibas]*mo.c[jbas];
          if (jbas!=ibas)
            matrix[k] += w*mo.c[ibas]*mo.c[jbas];
          k++;
        }
      }
    }
  }
}



PlotManager::PlotManager( Parakata *is , Viewer* view)
           : QWidget( 0 )
{
  prk = is;
  viewer = view;
  gi = new Integrator();

  tabWidget = new QTabWidget;


  // Initialize
  plane_axis = 2;
  lastplaneaxis = 2;
  plane_point = 0;

  QGroupBox *settBox = new QGroupBox(tr("Settings"));
  QHBoxLayout *settLayout = new QHBoxLayout;
  settLayout->addStretch(1);

  type = 0;
  QComboBox* typeC = new QComboBox( settBox );
  for ( int i = 0 ; i < NPLOT_TYPES ; i++ )
  {
    typeC->addItem( typenames[i] );
  }; 
  connect( typeC , SIGNAL( activated( const QString & ) ) , 
           this , SLOT( ChangeType( const QString & ) ) );
  settLayout->addWidget( typeC );

  style = 0;
  QComboBox* styleC = new QComboBox( settBox );
  for ( int i = 0 ; i < NPLOT_STYLES ; i++ )
  {
    styleC->addItem( stylenames[i] );
  }; 
  connect( styleC , SIGNAL( activated( const QString & ) ) , 
           this , SLOT( ChangeStyle( const QString & ) ) );
  settLayout->addWidget( styleC );

  QPushButton *freshB = new QPushButton(tr("Refresh"));
  connect( freshB, SIGNAL(clicked()), this, SLOT(Refresh()));
  settLayout->addWidget( freshB );

  QPushButton *plotB = new QPushButton(tr("Build"));
  connect( plotB, SIGNAL(clicked()), this, SLOT(Plot()));
  settLayout->addWidget( plotB );

  QPushButton *doneB = new QPushButton(tr("Done"));
  connect( doneB, SIGNAL(clicked()), this, SLOT(close()));
  settLayout->addWidget( doneB );

  settBox->setLayout( settLayout );

  QGroupBox *valBox = new QGroupBox(tr("Values: ISO/MIN/MAX/STEP"));
  QHBoxLayout *valLayout = new QHBoxLayout;
  valLayout->addStretch(1);

  double valrange = 10000.0;
  isovalue = 0.02;
  QDoubleSpinBox *isoS = new QDoubleSpinBox( valBox );
  isoS->setRange( -valrange , valrange );
  isoS->setSingleStep( 0.00001 );
  isoS->setDecimals( 10 );
  isoS->setValue( isovalue );
  connect( isoS , SIGNAL(valueChanged( double )) ,
           this , SLOT( ChangeIso( double ) ) );
  valLayout->addWidget( isoS );

  fmin = -1.0;
  QDoubleSpinBox *cminS = new QDoubleSpinBox( valBox );
  cminS->setRange( -valrange , valrange );
  cminS->setSingleStep( 0.01 );
  cminS->setDecimals( 10 );
  cminS->setValue( fmin );
  valLayout->addWidget( cminS );
  connect( cminS , SIGNAL(valueChanged( double )) ,
           this , SLOT( ChangeContourMin( double ) ) );

  fmax = 1.0;
  QDoubleSpinBox *cmaxS = new QDoubleSpinBox( valBox );
  cmaxS->setRange( -valrange , valrange );
  cmaxS->setSingleStep( 0.01 );
  cmaxS->setValue( fmax );
  cmaxS->setDecimals( 10 );
  valLayout->addWidget( cmaxS );
  connect( cmaxS , SIGNAL(valueChanged( double )) ,
           this , SLOT( ChangeContourMax( double ) ) );

  fstep = 0.1;
  QDoubleSpinBox *cstepS = new QDoubleSpinBox( valBox );
  cstepS->setRange( 0.00000001 , valrange );
  cstepS->setSingleStep( 0.001 );
  cstepS->setDecimals( 10 );
  cstepS->setValue( fstep );
  valLayout->addWidget( cstepS );
  connect( cstepS , SIGNAL(valueChanged( double )) ,
           this , SLOT( ChangeContourStep( double ) ) );

  valBox->setLayout(valLayout);

  QGroupBox *drawBox = new QGroupBox(tr("Drawing: Grid/Line/Point"));
  QHBoxLayout *drawLayout = new QHBoxLayout;
  drawLayout->addStretch(1);

  boxmesh = 0.6;
  QDoubleSpinBox *meshS = new QDoubleSpinBox( drawBox );
  meshS->setRange( 0.01 , 10.0 );
  meshS->setSingleStep( 0.01 );
  meshS->setDecimals( 2 );
  meshS->setValue( boxmesh );
  drawLayout->addWidget( meshS );
  connect( meshS , SIGNAL(valueChanged( double )) ,
           this , SLOT( ChangeMesh( double ) ) );

  line_width = 2;
  QSpinBox *lineS = new QSpinBox( drawBox );
  lineS->setRange( 1 , 10 );
  lineS->setSingleStep( 1 );
  lineS->setValue( line_width );
  drawLayout->addWidget( lineS );
  connect( lineS , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangeLineWidth( int ) ) );

  point_size = 3;
  QSpinBox *psS = new QSpinBox( drawBox );
  psS->setRange( 1 , 10 );
  psS->setSingleStep( 1 );
  psS->setValue( point_size );
  drawLayout->addWidget( psS );
  connect( lineS , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangePointSize( int ) ) );

  drawBox->setLayout( drawLayout );

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget( tabWidget );
  mainLayout->addWidget( drawBox );
  mainLayout->addWidget( valBox );
  mainLayout->addWidget( settBox );
  SetupBox( mainLayout );
  setLayout( mainLayout );
  setWindowTitle(tr( "Plot Manager" ));
  tabs.clear();
  listados.clear();
  plot_entries.clear();
  use_lists = true;
  cp.clear();
  bondpath.clear();
  grdvec.clear();
}



PlotManager::~PlotManager()
{
  tabs.clear();
  listados.clear();
  plot_entries.clear();
  cp.clear();
  bondpath.clear();
}

int PlotManager::NSurf()
{
  return (int) plot_entries.size();
}

void PlotManager::PlotISO( int n )
{
  EUCLID = 3;
  ChangeBox(0.0);

  if ( plot_entries[n].style == SOLID_STYLE )
  {
    glBegin( GL_TRIANGLES );
      DoPlot(n); 
    glEnd();
  }
  else if ( plot_entries[n].style == TRANSLUCENT_STYLE )
  {
    glEnable( GL_CULL_FACE );
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE);
    glBegin( GL_TRIANGLES );
      DoPlot(n); 
    glEnd();
    glDisable(GL_BLEND);
    glDisable( GL_CULL_FACE );
  }
  else if ( plot_entries[n].style == MESH_STYLE )
  {
    glLineWidth( plot_entries[n].line_width );
    glBegin( GL_LINES );
      DoPlot(n); 
    glEnd();
  }
  else if ( plot_entries[n].style == POINTS_STYLE )
  {
    glPointSize( plot_entries[n].point_size );
    glBegin( GL_POINTS );
      DoPlot(n); 
    glEnd();
  };
}

void PlotManager::PlotCURVE( int n )
{
  EUCLID = 1;
  ChangeBox(0.0);

  glLineWidth( plot_entries[n].line_width );
  glDisable( GL_LIGHTING );
  glBegin( GL_LINES );
    DoPlot(n); 
  glEnd();
  glEnable( GL_LIGHTING );
}

void PlotManager::PlotCONTOUR( int n )
{
  double ccolor[3];

  EUCLID = 2;
  ChangeBox(0.0);

  if ( strstr( plot_entries[n].name , "Neg" ) != NULL ) return;
  glDisable( GL_LIGHTING );
  glLineWidth( plot_entries[n].line_width );
  plot_entries[n].iso = fmin;
  while ( plot_entries[n].iso < fmax )
  {
    RealToRGB( plot_entries[n].iso, ccolor );
    glColor3dv( ccolor );
    glBegin( GL_LINES );
      DoPlot(n); 
    glEnd();
    plot_entries[n].iso += fstep;
  };
  glEnable( GL_LIGHTING );
}

void PlotManager::PlotCOLOURED_ISO( int n )
{
  EUCLID = 3;
  ChangeBox(0.0);

  plot_entries[n].style = SOLID_STYLE;
  glEnable( GL_LIGHTING );
  glBegin( GL_TRIANGLES );
    DoPlot(n); 
  glEnd();
}

void PlotManager::PlotCOLOURED_PLANE( int n )
{
  EUCLID = 2;
  ChangeBox(0.0);

  if ( strstr( plot_entries[n].name , "Neg" ) != NULL ) return;
  plot_entries[n].style = SOLID_STYLE;
  glDisable( GL_LIGHTING );
  glBegin( GL_TRIANGLES );
    DoPlot(n); 
  glEnd();
  glEnable( GL_LIGHTING );
}

void PlotManager::PlotSurface( int n )
{
  if ( (n < NSurf()) && !hide_surfaces )
  {
    GLfloat noneDir[4] = {0.2, 0.2, 0.2, 0.0};
    GLfloat whiteDir[4] = {2.0, 2.0, 2.0, 1.0};
    glColor3dv( plot_entries[n].color );
    if ( ! plot_entries[n].shiny )
    {
      glMaterialfv(GL_FRONT, GL_SPECULAR, noneDir);
    };
    if ( (! plot_entries[n].list ) || ( ! use_lists ) )
    {
      prk->Report( "Generating plot %s ...",
                             plot_entries[n].name );
      PLOT_TYPE = plot_entries[n].type;
      if (use_lists)
      {
        plot_entries[n].list = glGenLists( 1 );
        glNewList( plot_entries[n].list , GL_COMPILE );
      }
      if ( PLOT_TYPE == ISO ) PlotISO( n );
      else if  ( PLOT_TYPE == CURVE ) PlotCURVE( n );
      else if  ( PLOT_TYPE == CONTOUR ) PlotCONTOUR( n );
      else if  ( PLOT_TYPE == COLOURED_ISO ) PlotCOLOURED_ISO( n );
      else if  ( PLOT_TYPE == COLOURED_PLANE ) PlotCOLOURED_PLANE( n );
      if (use_lists) glEndList();
      prk->Report( "Generating plot %s ... DONE",
                             plot_entries[n].name );
    };
    if ( use_lists ) glCallList( plot_entries[n].list );
    if ( ! plot_entries[n].shiny )
      glMaterialfv(GL_FRONT, GL_SPECULAR, whiteDir);
  }
  else
  {
    if ( drawref )
    {
      Vector vx(3.0,0.0,0.0);
      Vector vy(0.0,3.0,0.0);
      Vector vz(0.0,0.0,3.0);
   // glDisable( GL_LIGHTING );
      glEnable (GL_LIGHTING);
      glColor4d( 0.5 , 0.5 , 0.5 , 0.0 );
      glLineWidth( line_width );
      glBegin( GL_LINE_STRIP );
        glVertex3d( 0.0 , 0.0 , 0.0 );
        glVertex3d( vx[0] , vx[1] , vx[2] );
        glVertex3d( 0.0 , 0.0 , 0.0 );
        glVertex3d( vy[0] , vy[1] , vy[2] );
        glVertex3d( 0.0 , 0.0 , 0.0 );
        glVertex3d( vz[0] , vz[1] , vz[2] );
      glEnd();
      glEnable( GL_LIGHTING );
      viewer->Write("X",vx,32);
      viewer->Write("Y",vy,32);
      viewer->Write("Z",vz,32);
    }
    if ( drawbox )
    {
    //  glDisable( GL_LIGHTING );
      glEnable(GL_LIGHTING);
      glColor4d( 0.5 , 0.5 , 0.5 , 0.0 );
      glLineWidth( line_width );
      glBegin( GL_LINE_STRIP );
      glVertex3d( LATTICE[0] , LATTICE[1] , LATTICE[2] );
      glVertex3d( LATTICE[3] , LATTICE[4] , LATTICE[5] );
      if ( EUCLID > 1 )
      {
        glVertex3d( LATTICE[3] + LATTICE[6] - LATTICE[0] , 
                    LATTICE[4] + LATTICE[7] - LATTICE[1] , 
                    LATTICE[5] + LATTICE[8] - LATTICE[2] );
        glVertex3d( LATTICE[6] , LATTICE[7] , LATTICE[8] );
        glVertex3d( LATTICE[0] , LATTICE[1] , LATTICE[2] );
        if ( EUCLID > 2 )
        {
          glVertex3d(LATTICE[9] , LATTICE[10], LATTICE[11] );
          glVertex3d(LATTICE[6] + LATTICE[9] - LATTICE[0] , 
                     LATTICE[7] + LATTICE[10]- LATTICE[1] , 
                     LATTICE[8] + LATTICE[11]- LATTICE[2] );
          glVertex3d(LATTICE[6] , LATTICE[7] , LATTICE[8] );
          glVertex3d(LATTICE[6] + LATTICE[9] - LATTICE[0] , 
                     LATTICE[7] + LATTICE[10]- LATTICE[1] , 
                     LATTICE[8] + LATTICE[11]- LATTICE[2] );
          glVertex3d(LATTICE[3] + LATTICE[6] +LATTICE[9]-2.0*LATTICE[0],
                     LATTICE[4] + LATTICE[7] +LATTICE[10]-2.0*LATTICE[1],
                     LATTICE[5] + LATTICE[8] +LATTICE[11]-2.0*LATTICE[2]);
          glVertex3d(LATTICE[3] + LATTICE[6] - LATTICE[0] , 
                     LATTICE[4] + LATTICE[7] - LATTICE[1] , 
                     LATTICE[5] + LATTICE[8] - LATTICE[2] );
          glVertex3d(LATTICE[3] + LATTICE[6] +LATTICE[9]-2.0*LATTICE[0] , 
                     LATTICE[4] + LATTICE[7] +LATTICE[10]-2.0*LATTICE[1], 
                     LATTICE[5] + LATTICE[8] +LATTICE[11]-2.0*LATTICE[2]);
          glVertex3d(LATTICE[3] + LATTICE[9] - LATTICE[0] , 
                     LATTICE[4] + LATTICE[10]- LATTICE[1] , 
                     LATTICE[5] + LATTICE[11]- LATTICE[2] );
          glVertex3d(LATTICE[3] , LATTICE[4] , LATTICE[5] );
          glVertex3d(LATTICE[3] + LATTICE[9] - LATTICE[0] , 
                     LATTICE[4] + LATTICE[10]- LATTICE[1] , 
                     LATTICE[5] + LATTICE[11]- LATTICE[2] );
          glVertex3d(LATTICE[9] , LATTICE[10], LATTICE[11] );
        };
      };
      glEnd();
      glEnable( GL_LIGHTING );
    };
  }
}

void PlotManager::GetColor( int n , double *color )
{
  if ( n < NSurf() )
  {
    for ( int i = 0 ; i < 4 ; i++ )
    {
      color[i] = plot_entries[n].color[i];
    };
  };
}

void PlotManager::SetColor( int n , double *color )
{
  if ( n < NSurf() )
  {
    for ( int i = 0 ; i < 4 ; i++ )
    {
      plot_entries[n].color[i] = color[i];
    };
  };
}

void PlotManager::NewList( char** entries , int nentries , char* name ) 
{
  tabs.push_back( new QWidget( 0 ) );
  QListWidget *slListBox = new QListWidget;
  QStringList sllist;

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

void PlotManager::Refresh()
{
  int n;

  n = NSurf()-1;
  if (n<0) return;

  if (strstr(plot_entries[n].name,"Neg")!=NULL)
  {
    plot_entries[n].iso = -isovalue;
    plot_entries[n-1].iso = isovalue;
    plot_entries[n-1].list = 0;
  }
  else plot_entries[n].iso = isovalue;
  plot_entries[n].list = 0;

  viewer->Redraw();
}

void PlotManager::Plot()
{
  if ( ! hold )
  {
    plot_entries.clear();
  };

  PlotEntry pe;
  int table = tabWidget->currentIndex();
  if ( listados[table]->currentItem() == NULL ) 
  {
    errmsg(__FILE__,__LINE__, 
           "You must choose a field first" , 0 );
    return;
  };
  QString st = (listados[table]->currentItem())->text();
  if ( st.isNull() )
  {
  }
  strcpy( pe.name , (st.toStdString()).c_str());
  pe.type = type; 
  pe.style = style; 
  if ( strstr(pe.name,"GRID") != 0 )
  {
    pe.style = POINTS_STYLE; 
  }
  pe.shiny = shiny;
  pe.line_width = line_width;
  pe.point_size = point_size;

  // Green for orbitals, yellow for others
  if ( strncmp(pe.name,"AMO",3) == 0 || strncmp(pe.name,"BMO",3) == 0 ||
       strncmp(pe.name,"ADY",3) == 0 || strncmp(pe.name,"BDY",3) == 0 || 
       strncmp(pe.name,"AKMOR",5) == 0 || strncmp(pe.name,"BKMOR",5) == 0 || 
       strncmp(pe.name,"AKMOI",5) == 0 || strncmp(pe.name,"BKMOI",5) == 0 || 
       strncmp(pe.name,"AO",2) == 0 ) 
  {
    pe.color[0] = 0.0;
    pe.color[1] = 0.0;
    pe.color[2] = 1.0;
    pe.color[3] = 0.0;
  }
  else
  {
    pe.color[0] = 1.0;
    pe.color[1] = 1.0;
    pe.color[2] = 0.0;
    pe.color[3] = 0.0;
  }

  pe.iso = isovalue;  
  pe.list = 0; // No list available yet
  plot_entries.push_back( pe );
  int n = plot_entries.size()-1;
  BuildCube(n);

  if ( strstr(pe.name,"DENSITY") == 0 && 
       strstr(pe.name,"REACTIVE") == 0 &&
       strstr(pe.name,"ELF") == 0 &&
       strstr(pe.name,"DLE") == 0 ) 
  {
    sprintf( pe.name , "%s Neg", (st.toStdString()).c_str());
    pe.iso = -pe.iso;
    // Red for orbitals, blue for others
    if ( strncmp(pe.name,"AMO",3) == 0 || strncmp(pe.name,"BMO",3) == 0 ||
         strncmp(pe.name,"ADY",3) == 0 || strncmp(pe.name,"BDY",3) == 0 ||
         strncmp(pe.name,"AKMOR",5) == 0 || strncmp(pe.name,"BKMOR",5) == 0 ||
         strncmp(pe.name,"AKMOI",5) == 0 || strncmp(pe.name,"BKMOI",5) == 0 ||
         strncmp(pe.name,"AO",2) == 0) 
    {
      pe.color[0] = 1.0;
      pe.color[1] = 0.0;
      pe.color[2] = 0.0;
      pe.color[3] = 0.0;
    }
    else
    {
      pe.color[0] = 0.0;
      pe.color[1] = 0.0;
      pe.color[2] = 1.0;
      pe.color[3] = 0.0;
    }
    plot_entries.push_back( pe );
    n = plot_entries.size()-1;
    BuildCube(n);
  };

  viewer->Redraw();
}

void PlotManager::Launch( void )
{
  if ( isVisible() )
  {
    hide();
    plot_entries.clear();
    return;
  }
  else
  {
    if (DENSITY_MATRIX) free(DENSITY_MATRIX);
    DENSITY_MATRIX = (double*) malloc(NBASIS*NBASIS*sizeof(double));
    build_density_matrix(DENSITY_MATRIX);
    if (AUX_MATRIX) free(AUX_MATRIX);
    AUX_MATRIX = (double*) malloc(NBASIS*NBASIS*sizeof(double));
  };
  int imo;
  int NMAX = MAX(NBASIS,MAX(NORB,50));
  char *names[NMAX];

  tabs.clear();
  listados.clear();
  tabWidget->clear();

  for ( imo = 0 ; imo  < NMAX ; imo++ )
  { 
    names[imo] = (char*) malloc(256*sizeof(char));
    if ( ! names[imo] )
    {
      errmsg(__FILE__,__LINE__, 
              "Memory allocation failed, aborting" , 2);
      exit( EXIT_FAILURE );
    }
  }

  int ns = 0; 
  //if (rhog.size()>0||rhod.size()>0) 
  //{ 
    sprintf( names[ns],"AZD FIELD");
    ns++;
    sprintf( names[ns],"LOADED FIELD");
    ns++;
  //}
  if (NORB>0) 
  {
    sprintf( names[ns],"DENSITY");
    ns++;
    sprintf( names[ns],"MEP");
    ns++;
    sprintf( names[ns],"LAPLACIAN");
    ns++;
    sprintf( names[ns],"KED");
    ns++;
    sprintf( names[ns],"ELF");
    ns++;
    sprintf( names[ns],"DLE");
    ns++;
  }
  if (FUKUIL_MATRIX) 
  {
    sprintf( names[ns],"FUKUI LEFT");  
    ns++;
  }
  if (FUKUIR_MATRIX) 
  {
    sprintf( names[ns],"FUKUI RIGHT");  
    ns++;
  }
  if (FUKUIL_MATRIX&&FUKUIR_MATRIX) 
  {
    sprintf( names[ns],"FUKUI AVERAGE");  
    ns++;
    sprintf( names[ns],"FUKUI DUAL");  
    ns++;
  }
  if (DMR1) 
  {
    sprintf( names[ns],"DR 1");  
    ns++;
  }
  if (DMR2) 
  {
    sprintf( names[ns],"DR 2");  
    ns++;
  }
  if (DMR3) 
  {
    sprintf( names[ns],"DR 3");  
    ns++;
  }
  if (NORB > 0)
  {
    sprintf( names[ns],"SPIN");
    ns++;
    sprintf( names[ns],"REACTIVE ORBITAL SPACE");
    ns++;
    sprintf( names[ns],"SHAPE ENTROPY");
    ns++;
    sprintf( names[ns],"SHAPE ENTROPY ALT");
    ns++;
    if (FUKUIL_MATRIX)
    {
      sprintf( names[ns],"IONIZATION RESPONSE LEFT");
      ns++;
    }
    if (FUKUIR_MATRIX)
    {
      sprintf( names[ns],"IONIZATION RESPONSE RIGHT");
      ns++;
    }
  }

  NewList( names , ns , "SCALAR" );

  for ( int iao = 0 ; iao < NBASIS ; iao++ )
  {
    sprintf( names[iao],"AO%d \n",iao+1);
  }
  NewList( names , NBASIS , "AOs" );

  if (NORB>0) 
  {
    int jmo;
    int fa = -1;
    int fb = -1;
    int fad = -1;
    int fbd = -1;
    for ( imo = 0 ; imo < NORB ; imo++ )
    {
      if (strcmp(MOS[imo].spin,"ALPHA")==0) 
      {
        fb = fb + 1;
        fad = fad + 1;
        fbd = fbd + 1;
        jmo = imo - fa;
        sprintf( names[imo],"AMO%d  %20.10f     %20.10f(%7.3f eV)\n",
               jmo,MOS[imo].occ,MOS[imo].energy,TOEV(MOS[imo].energy));
      }
      else if (strcmp(MOS[imo].spin,"BETA")==0) 
      {
        fad = fad + 1;
        fbd = fbd + 1;
        jmo = imo - fb;
        sprintf( names[imo],"BMO%d  %20.10f     %20.10f(%7.3f eV)\n",
               jmo,MOS[imo].occ,MOS[imo].energy,TOEV(MOS[imo].energy));
      }
      else if (strcmp(MOS[imo].spin,"ALPHAD")==0) 
      {
        fbd = fbd + 1;
        jmo = imo - fad;
        sprintf( names[imo],"ADY%d  ps(%5.3f)  pe(%7.3f eV)\n",
               jmo,MOS[imo].occ,MOS[imo].energy);
      }
      else if (strcmp(MOS[imo].spin,"BETAD")==0) 
      {
        jmo = imo - fbd;
        sprintf( names[imo],"BDY%d  ps(%5.3f)  pe(%7.3f eV)\n",
               jmo,MOS[imo].occ,MOS[imo].energy);
      }
      else if (strcmp(MOS[imo].spin,"UNDEF")==0) 
      {
        sprintf( names[imo],"UMO%d  %20.10f     %20.10f(kt=%7.3f eV)\n",
        imo+1,MOS[imo].occ,MOS[imo].energy,TOEV(MOS[imo].energy));
      }
    }
    NewList( names , NORB , "MOs" );
  }

  if (KNORB>0) 
  {
    int jmo;
    int fa = -1;
    int fb = -1;
    int fad = -1;
    int fbd = -1;
    for ( imo = 0 ; imo < KNORB ; imo++ )
    {
      if (strcmp(KMOSR[imo].spin,"ALPHA")==0) 
      {
        fb = fb + 1;
        fad = fad + 1;
        fbd = fbd + 1;
        jmo = imo - fa;
        sprintf( names[imo],"AKMOR%d  %20.10f     %20.10f(%7.3f eV)\n",
               jmo,KMOSR[imo].occ,KMOSR[imo].energy,TOEV(KMOSR[imo].energy));
      }
      else if (strcmp(KMOSR[imo].spin,"BETA")==0) 
      {
        fad = fad + 1;
        fbd = fbd + 1;
        jmo = imo - fb;
        sprintf( names[imo],"BKMOR%d  %20.10f     %20.10f(%7.3f eV)\n",
               jmo,KMOSR[imo].occ,KMOSR[imo].energy,TOEV(KMOSR[imo].energy));
      }
    }
    NewList( names , KNORB , "RKMOs" );
    for ( imo = 0 ; imo < KNORB ; imo++ )
    {
      if (strcmp(KMOSI[imo].spin,"ALPHA")==0) 
      {
        fb = fb + 1;
        fad = fad + 1;
        fbd = fbd + 1;
        jmo = imo - fa;
        sprintf( names[imo],"AKMOI%d  %20.10f     %20.10f(%7.3f eV)\n",
               jmo,KMOSI[imo].occ,KMOSI[imo].energy,TOEV(KMOSI[imo].energy));
      }
      else if (strcmp(KMOSI[imo].spin,"BETA")==0) 
      {
        fad = fad + 1;
        fbd = fbd + 1;
        jmo = imo - fb;
        sprintf( names[imo],"BKMOI%d  %20.10f     %20.10f(%7.3f eV)\n",
               jmo,KMOSI[imo].occ,KMOSI[imo].energy,TOEV(KMOSI[imo].energy));
      }
    }
    NewList( names , KNORB , "IKMOs" );
  }

/*
  if (NDORB>0) 
  {
    int jmo;
    int lastalpha = -1;
    for ( imo = 0 ; imo < NDORB ; imo++ )
    {
      if (toupper(DMOS[imo].spin[0])=='A') 
      {
        lastalpha++;
        jmo = imo + 1;
      }
      else jmo = imo - lastalpha;
      sprintf( names[imo],"%cDY%d ps(%5.3f) pe(%7.3f eV)\n",
               DMOS[imo].spin[0],jmo,DMOS[imo].occ,TOEV(DMOS[imo].energy));
    }
    NewList( names , NDORB , "Dyson MOs" );
  }
*/

  for ( imo = 0 ; imo  < NMAX ; imo++ )
  { 
    free ( names[imo] ); 
  };

  // Adapt Box
  //KPU
  FitBox();

  show();
}

void PlotManager::ChangeMesh( double newval ) 
{ 
  boxmesh = newval; 
  ChangeBox( 0.0 );
}

void PlotManager::ChangeIso( double newiso ) 
{ isovalue = newiso; }

void PlotManager::ChangeContourMin( double cmin ) 
{ fmin = cmin; }

void PlotManager::ChangeContourMax( double cmax ) 
{ fmax = cmax; }

void PlotManager::ChangeContourStep( double cstep ) 
{ fstep = cstep; }

void PlotManager::ChangeLineWidth( int lw ) 
{ line_width = lw; }

void PlotManager::ChangePointSize( int ps ) 
{ point_size = ps; }

void PlotManager::ChangeType( const QString &newtype )
{
  for ( int i = 0 ; i < NPLOT_TYPES ; i++ )
  {
    if ( newtype == typenames[i] )
    {
      type = i;
      break;
    };
  };
}

void PlotManager::ChangeStyle( const QString &newstyle )
{
  for ( int i = 0 ; i < NPLOT_STYLES ; i++ )
  {
    if ( newstyle == stylenames[i] )
    {
      style = i;
      break;
    };
  };
}

void PlotManager::ChangePlaneAxis( int npa )
{
  plane_axis = npa-1;
 // if (plane_axis<0) plane_axis = 0;
 // if (plane_axis>2) plane_axis = 2;
  {
      if (Gaussbool == true)
      {
          switch (lastplaneaxis)
          {
          case 0:
              if (plane_axis == 1)
              {
                  lastplaneaxis = 1;
                  viewer->Redraw();
                  Plot();
              }
              if (plane_axis == 3)
              {
                  lastplaneaxis = 3;
                  viewer->Redraw();
                  Plot();
              }
              break;
          case 1:
              if (plane_axis == 0)
              {
                  lastplaneaxis = 0;
                  viewer->Redraw();
                  Plot();
              }
              if (plane_axis == 2)
              {
                  lastplaneaxis = 2;
                  viewer->Redraw();
                  Plot();
              }
              break;
              case 2:
              if (plane_axis == 1)
              {
                  lastplaneaxis = 1;
                  viewer->Redraw();
                  Plot();
              }
              break;
          }

      }
  }
}

void PlotManager::ChangePlanePoint( int npp )
{
  plane_point = npp-1;
  viewer->Redraw();
  Plot();
  }

void PlotManager::SetUsageOfLists( bool use_them )
{
  use_lists = use_them;
}

void PlotManager::SetupBox( QVBoxLayout *mainLayout )
{
  int i,j;

  QTabWidget *boxtab = new QTabWidget;
  QWidget *cube = new QWidget( 0 );
  QVBoxLayout *cubeLayout = new QVBoxLayout;

  double axisval[PLOT_TAB_NROW][PLOT_TAB_NCOL] = {{-5.0,-5.0,-5.0},
                                                  { 5.0,-5.0,-5.0},
                                                  {-5.0, 5.0,-5.0},
                                                  {-5.0,-5.0, 5.0}};


  char *cube_names[PLOT_TAB_NROW] = {"Origin", "Axis 1", 
                                     "Axis 2", "Axis 3"};
  QGroupBox* cube_g[PLOT_TAB_NROW];
  QHBoxLayout* cube_l[PLOT_TAB_NROW]; 

  for ( i = 0 ; i < PLOT_TAB_NROW ; i++ )
  {
    cube_g[i] = new QGroupBox(tr(cube_names[i]));
    cube_l[i] = new QHBoxLayout;
    cube_l[i]->addStretch(1);
    for ( j = 0 ; j < PLOT_TAB_NCOL ; j++ )
    {
      axis[i][j] = new QDoubleSpinBox( cube_g[i] );
      axis[i][j]->setRange( -100000.0 , 100000.0 );
      axis[i][j]->setSingleStep( 1.0 );
      axis[i][j]->setValue( axisval[i][j] );
      axis[i][j]->setDecimals( 3 );
      connect( axis[i][j] , SIGNAL(valueChanged( double )) ,
               this , SLOT( ChangeBox( double ) ) );
      cube_l[i]->addWidget( axis[i][j] );
    };
    cube_g[i]->setLayout( cube_l[i] );
    cubeLayout->addWidget( cube_g[i] );
  };
  cube->setLayout( cubeLayout );
  boxtab->addTab( cube , tr("Box") );

  QWidget *flags = new QWidget( 0 );
  QVBoxLayout *flagsLayout = new QVBoxLayout;

  char *flags_names[PLOT_TAB_NROW] = {"General","Visible","Normals","Nothing"};
  QGroupBox* flags_g[PLOT_TAB_NROW];
  QHBoxLayout* flags_l[PLOT_TAB_NROW]; 

  for ( i = 0 ; i < PLOT_TAB_NROW ; i++ )
  {
    flags_g[i] = new QGroupBox(tr( flags_names[i]));
    flags_l[i] = new QHBoxLayout;
    flags_l[i]->addStretch(1);
    for ( j = 0 ; j < PLOT_TAB_NCOL ; j++ )
    {
      flagbox[i][j] = new QCheckBox( flags_g[i] );
      flagbox[i][j]->setText( flagbox_names[i][j] ); 
      connect( flagbox[i][j] , SIGNAL( clicked() ) , 
               this , SLOT( ChangeFlags() ) );
      flags_l[i]->addWidget( flagbox[i][j] );
    };
    flags_g[i]->setLayout( flags_l[i] );
    flagsLayout->addWidget( flags_g[i] );
  };
  flags->setLayout( flagsLayout );
  boxtab->addTab( flags , tr("Flags") );

  QWidget *plane = new QWidget( 0 );
  QVBoxLayout *planeLayout = new QVBoxLayout;

  QGroupBox *axisBox = new QGroupBox(tr("Axis"));
  QHBoxLayout *axisLayout = new QHBoxLayout;
axisLayout->addStretch(1);

   QSpinBox *axisS = new QSpinBox( axisBox );
   axisS->setRange( 1 , 3 );
   axisS->setSingleStep( 1 );
   axisS->setValue( plane_axis+1 );
   connect( axisS , SIGNAL(valueChanged( int )) ,
            this , SLOT( ChangePlaneAxis( int ) ) );


  /*
  QRadioButton *radiok = new QRadioButton(tr("axis k"));
  QRadioButton *radioj = new QRadioButton(tr("axis j"));
  QRadioButton *radioi = new QRadioButton(tr("axis i"));
  radiok->setChecked(true);



  axisLayout->addWidget(radioi);
  axisLayout->addWidget(radioj);
  axisLayout->addWidget(radiok);

  connect(radiok,SIGNAL(clicked()),this, SLOT(radioplanechange()));
  connect(radioj,SIGNAL(clicked()),this, SLOT(radioplanechange()));
  connect(radioi,SIGNAL(clicked()),this, SLOT(radioplanechange()));
*/


  axisLayout->addWidget( axisS );
  axisBox->setLayout( axisLayout );

  QGroupBox *pointBox = new QGroupBox(tr("Point"));
  QHBoxLayout *pointLayout = new QHBoxLayout;
  pointLayout->addStretch(1);

  QSpinBox *pointS = new QSpinBox( pointBox );
  pointS->setRange( 1 , 100000 );
  pointS->setSingleStep( 1 );
  pointS->setValue( plane_point+1 );
  connect( pointS , SIGNAL(valueChanged( int )) ,
           this , SLOT( ChangePlanePoint( int ) ) );

  pointLayout->addWidget( pointS );
  pointBox->setLayout( pointLayout );

  planeLayout->addWidget( axisBox );
  planeLayout->addWidget( pointBox );
  plane->setLayout( planeLayout );
  boxtab->addTab( plane , tr("Plane") );

  mainLayout->addWidget( boxtab );
  hold = 0;
  drawbox = 0;
  drawref = 0;
  approximate_normals = 0;
  average_normals = 0;
  hide_surfaces = 0;
  shiny = 0;

  ChangeBox( 0.0 );
}

void PlotManager::ChangeBox( double )
{
  int i,j;
  Vector va,vb;

  if ((Gaussbool == true) || (Demonbool == true))
  {
  /*
  for ( i = 0 ; i < PLOT_TAB_NROW ; i++ )
  {
    for ( j = 0 ; j < PLOT_TAB_NCOL ; j++ )
    {
      LATTICE[3*i+j] = axis[i][j]->value();
    };
  };
 */
  va = Vector(&LATTICE[0]);
  vb = Vector(&LATTICE[3]);
  vb -= va;
//NLATTICE[1] = (int)(vb.Norm()/boxmesh) + 1;
  if (EUCLID>1)
  {
    vb = Vector(&LATTICE[6]);
    vb -= va;
  //  NLATTICE[2] = (int)(vb.Norm()/boxmesh)+ 1;
  }
 // else NLATTICE[2] = 1;
  if (EUCLID>2)
  {
    vb = Vector(&LATTICE[9]);
    vb -= va;
   // NLATTICE[3] = (int)(vb.Norm()/boxmesh)+ 1;
  }
 // else NLATTICE[3] = 1;
  NLATTICE[0] = NLATTICE[1]*NLATTICE[2]*NLATTICE[3];

  origin = Vector(&LATTICE[0]);
  if (NLATTICE[1]>1)
    dstep[0] = (Vector(&LATTICE[3])-origin)/(NLATTICE[1]-1);
  else
    dstep[0] = (Vector(&LATTICE[3])-origin);

  if (NLATTICE[2]>1)
    dstep[1] = (Vector(&LATTICE[6])-origin)/(NLATTICE[2]-1);
  else
    dstep[1] = (Vector(&LATTICE[6])-origin);

  if (NLATTICE[3]>1)
    dstep[2] = (Vector(&LATTICE[9])-origin)/(NLATTICE[3]-1);
  else
      dstep[2] = (Vector(&LATTICE[9])-origin);
} else {

    for ( i = 0 ; i < PLOT_TAB_NROW ; i++ )
    {
      for ( j = 0 ; j < PLOT_TAB_NCOL ; j++ )
      {
        LATTICE[3*i+j] = axis[i][j]->value();
      };
    };

    va = Vector(&LATTICE[0]);
    vb = Vector(&LATTICE[3]);
    vb -= va;
  NLATTICE[1] = (int)(vb.Norm()/boxmesh) + 1;
    if (EUCLID>1)
    {
      vb = Vector(&LATTICE[6]);
      vb -= va;
      NLATTICE[2] = (int)(vb.Norm()/boxmesh)+ 1;
    }
    else NLATTICE[2] = 1;
    if (EUCLID>2)
    {
      vb = Vector(&LATTICE[9]);
      vb -= va;
      NLATTICE[3] = (int)(vb.Norm()/boxmesh)+ 1;
    }
    else NLATTICE[3] = 1;
    NLATTICE[0] = NLATTICE[1]*NLATTICE[2]*NLATTICE[3];

    origin = Vector(&LATTICE[0]);
    if (NLATTICE[1]>1)
      dstep[0] = (Vector(&LATTICE[3])-origin)/(NLATTICE[1]-1);
    else
      dstep[0] = (Vector(&LATTICE[3])-origin);

    if (NLATTICE[2]>1)
      dstep[1] = (Vector(&LATTICE[6])-origin)/(NLATTICE[2]-1);
    else
      dstep[1] = (Vector(&LATTICE[6])-origin);

    if (NLATTICE[3]>1)
      dstep[2] = (Vector(&LATTICE[9])-origin)/(NLATTICE[3]-1);
    else
        dstep[2] = (Vector(&LATTICE[9])-origin);
}
}

void PlotManager::ChangeFlags()
{
  int i,j;
  for ( i = 0 ; i < PLOT_TAB_NROW ; i++ )
  {
    for ( j = 0 ; j < PLOT_TAB_NCOL ; j++ )
    {
      if ( flagbox[i][j]->text() == "Hold Previous" )
      {
        if ( hold != flagbox[i][j]->checkState() )
        {
          hold = flagbox[i][j]->checkState();
        };
      }
      if ( flagbox[i][j]->text() == "Hide all" )
      {
        if ( hide_surfaces != flagbox[i][j]->checkState() )
        {
          hide_surfaces = flagbox[i][j]->checkState();
          viewer->Redraw();
        };
      }
      if ( flagbox[i][j]->text() == "Reference" )
      {
        if ( drawref != flagbox[i][j]->checkState() )
        {
          drawref = flagbox[i][j]->checkState();
          viewer->Redraw();
        };
      }
      if ( flagbox[i][j]->text() == "Approximate" )
      {
        if ( approximate_normals != flagbox[i][j]->checkState() )
        {
          approximate_normals = flagbox[i][j]->checkState();
          viewer->Redraw();
        };
      }
      if ( flagbox[i][j]->text() == "Average" )
      {
        if ( average_normals != flagbox[i][j]->checkState() )
        {
          average_normals = flagbox[i][j]->checkState();
          viewer->Redraw();
        };
      }
      if ( flagbox[i][j]->text() == "Draw Box" )
      {
        if ( drawbox != flagbox[i][j]->checkState() )
        {
          drawbox = flagbox[i][j]->checkState();
          viewer->Redraw();
        };
      }
      if ( flagbox[i][j]->text() == "Shiny" )
      {
        if ( shiny != flagbox[i][j]->checkState() )
        {
          shiny = flagbox[i][j]->checkState();
          if ( NSurf() > 0 )
          {
            plot_entries[NSurf()-1].shiny = shiny;
            if ( strstr( plot_entries[NSurf()-1].name , "Neg" ) != NULL )
            {
              plot_entries[NSurf()-2].shiny = shiny;
            };
            viewer->Redraw();
          };
        };
      };
    };
  };
}

int get_orbital_number(char* plot_name)
{
  int orbital_number = 0;
  int i,n;
  char str[MAX_STR_SIZE];

  // Bare number
  if (strncmp(plot_name,"AO",2)==0) 
    strcpy( str , &plot_name[2] );
  else if (strncmp(plot_name,"AKMO",4)==0||strncmp(plot_name,"BKMO",4)==0) 
    strcpy( str , &plot_name[5] );
  else
    strcpy( str , &plot_name[3] );
  i = 0;
  while ( str[i] != ' ' ) i++;
  str[i] = '\0';
  n = atoi(str) - 1;

  i = 0;
  if (strncmp(plot_name,"AMO",3)==0) { }
  else if (strncmp(plot_name,"AKMOR",5)==0) { }
  else if (strncmp(plot_name,"AKMOI",5)==0) { }
  else if (strncmp(plot_name,"UMO",3)==0) { }
  else if (strncmp(plot_name,"AO",2)==0) { }
  else if (strncmp(plot_name,"BMO",3)==0)
    while (!((MOS[i].spin[0]=='B')&&(strcmp(MOS[i].spin,"BETA")==0))) i++;
  else if (strncmp(plot_name,"ADY",3)==0)
    while (!((MOS[i].spin[0]=='A')&&(strcmp(MOS[i].spin,"ALPHAD")==0))) i++;
  else if (strncmp(plot_name,"BDY",3)==0)
    while (!((MOS[i].spin[0]=='B')&&(strcmp(MOS[i].spin,"BETAD")==0))) i++;
  orbital_number = n + i;

  return orbital_number;
}

void PlotManager::BuildCube(int n)
{
  int i,j,k;

  orbital_number = get_orbital_number(plot_entries[n].name);

  // Build our cube file
  int pctg;
  char filename[MAX_STR_SIZE];
  sprintf(filename,"cube%d.prk",n);

  ofstream cube(filename);
  cube << "This is a Parakata cube file generated by " <<PROGRAM_NAME<< endl;

  // Name of the field
  cube << "@name{" << endl;
  cube << plot_entries[n].name << endl;
  cube << "}" << endl;

  // Box specifications
  cube << "@box{" << endl;
  for (i=0;i<4;i++)
  {
    for (j=0;j<3;j++)
      cube << fixed <<setw(10) <<setprecision(5)<<right
           << LATTICE[3*i+j]<<" ";
           //rfm << BohrToAngstrom(LATTICE[3*i+j])<<" ";
    cube <<endl;
  }
  for (j=0;j<3;j++)
    cube <<NLATTICE[j+1]<<" ";
  cube <<endl;
  cube << "}" << endl;

  // Molecule
  cube << "@molecule{" << endl;
  for ( k = 0 ; k < NATOM ; k++ )
  {
    cube << setw(2) << left << ELEMENT_SYMBOL[ATOMNO[k]] << " ";
    for ( i = 0 ; i < 3 ; i++ )
    {
      cube << fixed << setw(10) << setprecision(5) << right
           << COORD[1][k][i] << " ";

//           rfm << BohrToAngstrom(COORD[1][k][i]) << " ";
    }
    cube << "\n";
  }
  cube << "}" << endl;

  if ( strstr( plot_entries[n].name , "REACTIVE ORBITAL SPACE" ) != NULL )
  {
    bool ok;
    // Esto esta comentado para lograr compatibilidad con
    // las librerias de la maquina de Gustavo Rodríguez de CULAGOS
    int id = 1;//QInputDialog::getInt(viewer->prk, 
               //                   tr("Input Dialog"),
               //                   tr("Orbital number: "), 
              //                    1, 1, NORB, 1, &ok);
    if (ok)
    {
      id = id - 1;
      build_rod_matrix(AUX_MATRIX,id);
    }
  }
  if ( strstr( plot_entries[n].name , "SPIN" ) != NULL )
  {
    build_spin_matrix(AUX_MATRIX);
  }

  // Field values
    // Field values
  if ( strstr( plot_entries[n].name , "AZD FIELD" ) != NULL )
  {
    EUCLID = 3;
    ChangeBox(0.0);
    ifstream f("azd.data");
    if ( f.fail() )
    {
      errmsg(__FILE__,__LINE__,"Error opening file",0);
      return;
    };
    char line[MAX_STR_SIZE];
    char str[MAX_STR_SIZE];
    while (f.getline(line,MAX_STR_SIZE))
    {
      string stringLine=line;
      int sindex=stringLine.find("AZD FIELD");
      if((sindex>=0)&&(sindex<(signed)stringLine.size()))
      {
        int i,j,k,ni,nj,nk;
        double step;
        f >> step;
        ni = 40+1;
        nj = 20+1;
        nk = 8+1;
        cout << ni << " III " << NLATTICE[1] << endl;
        cout << nj << " JJJ " << NLATTICE[2] << endl;
        cout << nk << " KKK " << NLATTICE[3] << endl;
        cout << NLATTICE[0] << endl;
        if (SCALAR_FIELD) free(SCALAR_FIELD);
        SCALAR_FIELD = (double*) malloc(ni*nj*nk*sizeof(double));
        k = 0;
        for (i=0;i<ni*nj*nk;i++)
        {
          f >> SCALAR_FIELD[k];
          k++;
        }
      }
    }
    f.close();
  }

  cube << "@scalar{" << endl;
  for ( k = 0 ; k < NLATTICE[3] ; k++ )
  {
    for ( j = 0 ; j < NLATTICE[2] ; j++ )
    {
      grdval.clear();
      if (PLOT_TYPE==COLOURED_ISO)
      {
        grdvec.clear();
        for ( i = 0 ; i < NLATTICE[1] ; i++ )
          grdvec.push_back(Vector(origin+dstep[0]*i+dstep[1]*j+dstep[2]*k));
        EvaluateSet( "DENSITY" );
      }
      else if ( strstr( plot_entries[n].name , "AZD FIELD" ) != NULL )
      {
        int l,ll;
        ll = NLATTICE[1]*NLATTICE[2]*k;
        for ( i = 0 ; i < NLATTICE[1] ; i++ )
        {
          l = ll + NLATTICE[1]*j + i;
          grdval.push_back( SCALAR_FIELD[l] );
        }
      }
      else if ( strstr( plot_entries[n].name , "LOADED FIELD" ) != NULL )
      {
        int l,ll;
        ll = NLATTICE[1]*NLATTICE[2]*k;
        for ( i = 0 ; i < NLATTICE[1] ; i++ )
        {
          l = ll + NLATTICE[1]*j + i;
          grdval.push_back( SCALAR_FIELD[l] );
        }
      }
      else 
      {
        grdvec.clear();
        for ( i = 0 ; i < NLATTICE[1] ; i++ )
          grdvec.push_back(Vector(origin+dstep[0]*i+dstep[1]*j+dstep[2]*k));
        EvaluateSet( plot_entries[n].name );
      }
      for ( i = 0 ; i < (signed) grdval.size() ; i++ )
        cube << grdval[i]<<" ";
      cube << endl;
      pctg = int(100*((k*NLATTICE[2]+j)/((double)(NLATTICE[2]*NLATTICE[3]))));
      prk->Report( "Evaluating cube: (%d)",pctg);
    }
  }
  cube << "}" << endl;
  cube.close();
}

void PlotManager::DoPlot(int n)
{
  int i,ip,it;
  int j,jgp,k,nt;
  Vector normal,pos;
  double rgb[3];
  double t[3][5][3],tn[3][5][3];

  orbital_number = get_orbital_number(plot_entries[n].name);

  char filename[MAX_STR_SIZE];
  sprintf(filename,"cube%d.prk",n);

  // Now, load the field to SCALAR_FIELD
  read_Parakata_cube(filename);

  // Active box
  int mini,minj,mink,maxi,maxj,maxk;
  mini = minj = mink = 0;
  maxi = NLATTICE[1];
  if (EUCLID>1) maxj = NLATTICE[2];
  else maxj = 2;
  if (EUCLID>2) maxk = NLATTICE[3];
  else maxk = 2;

  if ( PLOT_TYPE == CONTOUR )
  {
    mink = min(plane_point,NLATTICE[3]-1);
    maxk = mink+2;
  }

  // Plotting
  for ( k = mink+1 ; k < maxk ; k++ )
  {
    for ( j = minj+1 ; j < maxj ; j++ )
    {
      for ( i = mini+1 ; i < maxi ; i++ )
      {
        // Get grid points for a cube, square or line
        grdvec.clear();
        grdval.clear();
        GetCubePoints(i,j,k);
        if (PLOT_TYPE==CONTOUR) 
        {
          int nl;
          double lines[2][2][3];
          LBox(plot_entries[n].iso,lines,&nl);
          for (int il = 0; il < nl ; il++)
          {
            glVertex3dv(&lines[il][0][0]);
            glVertex3dv(&lines[il][1][0]);
          }
          continue;
        }
        if (PLOT_TYPE==CURVE) 
        {
          glColor3d(0.5,0.5,0.5);
          for (ip=0;ip<2;ip++)
          {
            pos = grdvec[ip];
            pos += dstep[1]*MAX(fmin,MIN(fmax,grdval[ip]));
            glVertex3d(pos[0],pos[1],pos[2]);
          }
          continue;
        }
        // Search for triangles
        if (EUCLID==3) TBox(plot_entries[n].iso,t,&nt);
        // Coloured plane triangles
        if (PLOT_TYPE==COLOURED_PLANE) 
        {
          glNormal3d(dstep[2][0],dstep[2][1],dstep[2][2]);
          RealToRGB(grdval[0],rgb);
          glColor3dv(rgb);
          glVertex3d(grdvec[0][0],grdvec[0][1],grdvec[0][2]);
          RealToRGB(grdval[1],rgb);
          glColor3dv(rgb);
          glVertex3d(grdvec[1][0],grdvec[1][1],grdvec[1][2]);
          RealToRGB(grdval[3],rgb);
          glColor3dv(rgb);
          glVertex3d(grdvec[3][0],grdvec[3][1],grdvec[3][2]);
          glVertex3d(grdvec[3][0],grdvec[3][1],grdvec[3][2]);
          RealToRGB(grdval[1],rgb);
          glColor3dv(rgb);
          glVertex3d(grdvec[1][0],grdvec[1][1],grdvec[1][2]);
          RealToRGB(grdval[2],rgb);
          glColor3dv(rgb);
          glVertex3d(grdvec[2][0],grdvec[2][1],grdvec[2][2]);
          continue;
        }
        if (nt>0&&PLOT_TYPE!=CONTOUR) 
        {
          // Got some triangles!
          if (PLOT_TYPE==COLOURED_ISO)
          {
            int ngpn;

            GetNormals("DENSITY",plot_entries[n].iso,t,tn,nt,i,j,k);
            ngpn = 0;
            grdvec.clear();
            grdval.clear();
            for ( it = 0 ; it < nt ; it++ )
            {
              grdvec.push_back( Vector(&t[0][it][0]) );
              grdvec.push_back( Vector(&t[1][it][0]) );
              grdvec.push_back( Vector(&t[2][it][0]) );
              ngpn += 3;
            }
            EvaluateSet(plot_entries[n].name);
          }
          else if (EUCLID==3) GetNormals(plot_entries[n].name,plot_entries[n].iso,t,tn,nt,i,j,k);
          jgp = 0;
          for ( it = 0 ; it < nt ; it++ )
          {
            for ( ip = 0 ; ip < 3 ; ip++ )
            {
              if (PLOT_TYPE==COLOURED_ISO||PLOT_TYPE==COLOURED_PLANE)
              {
                RealToRGB(grdval[jgp+ip],rgb);
                glColor3dv(rgb);
              }
              glNormal3dv(&tn[ip][it][0]);
              glVertex3dv(&t[ip][it][0]);
            }
            jgp = jgp + 3;
          }
        }
      }
    }
    if (EUCLID==2) break;
  }

}

void PlotManager::GetCubePoints(int i,int j,int k)
{
  // Out of range
  if (i<1||j<1||k<1) return;
  if (i>=NLATTICE[1]) return;
  if (j>=NLATTICE[2]&&EUCLID>1) return;
  if (k>=NLATTICE[3]&&EUCLID>2) return;

  if (EUCLID>1) 
  {
    grdvec.push_back(Vector(origin+dstep[0]*(i-1)+dstep[1]*j    +dstep[2]*(k-1)));
    grdvec.push_back(Vector(origin+dstep[0]*i    +dstep[1]*j    +dstep[2]*(k-1)));
    grdvec.push_back(Vector(origin+dstep[0]*i    +dstep[1]*(j-1)+dstep[2]*(k-1)));
    grdvec.push_back(Vector(origin+dstep[0]*(i-1)+dstep[1]*(j-1)+dstep[2]*(k-1)));
    if (EUCLID>2)
    {
      grdvec.push_back( Vector(origin + dstep[0]*(i-1) + dstep[1]*j     + dstep[2]*k));
      grdvec.push_back( Vector(origin + dstep[0]*i     + dstep[1]*j     + dstep[2]*k));
      grdvec.push_back( Vector(origin + dstep[0]*i     + dstep[1]*(j-1) + dstep[2]*k));
      grdvec.push_back( Vector(origin + dstep[0]*(i-1) + dstep[1]*(j-1) + dstep[2]*k));
    }
  }
  else
  {
    grdvec.push_back( Vector(origin + dstep[0]*i     + dstep[1]*(j-1) + dstep[2]*(k-1)));
    grdvec.push_back( Vector(origin + dstep[0]*(i-1) + dstep[1]*(j-1) + dstep[2]*(k-1)));
  }

  // Pick up the values
  int l,ll;
  ll = NLATTICE[1]*NLATTICE[2]*(k-1);
  if (EUCLID>1)
  {
    l = ll + NLATTICE[1]*j + i-1;
    grdval.push_back( SCALAR_FIELD[l] );
    grdval.push_back( SCALAR_FIELD[l+1] );
    l = ll + NLATTICE[1]*(j-1) + i;
    grdval.push_back( SCALAR_FIELD[l] );
    grdval.push_back( SCALAR_FIELD[l-1]);
    if (EUCLID>2)
    {
      ll = NLATTICE[1]*NLATTICE[2]*k;
      l = ll + NLATTICE[1]*j + i-1;
      grdval.push_back( SCALAR_FIELD[l] );
      grdval.push_back( SCALAR_FIELD[l+1] );
      l = ll + NLATTICE[1]*(j-1) + i;
      grdval.push_back( SCALAR_FIELD[l] );
      grdval.push_back( SCALAR_FIELD[l-1] );
    }
  }
  else
  {
    l = ll + NLATTICE[1]*(j-1) + i;
    grdval.push_back( SCALAR_FIELD[l] );
    grdval.push_back( SCALAR_FIELD[l-1] );
  }
}

void PlotManager::EvaluateSet(char* name)
{
  if ( strstr( name , "DENSITY" ) != NULL )
    EvaluateRHO(DENSITY_MATRIX);
  else if ( strstr( name , "MEP" ) != NULL )
    EvaluateMEP(DENSITY_MATRIX);
  else if ( strstr( name , "LAPLACIAN" ) != NULL )
    EvaluateLAPRHO(DENSITY_MATRIX);
  else if ( strstr( name , "KED" ) != NULL )
    EvaluateKED(DENSITY_MATRIX);
  else if ( strcmp( name , "ELF") == 0 )
    EvaluateELF(DENSITY_MATRIX,false);
  else if ( strcmp( name , "DLE") == 0 )
    EvaluateELF(DENSITY_MATRIX,true);
  else if ( strstr( name , "SPIN" ) != NULL )
    EvaluateRHO(AUX_MATRIX);
  else if ( strstr( name , "REACTIVE ORBITAL SPACE" ) != NULL )
    EvaluateRHO(AUX_MATRIX);
  else if ( strstr( name , "DR 1" ) != NULL )
    EvaluateRHO(DMR1);
  else if ( strstr( name , "DR 2" ) != NULL )
    EvaluateRHO(DMR2);
  else if ( strstr( name , "DR 3" ) != NULL )
    EvaluateRHO(DMR3);
  else if ( strstr( name , "FUKUI LEFT" ) != NULL )
    EvaluateRHO(FUKUIL_MATRIX);
  else if ( strstr( name , "FUKUI RIGHT" ) != NULL )
    EvaluateRHO(FUKUIR_MATRIX);
  else if ( strstr( name , "FUKUI AVERAGE" ) != NULL )
  {
    vector<double> homo,lumo;

    EvaluateRHO(FUKUIL_MATRIX);
    homo.clear();
    for (size_t i=0;i<grdval.size();i++)
      homo.push_back( grdval[i] );

    EvaluateRHO(FUKUIR_MATRIX);
    lumo.clear();
    for (size_t i=0;i<grdval.size();i++)
      lumo.push_back( grdval[i] );

    grdval.clear();
    for (size_t i=0;i<homo.size();i++)
      grdval.push_back( 0.5*(homo[i]+lumo[i]) );
  }
  else if ( strstr( name , "FUKUI DUAL" ) != NULL )
  {
    vector<double> homo,lumo;

    EvaluateRHO(FUKUIL_MATRIX);
    homo.clear();
    for (size_t i=0;i<grdval.size();i++)
      homo.push_back( grdval[i] );

    EvaluateRHO(FUKUIR_MATRIX);
    lumo.clear();
    for (size_t i=0;i<grdval.size();i++)
      lumo.push_back( grdval[i] );

    grdval.clear();
    for (size_t i=0;i<homo.size();i++)
      grdval.push_back( (lumo[i]-homo[i]) );
  }
  else if ( strstr( name , "SHAPE ENTROPY" ) != NULL )
    EvaluateShapeEntropy(DENSITY_MATRIX); 
  else if ( strstr( name , "SHAPE ENTROPY ALT" ) != NULL )
  {
    EvaluateRHO(DENSITY_MATRIX);
    for (size_t i=0;i<grdval.size();i++)
    {
      grdval[i] = grdval[i]/NELECTRON;
      if (grdval[i]>0.0) 
        grdval[i] = -grdval[i]*log(grdval[i]);
    }
  }
  else if ( strstr( name , "IONIZATION RESPONSE LEFT" ) != NULL )
  {
    vector<double> shape;
    vector<double> fukui;

    EvaluateRHO(DENSITY_MATRIX);
    shape.clear();
    for (size_t i=0;i<grdval.size();i++)
      shape.push_back( grdval[i]/(double)NELECTRON );

    EvaluateRHO(FUKUIL_MATRIX);
    fukui.clear();
    for (size_t i=0;i<grdval.size();i++)
      fukui.push_back( grdval[i] );

    grdval.clear();
    for (size_t i=0;i<shape.size();i++)
	{
	 if (shape[i]>0.0)
         grdval.push_back( -(fukui[i]-shape[i])*(log(shape[i])+1.0)/(double)NELECTRON); else grdval.push_back(0.0);
	} 
  }
  else if ( strstr( name , "IONIZATION RESPONSE RIGHT" ) != NULL )
  {
    vector<double> shape;
    vector<double> fukui;

    EvaluateRHO(DENSITY_MATRIX);
    shape.clear();
    for (size_t i=0;i<grdval.size();i++)
      shape.push_back( grdval[i]/(double)NELECTRON );

    EvaluateRHO(FUKUIR_MATRIX);
    fukui.clear();
    for (size_t i=0;i<grdval.size();i++)
      fukui.push_back( grdval[i] );

    grdval.clear();
    for (size_t i=0;i<shape.size();i++)
       if (shape[i]>0.0)
        grdval.push_back( -(fukui[i]-shape[i])*(log(shape[i])+1.0)/(double)NELECTRON); else grdval.push_back(0.0);
  }
  else if ( strstr( name , "KMOR" ) != NULL )
    EvaluateKMO(KMOSR[orbital_number],KMOSI[orbital_number],true);
  else if ( strstr( name , "KMOI" ) != NULL )
    EvaluateKMO(KMOSR[orbital_number],KMOSI[orbital_number],false);
  else if ( strstr( name , "MO" ) != NULL || strstr( name , "DY" ) != NULL )
    EvaluateMO(MOS[orbital_number]);
  else if ( strstr( name , "AO" ) != NULL )
    EvaluateAOBas(orbital_number);
}

int PlotManager::EvaluateSetD(char* name,double dval[15][3])
{
  int ret = 1;
  if ( strstr( name , "DENSITY" ) != NULL )
    EvaluateDRHO(dval,DENSITY_MATRIX);
  else if ( strstr( name , "LAPLACIAN" ) != NULL )
    ret = 0;
  else if ( strstr( name , "KED" ) != NULL )
    ret = 0;
  else if ( strcmp( name , "ELF") == 0 )
    ret = 0;
  else if ( strcmp( name , "DLE") == 0 )
    ret = 0;
  else if ( strstr( name , "SPIN" ) != NULL )
    EvaluateDRHO(dval,AUX_MATRIX);
  else if ( strstr( name , "REACTIVE ORBITAL SPACE" ) != NULL )
    EvaluateDRHO(dval,AUX_MATRIX);
  else if ( strstr( name , "DR 1" ) != NULL )
    EvaluateDRHO(dval,DMR1);
  else if ( strstr( name , "DR 2" ) != NULL )
    EvaluateDRHO(dval,DMR2);
  else if ( strstr( name , "DR 3" ) != NULL )
    EvaluateDRHO(dval,DMR3);
  else if ( strstr( name , "FUKUI LEFT" ) != NULL )
    EvaluateDRHO(dval,FUKUIL_MATRIX);
  else if ( strstr( name , "FUKUI RIGHT" ) != NULL )
    EvaluateDRHO(dval,FUKUIR_MATRIX);
  else if ( strstr( name , "FUKUI AVERAGE" ) != NULL )
  {
    double dhomo[15][3];
    double dlumo[15][3];

    EvaluateDRHO(dhomo,FUKUIL_MATRIX);
    EvaluateDRHO(dlumo,FUKUIR_MATRIX);
    for (size_t i=0;i<grdvec.size();i++)
      for (size_t j=0;j<3;j++)
        dval[i][j] = 0.5*(dlumo[i][j]+dhomo[i][j]);
  }
  else if ( strstr( name , "FUKUI DUAL" ) != NULL )
  {
    double dhomo[15][3];
    double dlumo[15][3];

    EvaluateDRHO(dhomo,FUKUIL_MATRIX);
    EvaluateDRHO(dlumo,FUKUIR_MATRIX);
    for (size_t i=0;i<grdvec.size();i++)
      for (size_t j=0;j<3;j++)
        dval[i][j] = (dlumo[i][j]-dhomo[i][j]);
  }
//  else if ( strstr( name , "KMOR" ) != NULL )
//    EvaluateDKMO(dval,KMOSR[orbital_number],KMOSI[orbital_number],true);
//  else if ( strstr( name , "KMOI" ) != NULL )
//    EvaluateDKMO(dval,KMOSR[orbital_number],KMOSI[orbital_number],false);
  else if ( strstr( name , "MO" ) != NULL || strstr( name , "DY" ) != NULL )
    EvaluateDMO(dval,MOS[orbital_number]);
  else if ( strstr( name , "AO" ) != NULL )
    EvaluateDAOBas(dval,orbital_number);

  return ret;
}
// Purpose: Get isosurface triangles from a cubic box.
//
// (c) This an adaption of Paul Bourke's code.
void PlotManager::TBox(double iso,
                       double t[3][5][3],
                       int* nt)
{
/*
C     Order of the cube vertices:
C
C          4----------5    
C         /|         /|   
C        / |        / |
C       7----------6  |  
C       |  |       |  | 
C       |  |       |  |   
C       |  0-------|--1  
C       | /        | /  
C       |/         |/  
C       3----------2
C
*/
  int i,j,k,indx; 
  double a,b;
  double vlist[12][3];

  int edgetable[256]={
  0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
  0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
  0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
  0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
  0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
  0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
  0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
  0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
  0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
  0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
  0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
  0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
  0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
  0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
  0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
  0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
  0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
  0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
  0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
  0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
  0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
  0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
  0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
  0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
  0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
  0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
  0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
  0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
  0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
  0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
  0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
  0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

  int tritable[256][16] =
  {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
  {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
  {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
  {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
  {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
  {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
  {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
  {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
  {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
  {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
  {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
  {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
  {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
  {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
  {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
  {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
  {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
  {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
  {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
  {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
  {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
  {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
  {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
  {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
  {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
  {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
  {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
  {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
  {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
  {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
  {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
  {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
  {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
  {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
  {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
  {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
  {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
  {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
  {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
  {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
  {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
  {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
  {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
  {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
  {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
  {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
  {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
  {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
  {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
  {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
  {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
  {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
  {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
  {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
  {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
  {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
  {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
  {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
  {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
  {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
  {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
  {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
  {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
  {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
  {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
  {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
  {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
  {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
  {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
  {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
  {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
  {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
  {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
  {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
  {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
  {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
  {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
  {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
  {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
  {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
  {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
  {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
  {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
  {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
  {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
  {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
  {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
  {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
  {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
  {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
  {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
  {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
  {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
  {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
  {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
  {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
  {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

  int ptr[2][12] = {{0,1,2,3,4,5,6,7,0,1,2,3},
                    {1,2,3,0,5,6,7,4,4,5,6,7}};

  *nt = 0;

  indx = 0;
  j = 1;
  for ( i = 0 ; i < 8 ; i++ )
  {
    if (grdval[i]<iso) indx = indx + j;
    j = j + j;
  }

  if (indx==0||indx==255) return;

  j = 1;
  for ( i = 0 ; i < 12 ; i++ )
  {
    if (edgetable[indx] & j)
    {
      a = grdval[ptr[0][i]];
      b = grdval[ptr[1][i]];
      // Get our approximate intersection ISO values
      if (ABS(iso-a)<0.00001) 
      {
        vlist[i][0] = grdvec[ptr[0][i]][0];
        vlist[i][1] = grdvec[ptr[0][i]][1];
        vlist[i][2] = grdvec[ptr[0][i]][2];
      }
      else if(ABS(iso-b)<0.00001) 
      {
        vlist[i][0] = grdvec[ptr[1][i]][0];
        vlist[i][1] = grdvec[ptr[1][i]][1];
        vlist[i][2] = grdvec[ptr[1][i]][2];
      }
      else if (ABS(b-a)<0.00001)
      {
        vlist[i][0] = grdvec[ptr[0][i]][0];
        vlist[i][1] = grdvec[ptr[0][i]][1];
        vlist[i][2] = grdvec[ptr[0][i]][2];
      }
      else
      {
        vlist[i][0] = grdvec[ptr[0][i]][0] + (iso-a)/(b-a)*
                     (grdvec[ptr[1][i]][0]-grdvec[ptr[0][i]][0]);
        vlist[i][1] = grdvec[ptr[0][i]][1] + (iso-a)/(b-a)*
                     (grdvec[ptr[1][i]][1]-grdvec[ptr[0][i]][1]);
        vlist[i][2] = grdvec[ptr[0][i]][2] + (iso-a)/(b-a)*
                     (grdvec[ptr[1][i]][2]-grdvec[ptr[0][i]][2]);
      }
    }
    j = j + j;
  }

  i = 0; 
  while ( tritable[indx][i] != -1)
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      for ( k = 0 ; k < 3 ; k++ )
      {
        t[j][*nt][k] = vlist[tritable[indx][i+j]][k];
      }
    }
    (*nt)++;
    i += 3; 
  }
}

static void estimate_normals(double t[3][5][3], double tn[3][5][3], int nt)
{
  int i,it,iv;
  double s1[3],s2[3],normal[3];

  for (it = 0; it < nt ; it++ )
  {
    for ( i = 0 ; i < 3 ; i++ )
    {
      s1[i] = t[1][it][i] - t[0][it][i];
      s2[i] = t[2][it][i] - t[0][it][i];
    }
    normal[0] = s1[1]*s2[2]-s1[2]*s2[1];
    normal[1] = s1[2]*s2[0]-s1[0]*s2[2];
    normal[2] = s1[0]*s2[1]-s1[1]*s2[0];
    for (iv=0;iv<3;iv++)
      for ( i = 0 ; i < 3 ; i++ )
        tn[iv][it][i] = normal[i];
  }

  double pad = 1000000.0;
  for (it = 0; it < nt ; it++ )
    for (iv=0;iv<3;iv++)
      for ( i = 0 ; i < 3 ; i++ )
        tn[iv][it][i] = double(int(pad*tn[iv][it][i])/pad);
}

void PlotManager::GetNormals(char* name,
                             double iso,
                             double rt[3][5][3],
                             double rtn[3][5][3],
                             int rnt,
                             int I,
                             int J,
                             int K)
{
  int i,ii,it,iv,rit,riv,jj,kk,nt;
  double norma;
  double t[3][5][3];
  double tn[3][5][3];
  double v[3],rv[3];

  double dval[15][3];
  if (approximate_normals) goto XX;
    grdvec.clear();
    for (rit=0;rit<rnt;rit++)
      for (riv=0;riv<3;riv++)
      {
        for (i=0;i<3;i++) 
          v[i] = rt[riv][rit][i];
        grdvec.push_back( Vector(v) );
      }
    EvaluateSetD(name,dval);
    for (rit=0;rit<rnt;rit++)
      for (riv=0;riv<3;riv++)
        for (i=0;i<3;i++) 
          rtn[riv][rit][i] = dval[rit*3+riv][i];
XX:

if (approximate_normals&&average_normals)
{
  // Initialize
  for (rit=0;rit<rnt;rit++)
    for (riv=0;riv<3;riv++)
      for (i=0;i<3;i++)
        rtn[riv][rit][i] = 0.0;

  // Average
  for (ii=-1;ii<2;ii++)
  {
    for (jj=-1;jj<2;jj++)
    {
      for (kk=-1;kk<2;kk++)
      {
        grdvec.clear();
        grdval.clear();
        GetCubePoints(I+ii,J+jj,K+kk);
        TBox(iso,t,&nt);
        if (nt>0) 
        {
          estimate_normals(t,tn,nt);
          for (it=0;it<nt;it++)
            for (iv=0;iv<3;iv++)
            {
              for (i=0;i<3;i++) v[i] = t[iv][it][i];
              for (rit=0;rit<rnt;rit++)
                for (riv=0;riv<3;riv++)
                {
                  for (i=0;i<3;i++) rv[i] = v[i] - rt[riv][rit][i]; 
                  norma = rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2];
                  if (norma<0.0001)
                    for (i=0;i<3;i++) rtn[riv][rit][i] += tn[iv][it][i];
                }
            }
        }
      }
    }
  }
}
else if (approximate_normals)
{
  estimate_normals(rt,rtn,rnt);
}

  // Orientation 
  if (iso>=0.0)
    for (rit=0;rit<rnt;rit++)
      for (riv=0;riv<3;riv++)
        for (i=0;i<3;i++)
          rtn[riv][rit][i] *= -1.0;

}

// Get contour lines in a square box.
void PlotManager::LBox(double iso,double l[2][2][3],int* nl)
{
/*
C  
C
C             0 ----------- 1
C            /             /
C           /             /
C          /             /
C         3 ----------- 2
C
*/
  int i,il,ip,j,indx;
  int couple[2][2][2];
  double EPS=0.00001;

  indx = 0;
  j = 1;
  for (i=0;i<4;i++)
  {
    if (grdval[i]<iso) indx += j;
    j *= 2;
  }
  if (indx>7) indx = 15 - indx;

  *nl = 1;
  if (indx==1) 
  {
    couple[0][0][0] = 0;
    couple[0][0][1] = 1;
    couple[0][1][0] = 0;
    couple[0][1][1] = 3;
  }
  else if (indx==2)
  {
    couple[0][0][0] = 0;
    couple[0][0][1] = 1;
    couple[0][1][0] = 1;
    couple[0][1][1] = 2;
  }
  else if (indx==3)
  {
    couple[0][0][0] = 0;
    couple[0][0][1] = 3;
    couple[0][1][0] = 1;
    couple[0][1][1] = 2;
  }
  else if (indx==4) 
  {
    couple[0][0][0] = 1;
    couple[0][0][1] = 2;
    couple[0][1][0] = 2;
    couple[0][1][1] = 3;
  }
  else if (indx==5)
  {
    *nl = 2;
    couple[0][0][0] = 0;
    couple[0][0][1] = 1;
    couple[0][1][0] = 1;
    couple[0][1][1] = 2;
    couple[1][0][0] = 2;
    couple[1][0][1] = 3;
    couple[1][1][0] = 3;
    couple[1][1][1] = 0;
  }
  else if (indx==6)
  {
    couple[0][0][0] = 0;
    couple[0][0][1] = 1;
    couple[0][1][0] = 2;
    couple[0][1][1] = 3;
  }
  else if (indx==7)
  {
    couple[0][0][0] = 2;
    couple[0][0][1] = 3;
    couple[0][1][0] = 3;
    couple[0][1][1] = 0;
  }
  else if (indx==8)
  {
    *nl = 2;
    couple[0][0][0] = 0;
    couple[0][0][1] = 1;
    couple[0][1][0] = 0;
    couple[0][1][1] = 3;
    couple[1][0][0] = 1;
    couple[1][0][1] = 2;
    couple[1][1][0] = 2;
    couple[1][1][1] = 3;
  }
  else *nl = 0;
        
  double factor;
  for ( il = 0 ; il < *nl ; il++ )
  {
    for ( ip = 0 ; ip < 2 ; ip++ )
    {
      i = couple[il][ip][0];
      j = couple[il][ip][1];
      if (ABS(iso-grdval[i])<EPS)
      {
        l[il][ip][0] = grdvec[i][0];
        l[il][ip][1] = grdvec[i][1];
        l[il][ip][2] = grdvec[i][2];
      }
      else if (ABS(iso-grdval[j])<EPS)
      {
        l[il][ip][0] = grdvec[j][0];
        l[il][ip][1] = grdvec[j][1];
        l[il][ip][2] = grdvec[j][2];
      }
      else if (ABS(grdval[i]-grdval[j])<EPS)
      {
        l[il][ip][0] = grdvec[i][0];
        l[il][ip][1] = grdvec[i][1];
        l[il][ip][2] = grdvec[i][2];
      }
      else 
      {
        factor = (iso-grdval[i])/(grdval[j]-grdval[i]);
        l[il][ip][0] = grdvec[i][0] + factor*(grdvec[j][0]-grdvec[i][0]);
        l[il][ip][1] = grdvec[i][1] + factor*(grdvec[j][1]-grdvec[i][1]);
        l[il][ip][2] = grdvec[i][2] + factor*(grdvec[j][2]-grdvec[i][2]);
      }
    }
  }
}

// Basis functions
void PlotManager::EvaluateAO(Vector r, double ao[MAXNBAS])
{
  int i,iatom,ibas,igto,ishl,jbas,lmax,lx,ly,lz,skipatom,skipshl;
  double factor,expf,r2,term;
  double *px,*py,*pz;
  Vector rpos;

  ibas = 0;
  for (iatom=0;iatom<NATOM;iatom++)
  {
    rpos = Vector(r);
    rpos -= Vector(&COORD[1][iatom][0]);
    r2 = pow(rpos.Norm(),2.0);
    AtomBasis ab = mb.ab[iatom];
    if (r2>ab.maxrad) skipatom = 1;
    else skipatom = 0;
    if (!skipatom)
    {
      lmax = 0;
      for (ishl=0;ishl<ab.nshl;ishl++)
        lmax = MAX(lmax,ab.sb[ishl].l);
      // Allocate memory for local vectors
      int ml = lmax + 1;
      px = (double*)malloc(ml*sizeof(double));
      py = (double*)malloc(ml*sizeof(double));
      pz = (double*)malloc(ml*sizeof(double));
      px[0] = 1.0;
      py[0] = 1.0;
      pz[0] = 1.0;
      for ( i = 1 ; i <= lmax ; i++ )
      {
        px[i] = px[i-1]*rpos[0];
        py[i] = py[i-1]*rpos[1];
        pz[i] = pz[i-1]*rpos[2];
      }
    }
    for (ishl=0;ishl<ab.nshl;ishl++)
    {
      ShellBasis sb = ab.sb[ishl];
      if (r2>sb.maxrad) skipshl = 1;
      else skipshl = 0;
      if (!skipshl)
      {
        expf = 0.0;
        for ( igto = 0 ; igto < sb.k ; igto++ )
        {
          term = sb.d[igto]*exp(-sb.z[igto]*r2);
          expf += term;
        }
        jbas = 0;
      }
      for ( lx = sb.l ; lx >= 0 ; lx-- )
      {
        for ( ly = sb.l - lx ; ly >= 0 ; ly-- )
        {
          lz = sb.l - lx - ly;
          if (!skipshl)
          {
            factor = sb.ncsto[jbas];
            ao[ibas] = factor*px[lx]*py[ly]*pz[lz]*expf;
            jbas++;
          }
          else ao[ibas] = 0.0;
          ibas++;
        }
      }
    }
    if (!skipatom)
    {
      // Free dynamic memory
      free(px);
      free(py);
      free(pz);
    }
  }
}

// Basis function gradient
void PlotManager::EvaluateAOGRAD(Vector r, double aograd[MAXNBAS][3])
{
  int i,iatom,ibas,igto,ishl,jbas,lmax,lx,ly,lz,skipatom,skipshl;
  double dexpf,factor,expf,r2,term;
  double *px,*py,*pz;
  Vector rpos;

  ibas = 0;
  for (iatom=0;iatom<NATOM;iatom++)
  {
    rpos = Vector(r);
    rpos -= Vector(&COORD[1][iatom][0]);
    r2 = pow(rpos.Norm(),2.0);
    AtomBasis ab = mb.ab[iatom];
    if (r2>ab.maxrad) skipatom = 1;
    else skipatom = 0;
    if (!skipatom) 
    {
      lmax = 0;
      for (ishl=0;ishl<ab.nshl;ishl++)
        lmax = MAX(lmax,ab.sb[ishl].l);
      // Allocate memory for local vectors
      int ml = lmax + 2;
      px = (double*)malloc(ml*sizeof(double));
      py = (double*)malloc(ml*sizeof(double));
      pz = (double*)malloc(ml*sizeof(double));
      px[0] = 1.0;
      py[0] = 1.0;
      pz[0] = 1.0;
      for ( i = 1 ; i <= lmax+1 ; i++ )
      {
        px[i] = px[i-1]*rpos[0];
        py[i] = py[i-1]*rpos[1];
        pz[i] = pz[i-1]*rpos[2];
      }
    }
    for (ishl=0;ishl<ab.nshl;ishl++)
    {
      ShellBasis sb = ab.sb[ishl];
      if (r2>sb.maxrad) skipshl = 1;
      else skipshl = 0;
      if (!skipshl)
      {
        expf = 0.0;
        dexpf = 0.0;
        for ( igto = 0 ; igto < sb.k ; igto++ )
        {
          term = sb.d[igto]*exp(-sb.z[igto]*r2);
          expf += term;
          term = 2.0*sb.z[igto]*term;
          dexpf -= term;
        }
        jbas = 0;
      }
      for ( lx = sb.l ; lx >= 0 ; lx-- )
      {
        for ( ly = sb.l - lx ; ly >= 0 ; ly-- )
        {
          lz = sb.l - lx - ly;
          if (!skipshl)
          {
            factor = sb.ncsto[jbas];
            aograd[ibas][0] = factor*px[lx+1]*py[ly]*pz[lz]*dexpf;
            if (lx>0) aograd[ibas][0] += factor*px[lx-1]*py[ly]*pz[lz]*lx*expf;
            aograd[ibas][1] = factor*px[lx]*py[ly+1]*pz[lz]*dexpf;
            if (ly>0) aograd[ibas][1] += factor*px[lx]*py[ly-1]*pz[lz]*ly*expf;
            aograd[ibas][2] = factor*px[lx]*py[ly]*pz[lz+1]*dexpf;
            if (lz>0) aograd[ibas][2] += factor*px[lx]*py[ly]*pz[lz-1]*lz*expf;
            jbas++;
          }
          else
          {
            aograd[ibas][0] = 0.0;
            aograd[ibas][1] = 0.0;
            aograd[ibas][2] = 0.0;
          }
          ibas++;
        }
      }
    }
    if (!skipatom)
    {
      // Free dynamic memory
      free(px);
      free(py);
      free(pz);
    }
  }
}

// Laplacian of basis functions
void PlotManager::EvaluateAOLAP(Vector r, double aolap[MAXNBAS])
{
  int i,iatom,ibas,igto,ishl,jbas,lmax,lx,ly,lz,skipatom,skipshl;
  double factor,expf,dexpf,d2expf,r2,term;
  double *px,*py,*pz;
  Vector rpos;

  ibas = 0;
  for (iatom=0;iatom<NATOM;iatom++)
  {
    rpos = Vector(r);
    rpos -= Vector(&COORD[1][iatom][0]);
    r2 = pow(rpos.Norm(),2.0);
    AtomBasis ab = mb.ab[iatom];
    if (r2>ab.maxrad) skipatom = 1;
    else skipatom = 0;
    if (!skipatom)
    {
      lmax = 0;
      for (ishl=0;ishl<ab.nshl;ishl++)
        lmax = MAX(lmax,ab.sb[ishl].l);
      // Allocate memory for local vectors
      int ml = lmax + 3;
      px = (double*)malloc(ml*sizeof(double));
      py = (double*)malloc(ml*sizeof(double));
      pz = (double*)malloc(ml*sizeof(double));
      px[0] = 1.0;
      py[0] = 1.0;
      pz[0] = 1.0;
      for ( i = 1 ; i <= lmax+2 ; i++ )
      {
        px[i] = px[i-1]*rpos[0];
        py[i] = py[i-1]*rpos[1];
        pz[i] = pz[i-1]*rpos[2];
      }
    }
    for (ishl=0;ishl<ab.nshl;ishl++)
    {
      ShellBasis sb = ab.sb[ishl];
      if (r2>sb.maxrad) skipshl = 1;
      else skipshl = 0;
      if (!skipshl)
      {
        expf = 0.0;
        dexpf = 0.0;
        d2expf = 0.0;
        for ( igto = 0 ; igto < sb.k ; igto++ )
        {
          term = sb.d[igto]*exp(-sb.z[igto]*r2);
          expf += term;
          term = 2.0*sb.z[igto]*term;
          dexpf -= term;
          term = 2.0*sb.z[igto]*term;
          d2expf += term;
        }
        jbas = 0;
      }
      for ( lx = sb.l ; lx >= 0 ; lx-- )
      {
        for ( ly = sb.l - lx ; ly >= 0 ; ly-- )
        {
          lz = sb.l - lx - ly;
          if (!skipshl)
          {
            factor = sb.ncsto[jbas]*d2expf;
            aolap[ibas] = factor*px[lx+2]*py[ly]*pz[lz];
            aolap[ibas] += factor*px[lx]*py[ly+2]*pz[lz];
            aolap[ibas] += factor*px[lx]*py[ly]*pz[lz+2];
            factor = sb.ncsto[jbas]*dexpf*px[lx]*py[ly]*pz[lz];
            aolap[ibas] += factor*(2*lx+1);
            aolap[ibas] += factor*(2*ly+1);
            aolap[ibas] += factor*(2*lz+1);
            factor = sb.ncsto[jbas]*expf;
            if (lx>1) aolap[ibas] += factor*px[lx-2]*py[ly]*pz[lz]*lx*(lx-1);
            if (ly>1) aolap[ibas] += factor*px[lx]*py[ly-2]*pz[lz]*ly*(ly-1);
            if (lz>1) aolap[ibas] += factor*px[lx]*py[ly]*pz[lz-2]*lz*(lz-1);
            jbas++;
          }
          else
          {
            aolap[ibas] = 0.0;
          }
          ibas++;
        }
      }
    }
    if (!skipatom)
    {
      // Free dynamic memory
      free(px);
      free(py);
      free(pz);
    }
  }
}

// TESTING BEGINS
void PlotManager::EvaluateDegenerateMOSRight(int N,int ONE)
{
grdvaldegen.clear();
double n = N;
for (int i = 0; i < N; i++)
	{
	EvaluateMO(MOS[ONE+i]);
	for (int q=0; q < grdval.size(); q++)
		{
		if(i == 0) grdvaldegen.push_back(pow(grdval[q],2)); else if (i == N-1) {
									       grdvaldegen[q] = (grdvaldegen[q] + pow(grdval[q],2))/n;
										}

							     else {
								grdvaldegen[q] += pow(grdval[q],2);
								  }
		if (grdvaldegen[q] == 0) grdval[q] = 0.0; else	grdval[q] = (grdvaldegen[q] * log(grdvaldegen[q]));
		}
	}
}

void PlotManager::EvaluateDegenerateMOSLeft(int N,int ONE)
{
grdvaldegen.clear();
double n = N;
for (int i = 0; i < N; i++)
        {
        EvaluateMO(MOS[ONE-i]);
        for (int q=0; q < grdval.size(); q++)
                {
                if(i == 0) grdvaldegen.push_back(pow(grdval[q],2)); else if (i == N-1) {
                                                                               grdvaldegen[q] = (grdvaldegen[q] + pow(grdval[q],2))/n;
                                                                                }

                                                             else {
                                                                grdvaldegen[q] += pow(grdval[q],2);
                                                                  }
                        grdval[q] = (grdvaldegen[q] * log(grdvaldegen[q]));
                }
        }
}

void PlotManager::EvaluateDegenerateDUAL(int N,int M,int ONE, int TWO)
{
grdval.clear();
grdvaldegen.clear();
double n = N;
double m = M;
for (int i = 0; i < M; i++)
        {
        EvaluateMO(MOS[ONE-i]);
        for (int q=0; q < grdval.size(); q++)
                {
                if(i == 0) grdvaldegen.push_back((pow(grdval[q],2))*-1); else if (i == M-1) {
                                                                               grdvaldegen[q] = (grdvaldegen[q] - pow(grdval[q],2))/m;
                                                                                }

                                                             else {
                                                                grdvaldegen[q] -= pow(grdval[q],2);
                                                                  }
                }
        }

for (int i = 0; i < N; i++)
        {
        EvaluateMO(MOS[TWO+i]);
        for (int q=0; q < grdval.size(); q++)
                {
                 if (i == N-1) {
                               grdvaldegen[q] = (grdvaldegen[q] + pow(grdval[q],2))/n;
                                     }

                                                             else {
                                                                grdvaldegen[q] += pow(grdval[q],2);
                                                                  }
                }
        }

for (int i = 0; i < grdvaldegen.size(); i++) grdval[i] = grdvaldegen[i];


}

// TESTING ENDS

void PlotManager::EvaluateAOBas(int id)
{
  int igp;
  Vector pos;
  double value;
  double ao[MAXNBAS];

  grdval.clear();
  for ( igp = 0 ; igp < (signed) grdvec.size() ; igp++ )
  {
    pos = Vector(grdvec[igp]);
    EvaluateAO(pos,ao);
    value = ao[id];
    grdval.push_back( value );
  }
}

void PlotManager::EvaluateMO(MO mo)
{
  int ibas,igp;
  Vector pos;
  double value;
  double ao[MAXNBAS];

  grdval.clear();
  for ( igp = 0 ; igp < (signed) grdvec.size() ; igp++ )
  {
    value = 0.0;
    pos = Vector(grdvec[igp]);
    EvaluateAO(pos,ao);
    for (ibas=0;ibas<NBASIS;ibas++)
    {
      value += mo.c[ibas]*ao[ibas];
    }

    grdval.push_back( value );
  }
}

void PlotManager::EvaluateKMO(MO kmor,MO kmoi,bool real)
{
  int ibas,igp,l,m,n;
  double coef,ckds,kds,skds,value;
  double ao[MAXNBAS];
  Vector kvec,pos,spos,shift,sl,sm,sn;

  kvec = Vector( &PBCK[0] );

  grdval.clear();
  for ( igp = 0 ; igp < (signed) grdvec.size() ; igp++ )
  {
    pos = Vector(grdvec[igp]);
    value = 0.0;
    for (l=-PBCN1;l<=PBCN1;l++)
    {
      sl = Vector( &PBCCD1[0] ); 
      sl *= double(l);
      for (m=-PBCN2;m<=PBCN2;m++)
      {
        sm = Vector( &PBCCD2[0] ); 
        sm *= double(m);
        for (n=-PBCN3;n<=PBCN3;n++)
        {
          sn = Vector( &PBCCD3[0] ); 
          sn *= double(n);
          shift = sl;
          shift += sm;
          shift += sn;
	  kds = shift.Dot( kvec );
	  ckds = cos(kds);
	  skds = sin(kds);
          spos = pos;
          spos -= shift;
          EvaluateAO(spos,ao);
          for (ibas=0;ibas<NBASIS;ibas++)
          {
	    if (real) coef = kmor.c[ibas]*ckds - kmoi.c[ibas]*skds;  // Real
	    else coef = kmoi.c[ibas]*ckds + kmor.c[ibas]*skds; // Imaginary
            value += coef*ao[ibas];
          }
	}
      }
    }
    grdval.push_back( value );
  }
}

void PlotManager::EvaluateRHO(double* P)
{
  int ibas,igp,jbas,k;
  double value;
  double ao[MAXNBAS];
  Vector pos;

  if (!P) return;

  grdval.clear();
  for ( igp = 0 ; igp < (signed)grdvec.size() ; igp++ )
  {
    value = 0.0;
    pos = Vector(grdvec[igp]);
    EvaluateAO(pos,ao);
    k = 0;
    for (ibas=0;ibas<NBASIS;ibas++)
      if (ao[ibas]!=0.0)
      {
        for (jbas=ibas;jbas<NBASIS;jbas++)
        {
          if (ao[jbas]!=0.0)
            value += P[k]*ao[ibas]*ao[jbas];
          k++;
        }
      }
      else k += (NBASIS-ibas);
    grdval.push_back( value );
  }
}

void PlotManager::EvaluateShapeEntropy(double* P)
{
  int ibas,igp,jbas,k;
  double value;
  double ao[MAXNBAS];
  Vector pos;

  if (!P) return;

  grdval.clear();
  for ( igp = 0 ; igp < (signed)grdvec.size() ; igp++ )
  {
    value = 0.0;
    double NE;
    NE = NELECTRON;
    pos = Vector(grdvec[igp]);
    EvaluateAO(pos,ao);
    k = 0;
    for (ibas=0;ibas<NBASIS;ibas++)
      if (ao[ibas]!=0.0)
      {
        for (jbas=ibas;jbas<NBASIS;jbas++)
        {
          if (ao[jbas]!=0.0)
            value += P[k]*ao[ibas]*ao[jbas];
          k++;
        }
      }
      else k += (NBASIS-ibas);
if (value <= 0.0) value = 0; else {
      value = value/NE;
      value = ((value) * log(value))*-1;
                  }
    grdval.push_back( value );
  }
}

void PlotManager::EvaluatedSsigma(double* F,double* R)
{
  int ibas,igp,jbas,k;
  double value,fvalue,rvalue,NE;
  double ao[MAXNBAS];
  Vector pos;
  NE = NELECTRON;
  if (!F) return;
  if (!R) return;
  grdval.clear();
  for ( igp = 0 ; igp < (signed)grdvec.size() ; igp++ )
  {
    value = 0.0;
    pos = Vector(grdvec[igp]);
    EvaluateAO(pos,ao);
    k = 0;
    for (ibas=0;ibas<NBASIS;ibas++)
      if (ao[ibas]!=0.0)
      {
        for (jbas=ibas;jbas<NBASIS;jbas++)
        {
          if (ao[jbas]!=0.0)
            fvalue += F[k]*ao[ibas]*ao[jbas];
	    rvalue += R[k]*ao[ibas]*ao[jbas];
          k++;
        }
      }
      else k += (NBASIS-ibas);
    value = (1/NE)*(fvalue-(rvalue/NE))*(log(rvalue/NE)+1);
    grdval.push_back( value );
  }
}

void PlotManager::EvaluateLAPRHO(double* P)
{
  int ibas,igp,jbas,k;
  double value;
  double ao[MAXNBAS];
  double aograd[MAXNBAS][3];
  double aolap[MAXNBAS];
  Vector pos;

  if (!P) return;

  grdval.clear();
  for ( igp = 0 ; igp < (signed)grdvec.size() ; igp++ )
  {
    value = 0.0;
    pos = Vector(grdvec[igp]);
    EvaluateAO(pos,ao);
    EvaluateAOGRAD(pos,aograd);
    EvaluateAOLAP(pos,aolap);
    k = 0;
    for (ibas=0;ibas<NBASIS;ibas++)
      if (ao[ibas]!=0.0)
      {
        for (jbas=ibas;jbas<NBASIS;jbas++)
        {
          if (ao[jbas]!=0.0)
            value += P[k]*(ao[ibas]*aolap[jbas]+ao[jbas]*aolap[ibas]+
            2.0*aograd[ibas][0]*aograd[jbas][0]+
            2.0*aograd[ibas][1]*aograd[jbas][1]+
            2.0*aograd[ibas][2]*aograd[jbas][2]);
          k++;
        }
      }
      else k += (NBASIS-ibas);
    grdval.push_back( value );
  }
}

void PlotManager::EvaluateKED(double* P)
{
  int i,ibas,igp,jbas,k;
  double factor,g,value;
  double aograd[MAXNBAS][3];
  Vector pos;

  if (!P) return;

  if (MOS[NORB].spin[0]=='A') factor = 0.5;
  else factor = 1.0;

  grdval.clear();
  for ( igp = 0 ; igp < (signed)grdvec.size() ; igp++ )
  {
    value = 0.0;
    pos = Vector(grdvec[igp]);
    EvaluateAOGRAD(pos,aograd);
    k = 0;
    for (ibas=0;ibas<NBASIS;ibas++)
      for (jbas=ibas;jbas<NBASIS;jbas++)
      {
        g = 0.0;
        for (i=0;i<3;i++)
          g += aograd[ibas][i]*aograd[jbas][i];
        value += factor*P[k]*g;
        k++;
      }
    grdval.push_back( value );
  }
}

// This function evaluates electron localization function
// of Becke and Edgecombe and the density of localized electrons.
//
// Lit: A. D. Becke and K. E. Edgecombe, 
//      J. Chem. Phys. 92, 5397 (1990).
//
void PlotManager::EvaluateELF(double* P, bool DLE)
{
  if (!P) return;

  int i,ibas,igp,jbas,k;
  double ao[MAXNBAS];
  double aograd[MAXNBAS][3];
  double factor,rho,ked,d0,chi;
  double grho[3];
  Vector pos;

  if (MOS[NORB-1].spin[0]=='A') factor = 0.5;
  else factor = 1.0;

  grdval.clear();
  for ( igp = 0 ; igp < (signed)grdvec.size() ; igp++ )
  {
    pos = Vector(grdvec[igp]);
    EvaluateAO    (pos,ao    );
    rho = 0.0;
    k = 0;
    for (ibas=0;ibas<NBASIS;ibas++)
      for (jbas=ibas;jbas<NBASIS;jbas++)
      {
        rho += factor*P[k]*ao[ibas]*ao[jbas];
        k++;
      }
    if (rho<0.00001) 
    {
      grdval.push_back( 0.0 );
      continue;
    }
    d0 = (3.0/5.0)*pow(6.0*M_PI*M_PI,2.0/3.0)*pow(rho,5.0/3.0);
    EvaluateAOGRAD(pos,aograd);
    ked = 0.0;
    k = 0;
    for (ibas=0;ibas<NBASIS;ibas++)
      for (jbas=ibas;jbas<NBASIS;jbas++)
      {
        chi = 0.0;
        for (i=0;i<3;i++)
          chi += aograd[ibas][i]*aograd[jbas][i];
        ked += factor*P[k]*chi;
        k++;
      }
    k = 0;
    for (i=0;i<3;i++)
      grho[i] = 0.0;
    for (ibas=0;ibas<NBASIS;ibas++)
      for (jbas=ibas;jbas<NBASIS;jbas++)
      {
        if (jbas==ibas)
          for (i=0;i<3;i++)
            grho[i] += factor*2.0*P[k]*ao[ibas]*aograd[jbas][i];
        else
          for (i=0;i<3;i++)
            grho[i] += factor*P[k]*(ao[ibas]*aograd[jbas][i]+ao[jbas]*aograd[ibas][i]);
        k++;
      }
    chi = grho[0]*grho[0] + grho[1]*grho[1] + grho[2]*grho[2];
    chi = chi/rho;
    chi = (ked - 0.25*chi)/d0;
    if (DLE) 
    {
      grdval.push_back( rho/(1.0+chi*chi) );
    }
    else 
    {
      grdval.push_back( 1.0/(1.0+chi*chi) );
    }
  }
}

void PlotManager::EvaluateDAOBas(double dval[15][3], int id)
{
  int i,igp;
  double factor;
  double dao[MAXNBAS][3];
  Vector pos;

  for ( igp = 0 ; igp < (signed) grdvec.size() ; igp++ )
  {
    pos = Vector(grdvec[igp]);
    EvaluateAOGRAD(pos,dao);
   for (i=0;i<3;i++)
     dval[igp][i] = dao[id][i];
 }
}

void PlotManager::EvaluateDMO(double dval[15][3], MO mo)
{
  int i,ibas,igp;
  double factor;
  double dao[MAXNBAS][3];
  Vector pos;

  for ( igp = 0 ; igp < (signed) grdvec.size() ; igp++ )
  {
    pos = Vector(grdvec[igp]);
    EvaluateAOGRAD(pos,dao);

    for (i=0;i<3;i++)
      dval[igp][i] = 0.0;
    for (ibas=0;ibas<NBASIS;ibas++)
    {
      factor = mo.c[ibas];
      for (i=0;i<3;i++)
        dval[igp][i] += factor*dao[ibas][i];
    }
  }
}

void PlotManager::EvaluateDRHO(double dval[15][3], double* P)
{
  int i,ibas,igp,jbas,k;
  double ao[MAXNBAS];
  double dao[MAXNBAS][3];
  Vector pos;

  if (!P) return;
  for ( igp = 0 ; igp < (signed)grdvec.size() ; igp++ )
  {
    pos = Vector(grdvec[igp]);
    EvaluateAOGRAD(pos,dao);
    EvaluateAO(pos,ao);

    for (i=0;i<3;i++)
      dval[igp][i] = 0.0;
    k = 0;
    for (ibas=0;ibas<NBASIS;ibas++)
      if (ao[ibas]!=0.0)
        for (jbas=ibas;jbas<NBASIS;jbas++)
        {
          if (ao[jbas]!=0.0)
          {
            pos = Vector(&dao[jbas][0]);
            pos *= P[k]*ao[ibas];
            for (i=0;i<3;i++)
              dval[igp][i] += pos[i];
            pos = Vector(&dao[ibas][0]);
            pos *= P[k]*ao[jbas];
            for (i=0;i<3;i++)
              dval[igp][i] += pos[i];
          }
          k++;
        }
      else k += NBASIS - ibas;

  }
}

void PlotManager::FitBox( void )
{
  int iatom;
  double xmin,ymin,zmin,xmax,ymax,zmax;
  double room = 5.0; // Extra bohrs
if ((Gaussbool==0) && (Demonbool == 0))
  {
  xmin = 0.0;
  ymin = 0.0;
  zmin = 0.0;
  xmax = 0.0;
  ymax = 0.0;
  zmax = 0.0;
  for (iatom=0;iatom<NATOM;iatom++)
  {
    xmin = MIN(xmin,COORD[1][iatom][0]);
    xmax = MAX(xmax,COORD[1][iatom][0]);
    ymin = MIN(ymin,COORD[1][iatom][1]);
    ymax = MAX(ymax,COORD[1][iatom][1]);
    zmin = MIN(zmin,COORD[1][iatom][2]);
    zmax = MAX(zmax,COORD[1][iatom][2]);
  }
  xmin -= room;
  ymin -= room;
  zmin -= room;
  xmax += room;
  ymax += room;
  zmax += room;

  axis[0][0]->setValue(xmin);
  axis[0][1]->setValue(ymin);
  axis[0][2]->setValue(zmin);

  axis[1][0]->setValue(xmax);
  axis[1][1]->setValue(ymin);
  axis[1][2]->setValue(zmin);

  axis[2][0]->setValue(xmin);
  axis[2][1]->setValue(ymax);
  axis[2][2]->setValue(zmin);

  axis[3][0]->setValue(xmin);
  axis[3][1]->setValue(ymin);
  axis[3][2]->setValue(zmax);
} else{
    axis[0][0]->setValue(LATTICE[0]);
    axis[0][1]->setValue(LATTICE[1]);
    axis[0][2]->setValue(LATTICE[2]);

    axis[1][0]->setValue(LATTICE[3]);
    axis[1][1]->setValue(LATTICE[4]);
    axis[1][2]->setValue(LATTICE[5]);

    axis[2][0]->setValue(LATTICE[6]);
    axis[2][1]->setValue(LATTICE[7]);
    axis[2][2]->setValue(LATTICE[8]);

    axis[3][0]->setValue(LATTICE[9]);
    axis[3][1]->setValue(LATTICE[10]);
    axis[3][2]->setValue(LATTICE[11]);
}
  ChangeBox( 0.0 );
}

// Translate real value (0,1) to RGB color.
void PlotManager::RealToRGB(double r,double *rgb)
{
  int n,irgb[3];
  double eps = 0.00000001;

  n = int(255.0*MIN(1.0,(r-fmin)/(fmax-fmin+eps)));
  
  if (n<128) irgb[0] = 0;
  else if (n<192) irgb[0] = (n-128)*4 + 1;
  else irgb[0] = 255;
 
  if (n<1) irgb[1] = 0;
  else if (n<64) irgb[1] = 4*n-1;
  else if (n<192) irgb[1] = 255;
  else irgb[1] = 255 - 4*(n-192) - MOD(n,2);

  if (n<64) irgb[2] = 255;
  else if (n<128) irgb[2] = 255-4*(n-64)-MOD(n,2);
  else irgb[2] = 0;

  rgb[0] = irgb[0]/255.0;
  rgb[1] = irgb[1]/255.0;
  rgb[2] = irgb[2]/255.0;

  return; // rfm
  if (abs(n-128)<32) 
  {
    rgb[0] = 1.0;
    rgb[1] = 1.0;
    rgb[2] = 1.0;
  }
}


void PlotManager::EvaluateMEP(double *P)
{
  if (!P||(NBASIS<1)) return;

  int fa,fb,i,iatom,ibas,igp,ishl,j,jatom,jbas,jshl,k,la,lb,na,nb;
  double value;
  double ra[3],rb[3],rc[3];
  double **ib,**PP;

  ib = new double*[MAXNCO];
  for (i=0;i<MAXNCO;i++)
    ib[i] = new double[MAXNCO];

  PP = new double*[NBASIS];
  for (i=0;i<NBASIS;i++)
    PP[i] = new double[NBASIS];

  k = 0;
  for (ibas=0;ibas<NBASIS;ibas++)
  {
    for (jbas=ibas;jbas<NBASIS;jbas++)
    {
      if (ibas==jbas) 
      {
        PP[ibas][jbas] = P[k];
      }
      else
      {
        PP[ibas][jbas] = 0.5*P[k];
        PP[jbas][ibas] = 0.5*P[k];
      }
      k++;
    }
  }

  grdval.clear();
  for ( igp = 0 ; igp < (signed)grdvec.size() ; igp++ )
  {
    value = 0.0;
    rc[0] = grdvec[igp][0];
    rc[1] = grdvec[igp][1];
    rc[2] = grdvec[igp][2];

    fa = 0;
    for (iatom=0;iatom<NATOM;iatom++)
    {
      ra[0] = COORD[1][iatom][0];
      ra[1] = COORD[1][iatom][1];
      ra[2] = COORD[1][iatom][2];

      // Electronic contribution
      AtomBasis abi = mb.ab[iatom];
      for (ishl=0;ishl<abi.nshl;ishl++)
      {
        ShellBasis sa = abi.sb[ishl];
        la = sa.l;
        na = ((la+1)*(la+2))/2;
        fb = 0;
        for (jatom=0;jatom<NATOM;jatom++)
        {
          rb[0] = COORD[1][jatom][0];
          rb[1] = COORD[1][jatom][1];
          rb[2] = COORD[1][jatom][2];
          AtomBasis abj = mb.ab[jatom];
          for (jshl=0;jshl<abj.nshl;jshl++)
          {
            ShellBasis sb = abj.sb[jshl];
            lb = sb.l;
            nb = ((lb+1)*(lb+2))/2;
            gi->Core(ra,rb,rc,&sa,&sb,ib);
            for (i=0;i<na;i++)
            {
	      ibas = fa + i;
              for (j=0;j<nb;j++)
              {
	        jbas = fb + j;
		//cout << ibas << " " << jbas << " " << NBASIS << endl;
	        value = value - PP[ibas][jbas]*ib[i][j];
              }
	    }
	    fb = fb + nb;
	  }
	}
	fa = fa + na;
      }
      // Nuclear part
      ra[0] = ra[0] - rc[0];
      ra[1] = ra[1] - rc[1];
      ra[2] = ra[2] - rc[2];
      value = value + ATOMNO[iatom]/sqrt(ra[0]*ra[0]+ra[1]*ra[1]+ra[2]*ra[2]);
    }
    cout << value << endl;
    grdval.push_back( value );
  }

  for (i=0;i<NBASIS;i++)
    delete[] PP[i];
  delete[] PP;

  for (i=0;i<MAXNCO;i++)
    delete[] ib[i];
  delete[] ib;


}

/*
void PlotManager::TBoxRFM(double iso, double t[3][5][3], int* nt)
{
C     Order of the cube vertices:
C
C          4----------5    
C         /|         /|   
C        / |        / |
C       7----------6  |  
C       |  |       |  | 
C       |  |       |  |   
C       |  0-------|--1  
C       | /        | /  
C       |/         |/  
C       3----------2
C
  int i,j,k,indx; 
  double a,b,o[3];
  double vlist[12][3];

  // Origin
  o = {0.0,0.0,0.0};
  for ( i = 0 ; i < 8 ; i++ )
    for ( j = 0 ; j < 3 ; j++ )
      o[j] += grdvec[i][j];
  for ( j = 0 ; j < 3 ; j++ )
    o[j] = o[j]/8.0;

  *nt = 0;

  indx = 0;
  j = 1;
  for ( i = 0 ; i < 8 ; i++ )
  {
    if (grdval[i]<iso) indx = indx + j;
    j = j + j;
  }

  if (indx==0||indx==255) return;

  j = 1;
  for ( i = 0 ; i < 12 ; i++ )
  {
    if (edgetable[indx] & j)
    {
      a = grdval[ptr[0][i]];
      b = grdval[ptr[1][i]];
      // Get our approximate intersection ISO values
      if (ABS(iso-a)<0.00001) 
      {
        vlist[i][0] = grdvec[ptr[0][i]][0];
        vlist[i][1] = grdvec[ptr[0][i]][1];
        vlist[i][2] = grdvec[ptr[0][i]][2];
      }
      else if(ABS(iso-b)<0.00001) 
      {
        vlist[i][0] = grdvec[ptr[1][i]][0];
        vlist[i][1] = grdvec[ptr[1][i]][1];
        vlist[i][2] = grdvec[ptr[1][i]][2];
      }
      else if (ABS(b-a)<0.00001)
      {
        vlist[i][0] = grdvec[ptr[0][i]][0];
        vlist[i][1] = grdvec[ptr[0][i]][1];
        vlist[i][2] = grdvec[ptr[0][i]][2];
      }
      else
      {
        vlist[i][0] = grdvec[ptr[0][i]][0] + (iso-a)/(b-a)*
                     (grdvec[ptr[1][i]][0]-grdvec[ptr[0][i]][0]);
        vlist[i][1] = grdvec[ptr[0][i]][1] + (iso-a)/(b-a)*
                     (grdvec[ptr[1][i]][1]-grdvec[ptr[0][i]][1]);
        vlist[i][2] = grdvec[ptr[0][i]][2] + (iso-a)/(b-a)*
                     (grdvec[ptr[1][i]][2]-grdvec[ptr[0][i]][2]);
      }
    }
    j = j + j;
  }

  i = 0; 
  while ( tritable[indx][i] != -1)
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      for ( k = 0 ; k < 3 ; k++ )
      {
        t[j][*nt][k] = vlist[tritable[indx][i+j]][k];
      }
    }
    (*nt)++;
    i += 3; 
  }
}
*/

