// ******************************************************************
// Purpose: PlotManager class.
//
// 2003-2012: Roberto Flores-Moreno
// ******************************************************************

#ifndef PLOTMANAGER_H
#define PLOTMANAGER_H

#include <vector>

#include <QWidget>

#include "types.h"

#include "Vector.h"

#define PLOT_TAB_NROW 4
#define PLOT_TAB_NCOL 3

using namespace std;

class Parakata;
class Viewer;
class Integrator;
class Vector;
class QListWidget;
class QTabWidget;
class QVBoxLayout;
class QDoubleSpinBox;
class QCheckBox;
class QRadioButton;

class PlotManager : public QWidget
{
  Q_OBJECT

  public:
    PlotManager( Parakata* , Viewer* );
   ~PlotManager( void );

    Parakata *prk;
    Viewer *viewer;
    Integrator *gi;

    int NSurf( void ); 
    void GetColor( int , double* );
    void PlotSurface( int );
    void SetColor( int , double* );
    void SetUsageOfLists( bool );

    vector<vector<double> > cp;
    vector<vector<double> > bondpath;
    vector<double> temp;

  public slots:

    void ChangeIso( double );
    void ChangeContourMin( double );
    void ChangeContourMax( double );
    void ChangeContourStep( double );
    void ChangeType( const QString & );
    void ChangeStyle( const QString & );
    void ChangeLineWidth( int );
    void ChangePointSize( int );
    void ChangePlaneAxis( int );
    void ChangePlanePoint( int );
    void ChangeBox( double );
    void ChangeMesh( double );
    void ChangeFlags( void );  
    void Launch( void );  
    void Plot( void );  
    void Refresh( void );  

  protected:

    void BuildCube( int );
    void DoPlot( int );
    void NewList( char** , int , char* ); 
    void PlotISO( int );
    void PlotCONTOUR( int );
    void PlotCOLOURED_ISO( int );
    void PlotCOLOURED_PLANE( int );
    void PlotCURVE( int );
    void SetupBox( QVBoxLayout* );
    void FitBox(void);
    void TBox(double iso,double t[3][5][3],int* nt);
    void LBox(double iso,double l[2][2][3],int* nl);
    void GetNormals(char* name, double iso, 
                    double t[3][5][3], double tn[3][5][3],
                    int nt, int I, int J, int K);
    void GetCubePoints(int,int,int);
    void RealToRGB(double,double*);
    void EvaluateSet(char*);
    int EvaluateSetD(char*,double dval[15][3]);
    void EvaluateAO(Vector,double*);
    void EvaluateAOGRAD(Vector,double aog[MAXNBAS][3]);
    void EvaluateAOLAP(Vector,double*);
    void EvaluateMO(MO);
    void EvaluateKMO(MO,MO,bool);
    void EvaluateAOBas(int);
    void EvaluateDegenerateMOSRight(int,int);
    void EvaluateDegenerateMOSLeft(int,int);
    void EvaluateDegenerateDUAL(int,int,int,int);
    void EvaluateRHO(double*);
    void EvaluateMEP(double*);
    void EvaluateShapeEntropy(double*);
    void EvaluatedSsigma(double*,double*);
    void EvaluateLAPRHO(double*);
    void EvaluateKED(double*);
    void EvaluateELF(double*,bool);
    void EvaluateDMO(double dval[15][3],MO);
    void EvaluateDAOBas(double dval[15][3],int);
    void EvaluateDRHO(double dval[15][3], double*);

    int hold;
    int hide_surfaces;
    int drawbox;
    int drawref;
    int approximate_normals;
    int average_normals;

    int type;
    int style;
    int shiny;
    int line_width;
    int orbital_number;
    int point_size;
    int plane_axis;
    int lastplaneaxis;
    int plane_point;
    bool use_lists;
    double isovalue;
    double fmin;
    double fmax;
    double fstep;
    double boxmesh;
    QDoubleSpinBox* axis[PLOT_TAB_NROW][PLOT_TAB_NCOL];
    QCheckBox* flagbox[PLOT_TAB_NROW][PLOT_TAB_NCOL];
    QTabWidget *tabWidget;
    vector<QWidget*> tabs;
    vector<QListWidget*> listados;
    vector<PlotEntry> plot_entries;
    vector<double> grdval;
    vector<double> grdvaldegen;
    vector<Vector> grdvec;
    vector<Vector> grdvecdegen;
    Vector origin,dstep[3];

};

#endif  // PLOTMANAGER_H
