// ******************************************************************
// Purpose: Declaration of class Viewer
//
// 2003-2012: Roberto Flores-Moreno
// ******************************************************************

#ifndef VIEWER_H
#define VIEWER_H

#include <GL/glu.h>

#include <vector>
#include <QWheelEvent>
#include <QtOpenGL/QtOpenGL>

#include "types.h"

using namespace std;

class Vector;
class Coral;
class ImageRender;
class PlotManager;
class wheelEvent;

class Viewer : public QGLWidget 
{
  Q_OBJECT

  public:
    Viewer( Coral* );
   ~Viewer();

    Coral *coral;

    double trans[3];
    double scale;
    GLdouble rot[16];
    ImageRender* image_render;
    PlotManager *pmgr;

    void Project( void );
    void PickAtom( void );
    void PickSurface( void );

    void DrawSphere( double* , double, int );
    void DrawCylinder( double* , double* , double, int );
    void DrawCone( double* , double* , double, int );
    void DrawMolecule( int, int, double, double);
    void SaveSelection( int );
    void DrawCPs( void );
    void DrawSurfaces( void );
    void DrawStrings( void );
    void Rotate( double , double , double , double );
    void ResetRotation(void);
    void Write(char*,Vector,size_t);

    GLUquadricObj* q;
    vector<vector<size_t> > monitored_bonds;
    vector<vector<size_t> > monitored_angles;
    vector<vector<size_t> > monitored_dihedrals;

  public slots:
 
    void Redraw( void );
    void ZoomIn( void );
    void ZoomOut( void );
    void ChangeBackgroundColor( void );
    void ChangeTags( const QString & );
    void ChangeTagColor( void );
    void ChangeArrows( const QString & );
    void ChangeArrowColor( void );
    void wheelEvent(QWheelEvent *event);

  protected:

    QString tag_type;
    QString arrow_type;
    int mouse[2];
    float view_width;
    float view_height;
    float view_far;
    float view_near;
    vector<ImageLabel> user_label;

    void initializeGL( void );
    void resizeGL( int , int );
    void paintGL( void );

    int SelectedAtom( size_t );
    void DrawAtom( int , double, int );
    void DrawBond( int , int , double, int );

    double tag_color[4];
    double arrow_color[4];
    double cp_color[4];
    double bg[4];
    vector<size_t> selected_atoms;

  protected slots:

    void mousePressEvent( QMouseEvent* );
    void mouseMoveEvent( QMouseEvent* );
    void mouseReleaseEvent( QMouseEvent* );
    void keyPressEvent( QKeyEvent* );

};

#endif  // VIEWER_H
