// ******************************************************************
//
// Purpose: Declaration of class CurvePlotter 
//
// 2003-2012: Roberto Flores-Moreno
//
// ******************************************************************

#ifndef CURVEPLOTTER_H
#define CURVEPLOTTER_H

#include <vector>

#include <QtWidgets>

using namespace std;

class CurvePlotter : public QWidget
{
  Q_OBJECT

  public:

    enum Shape { Line, Points, Polyline, Polygon, 
                 Rect, RoundRect, Ellipse, Arc,
                 Chord, Pie, Path, Text };

    CurvePlotter(QWidget *parent = 0);
   ~CurvePlotter();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    void SetData( double* , double* , int );
    void SetBounds( double , double , double, double );

  public slots:

    void setShape(Shape shape);
    void setPen(const QPen &pen);
    void setBrush(const QBrush &brush);
    void setAntialiased(bool antialiased);
    void setTransformed(bool transformed);
    void SaveImage( void );
    void HighlightPoint( int );

  protected:

    void paintEvent(QPaintEvent *event);
    void DrawStrings( void );

    vector<double> x;
    vector<double> y;

  private:

    Shape shape;
    QPen pen;
    QBrush brush;
    bool antialiased;
    bool transformed;
    double xmin;
    double ymin;
    double xmax;
    double ymax;
    int sel;
};

#endif  // CURVEPLOTTER_H
