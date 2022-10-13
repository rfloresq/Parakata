// ******************************************************************
// Purpose: CurvePlotter object constructor/destructor.
//
// Roberto Flores-Moreno
// 2003-2005: Cinvestav-IPN
// 2007-2008: Universidad de Guanajuato
// 2008-2009: Cinvestav-IPN
// 2009-2010: Universidad de Guadalajara
// ******************************************************************

#include <iostream>

#include <CurvePlotter.h>

static int offset = 2;

CurvePlotter::CurvePlotter( QWidget* wdg )
            : QWidget( wdg )
{
  shape = Polygon;
  antialiased = false;

  sel = -1;

  x.clear();
  y.clear();
}

CurvePlotter::~CurvePlotter()
{
  x.clear();
  y.clear();
}

QSize CurvePlotter::minimumSizeHint() const
{
  return QSize(100, 100);
}

QSize CurvePlotter::sizeHint() const
{
  return QSize(400, 200);
}

void CurvePlotter::setShape(Shape shape)
{
  this->shape = shape;
  update();
}

void CurvePlotter::setPen(const QPen &pen)
{
    this->pen = pen;
    update();
}

void CurvePlotter::setBrush(const QBrush &brush)
{
    this->brush = brush;
    update();
}

void CurvePlotter::setAntialiased(bool antialiased)
{
    this->antialiased = antialiased;
    update();
}

void CurvePlotter::setTransformed(bool transformed)
{
    this->transformed = transformed;
    update();
}

void CurvePlotter::paintEvent(QPaintEvent *)
{
  int x1,x2,y1,y2,h,w;
  QPainter painter(this);
  painter.setPen( QPen(QColor( 255, 255, 255)));
  painter.setBrush(brush);
  if (antialiased)
    painter.setRenderHint(QPainter::Antialiasing);
  painter.setBackground( QBrush( QColor( 255, 255, 255, 0) ) );

  w = width();
  h = height();
  painter.fillRect( 0, 0 , w , h , QBrush(QColor( 85, 0, 127)));

  painter.drawLine( offset , h-offset , offset, offset );
  painter.drawLine( offset , h-offset , w-offset, h-offset );
  for ( size_t i = 1 ; i < x.size() ; i++ )
  {
    x1 = (int)( ((x[i-1]-xmin)/(xmax-xmin))*w);
    x2 = (int)( ((x[i]-xmin)/(xmax-xmin))*w);
    y1 = (int)( ((y[i-1]-ymin)/(ymax-ymin))*h);
    y2 = (int)( ((y[i]-ymin)/(ymax-ymin))*h);
    painter.drawLine( x1+2*offset , h-y1-2*offset , x2+2*offset, h-y2-2*offset );
  }
  for ( size_t i = 0 ; i < x.size() ; i++ )
  {
    x2 = (int)( ((x[i]-xmin)/(xmax-xmin))*w);
    y2 = (int)( ((y[i]-ymin)/(ymax-ymin))*h) + 1;
    painter.drawRect( x2+2*offset, h-y2-2*offset , 2 , 2 );
  }
  if ( sel >= 0 ) 
  {
    x2 = (int)( ((x[sel]-xmin)/(xmax-xmin))*w);
    y2 = (int)( ((y[sel]-ymin)/(ymax-ymin))*h) + 1;
    painter.fillRect( x2+2*offset-4, h-y2-2*offset-4, 12 , 12 , 
                      QBrush(QColor( 255, 170, 128)));
  }
}

void CurvePlotter::SetData( double* xx, double *yy, int n)
{
  x.clear();
  y.clear();
  for ( int i = 1 ; i <= n ; i++ )
  {
    x.push_back( xx[i-1] );
    y.push_back( yy[i-1] );
  }
  if ( x.size() == 0 || y.size() == 0) return;

  xmin = x[0];
  ymin = y[0];
  xmax = x[0];
  ymax = y[0];
  for ( int i = 0 ; i < n ; i++ )
  {
    if ( x[i] < xmin ) xmin = x[i];
    if ( x[i] > xmax ) xmax = x[i];
    if ( y[i] < ymin ) ymin = y[i];
    if ( y[i] > ymax ) ymax = y[i];
  }
}

void CurvePlotter::SetBounds( double xxmin, double xxmax ,
                              double yymin, double yymax )
{
  xmin = xxmin;
  xmax = xxmax;
  ymin = yymin;
  ymax = yymax;
}

void CurvePlotter::SaveImage()
{
}

void CurvePlotter::HighlightPoint(int n)
{
  if ( n < 0 || n >= x.size()) return;
  sel = n;
  update();
}
