// ******************************************************************
// Purpose: Declaration of class ImageRender
//
// 2003-2012: Roberto Flores-Moreno
// ******************************************************************

#ifndef IMAGERENDER_H
#define IMAGERENDER_H

#include <QWidget>

class QComboBox;
class QLineEdit;
class Viewer;

class ImageRender : public QWidget
{
  Q_OBJECT

  public:

    ImageRender( Viewer* );
   ~ImageRender( void );

    QLineEdit *label;

  public slots:

    void ChangeFormat( const QString & );
    void Write( void );

  protected:

    void DoAnimatedGIF( int );
    bool GetFileName( char* );

    Viewer* viewer;
    char filename[1024];
    QComboBox *formatCombo;
    QLineEdit *axised[3];

};

#endif  // IMAGERENDER_H
