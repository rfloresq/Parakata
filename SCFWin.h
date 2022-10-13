// ******************************************************************
// Purpose: Declaration of class SCFWin
//
// 2003-2012: Roberto Flores Moreno
// ******************************************************************

#ifndef SCFWIN_H
#define SCFWIN_H

#include <QtWidgets>

class Parakata;
class CurvePlotter;

class SCFWin : public QWidget
{
  Q_OBJECT

  public:
    SCFWin( Parakata* );
   ~SCFWin();

    Parakata *prk;

    void Enabled( bool );

  public slots:

    void Launch( void );

  protected:

    void Load( void );

    QTabWidget *tabWidget;
    QWidget* frame;
    CurvePlotter* SCFPlotter;
 
};

#endif // SCFWIN_H
