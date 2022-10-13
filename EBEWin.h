// ******************************************************************
// Purpose: EBEWin class.
//
// 2003-2012: Roberto Flores-Moreno
// ******************************************************************

#ifndef EBEWIN_H
#define EBEWIN_H

#include <vector>

#include <QWidget>

#include "types.h"

using namespace std;

class Parakata;
class QListWidget;
class QTabWidget;
class QToolButton;

class EBEWin : public QWidget
{
  Q_OBJECT

  public:
    EBEWin( Parakata* );
   ~EBEWin( void );

    Parakata *prk;

  public slots:

    void Launch( void );  

  protected:

    void NewList( char** , char**, int , char*); 

    QTabWidget *tabWidget;
    vector<QWidget*> tabs;
    vector<QListWidget*> listados;
};

#endif  // OPTWIN_H
