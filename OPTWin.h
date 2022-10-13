// ******************************************************************
// Purpose: OPTWin class.
//
// 2003-2012: Roberto Flores-Moreno
// ******************************************************************

#ifndef OPTWIN_H
#define OPTWIN_H

#include <vector>

#include <QWidget>

#include "types.h"

using namespace std;

class Parakata;
class QListWidget;
class QTabWidget;
class QToolButton;

class OPTWin : public QWidget
{
  Q_OBJECT

  public:
    OPTWin( Parakata* );
   ~OPTWin( void );

    Parakata *prk;

  public slots:

    void Launch( void );  
    void ChangeStructure( void );
    void Advance(void);
    void SetAnimation(void);
    void ChangeAnimSpeed(int);

  protected:

    void NewList( char** , int , char*); 

    QTabWidget *tabWidget;
    vector<QWidget*> tabs;
    vector<QListWidget*> listados;
    QToolButton *animButton;
    QTimer* timer;
    int anim_delay;
};

#endif  // OPTWIN_H
