// ******************************************************************
// Purpose: Declaration of class ListedCurve
//
// 2003-2012: Roberto Flores-Moreno
// ******************************************************************

#ifndef LISTEDCURVE_H
#define LISTEDCURVE_H

#include <vector>

#include <QtWidgets>

using namespace std;

class Parakata;
class CurvePlotter;

class ListedCurve : public QWidget
{
  Q_OBJECT

  public:
    ListedCurve( Parakata*, const char* );
   ~ListedCurve();

    Parakata *prk;

    void Enabled( bool );

  public slots:

    void Launch( void );
    void Load(QStringList,double*,double*,int);
    void ChangeSelection(int);
    void Advance(void);
    void SetAnimation(void);
    void ChangeAnimSpeed(int);

  Q_SIGNALS:

    void selectionChanged(int);

  protected:

    QListWidget *theListBox;
    QWidget* frame;
    CurvePlotter* plotter;
    QTimer* timer;
    QToolButton *animButton;
    int anim_delay;
 
};

#endif // LISTEDCURVE_H
