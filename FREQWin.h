// ******************************************************************
// Purpose: Declaration of class FREQWin
//
// 2003-2012: Roberto Flores Moreno
// ******************************************************************

#ifndef FREQWIN_H
#define FREQWIN_H

#include <vector>

#include <QtWidgets>

using namespace std;

class Parakata;
class CurvePlotter;

class FREQWin : public QWidget
{
  Q_OBJECT

  public:
    FREQWin( Parakata* );
   ~FREQWin();

    Parakata *prk;

    void Enabled( bool );

  public slots:

    void Launch(void);
    void Update(void);
    void ChangeSelection(int);
    void Advance(void);
    void SetAnimation(void);
    void ChangeAnimSpeed(int);

  protected:

    QListWidget *theListBox;
    QWidget* frame;
    CurvePlotter* plotter;
    QTimer* timer;
    QToolButton *animButton;
    int anim_delay;
    int anim_pos;
    int num_anim_steps;
    bool forward;
 
};

#endif // FREQWIN_H
