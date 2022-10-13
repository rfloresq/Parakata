// ******************************************************************
// Purpose: Declaration of class Runner
//
// 2012: Roberto Flores-Moreno
// ******************************************************************

#ifndef RUNNER_H
#define RUNNER_H

#include <QtWidgets>

class QComboBox;
class QCheckBox;
class QToolButton;
class QDoubleSpinBox;

class Runner : public QWidget 
{
  Q_OBJECT

  public:

    Runner(void);
   ~Runner();

  public slots:
 
    void Run(void);

  protected:

    QComboBox *programBox;
    QComboBox *theoryBox;
    QComboBox *basisBox;

    QCheckBox *optCheck;
    QCheckBox *freqCheck;
    QCheckBox *fukuiCheck;

    QToolButton *runButton;

    QDoubleSpinBox *chargeSpin;

  protected slots:

    void Prepare_deMon(void);
};

#endif // RUNNER_H

