// ******************************************************************
// Purpose: Declaration of class Parakata
//
// 2003-2022: Roberto Flores-Moreno
// ******************************************************************

#ifndef PARAKATA_H
#define PARAKATA_H

#include <vector>

#include <QtWidgets>

#define SFB    0x00000010L           // Selecting for bond
#define SFA    0x00000020L           // Selecting for angle
#define SFD    0x00000040L           // Selecting for dihedral
#define SFRV   0x00000080L           // Selecting for distance value
#define SFAV   0x00000100L           // Selecting for angle value
#define SFDV   0x00000200L           // Selecting for diheral value
#define SFSURF 0x00000400L           // Selecting for surface
#define SFPG   0x00000800L           // Selecting for point group
#define SFFV   0x00001000L           // Selecting for fitting value
#define SFPB   0x00002000L           // rfm: No se que es esto, Tu sabes Kayim? 
 

using namespace std;

class QSpinBox;
class QComboBox;
class QToolButton;
class QDoubleSpinBox;
class QResizeEvent;
class QCheckBox;

class ElementSelector;
class GeometryEditor;
class Viewer;
class EBEWin;
class FREQWin;
class SCFWin;
class OPTWin;
class GAFWin;
class Runner;

class Parakata : public QMainWindow
{
  Q_OBJECT

  public:

    Parakata( QApplication*, int argc, char** argv );
   ~Parakata();

    void Report( const char*, ... );

    QApplication* app;
    GeometryEditor *geoeditor;
    Viewer *viewer;
    FREQWin *freqwin;
    Runner *runner;
    QMenuBar* menu;

    EBEWin *ebewin;
    SCFWin *scfwin;
    OPTWin *optwin;
    GAFWin *gafwin;
    ElementSelector *es;

    long status;

    void SetCursor( const char* );

  public slots:
 
    void ChangeModel( const QString & );
    void About( void );
    void MakeStyle( const QString & );
    void DoReadFile(char*,char*);
    void ReadFile( void );
    void WriteFile( void );
    void LaunchGeoEditor( void );
    void LaunchSCFWin( void );
    void LaunchFREQWin( void );

  protected:

    void CreateActions();
    void CreateMenus();

    QStatusBar *statusbar;

    QMenu *fileMenu;
    QMenu *editMenu;
    QMenu *computeMenu;
    QMenu *colorMenu;
    QMenu *readMenu;
    QMenu *writeMenu;
    QMenu *helpMenu;

    vector<QAction*> RFMTA;
    vector<QAction*> WFMTA;
    QAction *quitAct;
    QAction *aboutAct;
    QAction *imgAct;
    QAction *coloratomAct;
    QAction *colorbgAct;
    QAction *colortagAct;
    QAction *colorarrowAct;
    QAction *colorsurfAct;
    QAction *ebeAct;
    QAction *optAct;
    QAction *plotAct;
    QAction *runAct;
    QAction *exitAction;
    QAction *buildAct;

    QToolButton *buildButton;
    QToolButton *geoButton;
    QToolButton *scfButton;
    QToolButton *freqButton;
    QToolButton *inButton;
    QToolButton *outButton;
    QToolButton *giraButton;
    QToolButton *numAtomButton;

    QComboBox *modelComboBox;
    QComboBox *styleComboBox;
    QComboBox *tagComboBox;
    QComboBox *arrowComboBox;
   
    QTimer* timer;
    QSpinBox *lineS;

    QTabWidget *tabWidget;

  protected slots:

    void NewMolecule(void);
    void ChangeAtomColor(void);
    void ChangeSurfaceColor(void);
    void WriteImage(void);
    void SetSpin(void);
    void AdvanceSpin(void);
    void SetupMeasurement(const QString &);

};

#endif // PARAKATA_H

