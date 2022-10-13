// *********************************************************
// Purpose: Declaration of class GeometryEditor
//
// 2005-2012, Roberto Flores Moreno
// *********************************************************

#ifndef GEOMETRYEDITOR_H
#define GEOMETRYEDITOR_H

#include <vector>

#include <QtWidgets>
#include <QStringList>

using namespace std;

class Parakata;

class GeometryEditor : public QWidget
{
  Q_OBJECT

  public:
    GeometryEditor( Parakata* );
   ~GeometryEditor();

    Parakata *prk;
    int default_element;
    bool use_default;
    char gname[80];
    QList<QString> elebas;
    QList<QString> nucbas;

    void Enabled( bool );
    void AddAtom( int, int, int );

  public slots:

    void Launch( void );
    void AddAtomStart( void );
    void DeleteAtom( void );
    void Edited(QTableWidgetItem*);
    void SetupMonitor( const QString & );
    void BuildGroup( const QString & );    
    void ChangeCoordinates( void );
    void Update( void );
    void FitDistance( void );
    void UseDefault( void );
    void GoRIS(int);
    void BuildTube(int);

  protected:

    QTableWidget *tw;
    QStringList *Headers;
};

#endif // GEOMETRYEDITOR_H
