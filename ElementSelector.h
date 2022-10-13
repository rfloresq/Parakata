// ******************************************************************
// Purpose: Declaration of class ElementSelector 
//
// 2003-2012: Roberto Flores-Moreno
// ******************************************************************

#ifndef ELEMENTSELECTOR_H
#define ELEMENTSELECTOR_H

#include "qdialog.h"

class QToolButton;

class ElementSelector : public QDialog 
{
  Q_OBJECT

  public:
    ElementSelector( QWidget *parent = 0 , int guess=0 );
   ~ElementSelector();

  public slots:
    void SetDefaultElement( void );
    int GetSelectedElement( void );

  protected:
    QToolButton *element[120]; 
    int selected_element;

};

#endif
