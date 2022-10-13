#include <math.h>

#include <QtGui>

#include <fcns.h>
#include <global.h>
#include <Parakata.h>
#include <Viewer.h>
#include <BasisSelector.h>

BasisSelector::BasisSelector( Parakata *is)
 : QWidget( 0 )
{

prk = is;

setWindowTitle( "Basis Selector" );

basSelectLayout = new QHBoxLayout(this);
encabezados << "No." << "Atom" << "Isotope" << "e-basis" << "nuc-basis";
basisSelector = new QTableWidget(this);
but1 = new QPushButton("Done",this);
basisSelector->setRowCount(5);
basisSelector->setHorizontalHeaderLabels(encabezados);
basSelectLayout->addWidget(basisSelector);
basSelectLayout->addWidget(but1);
}

BasisSelector::~BasisSelector()
{
}
