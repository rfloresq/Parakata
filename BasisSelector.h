#ifndef BASISSELECTOR_H
#define BASISSELECTOR_H
 
#include <QtGui>
#include<QStringList>
#include<QTableWidget>
#include<QHBoxLayout>
#include<QPushButton>
#include<QLabel>
#include<vector>


using namespace std;

class Parakata;

class BasisSelector : public QWidget
{
    Q_OBJECT
 
public:
     BasisSelector( Parakata* );
    ~BasisSelector();
 
private:
   QTableWidget *basisSelector;
   QHBoxLayout *basSelectLayout;
   QStringList encabezados;
   QPushButton *but1;
   QLabel *lbl;
};
 
#endif // BASISSELECTOR_H
