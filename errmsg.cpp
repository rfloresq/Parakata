// *********************************************************
// Purpose: Error Message.
//
// Roberto Flores-Moreno
// *********************************************************

#include <QtWidgets>

#include <global.h>

void errmsg( const char* filename , int linenumber, const char*msg , int el )
{
  QString title = "Attention";
  QString text;
  if ( el )
  {
    text = "Error occurred in ";
  }
  else
  {
    text = "Warning received from "; 
  };
  text += filename;
  text += ", line number ";
  text += QString::number( linenumber );
  text += "\n";
  text += msg;
  text += "\n";

  if ( el )
  {
    (void)QMessageBox::critical( 0, title, text, QMessageBox::Ok,
                                 QMessageBox::NoButton, QMessageBox::NoButton);
  }
  else
  {
    (void)QMessageBox::warning( 0, title, text, QMessageBox::Ok,
                                 QMessageBox::NoButton, QMessageBox::NoButton);
  };

  if ( el )
  { 
    if ( prk->app ) 
    {
      prk->app->closeAllWindows();
      delete prk;
    }
  };
}
