// ******************************************************************
// Purpose: Coral Program
//
// Roberto Flores Moreno
//
// ******************************************************************

#include <fcns.h>
#include <global.h>
#include <iostream>

#include <Coral.h>

using namespace std;

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("Coral developer");
  QCoreApplication::setApplicationName("Coral");
  QCoreApplication::setApplicationVersion(QT_VERSION_STR);
  QCommandLineParser parser;

  parser.setApplicationDescription(QCoreApplication::applicationName());
  parser.addHelpOption();
  parser.addVersionOption();
  parser.addPositionalArgument("file", "The file to open.");
  parser.process(app);

  Coral coral(&app,argc,argv);
  coral.show();
  return app.exec();
}
