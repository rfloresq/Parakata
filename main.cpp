// ******************************************************************
// Purpose: Parakata Program
//
// Roberto Flores Moreno
//
// ******************************************************************

#include <fcns.h>
#include <global.h>
#include <iostream>

#include <Parakata.h>

using namespace std;

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("Parakata developer");
  QCoreApplication::setApplicationName("Parakata");
  QCoreApplication::setApplicationVersion(QT_VERSION_STR);
  QCommandLineParser parser;

  parser.setApplicationDescription(QCoreApplication::applicationName());
  parser.addHelpOption();
  parser.addVersionOption();
  parser.addPositionalArgument("file", "The file to open.");
  parser.process(app);

  Parakata prk(&app,argc,argv);
  prk.show();
  return app.exec();
}
