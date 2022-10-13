TEMPLATE = app
TARGET = parakata
DEPENDPATH += .
INCLUDEPATH += .
LIBS += -lGLU

# Input
QT += opengl 
HEADERS += CurvePlotter.h \
           ElementSelector.h \
           fcns.h \
           FREQWin.h \
           global.h \
           ImageRender.h \
           ListedCurve.h \
           macros.h \
           EBEWin.h \
           OPTWin.h \
           param.h \
           PlotManager.h \
           SCFWin.h \
           GeometryEditor.h \ 
           Parakata.h \
           types.h \
           Vector.h \
           Runner.h \
           Math.h \
           Integrator.h \
           Viewer.h 
SOURCES += basis.cpp \
           center.cpp \
           CurvePlotter.cpp \
           distance.cpp \
           ElementSelector.cpp \
           errmsg.cpp \
           factorial.cpp \
           GeometryEditor.cpp \
           FREQWin.cpp \
           global.cpp \
           ImageRender.cpp \
           ListedCurve.cpp \
           main.cpp \
           EBEWin.cpp \
           OPTWin.cpp \
           PlotManager.cpp \
           read_deMon_inp.cpp \
           read_deMon_out.cpp \
           read_Nagual_out.cpp \
           read_Parakata_cube.cpp \
           readcif.cpp \
           reader.cpp \
           readxyz.cpp \
           SCFWin.cpp \
           Parakata.cpp \
           structure.cpp \
           symtoan.cpp \
           Vector.cpp \
           Viewer.cpp \
           writer.cpp \
           writexyz.cpp \
    zctocc.cpp \
    sgdrv.cpp \
    sgbdrv.cpp \
    pgbdrv.cpp \
    apptrans.cpp \
    appcgroup.cpp \
    appdgroup.cpp \
    appsgroup.cpp \
    apptgroup.cpp \
    appogroup.cpp \
    appigroup.cpp \
    clonatom.cpp \
    appcaxis.cpp \
    appsaxis.cpp \
    appsigma.cpp \
    Math.cpp \
    Integrator.cpp \
    Runner.cpp 
