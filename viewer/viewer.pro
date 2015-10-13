VPATH += ../shared
INCLUDEPATH += ../shared

HEADERS       = glwidget.h \
                window.h \
                qtlogo.h \
    datareader.h \
    mainwindow.h \
    formleadparams.h
SOURCES       = glwidget.cpp \
                main.cpp \
                window.cpp \
                qtlogo.cpp \
    datareader.cpp \
    mainwindow.cpp \
    formleadparams.cpp
QT           += opengl
QT           += xml
LIBS += -lGLU
# install
INSTALLS += target sources

FORMS += \
    mainwindow.ui \
    formleadparams.ui
