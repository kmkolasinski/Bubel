VPATH += ../shared
INCLUDEPATH += ../shared

HEADERS       = glwidget.h \
                window.h \
    datareader.h \
    mainwindow.h \
    formleadparams.h
SOURCES       = glwidget.cpp \
                main.cpp \
                window.cpp \
    datareader.cpp \
    mainwindow.cpp \
    formleadparams.cpp
QT           += opengl widgets
QT           += xml
LIBS += -lGLU
# install
INSTALLS += target sources

FORMS += \
    mainwindow.ui \
    formleadparams.ui
