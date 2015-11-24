TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp
INCLUDEPATH += /usr/include/eigen3 "/media/jaineel/Macintosh HD/Users/jaineeldalal/Robotics/Statistical Techniques in Robotics/Project Platypus/"

unix:!macx: LIBS += -L$$PWD/../../libgp/ -lgp

INCLUDEPATH += $$PWD/../../libgp
DEPENDPATH += $$PWD/../../libgp

unix:!macx: PRE_TARGETDEPS += $$PWD/../../libgp/libgp.a
