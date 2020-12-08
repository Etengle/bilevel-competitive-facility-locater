#-------------------------------------------------
#
# Project created by QtCreator 2019-07-23T20:00:33
#
#-------------------------------------------------

QT       += core gui

DEFINES += USE_QT

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

QMAKE_CXXFLAGS += -std=c++11 -O3 -fopenmp

QMAKE_LFLAGS += -lm   -lboost_system -lboost_filesystem
#-lgurobi_c++ -lgurobi81 -m64

LIBS += -L/home/jdp17/dcmaster/gurobi811/linux64/lib -L../lib


TARGET = BCFL
TEMPLATE = app


SOURCES += main.cpp\
		mainwindow.cpp \
		layer.cpp \
		parser.cpp

HEADERS  += mainwindow.h \
		layer.h \
		parser.h \
		util.h \
    location.h \
    family.h \
    facility.h

FORMS    += mainwindow.ui

DISTFILES += \
    parameter.txt \
    parameter_out.txt \
    parameterB.txt \
    parameterB2.txt \
    initialSol.txt \
    parameterSSSS.txt \
    parameter1111.txt \
    parameter50%.txt
