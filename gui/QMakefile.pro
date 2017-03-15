######################################################################
# PRad Event Viewer, project file for qmake
# It provides several optional components, they can be disabled by
# comment out the corresponding line
# Chao Peng
# 10/07/2016
######################################################################

greaterThan(QT_MAJOR_VERSION, 4) {
    QT += widgets concurrent
}

# enable multi threading
DEFINES += MULTI_THREAD

######################################################################
# optional components
######################################################################

# empty components
COMPONENTS = 

# enable online mode, it requires Event Transfer,
# it is the monitoring process from CODA group
COMPONENTS += ONLINE_MODE

# enable high voltage control, it requires CAENHVWrapper library
COMPONENTS += HV_CONTROL

# use standard evio libraries instead of self-defined function to read
# evio data files
#COMPONENTS += STANDARD_EVIO

# enable the reconstruction display in GUI
COMPONENTS += RECON_DISPLAY

######################################################################
# optional components end
######################################################################

######################################################################
# general config for qmake
######################################################################

CONFIG += c++11
QMAKE_CXXFLAGS += -std=c++11

OBJECTS_DIR = obj
MOC_DIR = qt_moc

TEMPLATE = app
TARGET = EventViewer
DEPENDPATH += 
INCLUDEPATH += include \
               ../libs/include \
               $$(ROOTSYS)/include

# Input
HEADERS += include/PRadEventViewer.h \
           include/HyCalModule.h \
           include/HyCalScene.h \
           include/HyCalView.h \
           include/Spectrum.h \
           include/SpectrumSettingPanel.h \
           include/HtmlDelegate.h \
           include/QRootCanvas.h \
           include/HistCanvas.h \
           include/LogsBox.h

SOURCES += src/main.cpp \
           src/PRadEventViewer.cpp \
           src/HyCalModule.cpp \
           src/HyCalScene.cpp \
           src/HyCalView.cpp \
           src/Spectrum.cpp \
           src/SpectrumSettingPanel.cpp \
           src/HtmlDelegate.cpp \
           src/QRootCanvas.cpp \
           src/HistCanvas.cpp \
           src/LogsBox.cpp

LIBS += -L../libs -lprana \
        -lexpat -lgfortran \
        -L$$(ROOTSYS)/lib -lCore -lRint -lRIO -lNet -lHist \
                          -lGraf -lGraf3d -lGpad -lTree \
                          -lPostscript -lMatrix -lPhysics \
                          -lMathCore -lThread -lGui -lSpectrum

######################################################################
# general config end
######################################################################

######################################################################
# implement self-defined components
######################################################################

contains(COMPONENTS, ONLINE_MODE) {
    DEFINES += USE_ONLINE_MODE
    HEADERS += include/PRadETChannel.h \
               include/PRadETStation.h \
               include/ETSettingPanel.h
    SOURCES += src/PRadETChannel.cpp \
               src/PRadETStation.cpp \
               src/ETSettingPanel.cpp
    INCLUDEPATH += $$(ET_INC)
    LIBS += -L$$(ET_LIB) -let
}

contains(COMPONENTS, HV_CONTROL) {
    DEFINES += USE_CAEN_HV
    HEADERS += include/PRadHVSystem.h \
               include/CAENHVSystem.h
    SOURCES += src/PRadHVSystem.cpp \
               src/CAENHVSystem.cpp
    INCLUDEPATH += ../thirdparty/include
    LIBS += -L$$(THIRD_LIB) -lcaenhvwrapper
}

contains(COMPONENTS, STANDARD_EVIO) {
    DEFINES += USE_EVIO_LIB
    !contains(INCLUDEPATH, ../thirdparty/include) {
        INCLUDEPATH += ../thirdparty/include
    }
    LIBS += -L$$(THIRD_LIB) -levio -levioxx
}

contains(COMPONENTS, RECON_DISPLAY) {
    DEFINES += RECON_DISPLAY
    HEADERS += include/ReconSettingPanel.h \
               include/MarkSettingWidget.h
    SOURCES += src/ReconSettingPanel.cpp \
               src/MarkSettingWidget.cpp
}

######################################################################
# self-defined components end
######################################################################

