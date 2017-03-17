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
# it can be changed by command line `qmake "COMPONENTS += something"`

# enable online mode, it requires Event Transfer,
# it is the monitoring process from CODA group
#COMPONENTS += ONLINE_MODE

# enable high voltage control, it requires CAENHVWrapper library
#COMPONENTS += HV_CONTROL

# use standard evio libraries instead of self-defined function to read
# evio data files
#COMPONENTS += STANDARD_EVIO

# enable the reconstruction display in GUI
#COMPONENTS += RECON_DISPLAY

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
    HEADERS += include/online_monitor/PRadETChannel.h \
               include/online_monitor/PRadETStation.h \
               include/online_monitor/ETSettingPanel.h
    SOURCES += src/online_monitor/PRadETChannel.cpp \
               src/online_monitor/PRadETStation.cpp \
               src/online_monitor/ETSettingPanel.cpp
    INCLUDEPATH += $$(ET_INC)
    LIBS += -L$$(ET_LIB) -let
    message("Online Monitor = Enabled")
} else {
    message("Online Monitor = Disabled")
}

contains(COMPONENTS, HV_CONTROL) {
    DEFINES += USE_CAEN_HV
    HEADERS += include/high_voltage/PRadHVSystem.h \
               include/high_voltage/CAENHVSystem.h
    SOURCES += src/high_voltage/PRadHVSystem.cpp \
               src/high_voltage/CAENHVSystem.cpp
    INCLUDEPATH += ../thirdparty/include
    LIBS += -L$$(THIRD_LIB) -lcaenhvwrapper
    message("High Voltage Control = Enalbed")
} else {
    message("High Voltage Control = Disabled")
}

contains(COMPONENTS, STANDARD_EVIO) {
    DEFINES += USE_EVIO_LIB
    !contains(INCLUDEPATH, ../thirdparty/include) {
        INCLUDEPATH += ../thirdparty/include
    }
    LIBS += -L$$(THIRD_LIB) -levio -levioxx
    message("EVIO Reading = Standard library")
} else {
    message("EVIO Reading = PRad specific")
}

contains(COMPONENTS, RECON_DISPLAY) {
    DEFINES += RECON_DISPLAY
    HEADERS += include/ReconSettingPanel.h \
               include/MarkSettingWidget.h
    SOURCES += src/ReconSettingPanel.cpp \
               src/MarkSettingWidget.cpp
    message("Reconstruct Events Display = Enabled")
} else {
    message("Reconstruct Events Display = Disabled")
}

######################################################################
# self-defined components end
######################################################################

