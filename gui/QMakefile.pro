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
# it can be changed by command line `qmake "GUI_OPTION += something"`

# enable online mode, it requires Event Transfer,
# it is the monitoring process from CODA group
#GUI_OPTION += ONLINE_MODE

# enable high voltage control, it requires CAENHVWrapper library
#GUI_OPTION += HV_CONTROL

# use standard evio libraries instead of self-defined function to read
# evio data files
#GUI_OPTION += STANDARD_EVIO

# enable the reconstruction display in GUI
#GUI_OPTION += RECON_DISPLAY

######################################################################
# optional components end
######################################################################

######################################################################
# general config for qmake
######################################################################

QMAKE_CXXFLAGS += $$system(root-config --cflags)

OBJECTS_DIR = obj
MOC_DIR = qt_moc

TEMPLATE = app
TARGET = EventViewer
DEPENDPATH += 
INCLUDEPATH += include \
               $$(PRAD_INC) \
	       $$system(root-config --incdir)

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

LIBS += -L$$(PRAD_LIB) -lprana -lprconf -lcana \
        -lgfortran \
	$$system(root-config --libs) \
	-L$$system(root-config --libdir) -lSpectrum


######################################################################
# general config end
######################################################################

######################################################################
# implement self-defined components
######################################################################

contains(GUI_OPTION, ONLINE_MODE) {
    DEFINES += USE_ONLINE_MODE
    HEADERS += include/online_monitor/PRadETChannel.h \
               include/online_monitor/PRadETStation.h \
               include/online_monitor/ETSettingPanel.h
    SOURCES += src/online_monitor/PRadETChannel.cpp \
               src/online_monitor/PRadETStation.cpp \
               src/online_monitor/ETSettingPanel.cpp
    INCLUDEPATH += $$(ET_INC)
    LIBS += -L$$(ET_LIB) -let -lexpat
    message("Online Monitor = Enabled")
} else {
    message("Online Monitor = Disabled")
}

contains(GUI_OPTION, HV_CONTROL) {
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

contains(GUI_OPTION, STANDARD_EVIO) {
    DEFINES += USE_EVIO_LIB
    !contains(INCLUDEPATH, ../thirdparty/include) {
        INCLUDEPATH += ../thirdparty/include
    }
    LIBS += -L$$(THIRD_LIB) -levio -levioxx
    message("EVIO Reading = Standard library")
} else {
    message("EVIO Reading = PRad structure")
}

contains(GUI_OPTION, RECON_DISPLAY) {
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

