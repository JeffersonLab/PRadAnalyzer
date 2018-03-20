################################################################################
# Makefile for building the core library for PRad analysis                     #
################################################################################


# passed by command line
#LIB_OPTION    = PRIMEX_METHOD, MULTI_THREAD

include ../general.mk
FFLAGS        = -fPIC -cpp -ffixed-line-length-none
CXXFLAGS      = -shared -std=c++11 -m64 -O2 -g -pipe -Wall $(DEFINES) \
				-fstack-protector-strong --param=ssp-buffer-size=4 \
				-grecord-gcc-switches -mtune=generic -fPIC
INCPATH       = -I$(ROOTSYS)/include -I../conf/include -I../cana/include
LIBS          = -lpthread -lgfortran \
				-L.. -lprconf -lcana \
                -L$(ROOTSYS)/lib -lCore -lRint -lRIO -lNet -lHist \
				-lGraf -lGraf3d -lGpad -lTree -lPostscript -lMatrix \
				-lPhysics -lMathCore -lThread -lGui -lSpectrum

####### Targets, add new objects here

# lib
TARGET_LIB    = libprana.so
OBJECTS_DIR   = obj

# c++ part
CXX_SOURCES   = PRadDAQChannel \
                PRadADCChannel \
                PRadTDCChannel \
                PRadCalibConst \
                PRadEvioParser \
                PRadDSTParser \
                PRadDataHandler \
                PRadException \
                PRadBenchMark \
                PRadDetector \
                PRadHyCalSystem \
                PRadHyCalDetector \
                PRadHyCalModule \
                PRadHyCalCluster \
                PRadHyCalReconstructor \
                PRadClusterProfile \
                PRadClusterDensity \
                PRadSquareCluster \
                PRadIslandCluster \
                PRadGEMSystem \
                PRadGEMDetector \
                PRadGEMPlane \
                PRadGEMFEC \
                PRadGEMAPV \
                PRadGEMCluster \
                PRadEventFilter \
                PRadCoordSystem \
                PRadDetMatch \
                PRadEPICSystem \
                PRadInfoCenter \
                PRadTaggerSystem \
                PRadTriggerConst \
                PRadMollerGen \
                PRadEpElasGen


# fortran part
FOR_SOURCES   = merad_util \
                merad_fsir \
                merad_sig \
                merad_sig2

HEADER_FILES  = include/*

###### Components related

# enable multi-threading in the code
ifneq (, $(findstring MULTI_THREAD,$(LIB_OPTION)))
	DEFINES     += -DMULTI_THREAD
endif

include ../rules.mk
