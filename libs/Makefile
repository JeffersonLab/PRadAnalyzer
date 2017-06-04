################################################################################
# Makefile for building PRad analysis libraries                                #
################################################################################

MAKEFILE      = Makefile
# passed by command line
#LIB_OPTION    = PRIMEX_METHOD, MULTI_THREAD

# neural network directory
NEURAL 	      = cneural

####### Compiler, tools and options
CC            = gcc
CXX           = g++
FORTRAN       = gfortran
FFLAGS        = -fPIC -cpp -ffixed-line-length-none
CXXFLAGS      = -shared -std=c++11 -m64 -O2 -g -pipe -Wall $(DEFINES) \
				-fstack-protector-strong --param=ssp-buffer-size=4 \
				-grecord-gcc-switches -mtune=generic -fPIC
INCPATH       = -I$(ROOTSYS)/include
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
LINK          = g++
LFLAGS        = -shared
LIBS          = -lpthread -lexpat -lgfortran \
                -L$(ROOTSYS)/lib -lCore -lRint -lRIO -lNet -lHist \
				-lGraf -lGraf3d -lGpad -lTree -lPostscript -lMatrix \
				-lPhysics -lMathCore -lThread -lGui -lSpectrum
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = 

####### Targets, add new objects here

# lib
TARGET_LIB    = libprana.so
OBJECTS_DIR   = obj

# c++ part
CXX_SRC_DIR   = src
CXX_INC_DIR   = include
CXX_SUFFIX    = cpp
CXX_SOURCES   = canalib \
                PRadDAQChannel \
                PRadADCChannel \
                PRadTDCChannel \
                PRadCalibConst \
                PRadEvioParser \
                PRadDSTParser \
                PRadDataHandler \
                PRadException \
                PRadBenchMark \
                ConfigParser \
                ConfigValue \
                ConfigObject \
                PRadDetector \
                PRadHyCalSystem \
                PRadHyCalDetector \
                PRadHyCalModule \
                PRadHyCalCluster \
                PRadClusterProfile \
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
				PRadMollerGen


# fortran part
FOR_SRC_DIR   = fortran
FOR_INC_DIR   = fortran/include
FOR_SOURCES   = merad_util \
                merad_fsir \
                merad_sig \
                merad_sig2
FOR_SUFFIX    = f


###### Components related

# use primex clustering method
ifneq (,$(findstring PRIMEX_METHOD,$(LIB_OPTION)))
	FOR_SOURCES += island
	CXX_SOURCES += PRadPrimexCluster
	DEFINES     += -DUSE_PRIMEX_METHOD
endif

# enable multi-threading in the code
ifneq (, $(findstring MULTI_THREAD,$(LIB_OPTION)))
	DEFINES     += -DMULTI_THREAD
endif


# find all objects to be made
CXX_OBJECTS   = $(addprefix $(OBJECTS_DIR)/, $(CXX_SOURCES:=.$(CXX_SUFFIX).o))
FOR_OBJECTS   = $(addprefix $(OBJECTS_DIR)/, $(FOR_SOURCES:=.$(FOR_SUFFIX).o))
OBJECTS       = $(CXX_OBJECTS) $(FOR_OBJECTS)

####### Build rules
first: all neural

.PHONY: neural cleanobj cleanlib

all:
	$(MAKE) $(TARGET_LIB)

neural:
	$(MAKE) -C $(NEURAL)

$(TARGET_LIB):  $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET_LIB) $(OBJECTS) $(LIBS)

# cpp objects
$(OBJECTS_DIR)/%.$(CXX_SUFFIX).o: $(CXX_SRC_DIR)/%.$(CXX_SUFFIX)
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -I$(CXX_INC_DIR) -o $@ $<

# fortran objects
$(OBJECTS_DIR)/%.$(FOR_SUFFIX).o: $(FOR_SRC_DIR)/%.$(FOR_SUFFIX)
	$(FORTRAN) -c $(FFLAGS) $(INCPATH) -I$(FOR_INC_DIR) -o $@ $<

####### Clean
clean: cleanobj cleanlib

cleanobj:
	$(DEL_FILE) $(OBJECTS_DIR)/*.o
	$(MAKE) -C $(NEURAL) cleanobj

cleanlib:
	$(DEL_FILE) $(TARGET_LIB)
	$(MAKE) -C $(NEURAL) cleanlib

