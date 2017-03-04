################################################################################
# Makefile for building PRad analysis libraries                                #
################################################################################

MAKEFILE      = Makefile
COMPONENTS    = PRIMEX_METHOD, MULTI_THREAD

# neural network directory
NEURAL 		  = cneural

####### Compiler, tools and options
CC            = gcc
CXX           = g++
FORTRAN       = gfortran
FFLAGS        = -fPIC
CXXFLAGS      = -shared -std=c++11 -m64 -O2 -g -pipe -Wall $(DEFINES) \
				-fstack-protector-strong --param=ssp-buffer-size=4 \
				-grecord-gcc-switches -mtune=generic -fPIC
INCPATH       = -Iinclude -Iinclude -I$(ROOTSYS)/include
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
LFLAGS        = -shared -Wl,-O1 -Wl,-z,relro
LIBS          = -lpthread -lexpat -lgfortran \
                -L$(ROOTSYS)/lib -lCore -lRint -lRIO -lNet -lHist \
				-lGraf -lGraf3d -lGpad -lTree -lPostscript -lMatrix \
				-lPhysics -lMathCore -lThread -lGui -lSpectrum
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = 

####### Targets, add new objects here

# lib PRadDecoder
LIB_DIR       = .
TARGET_LIB    = libprana.so
LIB_SRC_DIR   = src
FORTRAN_SRC   = fortran
FORTRAN_INC   = fortran
LIB_OBJ_DIR   = $(LIB_DIR)/obj
LIB_CLASSES   = PRadDAQChannel \
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
                PRadTaggerSystem


LIB_OBJECTS   = $(addprefix $(LIB_OBJ_DIR)/, $(LIB_CLASSES:=.o))


###### Components related

# use primex clustering method
ifneq (,$(findstring PRIMEX_METHOD,$(COMPONENTS)))
	LIB_CLASSES +=  island PRadPrimexCluster
	DEFINES     +=  -DUSE_PRIMEX_METHOD
endif

# enable multi-threading in the code
ifneq (, $(findstring MULTI_THREAD,$(COMPONENTS)))
	DEFINES     += -DMULTI_THREAD
endif

####### Build rules
first: all neural

.PHONY: neural cleanobj cleanlib

all:
	$(MAKE) $(TARGET_LIB)

neural:
	$(MAKE) -C $(NEURAL)

$(TARGET_LIB):  $(LIB_OBJECTS)
	$(LINK) $(LFLAGS) -o $(LIB_DIR)/$(TARGET_LIB) $(LIB_OBJECTS) $(LIBS)

$(LIB_OBJ_DIR)/%.o: $(LIB_SRC_DIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

$(LIB_OBJ_DIR)/island.o: $(FORTRAN_SRC)/island.F
	$(FORTRAN) -c $(FFLAGS) $(INCPATH) -I$(FORTRAN_INC) -o $@ $<

####### Clean
clean: cleanobj cleanlib

cleanobj:
	$(DEL_FILE) $(LIB_OBJ_DIR)/*.o $(EXE_OBJ_DIR)/*.o
	$(MAKE) -C $(NEURAL) cleanobj

cleanlib:
	$(DEL_FILE) $(LIB_DIR)/$(TARGET_LIB)
	$(MAKE) -C $(NEURAL) cleanlib

