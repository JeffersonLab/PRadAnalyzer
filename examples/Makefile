################################################################################
# Makefile for building some example programs using PRad analysis libraries    #
################################################################################

MAKEFILE      = Makefile
PRAD_LIBS     = ../libs

####### Compiler, tools and options
CC            = gcc
CXX           = g++
CXXFLAGS      = -std=c++11 -O1 -g -pipe -Wall -m64 -mtune=generic -fPIC \
				-fstack-protector-strong --param=ssp-buffer-size=4 \
				-grecord-gcc-switches
INCPATH       = -Iinclude -I$(PRAD_LIBS)/include -I$(ROOTSYS)/include
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
LFLAGS        = -Wl,-O1 -Wl,-z,relro
LIBS          = -lpthread -lexpat \
                -L$(PRAD_LIBS) -lprana -lcneural \
				-L$(ROOTSYS)/lib -lCore -lRint -lRIO -lNet -lHist \
				-lGraf -lGraf3d -lGpad -lTree -lPostscript -lMatrix \
				-lPhysics -lMathCore -lThread -lGui -lSpectrum
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = 

####### Targets, add new objects here

EXAMPLE_SRC   = src
EXAMPLE_DIR   = bin
EXAMPLE       = testGEM \
                testHyCal \
				testDST \
                testLMS \
                testEPICS \
                testMatch \
                testSim \
				testPerform \
                getAvgGain \
                replay \
                eventSelect \
                beamChargeCount \
				messReject \
				mollerCenter \
				physCalib \
				reconstruct \
				testNeural \
				neuralTrain \
				neuralReject

EXAMPLE_BINS  = $(addprefix $(EXAMPLE_DIR)/, $(EXAMPLE))


####### Build rules
first: all

all: $(EXAMPLE_BINS)

$(EXAMPLE_DIR)/%: $(EXAMPLE_SRC)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCPATH) $(LIBS)

####### Clean
clean:
	$(DEL_FILE) $(EXAMPLE_BINS)

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

