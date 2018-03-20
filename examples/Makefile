################################################################################
# Makefile for building some example programs using PRad analysis libraries    #
################################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options
CC            = gcc
CXX           = g++
CXXFLAGS      = -std=c++11 -O1 -g -pipe -Wall -m64 -mtune=generic -fPIC \
				-fstack-protector-strong --param=ssp-buffer-size=4 \
				-grecord-gcc-switches
INCPATH       = -Iinclude -I$(PRAD_INC) -I$(ROOTSYS)/include
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
LIBS          = -lpthread \
                -L$(PRAD_LIB) -lprana -lprconf -lcana -lcneural \
				-L$(ROOTSYS)/lib -lCore -lRint -lRIO -lNet -lHist \
				-lGraf -lGraf3d -lGpad -lTree -lPostscript -lMatrix \
				-lPhysics -lMathCore -lThread -lGui -lSpectrum
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = 

####### Targets, add new objects here

SOURCE_DIR    = src
BUILD_DIR     = bin
SOURCES      := $(shell find $(SOURCE_DIR) -name '*.cpp')
TARGETS      := $(addprefix $(BUILD_DIR)/, $(SOURCES:$(SOURCE_DIR)/%.cpp=%))


####### Build rules
first: all

all: $(TARGETS)

$(BUILD_DIR)/%: $(SOURCE_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCPATH) $(LIBS)

####### Clean
clean:
	$(DEL_FILE) $(TARGETS)

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

