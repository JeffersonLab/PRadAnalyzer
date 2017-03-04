################################################################################
# Makefile for building: A simple neural network for pattern recognition       #
################################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options
CC            = gcc
CXX           = g++
FORTRAN       = gfortran
FFLAGS        = -fPIC
CXXFLAGS      = -shared -std=c++11 -O2 -g -pipe -Wall -m64 -fPIC
INCPATH       = -Iinclude
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
LIBS          = 
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = 

####### Targets, add new objects here
LIB_DIR       = ..
TARGET_LIB    = libcneural.so
LIB_OBJ_DIR   = obj
LIB_CLASSES   = CNeuron \
				CNeuronLayer \
				CNeuralNetwork

LIB_OBJECTS   = $(addprefix $(LIB_OBJ_DIR)/, $(LIB_CLASSES:=.o))


####### Build rules
first: all

all: lib

lib: Makefile $(TARGET_LIB)

$(LIB_OBJ_DIR)/%.o: src/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

$(TARGET_LIB):  $(LIB_OBJECTS)
	$(LINK) $(LFLAGS) -o $(LIB_DIR)/$(TARGET_LIB) $(LIB_OBJECTS) $(LIBS)

####### Clean
clean: cleanobj cleanlib

cleanobj:
	$(DEL_FILE) $(LIB_OBJ_DIR)/*.o

cleanlib:
	$(DEL_FILE) $(LIB_DIR)/$(TARGET_LIB)

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

