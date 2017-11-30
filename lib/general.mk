################################################################################
# Makefile for building: A simple neural network for pattern recognition       #
################################################################################

MAKEFILE      = Makefile
TARGET_DIR    = $(PRAD_LIB)
TARGET_INC    = $(PRAD_INC)

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
LFLAGS        = -shared
LIBS          = 
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = 

