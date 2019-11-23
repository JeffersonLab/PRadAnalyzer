################################################################################
# Makefile for building: A simple neural network for pattern recognition       #
################################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options
CC            = gcc
CXX           = g++
FORTRAN       = gfortran
FFLAGS        = -fPIC
CXXFLAGS      = -O2 -g -pipe -Wall -fPIC
CXXFLAGS     += $(shell root-config --cflags)
INCPATH       = -I. -Iinclude
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

# Cpp source files and fortran source files
CXX_SRC_DIR   = src
CXX_INC_DIR   = include
CXX_SUFFIX    = cpp

FOR_SRC_DIR   = fortran
FOR_INC_DIR   = fortran/include
FOR_SUFFIX    = f

