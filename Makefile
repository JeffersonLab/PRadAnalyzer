################################################################################
# Makefile for building PRad analysis libraries, gui and examples              #
################################################################################

# directories
LIB_DIR = libs
GUI_DIR = gui
EXE_DIR = examples

# QT binary
QT_MAKE = qmake-qt4

####### Build rules
first: all

.PHONY: lib gui exe

all: lib exe gui

lib:
	$(MAKE) -C $(LIB_DIR)

gui:
	$(QT_MAKE) $(GUI_DIR) -o $(GUI_DIR)
	$(MAKE) -C $(GUI_DIR)

exe:
	$(MAKE) -C $(EXE_DIR)

####### Clean
clean: cleanlib cleanexe cleangui

.PHONY: cleanlib cleangui cleanexe

cleanlib:
	$(MAKE) -C $(LIB_DIR) clean

cleangui:
	$(MAKE) -C $(GUI_DIR) clean

cleanexe:
	$(MAKE) -C $(EXE_DIR) clean
