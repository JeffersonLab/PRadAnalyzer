################################################################################
# Makefile for building PRad analysis libraries                                #
################################################################################

LIB_DIR = libs
GUI_DIR = gui
EXE_DIR = examples
QT_MAKE = qmake-qt4

####### Build rules
first: all

.PHONY: lib gui exe

all: lib exe

lib:
	$(MAKE) -C $(LIB_DIR)

gui:
	$(QT_MAKE) $(GUI_DIR) -o $(GUI_DIR)
	$(MAKE) -C $(GUI_DIR)

exe:
	$(MAKE) -C $(EXE_DIR)

####### Clean
clean:
	$(MAKE) -C $(LIB_DIR) clean
	$(MAKE) -C $(GUI_DIR) clean
	$(MAKE) -C $(EXE_DIR) clean
