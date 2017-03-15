################################################################################
# Makefile for building PRad analysis libraries, gui and examples              #
################################################################################

# directories
LIB_DIR = libs
GUI_DIR = gui
EXE_DIR = examples

# makefile name
MAKE_FILE = Makefile

# QT binary
QT_MAKE = qmake

####### Build rules
first: all

.PHONY: lib gui exe

all: lib exe gui

lib:
	$(MAKE) -C $(LIB_DIR) -f $(MAKE_FILE)

exe:
	$(MAKE) -C $(EXE_DIR) -f $(MAKE_FILE)

gui:
	$(QT_MAKE) $(GUI_DIR) -o $(GUI_DIR)/$(MAKE_FILE)
	$(MAKE) -C $(GUI_DIR) -f $(MAKE_FILE)

####### Clean
clean: cleanlib cleanexe cleangui

.PHONY: cleanlib cleanexe cleangui

cleanlib:
	$(MAKE) -C $(LIB_DIR) -f $(MAKE_FILE) clean

cleanexe:
	$(MAKE) -C $(EXE_DIR) -f $(MAKE_FILE) clean

# distclean is set by Qt, will remove the Makefile and target
# in addition to the objects
cleangui:
	$(MAKE) -C $(GUI_DIR) -f $(MAKE_FILE) distclean
