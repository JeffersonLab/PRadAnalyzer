################################################################################
# Makefile for building PRad analysis libraries, gui and examples              #
################################################################################

# directories
LIB_DIR = lib
GUI_DIR = gui
EXE_DIR = examples

# makefile name
MAKE_FILE = Makefile

# Qt binary
QT_MAKE = qmake

# LIB setting, set options for the library
# MULTI_THREAD    multi threading support on decoding raw data file
LIB_OPTION = MULTI_THREAD

# GUI Setting, enable optional components for GUI here, separated by spaces
# HV_CONTROL      high voltage monitor and control
# ONLINE_MODE     online events monitor
# STANDARD_EVIO   using standard evio library instead of our specific code to read evio file
# RECON_DISPLAY   add reconstruction display setting panel, show reconstructed hits
GUI_OPTION = RECON_DISPLAY

####### Build rules
first: all

.PHONY: lib gui exe

all: lib exe

lib:
	$(MAKE) -C $(LIB_DIR) -f $(MAKE_FILE) "LIB_OPTION = $(LIB_OPTION)"

exe: lib
	$(MAKE) -C $(EXE_DIR) -f $(MAKE_FILE)

gui: lib
	$(QT_MAKE) $(GUI_DIR) -o $(GUI_DIR)/$(MAKE_FILE) "GUI_OPTION = $(GUI_OPTION)"
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
ifneq ("$(wildcard $(GUI_DIR)/$(MAKE_FILE))","")
	$(MAKE) -C $(GUI_DIR) -f $(MAKE_FILE) distclean
endif
# do nothing if the make file doesn't exist
