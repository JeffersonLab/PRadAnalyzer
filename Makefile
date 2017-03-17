################################################################################
# Makefile for building PRad analysis libraries, gui and examples              #
################################################################################

# directories
LIB_DIR = libs
GUI_DIR = gui
EXE_DIR = examples

# makefile name
MAKE_FILE = Makefile

# Qt binary
QT_MAKE = qmake

# LIB setting, set options for the library
# MULTI_THREAD    multi threading support on decoding raw data file
# PRIMEX_METHOD   add original PrimEx island reconstruction method
LIB_OPTION = "OPTION = MULTI_THREAD PRIMEX_METHOD"

# GUI Setting, enable optional components for GUI here
# HV_CONTROL      high voltage monitor and control
# ONLINE_MODE     online events monitor
# STANDARD_EVIO   using standard evio library instead of our specific code to read evio file
# RECON_DISPLAY   add reconstruction display setting panel, show reconstructed hits
GUI_OPTION = "COMPONENTS = RECON_DISPLAY"

####### Build rules
first: all

.PHONY: lib gui exe

all: lib exe gui

lib:
	$(MAKE) -C $(LIB_DIR) -f $(MAKE_FILE) $(LIB_OPTION)

exe:
	$(MAKE) -C $(EXE_DIR) -f $(MAKE_FILE)

gui:
	$(QT_MAKE) $(GUI_DIR) -o $(GUI_DIR)/$(MAKE_FILE) $(GUI_OPTION)
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
