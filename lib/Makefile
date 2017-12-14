################################################################################
# Makefile for building PRad analysis libraries                                #
################################################################################

include general.mk

####### Build rules
first: all

.PHONY: cana conf cneural prana

all: prana

prana: cana conf cneural
	$(MAKE) -C prana -f $(MAKEFILE) "LIB_OPTION = $(LIB_OPTION)"
	$(MAKE) -C prana -f $(MAKEFILE) install "INSTALL_DIR = $(PRAD_PATH)"

cana:
	$(MAKE) -C cana -f $(MAKEFILE)
	$(MAKE) -C cana -f $(MAKEFILE) install "INSTALL_DIR = $(PRAD_PATH)"

conf:
	$(MAKE) -C conf -f $(MAKEFILE)
	$(MAKE) -C conf -f $(MAKEFILE) install "INSTALL_DIR = $(PRAD_PATH)"

cneural:
	$(MAKE) -C cneural -f $(MAKEFILE)
	$(MAKE) -C cneural -f $(MAKEFILE) install "INSTALL_DIR = $(PRAD_PATH)"

####### Clean
clean:
	$(MAKE) -C cana -f $(MAKEFILE) clean
	$(MAKE) -C conf -f $(MAKEFILE) clean
	$(MAKE) -C cneural -f $(MAKEFILE) clean
	$(MAKE) -C prana -f $(MAKEFILE) clean
	$(MAKE) -C cana -f $(MAKEFILE) uninstall "INSTALL_DIR = $(PRAD_PATH)"
	$(MAKE) -C conf -f $(MAKEFILE) uninstall "INSTALL_DIR = $(PRAD_PATH)"
	$(MAKE) -C cneural -f $(MAKEFILE) uninstall "INSTALL_DIR = $(PRAD_PATH)"
	$(MAKE) -C prana -f $(MAKEFILE) uninstall "INSTALL_DIR = $(PRAD_PATH)"

