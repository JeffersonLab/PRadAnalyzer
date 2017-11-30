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

cana:
	$(MAKE) -C cana -f $(MAKEFILE)

conf:
	$(MAKE) -C conf -f $(MAKEFILE)

cneural:
	$(MAKE) -C cneural -f $(MAKEFILE)

####### Clean
clean:
	$(MAKE) -C cana -f $(MAKEFILE) clean
	$(MAKE) -C conf -f $(MAKEFILE) clean
	$(MAKE) -C cneural -f $(MAKEFILE) clean
	$(MAKE) -C prana -f $(MAKEFILE) clean

