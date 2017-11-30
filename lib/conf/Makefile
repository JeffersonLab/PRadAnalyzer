################################################################################
# Makefile for building: A configuration library for parsing strings           #
################################################################################

include ../general.mk

####### Targets, add new objects here
TARGET_LIB    = libprconf.so
LIB_OBJ_DIR   = obj
LIB_CLASSES   = ConfigValue \
                ConfigParser \
                ConfigObject \
                ConfigOption

LIB_OBJECTS   = $(addprefix $(LIB_OBJ_DIR)/, $(LIB_CLASSES:=.o))

include ../rules.mk
