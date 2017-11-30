################################################################################
# Makefile for building: A simple neural network for pattern recognition       #
################################################################################

include ../general.mk

####### Targets, add new objects here
TARGET_LIB    = libcneural.so
LIB_OBJ_DIR   = obj
LIB_CLASSES   = CNeuralNetwork \
                CNeuronLayer \
                CNeuron

LIB_OBJECTS   = $(addprefix $(LIB_OBJ_DIR)/, $(LIB_CLASSES:=.o))

include ../rules.mk
