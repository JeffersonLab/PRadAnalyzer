####### Build rules
first: all

all: dir lib

dir:
	$(MKDIR) $(LIB_OBJ_DIR)
	$(MKDIR) $(TARGET_INC)

lib: $(MAKEFILE) $(TARGET_LIB)
	$(MOVE) $(TARGET_LIB) $(TARGET_DIR)

$(LIB_OBJ_DIR)/%.o: src/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

$(TARGET_LIB):  $(LIB_OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET_LIB) $(LIB_OBJECTS) $(LIBS)
	$(COPY) include/* $(TARGET_INC)

####### Clean
clean: cleanobj cleanlib

cleanobj:
	$(DEL_FILE) $(LIB_OBJ_DIR)/*.o

cleanlib:
	$(DEL_FILE) $(TARGET_DIR)/$(TARGET_LIB)

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:
