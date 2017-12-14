
# firstly find all objects to be made
CXX_OBJECTS   = $(addprefix $(OBJECTS_DIR)/, $(CXX_SOURCES:=.$(CXX_SUFFIX).o))
FOR_OBJECTS   = $(addprefix $(OBJECTS_DIR)/, $(FOR_SOURCES:=.$(FOR_SUFFIX).o))
OBJECTS       = $(CXX_OBJECTS) $(FOR_OBJECTS)

####### Build rules
first: all

all: dir lib

dir:
	$(MKDIR) $(OBJECTS_DIR)

lib: $(MAKEFILE) $(TARGET_LIB)

$(TARGET_LIB):  $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET_LIB) $(OBJECTS) $(LIBS)

# cpp objects
$(OBJECTS_DIR)/%.$(CXX_SUFFIX).o: $(CXX_SRC_DIR)/%.$(CXX_SUFFIX)
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -I$(CXX_INC_DIR) -o $@ $<

# fortran objects
$(OBJECTS_DIR)/%.$(FOR_SUFFIX).o: $(FOR_SRC_DIR)/%.$(FOR_SUFFIX)
	$(FORTRAN) -c $(FFLAGS) $(INCPATH) -I$(FOR_INC_DIR) -o $@ $<


####### Clean
clean:
	$(DEL_FILE) $(OBJECTS_DIR)/*.o

####### Install

install:
	$(MOVE) $(TARGET_LIB) $(INSTALL_DIR)/lib
	$(SYMLINK) -r -t $(INSTALL_DIR)/include $(HEADER_FILES)

uninstall:
	$(DEL_FILE) $(INSTALL_DIR)/lib/$(TARGET_LIB)

FORCE:
