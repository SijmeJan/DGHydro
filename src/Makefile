# Directory to put binaries in
BINDIR = /Users/sjp/Codes/DGHydro/bin

# Path to Boost library
BOOST_PATH = /Users/sjp/Codes/DGHydro/libraries/boost_1_74_0

# In case SHELL is inherited from OS
SHELL = /bin/sh
# No implicit rules needed
.SUFFIXES:

################################################################################
# Create list of source, object and dependency files
################################################################################

# List of modules (must be directories in src/astrix)
MODULES := mesh config_file simulation state gauss_quad flux

# Create list of source files in module directories: list all .cpp and .cu files
SRC := $(wildcard *.cpp) $(foreach sdir,$(MODULES),$(wildcard $(sdir)/*.cpp))

# Create list of objects from .cu and .cpp in module directories
OBJ = $(patsubst %.cpp,%.o,$(SRC))

# Create list of dependency files from .cu and .cpp in module directories
DEP = $(patsubst %.cpp,%.d,$(SRC))

################################################################################
# Compiler and linker flags
################################################################################

CC=mpicxx

# Standard compiler and linker flags
CCFLAGS     := -O3 -std=c++17 -Wall
#-Wno-unused-private-field -Wno-unused-function
LDFLAGS     :=
#-lboost_program_options

# Common includes and paths for CUDA
INCLUDES  := -I $(BOOST_PATH)
LIBRARIES :=
#-L $(BOOST_PATH)/stage/lib
#LIBRARIES := -$(BOOST_PATH)/stage/lib/libboost_program_options.a

################################################################################
# Target rules
################################################################################

all: $(BINDIR)/dghydro

# Build executable
$(BINDIR)/dghydro: $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $+ $(LIBRARIES)
#$(BOOST_PATH)/stage/lib/libboost_program_options.a

# Clean up
clean:
	$(foreach sdir,$(MODULES),rm -f $(sdir)/*.o $(sdir)/*.d $(sdir)/*~)
	rm -f *.o *.d *~
	rm -f $(BINDIR)/dghydro

################################################################################
# Rules to convert .cpp into .o. Generate new dependency file .d
################################################################################

%.o:%.cpp
	$(CC) $(INCLUDES) $(CCFLAGS) -o $@ -c $<
	$(CC) $(INCLUDES) -E -MT $@ -MM -o $*.d $<

# Keep .d files
.PRECIOUS: %.d

##############################################################################
# Include dependency files
##############################################################################

-include $(DEP)
