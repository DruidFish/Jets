SHELL = /bin/sh
UNAME = $(shell uname)

# Root variables
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

################
##linux
CXX          = g++
RM           = rm -f
AR           = ar cru

##Flags
CXXFLAGS     = -O3 -g -fPIC -funroll-loops -Wall -std=c++11


EXENAME		= benJet
SRCEXT   	= cpp
SRCDIR  	= src
INCDIR   	= include
OBJDIR   	= build
EXEDIR  	= bin
SRCS    	:= $(shell find $(SRCDIR) -name '*.$(SRCEXT)')
OBJS    	:= $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))

GARBAGE  = $(OBJDIR)/*.o $(EXEDIR)/$(EXENAME)

#################
##Dependencies
# Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) $(ROOTCFLAGS) #-I$(GSLINC)
LINKFLAGS    = -g $(shell root-config --nonew) $(shell root-config --ldflags) -Wl,--no-as-needed
endif

# OS X
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) $(ROOTCFLAGS) #-I$(GSLINC)
LINKFLAGS    =
endif

##Libraries
LIBS       += $(ROOTLIBS) -ltbb

##Targets
all : $(EXEDIR)/$(EXENAME)

$(EXEDIR)/$(EXENAME) : $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LINKFLAGS) $(LIBS)

$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean   :
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE)
