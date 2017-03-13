# Package information
LIBRARY = SubstructureAnalysis
OBJDIR  = obj
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include

# configure fastjet include:
INCLUDES += -I$(FASTJETDIR)/../include
INCLUDES += -I$(FASTJETDIR)/../include/contribs/RecursiveTools

#INCLUDES += -I$(LHAPDFDIR)/include
#INCLUDES += -I/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.8/x86_64-slc5-gcc46-opt/include
INCLUDES += -I/nfs/dust/cms/user/dreyert/LHAPDF-install2/include/

ifneq ($(BOOSTDIR),)
INCLUDES += -I$(BOOSTDIR)
endif

# in case you need to link to an external library, use USERLDFLAGS, e.g. like that:
#USERLDFLAGS += $(ROOTLIBS) -lMinuit -L /nfs/dust/cms/user/dreyert/LHAPDF-install/lib/libLHAPDF.so
USERLDFLAGS += $(ROOTLIBS) -lMinuit
USERLDFLAGS +=  /nfs/dust/cms/user/dreyert/fastjet-install/lib/libRecursiveTools.a
USERLDFLAGS +=  /nfs/dust/cms/user/dreyert/fastjet-install/lib/libJetsWithoutJets.a
USERLDFLAGS +=  /nfs/dust/cms/user/dreyert/fastjet-install/lib/libEnergyCorrelator.a
USERLDFLAGS +=  /nfs/dust/cms/user/dreyert/fastjet-install/lib/libVariableR.a	
USERLDFLAGS +=  /nfs/dust/cms/user/dreyert/fastjet-install/lib/libClusteringVetoPlugin.a
USERLDFLAGS +=  /nfs/dust/cms/user/dreyert/fastjet-install/lib/libfastjetplugins.a
#USERLDFLAGS +=  /nfs/dust/cms/user/dreyert//LHAPDF-install2/lib/libLHAPDF.a

# to pass additional compiler flags, set USERCXXFLAGS. Add debugging info here:
USERCXXFLAGS := -g 

# Include the generic compilation rules
#include $(SFRAME_DIR)/Makefile.common
include $(SFRAME_DIR)/SFrameTools/Makefile.defs
