# $Id: GNUmakefile,v 1.1 1999-01-07 16:05:42 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := hdgeant4
G4TARGET := $(name)
G4EXLIB := true

CPPFLAGS += -I$(HDDS_HOME) -I./src 
CPPFLAGS += -I$(HALLD_HOME)/$(BMS_OSNAME)/include
CPPFLAGS += -I$(JANA_HOME)/include
CPPFLAGS += -I/usr/include/Qt
CPPFLAGS += -DUSE_SSE2
#CPPFLAGS += -I/usr/include/Qt
#CPPFLAGS += -DMOD_SPONCE
#CPPFLAGS += -DLINUX_CPUTIME_PROFILING=1
#CPPFLAGS += -DDEBUG_SECTIONPLANE
#CPPFLAGS += -DDEBUG_SECTIONPLANE_ZAVE
#CPPFLAGS += -DCHECK_OVERLAPS_MM=1e-4
CPPFLAGS += -DBYPASS_DRAWING_CLIPPED_VOLUMES
#CPPFLAGS += -DG4UI_USE_EXECUTIVE
#CPPFLAGS += -DDEBUG_PLACEMENT
CPPFLAGS += -DBP_DEBUG
CPPFLAGS += -DG4VIS_BUILD_OPENGL_DRIVER
CPPFLAGS += -DG4VIS_BUILD_OPENGLX_DRIVER

G4LIB_USE_GDML = 1
CPPVERBOSE = 1
#G4DEBUG = 1

EXTRALIBS = $(HDDS_HOME)/hddsCommon.cpp \
            $(HDDS_HOME)/XParsers.cpp \
            $(HDDS_HOME)/XString.cpp \
            $(HDDS_HOME)/md5.c

LDLIBS2EXTRA = -L$(HALLD_HOME)/$(BMS_OSNAME)/lib -lHDGEOMETRY -lDANA \
               -lANALYSIS -lBCAL -lCCAL -lCDC -lCERE -lDIRC -lFCAL \
               -lFDC -lFMWPC -lHDDM -lPAIR_SPECTROMETER -lPID -lRF \
               -lSTART_COUNTER -lTAGGER -lTOF -lTPOL -lTRACKING \
               -lTRIGGER -lDAQ -lTTAB -lHDGEOMETRY \
               -lxstream -lbz2 -lz \
               -L$(JANA_HOME)/lib -lJANA \
               -L$(CCDB_HOME)/lib -lccdb \
               -L$(EVIOROOT)/lib -levioxx -levio \
               -L$(ROOTSYS)/lib -lCore -lPhysics -lTree -lHist -lGeom \
               -L$(CLHEP_LIB_DIR) -lCLHEP

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

CXXFLAGS = -fPIC -W -Wall -pedantic -Wno-non-virtual-dtor -Wno-long-long

test-shared-lib:
	@echo $(CXX) -Wl,-soname,$(@F) -shared -o $$libdir/$(@F) $(INTYLIBS) *.o
