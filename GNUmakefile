# $Id: GNUmakefile,v 1.1 1999-01-07 16:05:42 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := hdgeant4
G4TARGET := $(name)
G4EXLIB := true
G4LIB_BUILD_SHARED := true

CPPFLAGS += -I$(HDDS_HOME) -I./src 
CPPFLAGS += -I$(HALLD_HOME)/$(BMS_OSNAME)/include
CPPFLAGS += -I$(JANA_HOME)/include
CPPFLAGS += -I/usr/include/Qt
CPPFLAGS += -I/usr/include/python2.6
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

CORELIBS = -L$(HALLD_HOME)/$(BMS_OSNAME)/lib -lHDGEOMETRY -lDANA \
           -lANALYSIS -lBCAL -lCCAL -lCDC -lCERE -lDIRC -lFCAL \
           -lFDC -lFMWPC -lHDDM -lPAIR_SPECTROMETER -lPID -lRF \
           -lSTART_COUNTER -lTAGGER -lTOF -lTPOL -lTRACKING \
           -lTRIGGER -lDAQ -lTTAB \
           -lxstream -lbz2 -lz \
           -L$(JANA_HOME)/lib -lJANA \
           -L$(CCDB_HOME)/lib -lccdb \
           -L$(EVIOROOT)/lib -levioxx -levio \
           -L$(ROOTSYS)/lib -lCore -lPhysics -lTree -lHist -lGeom \
           -L$(CLHEP_LIB_DIR) -lCLHEP

.PHONY: all
all: hdds fixes exe lib bin g4py

include $(G4INSTALL)/config/binmake.gmk

CXXFLAGS = -g -fPIC -W -Wall -pedantic -Wno-non-virtual-dtor -Wno-long-long
INTYLIBS += -Wl,--whole-archive $(CORELIBS) -Wl,--no-whole-archive
#INTYLIBS += $(wildcard $(HALLD_HOME)/src/.$(BMS_OSNAME)/libraries/HDGEOMETRY/*.o)
INTYLIBS += -fPIC -I$(HDDS_HOME) -I/usr/local/xerces/include $(HDDSLIBS)
INTYLIBS += -L/usr/local/xerces/lib -lxerces-c
INTYLIBS += -L$(G4TMPDIR) -lhdds
INTYLIBS += -lboost_python

G4shared_libs := $(wildcard $(G4ROOT)/lib64/*.so)
INTYLIBS += -L $(G4ROOT)/lib64 $(patsubst $(G4ROOT)/lib64/lib%.so, -l%, $(G4shared_libs))

G4fix_sources := $(wildcard src/G4fixes/*.cc)
G4fix_objects := $(patsubst src/G4fixes/%.cc, $(G4TMPDIR)/%.o, $(G4fix_sources))

fixes: $(G4fix_objects)
hdds:  $(G4TMPDIR)/libhdds.so

HDDSDIR := $(G4TMPDIR)/hdds
HDDS_sources := $(HDDS_HOME)/XString.cpp $(HDDS_HOME)/XParsers.cpp $(HDDS_HOME)/hddsCommon.cpp
HDDS_objects := $(patsubst $(HDDS_HOME)/%.cpp, $(HDDSDIR)/%.o, $(HDDS_sources))
$(G4TMPDIR)/libhdds.so: $(HDDS_objects)
	@$(CXX) -Wl,-soname,$@ -shared -o $@ $^

$(G4TMPDIR)/%.o: src/G4fixes/%.cc
ifdef CPPVERBOSE
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
else
	@echo Compiling $*.cc ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
endif

$(HDDSDIR)/%.o: $(HDDS_HOME)/%.cpp
	@mkdir -p $(HDDSDIR)
ifdef CPPVERBOSE
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
else
	@echo Compiling $*.cc ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
endif

exe:
	@mkdir -p $(G4LIBDIR)/$@

g4py: $(G4LIBDIR)/../../../g4py/HDGeant4/libhdgeant4.so

$(G4LIBDIR)/../../../g4py/HDGeant4/libhdgeant4.so: $(G4LIBDIR)/libhdgeant4.so
	@rm -f $@
	@cd g4py/HDGeant4 && ln -s ../../tmp/*/hdgeant4/libhdgeant4.so .
