#
# $Id: GNUmakefile,v 1.8 2009/11/20 19:04:17 soderber Exp $
#
# brebel@fnal.gov
#
#
ROOTCINT    := YES
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

ROOTLIBS  = $(shell root-config --evelibs) -lFFTW -lReflex -lEGPythia6

override LDFLAGS += $(shell root-config --evelibs) 

override LOADLIBES += -L$(ART_LIB) -lart_Framework_IO_Root_RootInput_source $(ROOTLIBS) $(ROOTGLIBS) -lRGL -lGed -lX3d

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)


########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_nutools.mk

override LIBLIBS += $(LOADLIBES) -L$(NUTOOLS_LIB) -lSimulationBase -lEventDisplayBase -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(LARSOFT)/lib/$(SRT_SUBDIR) -lFilters -lGeometry -lSimulation -lRawData -lRecoBase -lRecoObjects -lUtilities -lAnalysisBase -lAnalysisAlg

