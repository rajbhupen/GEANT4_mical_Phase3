# $Id: GNUmakefile 68058 2013-03-13 14:47:43Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for inoical_v20161223.  
# Apoorva Bhatt, (apoorva.bhatt@tifr.res.in)
# V. Pavan Kumar, (pavan.v@tifr.res.in)
# 23rd Dec, 2016
# --------------------------------------------------------------

name := mICAL
G4TARGET := $(name)
G4EXLIB := true

ifndef INOICAL
  INOICAL := true
endif

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

CPPFLAGS += $(shell root-config --cflags)
EXTRALIBS := $(shell root-config --libs) $(shell root-config --glibs) -lconfig++ -lpq  -ltbb -lGeom 

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

