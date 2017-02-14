//
// GlueXPathFinder class header
//
// author: richard.t.jones at uconn.edu
// version: aug 12, 2015
//
// This is a wrapper around the standard G4PathFinder
// singleton class that adds some functionality needed
// for efficient lookup of volume information in a
// multi-layer mass geometry.
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.

#ifndef GlueXPathFinder_h
#define GlueXPathFinder_h 1

#include "G4VPhysicalVolume.hh"
#include "G4TouchableHandle.hh"

class GlueXPathFinder
{
 public:
   static G4VPhysicalVolume* GetLocatedVolume();
   static G4TouchableHandle CreateTouchableHandle();

 private:
   GlueXPathFinder() {}
   GlueXPathFinder(GlueXPathFinder&) {}
   GlueXPathFinder &operator=(GlueXPathFinder&) { return *this; }
};

#endif
