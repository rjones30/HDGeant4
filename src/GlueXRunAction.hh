//
// GlueXRunAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.
 
#ifndef GlueXRunAction_h
#define GlueXRunAction_h 1

#include "G4UserRunAction.hh"
#include "GlueXPhysicsList.hh"
#include "G4Run.hh"

class G4Run;

class GlueXRunAction : public G4UserRunAction
{
 public:
   GlueXRunAction(GlueXPhysicsList *plist);
   void BeginOfRunAction(const G4Run* run);
   void EndOfRunAction(const G4Run* run);

 protected:
   GlueXPhysicsList *fPhysicsList;
};

#endif
