//
// GlueXRunAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread,
// and one for the master thread as well.

#ifndef GlueXRunAction_h
#define GlueXRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4RotationMatrix.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "globals.hh"

class G4Run;

class GlueXRunAction : public G4UserRunAction
{
 public:
   void BeginOfRunAction(const G4Run*);
   void EndOfRunAction(const G4Run*);

 private:
   static G4Mutex fMutex;
};

#endif
