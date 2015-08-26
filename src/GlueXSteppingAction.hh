//
// GlueXSteppingAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.

#ifndef GlueXSteppingAction_h
#define GlueXSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class GlueXSteppingAction : public G4UserSteppingAction
{
  public:
    void UserSteppingAction(const G4Step*);
};

#endif
