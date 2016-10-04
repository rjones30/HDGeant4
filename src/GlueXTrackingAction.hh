//
// GlueXTrackingAction class header
//
// author: richard.t.jones at uconn.edu
// version: september 23, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.

#ifndef _GLUEXTRACKINGACTION_
#define _GLUEXTRACKINGACTION_

#include "G4UserTrackingAction.hh"
#include "GlueXUserTrackInformation.hh"

class GlueXTrackingAction : public G4UserTrackingAction
{
 public:  
   GlueXTrackingAction();
   ~GlueXTrackingAction();
   
   virtual void PreUserTrackingAction(const G4Track*);
   virtual void PostUserTrackingAction(const G4Track*);
};

#endif // _GLUEXTRACKINGACTION_
