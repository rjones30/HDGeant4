//
// GlueXSteppingVerbose class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.
// However its output actions need to be serialized, so it keeps
// its own set of interlocks for this purpose.

class GlueXSteppingVerbose;

#ifndef GlueXSteppingVerbose_h
#define GlueXSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

class GlueXSteppingVerbose : public G4SteppingVerbose 
{
 public:
   GlueXSteppingVerbose();
   GlueXSteppingVerbose(const GlueXSteppingVerbose &src);
   ~GlueXSteppingVerbose();
   void StepInfo();
   void TrackingStarted();

 private:
   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
