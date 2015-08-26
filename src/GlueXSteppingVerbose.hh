//
// GlueXSteppingVerbose class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

class GlueXSteppingVerbose;

#ifndef GlueXSteppingVerbose_h
#define GlueXSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"
#include <pthread.h>

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
   static pthread_mutex_t *fMutex;
};

#endif
