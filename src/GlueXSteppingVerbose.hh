//
// GlueXSteppingVerbose class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

class GlueXSteppingVerbose;

#ifndef GlueXSteppingVerbose_h
#define GlueXSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class GlueXSteppingVerbose : public G4SteppingVerbose 
{
 public:
   
  GlueXSteppingVerbose();
 ~GlueXSteppingVerbose();

  void StepInfo();
  void TrackingStarted();

};

#endif
