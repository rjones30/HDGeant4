//
// GlueXSteppingAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#ifndef GlueXSteppingAction_h
#define GlueXSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class GlueXSteppingAction : public G4UserSteppingAction
{
  public:
    GlueXSteppingAction();
   ~GlueXSteppingAction(){};

    void UserSteppingAction(const G4Step*);
};

#endif
