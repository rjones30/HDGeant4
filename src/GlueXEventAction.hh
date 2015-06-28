//
// GlueXEventAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
 
#ifndef GlueXEventAction_h
#define GlueXEventAction_h 1

#include "G4UserEventAction.hh"

class G4Event;

class GlueXEventAction : public G4UserEventAction
{
  public:
    GlueXEventAction();
   ~GlueXEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
};

#endif
