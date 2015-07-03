//
// GlueXRunAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#ifndef GlueXRunAction_h
#define GlueXRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class GlueXRunAction : public G4UserRunAction
{
  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
};

#endif
