//
// GlueXPhysicsList class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#ifndef GlueXPhysicsList_h
#define GlueXPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class GlueXPhysicsList: public G4VUserPhysicsList
{
  public:
    GlueXPhysicsList();

  protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
 
    void SetCuts();

    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();

  // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    void AddStepMax();
};

#endif
