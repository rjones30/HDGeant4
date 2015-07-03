//
// GlueXPrimaryGeneratorAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
 
#ifndef GlueXPrimaryGeneratorAction_h
#define GlueXPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class GlueXDetectorConstruction;
class G4ParticleGun;
class G4Event;

class GlueXPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    GlueXPrimaryGeneratorAction(GlueXDetectorConstruction*);    
    GlueXPrimaryGeneratorAction(const GlueXPrimaryGeneratorAction &src);
    GlueXPrimaryGeneratorAction &operator=(const GlueXPrimaryGeneratorAction &src);
   ~GlueXPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun* particleGun;
    GlueXDetectorConstruction* myDetector;
};

#endif
