//
// GlueXPrimaryGeneratorAction class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

GlueXPrimaryGeneratorAction::GlueXPrimaryGeneratorAction(
                                               GlueXDetectorConstruction* myDC)
:myDetector(myDC)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

// default particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
  
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.1,0.3));
  particleGun->SetParticleEnergy(0.3*GeV);
}

GlueXPrimaryGeneratorAction::~GlueXPrimaryGeneratorAction()
{
  delete particleGun;
}

void GlueXPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  G4double position = 0.65*m;
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));
  
  particleGun->GeneratePrimaryVertex(anEvent);
}
