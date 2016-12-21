//
// GlueXSteppingAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXPathFinder.hh"

G4Mutex GlueXSteppingAction::fMutex = G4MUTEX_INITIALIZER;

void GlueXSteppingAction::UserSteppingAction(const G4Step* step)
{ 
   G4Track *track = step->GetTrack();
   G4int id = track->GetTrackID();
   G4int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   G4VPhysicalVolume *pvol = GlueXPathFinder::GetLocatedVolume();
   if (pvol) {
      if (pvol->GetName() == "PTAR" && id == 1 && pdgtype == 22) {
         const GlueXPrimaryGeneratorAction *generator;
         generator = GlueXPrimaryGeneratorAction::GetInstance();
         generator->GenerateBeamPairConversion(step);
         track->SetTrackStatus(fStopAndKill);
      }
   }
}
