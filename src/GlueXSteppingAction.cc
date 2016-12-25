//
// GlueXSteppingAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserOptions.hh"
#include "GlueXPathFinder.hh"

// Jack up this threshold to 20GeV to disable this feature
// and let your application run without Dirac++ support.
#define FORCED_PTAR_PAIR_CONVERSION_THRESHOLD 3*GeV

G4Mutex GlueXSteppingAction::fMutex = G4MUTEX_INITIALIZER;

GlueXSteppingAction::GlueXSteppingAction()
{
   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXSteppingAction constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   std::map<int,double> showersInCol;
   if (user_opts->Find("SHOWERSINCOL", showersInCol)) {
      fShowersInCollimator = (showersInCol[1] != 0);
   }
}

GlueXSteppingAction::GlueXSteppingAction(const GlueXSteppingAction &src)
{
   fShowersInCollimator = src.fShowersInCollimator;
}

GlueXSteppingAction::~GlueXSteppingAction()
{}

void GlueXSteppingAction::UserSteppingAction(const G4Step* step)
{ 
   G4Track *track = step->GetTrack();
   G4int id = track->GetTrackID();
   G4int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   G4VPhysicalVolume *pvol = GlueXPathFinder::GetLocatedVolume();
   if (pvol && pvol->GetName() == "PTAR" && id == 1 && pdgtype == 22 &&
       track->GetKineticEnergy() > FORCED_PTAR_PAIR_CONVERSION_THRESHOLD)
   {
      const GlueXPrimaryGeneratorAction *generator;
      generator = GlueXPrimaryGeneratorAction::GetInstance();
      ((GlueXPrimaryGeneratorAction*)generator)->GenerateBeamPairConversion(step);
      track->SetTrackStatus(fStopAndKill);
   }
   else if (pvol && (pvol->GetName() == "INSU" ||
                     pvol->GetName() == "PCTT" ||
                     pvol->GetName() == "PCPB") )
   {
      track->SetTrackStatus(fStopAndKill);
   }
}
