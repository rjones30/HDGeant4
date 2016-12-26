//
// GlueXSteppingAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "GlueXUserOptions.hh"
#include "GlueXPathFinder.hh"

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
      fStopTracksInCollimator = (showersInCol[1] == 0);
   }
}

GlueXSteppingAction::GlueXSteppingAction(const GlueXSteppingAction &src)
{
   fStopTracksInCollimator = src.fStopTracksInCollimator;
}

GlueXSteppingAction::~GlueXSteppingAction()
{}

void GlueXSteppingAction::UserSteppingAction(const G4Step* step)
{ 
   if (fStopTracksInCollimator) {
      G4VPhysicalVolume *pvol = GlueXPathFinder::GetLocatedVolume();
      if (pvol && (pvol->GetName() == "INSU" ||
                   pvol->GetName() == "PCTT" ||
                   pvol->GetName() == "PCPB") )
      {
         step->GetTrack()->SetTrackStatus(fStopAndKill);
      }
   }
}
