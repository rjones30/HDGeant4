//
// GlueXSteppingAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXSteppingAction.hh"
#include "GlueXUserOptions.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXPathFinder.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Decay.hh"

#include <stdio.h>

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
   G4Track *track = (G4Track*)step->GetTrack();

   // Kill tracks at entry to primary collimator or active collimator
   // if this was asked for in the control file.
   if (fStopTracksInCollimator) {
      G4VPhysicalVolume *pvol = GlueXPathFinder::GetLocatedVolume();
      if (pvol && (pvol->GetName() == "INSU" ||
                   pvol->GetName() == "PCTT" ||
                   pvol->GetName() == "PCPB") )
      {
         track->SetTrackStatus(fStopAndKill);
      }
   }

   // Post new vertices to the MC record for primary particle decays
   GlueXUserTrackInformation *trackinfo;
   trackinfo = (GlueXUserTrackInformation*)track->GetUserInformation();
   if (trackinfo) {
      int primeID = trackinfo->GetGlueXTrackID();
      if (primeID > 0) {
         const G4VProcess* process;
         process = step->GetPostStepPoint()->GetProcessDefinedStep();
         if (dynamic_cast<const G4Decay*>(process)) {
            const G4Event *event = G4RunManager::GetRunManager()
                                                 ->GetCurrentEvent();
            GlueXUserEventInformation *eventinfo = (GlueXUserEventInformation*)
                                                   event->GetUserInformation();
            G4TrackVector &secondary = *(G4TrackVector*)
                                        step->GetSecondaryInCurrentStep();
            G4cout << "got a decay of track " << track->GetTrackID()
                   << "(gluex ID " << primeID << ")" << G4endl;
            G4TrackVector::iterator iter;
            for (iter = secondary.begin(); iter != secondary.end(); ++iter) {
               int newID = eventinfo->AssignNextGlueXTrackID();
               trackinfo = new GlueXUserTrackInformation();
               if (eventinfo) {
                  trackinfo->SetGlueXTrackID(newID);
               }
               (*iter)->SetUserInformation(trackinfo);
            }
            if (eventinfo) {
               char mech[5];
               sprintf(mech, "%c%c%c%c", 'D', 'C', 'A', 'Y');
               eventinfo->AddSecondaryVertex(secondary, primeID, *(int*)mech);
            }
            G4cout << "(total of " << secondary.size() << " daughters)" << G4endl;
         }
      }
   }

   // Kill neutrinos as soon as they start
   int PDGtype = track->GetDefinition()->GetPDGEncoding();
   switch (PDGtype) {
      case 12:
      case -12:
      case 14:
      case -14:
      case 16:
      case -16:
         G4cout << "Die, neutrino, die!" << G4endl;
         track->SetTrackStatus(fStopAndKill);
   }
}
