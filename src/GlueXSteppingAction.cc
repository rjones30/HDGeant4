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
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Decay.hh"

#include <stdio.h>

#define NEUTRON_KINETIC_ENERGY_CUTOFF_MEV 0.1

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

   std::map<int,int> trajectories;
   if (user_opts->Find("TRAJECTORIES", trajectories)) {
      fSaveTrajectories = trajectories[1];
   }
   else {
      fSaveTrajectories = 0;
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
   const G4Event *event = G4RunManager::GetRunManager()
                                      ->GetCurrentEvent();
   GlueXUserEventInformation *eventinfo = (GlueXUserEventInformation*)
                                          event->GetUserInformation();
   GlueXUserTrackInformation *trackinfo;
   trackinfo = (GlueXUserTrackInformation*)track->GetUserInformation();
   G4VPhysicalVolume *pvol = GlueXPathFinder::GetLocatedVolume();

   // Kill tracks at entry to primary collimator or active collimator
   // if this was asked for in the control file.
   if (fStopTracksInCollimator) {
      if (pvol && (pvol->GetName() == "INSU" ||
                   pvol->GetName() == "PCTT" ||
                   pvol->GetName() == "PCPB") )
      {
         track->SetTrackStatus(fStopAndKill);
      }
   }

   // Post new vertices to the MC record for primary particle decays
   if (trackinfo) {
      int primeID = trackinfo->GetGlueXTrackID();
      if (primeID > 0) {
         const G4VProcess* process;
         process = step->GetPostStepPoint()->GetProcessDefinedStep();
         if (dynamic_cast<const G4Decay*>(process)) {
            G4TrackVector &secondary = *(G4TrackVector*)
                                        step->GetSecondaryInCurrentStep();
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
         }
      }
   }

   // Kill neutrinos as soon as they start,
   // neutrons below some energy cutoff
   int PDGtype = track->GetDefinition()->GetPDGEncoding();
   switch (PDGtype) {
      case 12:
      case -12:
      case 14:
      case -14:
      case 16:
      case -16:
         track->SetTrackStatus(fStopAndKill);
         break;
      case 2112:
         double Ekin = step->GetPreStepPoint()->GetKineticEnergy();
         if (Ekin < NEUTRON_KINETIC_ENERGY_CUTOFF_MEV*MeV)
            track->SetTrackStatus(fStopAndKill);
         break;
   }

   // Save mc trajectory information if requested
   if (fSaveTrajectories) {
      eventinfo->AddMCtrajectoryPoint(*step, fSaveTrajectories);
   }
}
