//
// GlueXTrackingAction class header
//
// author: richard.t.jones at uconn.edu
// version: september 23, 2016

#include "GlueXTrackingAction.hh"
#include "G4RunManager.hh"
#include "G4TrackVector.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserEventInformation.hh"

GlueXTrackingAction::GlueXTrackingAction()
{;}

GlueXTrackingAction::~GlueXTrackingAction()
{;}

void GlueXTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
   if (aTrack->GetParentID() == 0 && aTrack->GetUserInformation() == 0) {
      GlueXUserTrackInformation* anInfo = new GlueXUserTrackInformation(aTrack);
      G4Track* theTrack = (G4Track*)aTrack;
      theTrack->SetUserInformation(anInfo);
      const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
      GlueXUserEventInformation *eventinfo = (GlueXUserEventInformation*)
                                             event->GetUserInformation();
      anInfo->SetGlueXTrackID(eventinfo->GetGlueXTrackID(aTrack));
      G4cout << "starting tracking of primary " << aTrack->GetTrackID()
             << " with gluexID " << eventinfo->GetGlueXTrackID(aTrack)
             << " and type " << aTrack->GetDefinition()->GetParticleName()
             << G4endl;
   }
   else if (aTrack->GetUserInformation() == 0) {
      G4cerr << "GlueXTrackingAction::PreUserTrackingAction error - "
                "secondary track found without UserTrackInformation, "
                "fatal error, aborting." << G4endl;
      exit(91);
   }
   else {
      GlueXUserTrackInformation* info = (GlueXUserTrackInformation*)
                                        aTrack->GetUserInformation();
      int gluexID = info->GetGlueXTrackID();
      if (gluexID > 0 && fTrackingCounter[gluexID]++ == 0) {
         G4cout << "starting tracking of secondary " << aTrack->GetTrackID()
                << " with gluexID " << gluexID
                << " and type " << aTrack->GetDefinition()->GetParticleName()
                << G4endl;
      }
   }
}

void GlueXTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
   G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
   if (secondaries) {
      GlueXUserTrackInformation* info = (GlueXUserTrackInformation*)
                                        aTrack->GetUserInformation();
      size_t nSeco = secondaries->size();
      for (size_t i=0; i < nSeco; i++) { 
         if ((*secondaries)[i]->GetUserInformation() == 0) {
            GlueXUserTrackInformation* infoNew = new GlueXUserTrackInformation(info);
            int gluexID = infoNew->GetGlueXTrackID();
            if (gluexID > 0)
               infoNew->SetGlueXTrackID(-gluexID);
            (*secondaries)[i]->SetUserInformation(infoNew);
         }
      }
   }
}
