//
// GlueXTrackingAction class header
//
// author: richard.t.jones at uconn.edu
// version: september 23, 2016

#include "GlueXTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "GlueXUserTrackInformation.hh"

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
   }
}

void GlueXTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
   G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
   if (secondaries) {
      GlueXUserTrackInformation* info = (GlueXUserTrackInformation*)
                                        aTrack->GetUserInformation();
      size_t nSeco = secondaries->size();
      if (nSeco > 0) {
         for(size_t i=0; i < nSeco; i++) { 
            GlueXUserTrackInformation* infoNew = new GlueXUserTrackInformation(info);
            (*secondaries)[i]->SetUserInformation(infoNew);
         }
      }
   }
}
