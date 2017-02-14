//
// GlueXUserTrackInformation class header
//
// author: richard.t.jones at uconn.edu
// version: aug 1, 2015
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.

#ifndef _GLUEXUSERTRACKINFORMATION_
#define _GLUEXUSERTRACKINFORMATION_

#include "G4VUserTrackInformation.hh"
#include "GlueXUserEventInformation.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"

class GlueXUserTrackInformation: public G4VUserTrackInformation
{
 public:
   GlueXUserTrackInformation() 
   : fGlueXTrackID(0), fGlueXHistory(0), fBCALincidentID(0)
   {}

   GlueXUserTrackInformation(const G4Track *trk)
   : fGlueXTrackID(trk->GetTrackID()), fGlueXHistory(0), fBCALincidentID(0)
   {}

   GlueXUserTrackInformation(GlueXUserTrackInformation *info)
   : fGlueXTrackID(info->GetGlueXTrackID()),
     fGlueXHistory(info->GetGlueXHistory()),
     fBCALincidentID(info->GetBCALincidentID())
   {}

   int GetGlueXTrackID() { return fGlueXTrackID; }
   int GetGlueXHistory() { return fGlueXHistory; }
   int GetBCALincidentID() { return fBCALincidentID; }

   void SetGlueXTrackID(int trackID) {
      fGlueXTrackID = trackID;
   }
   void SetGlueXHistory(int history) {
      fGlueXHistory = history;
   }
   void AssignBCALincidentID(G4Track *track) {
      const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
      GlueXUserEventInformation *eventinfo = (GlueXUserEventInformation*)
                                             event->GetUserInformation();
      fBCALincidentID = eventinfo->AssignBCALincidentID(track);
   }

   void Print() const {
      G4cout << "GlueXUserTrackInformation: id=" << fGlueXTrackID
             << ", history=" << fGlueXHistory << G4endl;
   }

 private:
   int fGlueXTrackID;
   int fGlueXHistory;
   int fBCALincidentID;
};

#endif // _GLUEXUSERTRACKINFORMATION_
