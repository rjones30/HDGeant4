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

class GlueXUserTrackInformation: public G4VUserTrackInformation
{
 public:
   int GetGlueXTrackID() { return fGlueXTrackID; }

   void SetGlueXTrackID(int trackID) {
      fGlueXTrackID = trackID;
   }

 private:
   int fGlueXTrackID;
};

#endif // _GLUEXUSERTRACKINFORMATION_
