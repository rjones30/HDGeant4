//
// GlueXSensitiveDetectorCERE - class header
//
// author: richard.t.jones at uconn.edu
// version: october 29, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorCERE_h
#define GlueXSensitiveDetectorCERE_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitCEREtube.hh"
#include "GlueXHitCEREpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorCERE : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorCERE(const G4String& name);
   GlueXSensitiveDetectorCERE(const GlueXSensitiveDetectorCERE &right);
   GlueXSensitiveDetectorCERE &operator=(const GlueXSensitiveDetectorCERE &right);
   virtual ~GlueXSensitiveDetectorCERE();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* ROhist);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapCEREtube* fTubeHitsMap;
   GlueXHitsMapCEREpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_PE;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
