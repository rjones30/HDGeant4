//
// GlueXSensitiveDetectorPSC - class header
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorPSC_h
#define GlueXSensitiveDetectorPSC_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitPSCpaddle.hh"
#include "GlueXHitPSCpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorPSC : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorPSC(const G4String& name);
   GlueXSensitiveDetectorPSC(const GlueXSensitiveDetectorPSC &right);
   GlueXSensitiveDetectorPSC &operator=(const GlueXSensitiveDetectorPSC &right);
   virtual ~GlueXSensitiveDetectorPSC();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapPSCpaddle* fCounterHitsMap;
   GlueXHitsMapPSCpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static int NUM_MODULES_PER_ARM;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
