//
// GlueXSensitiveDetectorUPV - class header
//
// author: richard.t.jones at uconn.edu
// version: october 29, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorUPV_h
#define GlueXSensitiveDetectorUPV_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitUPVbar.hh"
#include "GlueXHitUPVpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorUPV : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorUPV(const G4String& name);
   GlueXSensitiveDetectorUPV(const GlueXSensitiveDetectorUPV &right);
   GlueXSensitiveDetectorUPV &operator=(const GlueXSensitiveDetectorUPV &right);
   virtual ~GlueXSensitiveDetectorUPV();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapUPVbar* fBarHitsMap;
   GlueXHitsMapUPVpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static double ATTENUATION_LENGTH;
   static double C_EFFECTIVE;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
