//
// GlueXSensitiveDetectorTPOL - class header
//
// author: richard.t.jones at uconn.edu
// version: december 16, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorTPOL_h
#define GlueXSensitiveDetectorTPOL_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitTPOLwedge.hh"
#include "GlueXHitTPOLpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorTPOL : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorTPOL(const G4String& name);
   GlueXSensitiveDetectorTPOL(const GlueXSensitiveDetectorTPOL &right);
   GlueXSensitiveDetectorTPOL &operator=(const GlueXSensitiveDetectorTPOL &right);
   virtual ~GlueXSensitiveDetectorTPOL();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapTPOLwedge* fHitsMap;
   GlueXHitsMapTPOLpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;
   const static int NSECTORS = 32;
   const static int NRINGS = 1;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
