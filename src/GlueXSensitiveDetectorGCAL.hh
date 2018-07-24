//
// GlueXSensitiveDetectorGCAL - class header
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorGCAL_h
#define GlueXSensitiveDetectorGCAL_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitGCALblock.hh"
#include "GlueXHitGCALpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorGCAL : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorGCAL(const G4String& name);
   GlueXSensitiveDetectorGCAL(const GlueXSensitiveDetectorGCAL &right);
   GlueXSensitiveDetectorGCAL &operator=(const GlueXSensitiveDetectorGCAL &right);
   virtual ~GlueXSensitiveDetectorGCAL();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapGCALblock* fBlockHitsMap;
   GlueXHitsMapGCALpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static double ATTENUATION_LENGTH;
   static double C_EFFECTIVE;
   static double LENGTH_OF_BLOCK;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
