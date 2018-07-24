//
// GlueXSensitiveDetectorCCAL - class header
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorCCAL_h
#define GlueXSensitiveDetectorCCAL_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitCCALblock.hh"
#include "GlueXHitCCALpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorCCAL : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorCCAL(const G4String& name);
   GlueXSensitiveDetectorCCAL(const GlueXSensitiveDetectorCCAL &right);
   GlueXSensitiveDetectorCCAL &operator=(const GlueXSensitiveDetectorCCAL &right);
   virtual ~GlueXSensitiveDetectorCCAL();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* ROhist);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapCCALblock* fBlocksMap;
   GlueXHitsMapCCALpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static int CENTRAL_COLUMN;
   static int CENTRAL_ROW;
   static double WIDTH_OF_BLOCK;
   static double LENGTH_OF_BLOCK;
   static double ATTENUATION_LENGTH;
   static double C_EFFECTIVE;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
