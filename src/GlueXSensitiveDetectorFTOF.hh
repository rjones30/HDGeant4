//
// GlueXSensitiveDetectorFTOF - class header
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorFTOF_h
#define GlueXSensitiveDetectorFTOF_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitFTOFbar.hh"
#include "GlueXHitFTOFpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorFTOF : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorFTOF(const G4String& name);
   GlueXSensitiveDetectorFTOF(const GlueXSensitiveDetectorFTOF &right);
   GlueXSensitiveDetectorFTOF &operator=(const GlueXSensitiveDetectorFTOF &right);
   virtual ~GlueXSensitiveDetectorFTOF();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapFTOFbar* fBarHitsMap;
   GlueXHitsMapFTOFpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static int MAX_HITS_PER_BAR;
   static double ATTENUATION_LENGTH;
   static double C_EFFECTIVE;
   static double FULL_BAR_LENGTH;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;
   static double MAX_TOF;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
