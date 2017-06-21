//
// GlueXSensitiveDetectorDIRC - class header
//
// author: richard.t.jones at uconn.edu
// version: november 28, 2016
//          may 7, 2017 R. Dzhygadlo - added adjustments for
//                      GlueXHitDIRCBar GlueXHitDIRCPmt    
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorDIRC_h
#define GlueXSensitiveDetectorDIRC_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitDIRCflash.hh"
#include "GlueXHitDIRCpoint.hh"
#include "GlueXHitDIRCPmt.hh"
#include "GlueXHitDIRCBar.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorDIRC : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorDIRC(const G4String& name);
   GlueXSensitiveDetectorDIRC(const GlueXSensitiveDetectorDIRC &right);
   GlueXSensitiveDetectorDIRC &operator=(const GlueXSensitiveDetectorDIRC &right);
   virtual ~GlueXSensitiveDetectorDIRC();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);
  
 private:
  std::vector<GlueXHitDIRCBar> fHitsBar;
  std::vector<GlueXHitDIRCPmt> fHitsPmt;
  
  static std::map<G4LogicalVolume*, int> fVolumeTable;
  
  static int MAX_HITS;
  // put all other detector response parameters here
  static double TWO_HIT_TIME_RESOL;
  
  static int instanceCount;
  static G4Mutex fMutex;
};

#endif
