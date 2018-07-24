//
// GlueXSensitiveDetectorSTC - class header
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorSTC_h
#define GlueXSensitiveDetectorSTC_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitSTCpaddle.hh"
#include "GlueXHitSTCpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorSTC : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorSTC(const G4String& name);
   GlueXSensitiveDetectorSTC(const GlueXSensitiveDetectorSTC &right);
   GlueXSensitiveDetectorSTC &operator=(const GlueXSensitiveDetectorSTC &right);
   virtual ~GlueXSensitiveDetectorSTC();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapSTCpaddle* fHitsMap;
   GlueXHitsMapSTCpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static double ATTENUATION_LENGTH;
   static double C_EFFECTIVE;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_MEV;
   static double LIGHT_GUIDE;
   static double ANGLE_COR;
   static double BENT_REGION;
   static double STRAIGHT_LENGTH;
   static double BEND_LENGTH;
   static double NOSE_LENGTH;
   const static int NCHANNELS = 30;
   static double STRAIGHT_ATTENUATION_A[NCHANNELS];
   static double STRAIGHT_ATTENUATION_B[NCHANNELS];
   static double STRAIGHT_ATTENUATION_C[NCHANNELS];
   static double BENDNOSE_ATTENUATION_A[NCHANNELS];
   static double BENDNOSE_ATTENUATION_B[NCHANNELS];
   static double BENDNOSE_ATTENUATION_C[NCHANNELS];
   static double STRAIGHT_PROPAGATION_A[NCHANNELS];
   static double STRAIGHT_PROPAGATION_B[NCHANNELS];
   static double BEND_PROPAGATION_A[NCHANNELS];
   static double BEND_PROPAGATION_B[NCHANNELS];
   static double NOSE_PROPAGATION_A[NCHANNELS];
   static double NOSE_PROPAGATION_B[NCHANNELS];

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
