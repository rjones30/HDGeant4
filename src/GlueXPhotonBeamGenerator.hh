//
// GlueXPhotonBeamGenerator class header
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state. It is
// invoked from the GlueXPhotonBeamGeneratorAction class, which is
// responsible for managing the interlocks to ensure thread safety.

#ifndef GlueXPhotonBeamGenerator_H
#define GlueXPhotonBeamGenerator_H

#include <G4VPrimaryGenerator.hh>
#include <G4GenericMessenger.hh>
#include <CobremsGeneration.hh>
#include <ImportanceSampler.hh>
#include <GlueXPseudoDetectorTAG.hh>
#include <G4Event.hh>

class GlueXPhotonBeamGenerator: public G4VPrimaryGenerator
{
 public:
   GlueXPhotonBeamGenerator(CobremsGeneration *gen);
   virtual ~GlueXPhotonBeamGenerator();

   virtual void GeneratePrimaryVertex(G4Event *event);
   virtual void GenerateBeamPhoton(G4Event *event, double t0);

   static double GenerateTriggerTime(const G4Event *event);
   static int GenerateTaggerHit(const G4Event *event, 
                                double energy, double time, int bg=0);
   static void GenerateRFsync(const G4Event *event);

 protected:
   CobremsGeneration *fCobrems;
   static GlueXPseudoDetectorTAG *fTagger;
   static int fGenerateNotSimulate;
   static int fBeamBackgroundTagOnly;

   static double fBeamBucketPeriod;
   static double fBeamStartZ;
   static double fBeamDiameter;
   static double fBeamVelocity;
   static double fBeamOffset[2];

   static ImportanceSampler fCoherentPDFx; 
   static ImportanceSampler fIncoherentPDFlogx;
   static ImportanceSampler fIncoherentPDFy;
   static double fIncoherentPDFtheta02;

   void prepareImportanceSamplingPDFs();

   static int fForceFixedPolarization;
   static double fFixedPolarization;
   static double fFixedPolarization_phi;

   G4GenericMessenger *fMessenger;

 public:
   static void setBeamDiameter(double D) {
      fBeamDiameter = D;
   }
   static double getBeamDiameter() {
      return fBeamDiameter;
   }
   static double getBeamVelocity() {
      return fBeamVelocity;
   }
   static void setBeamOffset(double x, double y) {
      fBeamOffset[0] = x;
      fBeamOffset[1] = y;
   }
   static double getBeamOffset(int i) {
      return fBeamOffset[i];
   }
   static double getRFreferencePlaneZ(int runno=0);
   static double getBeamBucketPeriod(int runno=0);
   static void setBeamBucketPeriod(double period) {
      fBeamBucketPeriod = period;
   }
   static void setBeamStartZ(double z) {
      fBeamStartZ = z;
   }
   static double getBeamStartZ() {
      return fBeamStartZ;
   }
   void enableFixedPolarization(double polar, double phi_deg) {
      fFixedPolarization = polar;
      fFixedPolarization_phi = phi_deg * M_PI/180;
      fForceFixedPolarization = true;
   }
   void disableFixedPolarization() {
      fForceFixedPolarization = false;
   }

 private:
   GlueXPhotonBeamGenerator(const GlueXPhotonBeamGenerator &src) {}
   GlueXPhotonBeamGenerator &operator=(const GlueXPhotonBeamGenerator &src) {
      return *this;
   }
};

#endif
