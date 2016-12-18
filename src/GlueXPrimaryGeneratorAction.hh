//
// GlueXPrimaryGeneratorAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread,
// but virtually all of its functions need to be serialized, so
// it maintains its own interlocks for this purpose. Resources
// are created once when the first object is instantiated, and
// destroyed once when the last object is destroyed.

#ifndef _GLUEXPRIMARYGENERATORACTION_H_
#define _GLUEXPRIMARYGENERATORACTION_H_

#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "CobremsGenerator.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"
#include "GlueXParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"

#include "globals.hh"

#include <HDDM/hddm_s.hpp>

#include <fstream>

class G4Event;

class GlueXPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
   
   enum source_type_t {
      SOURCE_TYPE_NONE,
      SOURCE_TYPE_PARTICLE_GUN,
      SOURCE_TYPE_COBREMS_GEN,
      SOURCE_TYPE_HDDM
   };

   GlueXPrimaryGeneratorAction();
   GlueXPrimaryGeneratorAction(const GlueXPrimaryGeneratorAction &src);
   GlueXPrimaryGeneratorAction &operator=(const GlueXPrimaryGeneratorAction &src);
   ~GlueXPrimaryGeneratorAction();
   
   virtual void GeneratePrimaries(G4Event* anEvent);
   void GeneratePrimariesHDDM(G4Event* anEvent);
   void GeneratePrimariesParticleGun(G4Event* anEvent);
   void GeneratePrimariesCobrems(G4Event* anEvent);
   void GenerateBeamPhoton(G4Event* anEvent, double t0);
   void GenerateBeamPairConversion(G4Step* step);

   static int ConvertGeant3ToPdg(int Geant3Type);
   static int ConvertPdgToGeant3(int PDGtype);
   static double GetMassPDG(int PDGtype);
   static double GetMass(int Geant3Type);
   static double GenerateTriggerTime();
 
 private:
   static int instanceCount;
   static source_type_t fSourceType;
   static std::ifstream *fHDDMinfile;
   static hddm_s::istream *fHDDMistream;
   static CobremsGenerator *fCobremsGenerator;
   static G4ParticleTable *fParticleTable;
   static GlueXParticleGun *fParticleGun;

 public:
   struct single_particle_gun_t {
      int geantType;
      int pdgType;
      G4ParticleDefinition *partDef;
      G4ThreeVector pos;
      double mom;
      double theta;
      double phi;
      double deltaR;
      double deltaZ;
      double deltaMom;
      double deltaTheta;
      double deltaPhi;
   };

 private:
   static single_particle_gun_t fGunParticle;

   static double fBeamBucketPeriod;
   static double fBeamBackgroundRate;
   static double fBeamBackgroundGateStart;
   static double fBeamBackgroundGateStop;
   static double fL1triggerTimeSigma;
   static double fBeamStartZ;
   static double fBeamVelocity;

   static int fEventCount;

   // The following parameters describe the dimensions of the target
   // that are used when generating the primary interaction vertex for
   // events from an external generator. An external event generator
   // knows nothing about the simulation geometry, so it makes sense
   // that this should be modeled in the simulation. They only apply
   // to the HDDM input source.  They are initialized to default values
   // for the GlueX liquid hydrogen target in the constructor, but
   // can be accessed/changed by the getter/setter methods below.
   static double fTargetCenterZ;
   static double fTargetLength;
   static double fBeamDiameter;

 public:
   static void setTargetCenterZ(double Z) {
      fTargetCenterZ = Z;
   }
   static void setTargetLength(double L) {
      fTargetLength = L;
   }
   static void setBeamDiameter(double D) {
      fBeamDiameter = D;
   }
   static double getTargetCenterZ() {
      return fTargetCenterZ;
   }
   static double getTargetLength() {
      return fTargetLength;
   }
   static double getBeamDiameter() {
      return fBeamDiameter;
   }
   static double getBeamVelocity() {
      return fBeamVelocity;
   }

   static double getBeamBucketPeriod(int runno=0);

   int getEventCount() {
      return fEventCount;
   }

   static void setBeamBucketPeriod(double period) {
      fBeamBucketPeriod = period;
   }
   static void setL1triggerTimeSigma(double sigma) {
      fL1triggerTimeSigma = sigma;
   }
   static void setBeamStartZ(double z) {
      fBeamStartZ = z;
   }
   static double getL1triggerTimeSigma() {
      return fL1triggerTimeSigma;
   }
   static double getBeamStartZ() {
      return fBeamStartZ;
   }

   // The following tables contain PDFs for importance-sampling the
   // kinematic variables in coherent bremsstrahlung beam generation.
 
   struct ImportanceSampler {
      std::vector<double> randvar;
      std::vector<double> density;
      std::vector<double> integral;
      double Psum;
      double Pcut;
      double Pmax;
      int Nfailed;
      int Npassed;

      ImportanceSampler()
       : Psum(0), Pcut(1), Pmax(0), Nfailed(0), Npassed(0) {}
   };

 private:
   static ImportanceSampler fCoherentPDFx; 
   static ImportanceSampler fIncoherentPDFlogx;
   static ImportanceSampler fIncoherentPDFy;

   static double fIncoherentPDFtheta02;

   void prepareCobremsImportanceSamplingPDFs();

 private:
   static G4Mutex fMutex;
};

#endif
