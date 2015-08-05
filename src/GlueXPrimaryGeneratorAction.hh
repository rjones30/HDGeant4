//
// GlueXPrimaryGeneratorAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#ifndef _GLUEXPRIMARYGENERATORACTION_H_
#define _GLUEXPRIMARYGENERATORACTION_H_

#include "CobremsGenerator.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"
#include "GlueXParticleGun.hh"
#include "G4SystemOfUnits.hh"

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
   ~GlueXPrimaryGeneratorAction();
   
   virtual void GeneratePrimaries(G4Event* anEvent);
   void GeneratePrimariesHDDM(G4Event* anEvent);
   void GeneratePrimariesParticleGun(G4Event* anEvent);
   void GeneratePrimariesCobrems(G4Event* anEvent);
   void GenerateBeamPhoton(G4Event* anEvent, double t0);

   int ConvertGeant3ToPdg(int Geant3number) const;
   
private:
   source_type_t fSourceType;
   std::ifstream *fHDDMinfile;
   hddm_s::istream *fHDDMistream;
   CobremsGenerator *fCobremsGenerator;
   GlueXParticleGun *fParticleGun;
   G4ParticleTable *fParticleTable;
   int fGunParticleGeantType;
   int fGunParticlePDGType;
   G4ParticleDefinition *fGunParticle;
   G4ThreeVector fGunParticlePos;
   double fGunParticleMom;
   double fGunParticleMomTheta;
   double fGunParticleMomPhi;
   double fGunParticleDeltaPosR;
   double fGunParticleDeltaPosZ;
   double fGunParticleDeltaMom;
   double fGunParticleDeltaTheta;
   double fGunParticleDeltaPhi;
   double fBeamBucketPeriod;
   double fBeamBackgroundRate;
   double fBeamBackgroundGateStart;
   double fBeamBackgroundGateStop;
   double fL1triggerTimeSigma;
   double fBeamStartZ;
   int fEventCount;

   // The following parameters describe the dimensions of the target
   // that are used when generating the primary interaction vertex for
   // events from an external generator. An external event generator
   // knows nothing about the simulation geometry, so it makes sense
   // that this should be modeled in the simulation. They only apply
   // to the HDDM input source.  They are initialized to default values
   // for the GlueX liquid hydrogen target in the constructor, but
   // can be accessed/changed by the getter/setter methods below.
   double fTargetCenterZ;
   double fTargetLength;
   double fBeamDiameter;

 public:
   void setTargetCenterZ(double Z_cm) {
      fTargetCenterZ = Z_cm * cm;
   }
   void setTargetLength(double L_cm) {
      fTargetLength = L_cm * cm;
   }
   void setBeamDiameter(double D_cm) {
      fBeamDiameter = D_cm * cm;
   }
   double getTargetCenterZ() {
      return fTargetCenterZ / cm;
   }
   double getTargetLength() {
      return fTargetLength / cm;
   }
   double getBeamDiameter() {
      return fBeamDiameter / cm;
   }

   double getBeamBucketPeriod(int runno=0);

   int getEventCount() {
      return fEventCount;
   }

   void setBeamBucketPeriod(double period_ns) {
      fBeamBucketPeriod = period_ns * ns;
   }
   void setL1triggerTimeSigma(double sigma_ns) {
      fL1triggerTimeSigma = sigma_ns;
   }
   void setBeamStartZ(double Z_cm) {
      fBeamStartZ = Z_cm * cm;
   }
   double getL1triggerTimeSigma() {
      return fL1triggerTimeSigma;
   }
   double getBeamStartZcm() {
      return fBeamStartZ / cm;
   }

 private:
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

   ImportanceSampler fCoherentPDFx; 
   ImportanceSampler fIncoherentPDFlogx;
   ImportanceSampler fIncoherentPDFy;
   double fIncoherentPDFtheta02;

   void prepareCobremsImportanceSamplingPDFs();
};

#endif
