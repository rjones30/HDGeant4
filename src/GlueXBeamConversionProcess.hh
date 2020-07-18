//
// GlueXBeamConversionProcess class header
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.
// Resources are created once when the first object is instantiated,
// and destroyed once when the last object is destroyed.

#ifndef _GLUEXBEAMCONVERSIONPROCESS_H_
#define _GLUEXBEAMCONVERSIONPROCESS_H_

#include <G4VDiscreteProcess.hh>
#include <G4VParticleChange.hh>
#include <ImportanceSampler.hh>
#include <PairConversionGeneration.hh>
#include <AdaptiveSampler.hh>
#include <G4AutoLock.hh>
#include <G4Step.hh>

class GlueXBeamConversionProcess: public G4VDiscreteProcess
{
 public:
   GlueXBeamConversionProcess(const G4String &name, 
                              G4ProcessType aType=fGeneral);
   GlueXBeamConversionProcess(GlueXBeamConversionProcess &src);
   virtual ~GlueXBeamConversionProcess();

   virtual G4double PostStepGetPhysicalInteractionLength(const G4Track &track,
                                                  G4double previousStepSize,
                                                  G4ForceCondition *condition);
   virtual G4VParticleChange *PostStepDoIt(const G4Track &track, 
                                           const G4Step &step);

 protected:
   void GenerateBeamPairConversion(const G4Step &step);
   virtual G4double GetMeanFreePath(const G4Track &track, 
                                    G4double previousStepSize,
                                    G4ForceCondition *condition);
   void GenerateBetheHeitlerProcess(const G4Step &step);

#ifdef USING_DIRACXX
   PairConversionGeneration *fPairsGeneration;
#endif

   static int fStopBeamBeforeConverter;
   static int fStopBeamAfterConverter;
   static int fStopBeamAfterTarget;
   static G4double fBHpair_mass_min;

   void prepareImportanceSamplingPDFs();

   ImportanceSampler *fPaircohPDF;
   ImportanceSampler *fTripletPDF;


   AdaptiveSampler *fAdaptiveSampler;
   void prepareAdaptiveSampler();

 private:
   GlueXBeamConversionProcess operator=(GlueXBeamConversionProcess &src);

   static G4Mutex fMutex;
   static int fInitialized;
   static std::vector<AdaptiveSampler*> fAdaptiveSamplerRegistry;
};

#endif
