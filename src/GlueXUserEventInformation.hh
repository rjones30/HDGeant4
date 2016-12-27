//
// GlueXUserEventInformation class header
//
// author: richard.t.jones at uconn.edu
// version: aug 1, 2015
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.

#ifndef _GLUEXUSEREVENTINFORMATION_
#define _GLUEXUSEREVENTINFORMATION_

#include <iostream>
#include <fstream>
#include <vector>

#include "G4UImanager.hh"
#include "G4VUserEventInformation.hh"
#include "G4PrimaryVertex.hh"
#include "G4ThreeVector.hh"
#include <G4Track.hh>
#include "Randomize.hh"

#include <HDDM/hddm_s.hpp>

class GlueXUserEventInformation: public G4VUserEventInformation
{
 public:
   GlueXUserEventInformation(hddm_s::HDDM *hddmevent=NULL);
   ~GlueXUserEventInformation();

   void AddBeamParticle(int geanttype, double t0, const G4ThreeVector &pos, 
                                                  const G4ThreeVector &mom);
   void AddBeamParticle(int geanttype, double t0, const G4ThreeVector &pos, 
                                                  const G4ThreeVector &mom,
                                                  const G4ThreeVector &pol);
   void AddTargetParticle(int geanttype, double t0, const G4ThreeVector &pos, 
                                                    const G4ThreeVector &mom);
   void AddTargetParticle(int geanttype, double t0, const G4ThreeVector &pos, 
                                                    const G4ThreeVector &mom,
                                                    const G4ThreeVector &pol);
   void AddPrimaryVertex(const G4PrimaryVertex &vertex);
   void AddSecondaryVertex(const std::vector<G4Track*> &secondaries, 
                           int parentID);

   void SetRandomSeeds();
   void Print() const;

   hddm_s::HDDM *getOutputRecord() {
      return fOutputRecord;
   }

 protected:
   hddm_s::HDDM *fOutputRecord;
   bool fKeepEvent;
   int fNprimaries;

 private:
   GlueXUserEventInformation(const GlueXUserEventInformation &src);
   GlueXUserEventInformation &operator=(const GlueXUserEventInformation &src);
};

#endif // _GLUEXUSEREVENTINFORMATION_
