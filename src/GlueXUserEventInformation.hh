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
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <HDDM/hddm_s.hpp>

class GlueXUserEventInformation: public G4VUserEventInformation
{
 public:
   GlueXUserEventInformation(hddm_s::HDDM *hddmevent=NULL);
   GlueXUserEventInformation(int geanttype, G4ThreeVector &pos, 
                                            G4ThreeVector &mom);
   ~GlueXUserEventInformation();

   void SetRandomSeeds();
   void Print() const;

   hddm_s::HDDM *getOutputRecord() {
      return fOutputRecord;
   }

 protected:
   hddm_s::HDDM *fOutputRecord;
   bool fKeepEvent;
   int fNprimaries;
};

#endif // _GLUEXUSEREVENTINFORMATION_
