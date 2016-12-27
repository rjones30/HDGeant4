//
// GlueXPhysicsList class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. The
// Geant4 runtime makes sure that its constructor runs in the master
// thread, and then its Construct methods are run again each time
// a new worker thread is initialized.

#ifndef GlueXPhysicsList_h
#define GlueXPhysicsList_h 1

#include <GlueXDetectorConstruction.hh>
#include <GlueXBeamConversionProcess.hh>

#include <G4VUserPhysicsList.hh>
#include <QGSP_FTFP_BERT.hh>
#include "globals.hh"

class GlueXPhysicsList: public QGSP_FTFP_BERT
{
 public:
   GlueXPhysicsList(const GlueXDetectorConstruction *geometry=0);
   virtual ~GlueXPhysicsList();

   virtual void ConstructProcess();

 protected:
#if USING_DIRACXX
   static G4ThreadLocal GlueXBeamConversionProcess *fBeamConversion;
#endif
};

#endif
