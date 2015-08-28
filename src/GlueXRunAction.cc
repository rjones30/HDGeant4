//
// GlueXRunAction class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXRunAction.hh"

#include "G4Run.hh"

#include "G4MTRunManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVDivision.hh"

#include "assert.h"

G4Mutex GlueXRunAction::fMutex = G4MUTEX_INITIALIZER;

void GlueXRunAction::BeginOfRunAction(const G4Run* aRun)
{
   G4AutoLock barrier(&fMutex);
   G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

   // This section was added to deal with a specific bug related to
   // GPVDivision volumes that need to have a thread-local instance
   // of the rotation matrix. As soon as the Geant4 team releases a
   // proper fix for this bug, this section can go away.

   G4RunManager::RMType rmtype = G4RunManager::GetRunManager()->
                                               GetRunManagerType();
   if (rmtype == G4RunManager::workerRM) {
      G4PhysicalVolumeStore* const physVolStore = 
                             G4PhysicalVolumeStore::GetInstance();
      assert (physVolStore != NULL);
      for (G4PhysicalVolumeStore::const_iterator iter = physVolStore->begin();
           iter != physVolStore->end();
           ++iter)
      {
         G4VPhysicalVolume *pvol = *iter;
         if (dynamic_cast<G4PVDivision*>(pvol)) {
            G4RotationMatrix *masterRmat = pvol->GetRotation();
            G4RotationMatrix *workerRmat;
            if (masterRmat)
               workerRmat = new G4RotationMatrix(*masterRmat);
            else
               workerRmat = new G4RotationMatrix();
            pvol->SetRotation(workerRmat);
            usedRmatrix.push_back(workerRmat);
         }
      }
   }
}

void GlueXRunAction::EndOfRunAction(const G4Run*)
{
   // Clean up the thread-local clones of G4RotationMatrix that were
   // created in the BeginOfRunAction. As soon as the Geant4 team 
   // releases a proper fix for this bug, this section can go away.
 
   std::list<G4RotationMatrix*>::iterator iter;
   for (iter = usedRmatrix.begin(); iter != usedRmatrix.end(); ++iter)
      delete *iter;
   usedRmatrix.clear();

}
