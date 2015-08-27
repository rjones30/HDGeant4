//
// GlueXPhysicsList class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "globals.hh"
#include "GlueXPhysicsList.hh"

#include <G4TransportationManager.hh>
#include <G4ParallelWorldPhysics.hh>

GlueXPhysicsList::GlueXPhysicsList(const GlueXDetectorConstruction *geometry)
 : QGSP_FTFP_BERT()
{
  if (geometry == 0) {
    geometry = GlueXDetectorConstruction::GetInstance();
    if (geometry == 0) {
      G4cerr << "GlueXPhysicsList constructor error - "
             << "cannot construct GlueXPhysicsList until the detector "
             << "geometry has been constructed, aborting."
             << G4endl;
      exit(1);
    }
  }
  for (int para=1; para <= geometry->GetParallelWorldCount(); ++para) {
    G4String name = geometry->GetParallelWorldName(para);
    RegisterPhysics(new G4ParallelWorldPhysics(name, true));
  }
}
