//
// GlueXPhysicsList class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "globals.hh"
#include "GlueXPhysicsList.hh"

#include <G4TransportationManager.hh>
#include <G4ParallelWorldPhysics.hh>

GlueXPhysicsList::GlueXPhysicsList(GlueXDetectorConstruction *geometry)
 : QGSP_FTFP_BERT()
{
  for (int para=1; para <= geometry->GetParallelWorldCount(); ++para) {
    G4String name = geometry->GetParallelWorldName(para);
    RegisterPhysics(new G4ParallelWorldPhysics(name, true));
  }
}
