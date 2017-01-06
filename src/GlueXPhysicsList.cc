//
// GlueXPhysicsList class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXPhysicsList.hh"

#include <G4ParallelWorldPhysics.hh>

#if USING_DIRACXX
G4ThreadLocal GlueXBeamConversionProcess *GlueXPhysicsList::fBeamConversion = 0;
#endif

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

GlueXPhysicsList::~GlueXPhysicsList()
{}

void GlueXPhysicsList::ConstructProcess()
{
   QGSP_FTFP_BERT::ConstructProcess();

#if USING_DIRACXX
   // Attach the TPOL beam pair conversion process to the photon
   fBeamConversion = new GlueXBeamConversionProcess("TPOL beam conversion");
   G4ParticleDefinition *gamma = G4Gamma::GammaDefinition();
   G4ProcessManager *mgr = gamma->GetProcessManager();
   if (mgr != 0) {
      mgr->AddDiscreteProcess(fBeamConversion);
   }
   else {
      G4cerr << "Error in GlueXPhotonBeamGenerator::GeneratePrimaryVertex"
             << " - gamma process manager is null, cannot continue."
             << G4endl;
      exit(1);
   }
#endif
}
