//
// GlueXPathFinder class implementation
//
// author: richard.t.jones at uconn.edu
// version: aug 12, 2015
//
// This is a wrapper around the standard G4PathFinder
// singleton class that adds some functionality needed
// for efficient lookup of volume information in a
// multi-layer mass geometry.
//

#include "GlueXPathFinder.hh"
#include "G4TransportationManager.hh"

G4VPhysicalVolume* GlueXPathFinder::GetLocatedVolume()
{
   // Search through the parallel worlds from highest to lowest
   // and return the first physical volume found that has a defined
   // material, or a null pointer if not found.

   G4PathFinder *pathfinder = G4PathFinder::GetInstance();
   int Nworlds = G4TransportationManager::GetTransportationManager()
                                          ->GetNoWorlds();
   for (int world = Nworlds - 1; world >= 0; --world) {
      G4VPhysicalVolume *pvol = pathfinder->GetLocatedVolume(world);
      G4LogicalVolume *lvol = (pvol)? pvol->GetLogicalVolume() : 0;
      G4Material *mat = (lvol)? lvol->GetMaterial() : 0;
      if (mat)
         return pvol;
   }
   return 0;
}

G4TouchableHandle GlueXPathFinder::CreateTouchableHandle()
{
   // Search through the parallel worlds from highest to lowest
   // and return the first touchable found that has a defined
   // material, or a null pointer if not found.

   G4PathFinder *pathfinder = G4PathFinder::GetInstance();
   int Nworlds = G4TransportationManager::GetTransportationManager()
                                          ->GetNoWorlds();
   for (int world = Nworlds - 1; world >= 0; --world) {
      G4VPhysicalVolume *pvol = pathfinder->GetLocatedVolume(world);
      G4LogicalVolume *lvol = (pvol)? pvol->GetLogicalVolume() : 0;
      G4Material *mat = (lvol)? lvol->GetMaterial() : 0;
      if (mat)
         return pathfinder->CreateTouchableHandle(world);
   }
   return G4TouchableHandle();
}
