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
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"
#include "GlueXUserEventInformation.hh"
#include "G4Material.hh"
#include <sstream>

class MyPathFinder : public G4PathFinder {
 public:
   const G4Navigator *getNavigator(int world) { return GetNavigator(world); }
};

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
      G4ThreeVector lpos = ((MyPathFinder*)pathfinder)->getNavigator(world)->GetCurrentLocalCoordinate();
      G4ThreeVector gpos(lpos);
      ((MyPathFinder*)pathfinder)->getNavigator(world)->GetLocalToGlobalTransform().ApplyPointTransform(gpos);
      if (mat) {
std::stringstream msg;
msg << " GlueXPathFinder::GetLocatedVolume finds volume " << pvol->GetName()
    << " with material " << mat->GetName()
    << " in world " << world << " ("
    << ((MyPathFinder*)pathfinder)->getNavigator(world)->GetWorldVolume()->GetName()
    << ") at (" << gpos.x() << "," << gpos.y() << "," << gpos.z() << ")"
    << ", local=(" << lpos.x() << "," << lpos.y() << "," << lpos.z() << ")";
GlueXUserEventInformation::Dlog(msg.str());
         return pvol;
      }
      else {
std::stringstream msg;
msg << " GlueXPathFinder::GetLocatedVolume finds volume " << pvol->GetName()
    << " with material NULL"
    << " in world " << world << " ("
    << ((MyPathFinder*)pathfinder)->getNavigator(world)->GetWorldVolume()->GetName()
    << ") at (" << gpos.x() << "," << gpos.y() << "," << gpos.z() << ")"
    << ", local=(" << lpos.x() << "," << lpos.y() << "," << lpos.z() << ")";
GlueXUserEventInformation::Dlog(msg.str());
      }
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
      G4ThreeVector lpos = ((MyPathFinder*)pathfinder)->getNavigator(world)->GetCurrentLocalCoordinate();
      G4ThreeVector gpos(lpos);
      ((MyPathFinder*)pathfinder)->getNavigator(world)->GetLocalToGlobalTransform().ApplyPointTransform(gpos);
      G4TouchableHandle touch = pathfinder->CreateTouchableHandle(world);
      int copyno = (touch)? touch->GetCopyNumber() : 0;
      if (mat) {
std::stringstream msg;
msg << " GlueXPathFinder::CreateTouchableHandle finds volume " << pvol->GetName() << ":" << copyno
    << " with material " << mat->GetName()
    << " in world " << world << " ("
    << ((MyPathFinder*)pathfinder)->getNavigator(world)->GetWorldVolume()->GetName()
    << ") at (" << gpos.x() << "," << gpos.y() << "," << gpos.z() << ")"
    << ", local=(" << lpos.x() << "," << lpos.y() << "," << lpos.z() << ")";
GlueXUserEventInformation::Dlog(msg.str());
         return touch;
      }
      else {
std::stringstream msg;
msg << " GlueXPathFinder::CreateTouchableHandle finds volume " << pvol->GetName() << ":" << copyno
    << " with material NULL"
    << " in world " << world << " ("
    << ((MyPathFinder*)pathfinder)->getNavigator(world)->GetWorldVolume()->GetName()
    << ") at (" << gpos.x() << "," << gpos.y() << "," << gpos.z() << ")"
    << ", local=(" << lpos.x() << "," << lpos.y() << "," << lpos.z() << ")";
GlueXUserEventInformation::Dlog(msg.str());
      }
   }
   return G4TouchableHandle();
}
