//
// MyPrimaryPionZapper class
//
// Provides G4WrapperProcsess around G4Decay process to enable
// user-defined zapping zone in the detector where primary
// pions are forced to decay.
//

#include "MyPrimaryPionZapper.hh"
#include "GlueXPathFinder.hh"

G4double MyPrimaryPionZapper::PostStepGetPhysicalInteractionLength(
                              const G4Track &track,
                              G4double previousStepSize,
                              G4ForceCondition *condition)
{
   G4TouchableHandle touch = GlueXPathFinder::CreateTouchableHandle();
   G4VPhysicalVolume *pvol = (touch)? touch->GetVolume() : 0;
   if (track.GetTrackID() < 3 && pvol && pvol->GetName() == "BCL0") {
      return 1e-6;
   }
   return pRegProcess->PostStepGetPhysicalInteractionLength(
                       track, previousStepSize, condition);
}
