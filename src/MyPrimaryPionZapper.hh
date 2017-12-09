//
// MyPrimaryPionZapper class header
//
// Provides G4WrapperProcsess around G4Decay process to enable
// user-defined zapping zone in the detector where primary
// pions are forced to decay.
//

#include <G4WrapperProcess.hh>

class MyPrimaryPionZapper : public G4WrapperProcess {
   virtual G4double PostStepGetPhysicalInteractionLength(const G4Track &track,
                    G4double previousStepSize, G4ForceCondition *condition);
};
