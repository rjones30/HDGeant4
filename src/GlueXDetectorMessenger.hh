//
// GlueXDetectorMessenger class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#ifndef GlueXDetectorMessenger_h
#define GlueXDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class GlueXDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class GlueXDetectorMessenger: public G4UImessenger
{
  public:
    GlueXDetectorMessenger(GlueXDetectorConstruction*);
    GlueXDetectorMessenger(const GlueXDetectorMessenger& src);
    GlueXDetectorMessenger &operator=(const GlueXDetectorMessenger& src);
   ~GlueXDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    GlueXDetectorConstruction* myDetector;
    
    G4UIdirectory*             hdgeant4Dir;
    G4UIdirectory*             detDir;
    G4UIcmdWithADoubleAndUnit* FieldCmd;
    G4UIcmdWithADoubleAndUnit* StepMaxCmd;    
    G4UIcmdWithoutParameter*   OpenGeomCmd;
    G4UIcmdWithoutParameter*   CloseGeomCmd;
};

#endif
