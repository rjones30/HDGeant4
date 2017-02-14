//
// GlueXDetectorMessenger class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXDetectorMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "GlueXDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4GeometryManager.hh"

GlueXDetectorMessenger::GlueXDetectorMessenger(GlueXDetectorConstruction* myDet)
:myDetector(myDet)
{ 
  hdgeant4Dir = new G4UIdirectory("/hdgeant4/");
  hdgeant4Dir->SetGuidance("UI commands specific to hdgeant4");
  
  detDir = new G4UIdirectory("/hdgeant4/det/");
  detDir->SetGuidance("detector control");
  
  FieldCmd = new G4UIcmdWithADoubleAndUnit("/hdgeant4/det/setField",this);
  FieldCmd->SetGuidance("Define magnetic field.");
  FieldCmd->SetGuidance("Magnetic field will be in z direction.");
  FieldCmd->SetParameterName("Bz",false);
  FieldCmd->SetUnitCategory("Magnetic flux density");
  FieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
        
  StepMaxCmd = new G4UIcmdWithADoubleAndUnit("/hdgeant4/det/stepMax",this);
  StepMaxCmd->SetGuidance("Define a step max");
  StepMaxCmd->SetParameterName("stepMax",false);
  StepMaxCmd->SetUnitCategory("Length");
  StepMaxCmd->AvailableForStates(G4State_Idle);

  OpenGeomCmd = new G4UIcmdWithoutParameter("/hdgeant4/openGeometry",this);
  OpenGeomCmd->SetGuidance("Re-open the geometry.");
  OpenGeomCmd->SetGuidance("This command is only needed "
                           "if the geometry has already been closed.");
  OpenGeomCmd->AvailableForStates(G4State_Idle);

  CloseGeomCmd = new G4UIcmdWithoutParameter("/hdgeant4/closeGeometry",this);
  CloseGeomCmd->SetGuidance("Explicitly close the geometry.");
  CloseGeomCmd->AvailableForStates(G4State_Idle);
}

GlueXDetectorMessenger::~GlueXDetectorMessenger()
{
  delete FieldCmd;
  delete StepMaxCmd;
  delete detDir;
  delete hdgeant4Dir;
}

void GlueXDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == FieldCmd) {
    myDetector->SetUniformField(FieldCmd->GetNewDoubleValue(newValue)/tesla);
  }
  else if (command == StepMaxCmd) {
    myDetector->SetMaxStep(StepMaxCmd->GetNewDoubleValue(newValue));
  }
  else if (command == OpenGeomCmd) {
    G4GeometryManager::GetInstance()->OpenGeometry();
  }
  else if (command == CloseGeomCmd) {
    G4GeometryManager::GetInstance()->CloseGeometry();
  }
}
