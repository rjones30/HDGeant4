//
// hdgeant4.cc : main program for hdgeant4 physics simulation
//
// author: richard.t.jones at uconn.edu
// version: august 15, 2013
//

#include <GlueXUserOptions.hh>
#include <GlueXDetectorConstruction.hh>
#include <GlueXPhysicsList.hh>
#include <GlueXPrimaryGeneratorAction.hh>
#include <GlueXRunAction.hh>
#include <GlueXEventAction.hh>
#include <GlueXSteppingAction.hh>
#include <GlueXSteppingVerbose.hh>

#include <DANA/DApplication.h>

#include <G4RunManager.hh>
#include <G4UImanager.hh>

#ifdef G4VIS_USE
#include <G4VisExecutive.hh>
#endif

#ifdef G4UI_USE
#include <G4UIExecutive.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>
#endif

int main(int argc,char** argv)
{
  // Initialize the jana framework
  DApplication dapp(argc, argv);
  dapp.Init();

  // Read user options from file control.in
  GlueXUserOptions opts;
  if (! opts.ReadControl_in("control.in")) {
    exit(3);
  }

  // User Verbose output class
  G4VSteppingVerbose* verbosity = new GlueXSteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // Run manager handles the rest of the initialization
  G4RunManager* runManager = new G4RunManager;

  // Geometry initialization
  GlueXDetectorConstruction* geometry = new GlueXDetectorConstruction;
  int Npara = geometry->GetParallelWorldCount();
  for (int para=1; para <= Npara; ++para) {
    G4String name = geometry->GetParallelWorldName(para);
    G4LogicalVolume *topvol = geometry->GetParallelWorldVolume(para);
    GlueXParallelWorld *parallelWorld = new GlueXParallelWorld(name,topvol);
    geometry->RegisterParallelWorld(parallelWorld);
  }
  runManager->SetUserInitialization(geometry);

  // Physics process initialization
  G4VUserPhysicsList* physics = new GlueXPhysicsList;
  runManager->SetUserInitialization(physics);
   
  // User actions initialization
  G4VUserPrimaryGeneratorAction* event_gen = 
                                 new GlueXPrimaryGeneratorAction(geometry);
  runManager->SetUserAction(event_gen);
  G4UserRunAction* run_action = new GlueXRunAction;
  runManager->SetUserAction(run_action);
  G4UserEventAction* event_action = new GlueXEventAction;
  runManager->SetUserAction(event_action);
  G4UserSteppingAction* stepping_action = new GlueXSteppingAction;
  runManager->SetUserAction(stepping_action);

  // Initialize G4 kernel
  runManager->Initialize();
      
  // Initialize graphics
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // Start the user interface
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  
  if (argc > 1) {   // batch mode  
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }
  else {           // interactive mode, define UI session
#ifdef G4UI_USE
 #ifdef G4UI_USE_EXECUTIVE
    G4UIExecutive* ui = new G4UIExecutive(argc,argv);
 #else
    G4UIterminal* ui = new G4UIterminal(new G4UItcsh);
 #endif
    ui->SessionStart();
 #ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute vis.mac");     
 #endif
    delete ui;
#endif
  }
 
  // Clean up and exit
#ifdef G4VIS_USE
  delete visManager;
#endif     

  delete runManager;
  delete verbosity;
  return 0;
}
