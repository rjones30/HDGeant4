//
// hdgeant4.cc : main program for hdgeant4 physics simulation
//
// author: richard.t.jones at uconn.edu
// version: august 15, 2013
//

#include <GlueXUserOptions.hh>
#include <GlueXDetectorConstruction.hh>
#include <GlueXUserActionInitialization.hh>
#include <GlueXPhysicsList.hh>

#include <DANA/DApplication.h>
#include <unistd.h>

#include <G4MTRunManager.hh>
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

void usage()
{
  std::cout << std::endl
            << "Usage: hdgeant4 [-v] [-t <THREAD_COUNT>] [<filename.mac>]"
            << std::endl;
  exit(9);
}

int main(int argc,char** argv)
{
  // Initialize the jana framework
  DApplication dapp(argc, argv);
  dapp.Init();

  // Interpret special command-line arguments
  int use_visualization = 0;
  int worker_threads = 1;
  int argi=1;
  int c;
  while ((c = getopt(argc, argv, "vt:")) != -1) {
    if (c == 'v') {
      use_visualization = 1;
      argi += 1;
    }
    else if (c == 't') {
      worker_threads = atoi(optarg);
      argi += 2;
    }
    else {
      usage();
    }
  }

  // Read user options from file control.in
  GlueXUserOptions opts;
  if (! opts.ReadControl_in("control.in")) {
    exit(3);
  }

  // Declare our G4VSteppingVerbose implementation
  G4VSteppingVerbose::SetInstance(new GlueXSteppingVerbose());

  // Run manager handles the rest of the initialization
#ifdef G4MULTITHREADED
  G4MTRunManager runManager;
  runManager.SetNumberOfThreads(worker_threads);
#else
  G4RunManager runManager;
#endif

  // Geometry initialization
  GlueXDetectorConstruction *geometry = new GlueXDetectorConstruction();
  int Npara = geometry->GetParallelWorldCount();
  for (int para=1; para <= Npara; ++para) {
    G4String name = geometry->GetParallelWorldName(para);
    G4LogicalVolume *topvol = geometry->GetParallelWorldVolume(para);
    GlueXParallelWorld *parallelWorld = new GlueXParallelWorld(name,topvol);
    geometry->RegisterParallelWorld(parallelWorld);
  }
  runManager.SetUserInitialization(geometry);

  // Physics process initialization
  runManager.SetUserInitialization(new GlueXPhysicsList());
   
  // User actions initialization
  runManager.SetUserInitialization(new GlueXUserActionInitialization());

  // Initialize G4 kernel
  runManager.Initialize();
      
  // Initialize graphics (option -v)
  G4VisManager* visManager = 0;
  if (use_visualization) {
#ifdef G4VIS_USE
    visManager = new G4VisExecutive();
    visManager->Initialize();
#else
    std::cerr << "Visualization system not available,"
              << " please rebuild hdgeant4 with visualization enabled."
              << std::endl;
     exit(1);
#endif
  }

  // Start the user interface
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  
  if (argc > argi) {   // batch mode  
    G4String command = "/control/execute ";
    G4String fileName = argv[argi];
    UImanager->ApplyCommand(command + fileName);
  }
  else if (visManager) {    // interactive mode with visualization
#ifdef G4UI_USE_EXECUTIVE
    G4UIExecutive UIexec(argc, argv, "qt");
    UImanager->ApplyCommand("/control/execute vis.mac");     
    UIexec.SessionStart();
#else
    G4UIterminal UIterm(new G4UItcsh());
    UImanager->ApplyCommand("/control/execute vis.mac");     
    UIterm.SessionStart();
#endif
  }
  else {    // interactive mode without visualization
    G4UIterminal UIterm(new G4UItcsh());
    UIterm.SessionStart();
  }
 
  // Clean up and exit
  if (visManager)
    delete visManager;
  return 0;
}
