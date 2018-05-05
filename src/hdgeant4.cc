//
// hdgeant4.cc : main program for hdgeant4 physics simulation
//
// author: richard.t.jones at uconn.edu
// version: august 15, 2013
//

#include <GlueXUserOptions.hh>
#include <GlueXDetectorConstruction.hh>
#include <GlueXUserActionInitialization.hh>
#include <GlueXPrimaryGeneratorAction.hh>
#include <GlueXPhysicsList.hh>
#include <GlueXTimer.hh>
#include <HddmOutput.hh>
#include <Randomize.hh>

#include <DANA/DApplication.h>
#include <unistd.h>

#include <G4MTRunManager.hh>
#include <G4RunManager.hh>
#include <G4UImanager.hh>
#include <G4Timer.hh>

#ifdef G4VIS_USE
#include <G4VisExecutive.hh>
#endif

#ifdef G4UI_USE
#include <G4UIExecutive.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>
#endif

int run_number = 0;

void usage()
{
   G4cout << G4endl
          << "Usage: hdgeant4 [options] [<batch.mac>]" << G4endl
          << " where options include:" << G4endl
          << "    -v : open a graphics window for visualization" << G4endl
          << "    -tN : start N worker threads, default 1" << G4endl
          << "    -rN : set run to N, default taken from control.in" << G4endl
          << G4endl;
   exit(9);
}

void OpenGLXpreload();

int main(int argc,char** argv)
{
   // Initialize the jana framework
   DApplication dapp(argc, argv);
   dapp.create_event_buffer_thread = false;
   dapp.Init();

   // Interpret special command-line arguments
   int use_visualization = 0;
   int worker_threads = 1;
   int c;
   while ((c = getopt(argc, argv, "vt:r:")) != -1) {
      if (c == 'v') {
         use_visualization = 1;
      }
      else if (c == 't') {
         worker_threads = atoi(optarg);
      }
      else if (c == 'r') {
         run_number = atoi(optarg);
      }
      else {
         usage();
      }
   }

   // Initialize the graphics subsystem, if requested
   if (use_visualization) {
#ifdef G4VIS_USE
      OpenGLXpreload();
#endif
   }

   // Read user options from file control.in
   GlueXUserOptions opts;
   if (! opts.ReadControl_in("control.in")) {
      exit(3);
   }
   if (run_number == 0) {
      std::map<int, int> runno_opts;
      if (opts.Find("RUNNO", runno_opts) || opts.Find("RUNG", runno_opts)) {
         run_number = runno_opts[1];
      }
      else {
         G4cerr << "Warning - "
                << "no run number specified in control.in, "
                << "default value of 0 assumed." << G4endl;
         run_number = 0;
      }
   }

   G4Timer runtimer;
   G4Timer simtimer;
   runtimer.Start();

   HddmOutput *hddmOut = 0;
   std::map<int, std::string> outfile_opts;
   if (opts.Find("OUTFILE", outfile_opts)) {
      hddmOut = new HddmOutput(outfile_opts[1]);
      hddmOut->setRunNo(run_number);
   }

   G4Random::setTheEngine(new CLHEP::RanecuEngine);
   std::map<int, int> rndm_opts;
   if (opts.Find("RNDM", rndm_opts)) {
      if (rndm_opts.size() > 1) {
         long int seed[2];
         seed[0] = rndm_opts[1];
         seed[1] = rndm_opts[2];
         G4Random::setTheSeeds(seed);
      }
      else if (rndm_opts.size() == 1) {
         long int seed[2];
         G4Random::getTheTableSeeds(seed, rndm_opts[1]);
         G4Random::setTheSeeds(seed);
      }
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

	// Let user turn off geometry optimization for faster startup
	// (and slower running)
	std::map<int, int> geomopt;
	if ( opts.Find("GEOMOPT", geomopt) ) {
         if( geomopt[1] == 0 ) runManager.SetGeometryToBeOptimized( false );
	}

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
   GlueXPhysicsList *physicslist = new GlueXPhysicsList();
   runManager.SetUserInitialization(physicslist);

   // User actions initialization
   GlueXUserActionInitialization *userinit;
   userinit = new GlueXUserActionInitialization(physicslist);
   runManager.SetUserInitialization(userinit);

   // Initialize G4 kernel
   std::cout << "Initializing the Geant4 kernel..." << std::endl;
   runManager.Initialize();
   
   // Initialize graphics (option -v)
   G4VisManager* visManager = 0;
   if (use_visualization) {
#ifdef G4VIS_USE
      visManager = new G4VisExecutive();
      visManager->Initialize();
#else
      G4cerr << "Visualization system not available,"
             << " please rebuild hdgeant4 with visualization enabled."
             << G4endl;
      exit(1);
#endif
   }

   simtimer.Start();

   // Start the user interface
   G4UImanager * UImanager = G4UImanager::GetUIpointer();  
   if (argc > optind) {   // batch mode  
      G4String command = "/control/execute ";
      G4String fileName = argv[optind];
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

   double nsim = GlueXPrimaryGeneratorAction::GetInstance()->getEventCount();
   double norm = (nsim == 0)? 1e9 : nsim;
   runtimer.Stop();
   double realrun = runtimer.GetRealElapsed();
   double sysrun = runtimer.GetSystemElapsed();
   double userrun = runtimer.GetUserElapsed();
   simtimer.Stop();
   double realsim = simtimer.GetRealElapsed();
   double syssim = simtimer.GetSystemElapsed();
   double usersim = simtimer.GetUserElapsed();
   char perevent[200];
   snprintf(perevent, 200, "%18f%18f%18f", realsim/norm, syssim/norm,
                                                         usersim/norm);
   char totalused[200];
   snprintf(totalused, 200, "%18f%18f%18f", realrun, sysrun, userrun);
   
   std::cout << nsim << " events generated." << std::endl
             << "Processing usage report:"
             << "     real time (s)   system time (s)     user time (s)" 
             << std::endl
             << "      per event average:" << perevent
             << std::endl
             << "          total for run:" << totalused
             << std::endl;
 
   // Close output file and clean up
   if (visManager)
      delete visManager;
   if (hddmOut)
      delete hddmOut;

   // Invoke mcsmear to smear the results, if requested
   if (nsim > 0 && hddmOut) {
      std::map<int, int> postsmear_opts;
      if (opts.Find("POSTSMEAR", postsmear_opts) &&
          postsmear_opts.size() > 0 && postsmear_opts[1] > 0)
      {
         std::stringstream cmd;
         cmd << "mcsmear ";
         std::map<int,std::string> mcsmear_opts;
         if (opts.Find("MCSMEAROPTS", mcsmear_opts)) {
            cmd << mcsmear_opts[1];
         }
         cmd << " " << outfile_opts[1];
         std::cout << "Smearing data with: " << std::endl 
                   << cmd.str() << std::endl;
         int res = system(cmd.str().c_str());
         if (opts.Find("DELETEUNSMEARED", postsmear_opts) && 
             postsmear_opts.size() > 0 && postsmear_opts[1] > 0)
         {
            if (res == 0) {
               cmd.str("");
               cmd << "rm -f " << outfile_opts[1];
               res = system(cmd.str().c_str());
            }
            else {
               std::cout << "Not deleting unsmeared file "
                         << "because of problems with the smearing."
                         << std::endl;
            }
         }
         return res;
      }
   }

   GlueXTimer::PrintAll();
   return 0;
}
