//
// GlueXEventAction class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
 
#include "GlueXEventAction.hh"
#include "GlueXUserOptions.hh"

#include "G4ios.hh"

GlueXEventAction::GlueXEventAction()
 : G4UserEventAction()
{
   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXEventAction constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   fDebugPrint[0] = 0;
   fDebugPrint[1] = 0;
   fDebugPrint[2] = 1;
   std::map<int,int> debugpars;
   if (user_opts->Find("DEBU", debugpars)) {
      if (debugpars.find(1) != debugpars.end())
         fDebugPrint[0] = debugpars[1];
      if (debugpars.find(2) != debugpars.end())
         fDebugPrint[1] = debugpars[2];
      if (debugpars.find(3) != debugpars.end())
         fDebugPrint[2] = debugpars[3];
   }
}

void GlueXEventAction::BeginOfEventAction(const G4Event*)
{}

void GlueXEventAction::EndOfEventAction(const G4Event* evt)
{
   G4int event_id = evt->GetEventID();

   // periodic printing

   if (event_id % fDebugPrint[2] == 1 || 
      (event_id >= fDebugPrint[0] && event_id < fDebugPrint[1]))
   {
      G4cout << ">>> Event " << event_id << G4endl;
   }
}
