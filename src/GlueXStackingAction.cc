//
// GlueXStackingAction class implementation
//
// author: richard.t.jones at uconn.edu
// version: january 29, 2017

#include "GlueXStackingAction.hh"
#include "GlueXUserOptions.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXSensitiveDetectorDIRC.hh"

#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
#include "TMath.h"

GlueXStackingAction::GlueXStackingAction()
{
   if (!(G4ParticleTable::GetParticleTable()->GetReadiness())) {
      G4String msg;
      msg =  " You are instantiating GlueXStackingAction BEFORE your\n";
      msg += "G4VUserPhysicsList is instantiated and assigned to ";
      msg += "G4RunManager.\nSuch an instantiation is prohibited by ";
      msg += "Geant4 version 8.0. To fix this problem,\n";
      msg += "please make sure that your main() instantiates ";
      msg += "G4VUserPhysicsList AND\nset it to G4RunManager before ";
      msg += "instantiating other user action classes\n";
      msg += "such as GlueXStackingAction.";
      G4Exception("GlueXStackingAction::GlueXStackingAction()",
                  "Event0031",FatalException,msg);
   }

   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXStackingAction constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   nosecondaries = 0;
   std::map<int,int> opt;
   if (user_opts->Find("NOSECONDARIES", opt)) {
      nosecondaries = (opt[1] != 0);
   }
}

GlueXStackingAction::~GlueXStackingAction()
{}

G4ClassificationOfNewTrack GlueXStackingAction::ClassifyNewTrack(
                                                const G4Track *aTrack)
{
   // Determine what to do with this secondary, which has just been
   // made into a G4Track object and needs to be classified as follows.
   //
   //    enum G4ClassificationOfNewTrack
   //    {
   //      fUrgent,    // put into the urgent stack
   //      fWaiting,   // put into the waiting stack
   //      fPostpone,  // postpone to the next event
   //      fKill       // kill without stacking
   //    };
   //
   //    The parent_ID of the track indicates the origin of it.
   //                
   //    G4int parent_ID = aTrack->get_parentID();
   //   
   //      parent_ID = 0 : primary particle
   //                > 0 : secondary particle
   //                < 0 : postponed from the previous event

   if (nosecondaries && aTrack->GetParentID() != 0)
      return fKill;

   // apply detection efficiency for the DIRC at production stage:
   G4String ParticleName = aTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
   if (aTrack->GetParentID() > 0) { // particle is secondary
      if (ParticleName == "opticalphoton") {
         Double_t Ephoton = aTrack->GetMomentum().mag();
         Double_t ra = G4UniformRand();
         if (ra > GlueXSensitiveDetectorDIRC::GetDetectionEfficiency(Ephoton))
	        return fKill;
       }
   }
   
   return fUrgent;
}

void GlueXStackingAction::NewStage()
{
   // This method is called by G4StackManager when the urgentStack
   // becomes empty and contents in the waitingStack are transtered
   // to the urgentStack. Note that this method is not called at the
   // begining of each event, see "PrepareNewEvent" for that.
   //
   // In case re-classification of the stacked tracks is needed,
   // use the following method to request to G4StackManager.
   //
   //    stackManager->ReClassify();
   //
   // All of the stacked tracks in the waitingStack will be re-classified 
   // by "ClassifyNewTrack" method. To abort current event, do
   //
   //    stackManager->clear();
   //
   // Note that this way is valid and safe only for the case it is called
   // from this user class. The more global way of event abortion is
   //
   //    G4UImanager * UImanager = G4UImanager::GetUIpointer();
   //    UImanager->ApplyCommand("/event/abort");
}

void GlueXStackingAction::PrepareNewEvent()
{
   // This method is called by G4StackManager at the begining of each event.
   // Be careful that the urgentStack and the waitingStack of G4StackManager
   // are empty at this moment, because this method is called before 
   // accepting primary particles. Also, note that the postponeStack of
   // G4StackManager may have some postponed tracks.
}
