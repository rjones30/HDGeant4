//
// GlueXSteppingVerbose class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXSteppingVerbose.hh"
#include "GlueXPathFinder.hh"

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

int GlueXSteppingVerbose::instanceCount = 0;
pthread_mutex_t *GlueXSteppingVerbose::fMutex = 0;

GlueXSteppingVerbose::GlueXSteppingVerbose()
 : G4SteppingVerbose()
{
   pthread_mutex_t *mutex = new pthread_mutex_t;
   if (fMutex == 0) {
      fMutex = mutex;
      pthread_mutex_init(mutex, 0);
   }
   if (fMutex != mutex) {
      pthread_mutex_destroy(mutex);
      delete mutex;
   }
   ++instanceCount;
}

GlueXSteppingVerbose::GlueXSteppingVerbose(const GlueXSteppingVerbose &src)
 : G4SteppingVerbose(src)
{
   ++instanceCount;
}

GlueXSteppingVerbose::~GlueXSteppingVerbose()
{
   if (--instanceCount == 0) {
      pthread_mutex_destroy(fMutex);
      delete fMutex;
      fMutex = 0;
   }
}

void GlueXSteppingVerbose::StepInfo()
{
   pthread_mutex_lock(fMutex);
   CopyState();
   G4int prec = G4cout.precision(5);
 
   if (verboseLevel >= 1) {
     if (verboseLevel >= 4 ) VerboseTrack();
     if (verboseLevel >= 3) {
       G4cout << G4endl;
       G4cout << std::setw( 5) << "#Step#"     << " "
 	     << std::setw( 8) << "X"          << "    "
 	     << std::setw( 8) << "Y"          << "    "
 	     << std::setw( 8) << "Z"          << "    "
 	     << std::setw(11) << "KineE"      << " "
 	     << std::setw(11) << "dEStep"     << " "
 	     << std::setw(12) << "StepLeng"
 	     << std::setw(12) << "TrakLeng"
 	     << std::setw(10) << "Volume"    << "  "
 	     << std::setw(10) << "Process"   << G4endl;
     }
 
     G4TouchableHandle touch = GlueXPathFinder::CreateTouchableHandle();
     G4VPhysicalVolume *pvol = (touch)? touch->GetVolume() : 0;
     G4String volname = (pvol)? pvol->GetName() : "NULL";
     int copyno = (touch)? touch->GetCopyNumber() : 0;
     G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
 	<< std::setw(8) << G4BestUnit(fTrack->GetPosition().x(),"Length")
 	<< std::setw(8) << G4BestUnit(fTrack->GetPosition().y(),"Length")
 	<< std::setw(8) << G4BestUnit(fTrack->GetPosition().z(),"Length")
 	<< std::setw(8) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
 	<< std::setw(8) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
 	<< std::setw(8) << G4BestUnit(fStep->GetStepLength(),"Length")
 	<< std::setw(8) << G4BestUnit(fTrack->GetTrackLength(),"Length")
 	<< std::setw(10) << volname << ":" << copyno;
 
     const G4VProcess* process
                       = fStep->GetPostStepPoint()->GetProcessDefinedStep();
     G4String procName = " UserLimit";
     if (process) procName = process->GetProcessName();
     if (fStepStatus == fWorldBoundary) procName = "OutOfWorld";
     G4cout << "   " << std::setw(10) << procName;
     G4cout << G4endl;
 
     if (verboseLevel == 2) {
       G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
 	                    fN2ndariesAlongStepDoIt +
 	                    fN2ndariesPostStepDoIt;
       if(tN2ndariesTot>0){
 	G4cout << "\n    :----- List of secondaries ----------------"
 	       << G4endl;
         G4cout.precision(4);
 	for(size_t lp1=(*fSecondary).size()-tN2ndariesTot;
                         lp1<(*fSecondary).size(); lp1++){
 	  G4cout << "   "
 		 << std::setw(13)
 		 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName()
 		 << ":  energy ="
 		 << std::setw(6)
 		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
 		 << "  time ="
 		 << std::setw(6)
 		 << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time");	 
 	  G4cout << G4endl;
 	}
 
 	G4cout << "    :------------------------------------------\n"
 	       << G4endl;
       }
     }
 
   }

   G4cout.precision(prec);
   pthread_mutex_unlock(fMutex);
}

void GlueXSteppingVerbose::TrackingStarted()
{
   pthread_mutex_lock(fMutex);
   CopyState();
   G4int prec = G4cout.precision(5);

   if (verboseLevel > 0) {
     G4cout << std::setw( 5) << "Step#"      << " "
            << std::setw( 8) << "X"          << "    "
 	   << std::setw( 8) << "Y"          << "    "
 	   << std::setw( 8) << "Z"          << "    "
 	   << std::setw(11) << "KineE"      << " "
 	   << std::setw(11) << "dEStep"     << " "
 	   << std::setw(12) << "StepLeng"
 	   << std::setw(12) << "TrakLeng"
 	   << std::setw(10) << "Volume"     << "  "
 	   << std::setw(10) << "Process"    << G4endl;
 
     // PathFinder is not yet initialized for this track, so drill down
     // into the list of navigators and find the layer where this point lies.
 
     G4TouchableHistory *touch = 0;
     G4TransportationManager *tmanager =
                         G4TransportationManager::GetTransportationManager();
     std::vector<G4VPhysicalVolume*>::iterator iter = 
                                               tmanager->GetWorldsIterator();
     for (int world = tmanager->GetNoWorlds() - 1; world >= 0; --world) {
       G4Navigator *navigator = tmanager->GetNavigator(iter[world]);
       G4VPhysicalVolume *pvol = navigator->
                  LocateGlobalPointAndSetup(fTrack->GetPosition(),0,false);
       G4LogicalVolume *lvol = (pvol)? pvol->GetLogicalVolume() : 0;
       if (lvol && lvol->GetMaterial()) {
         touch = navigator->CreateTouchableHistory();
         break;
       }
     }
     G4VPhysicalVolume *pvol = (touch)? touch->GetVolume() : 0;
     G4String volname = (pvol)? pvol->GetName() : "NULL";
     int copyno = (touch)? touch->GetCopyNumber() : 0;
     G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
 	<< std::setw(8) << G4BestUnit(fTrack->GetPosition().x(),"Length")
 	<< std::setw(8) << G4BestUnit(fTrack->GetPosition().y(),"Length")
 	<< std::setw(8) << G4BestUnit(fTrack->GetPosition().z(),"Length")
 	<< std::setw(8) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
 	<< std::setw(8) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
 	<< std::setw(8) << G4BestUnit(fStep->GetStepLength(),"Length")
 	<< std::setw(8) << G4BestUnit(fTrack->GetTrackLength(),"Length")
 	<< std::setw(10) << volname << ":" << copyno
         << "   initStep" << G4endl;
   }

   G4cout.precision(prec);
   pthread_mutex_unlock(fMutex);
}
