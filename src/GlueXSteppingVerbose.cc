//
// GlueXSteppingVerbose class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXSteppingVerbose.hh"
#include "GlueXPathFinder.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

void GlueXSteppingVerbose::StepInfo()
{
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

    G4TouchableHandle phand = GlueXPathFinder::CreateTouchableHandle();
    G4VPhysicalVolume *pvol = (phand)? phand->GetVolume() : 0;
    G4String volname = (pvol)? pvol->GetName() : "NULL";
    int copyno = (phand)? phand->GetCopyNumber() : 0;
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
}

void GlueXSteppingVerbose::TrackingStarted()
{

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

    G4TouchableHandle phand = GlueXPathFinder::CreateTouchableHandle();
    G4VPhysicalVolume *pvol = (phand)? phand->GetVolume() : 0;
    G4String volname = (pvol)? pvol->GetName() : "NULL";
    int copyno = (phand)? phand->GetCopyNumber() : 0;
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
}
