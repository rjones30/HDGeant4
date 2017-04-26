//
// GlueXSensitiveDetectorDIRC - class implementation
//
// author: richard.t.jones at uconn.edu
// version: november 28, 2016

#include "GlueXSensitiveDetectorDIRC.hh"
#include "GlueXDetectorConstruction.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include <JANA/JApplication.h>

// Geant3-style particle index for optical photon
#define OPTICAL_PHOTON 50

// Cutoff on the total number of allowed hits
int GlueXSensitiveDetectorDIRC::MAX_HITS = 500;

// Minimum hit time difference for two hits on the same tube
double GlueXSensitiveDetectorDIRC::TWO_HIT_TIME_RESOL = 50*ns;

int GlueXSensitiveDetectorDIRC::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorDIRC::fMutex = G4MUTEX_INITIALIZER;

std::map<G4LogicalVolume*, int> GlueXSensitiveDetectorDIRC::fVolumeTable;

GlueXSensitiveDetectorDIRC::GlueXSensitiveDetectorDIRC(const G4String& name)
  : G4VSensitiveDetector(name)
{
  // The rest of this only needs to happen once, the first time an object
  // of this type is instantiated for this configuration of geometry and
  // fields. If the geometry or fields change in such a way as to modify
  // the drift-time properties of hits in the DIRC, you must delete all old
  // objects of this class and create new ones.

  G4AutoLock tuberier(&fMutex);
  if (instanceCount++ == 0) {
    extern int run_number;
    extern jana::JApplication *japp;
    if (japp == 0) {
      G4cerr << "Error in GlueXSensitiveDetector constructor - "
	     << "jana global DApplication object not set, "
	     << "cannot continue." << G4endl;
      exit(-1);
    }
    jana::JCalibration *jcalib = japp->GetJCalibration(run_number);
    if (japp == 0) {   // dummy
      jcalib = 0;
      G4cout << "DIRC: ALL parameters loaded from ccdb" << G4endl;
    }
  }
}

    GlueXSensitiveDetectorDIRC::GlueXSensitiveDetectorDIRC(const GlueXSensitiveDetectorDIRC &src)
      : G4VSensitiveDetector(src)
{
  ++instanceCount;
}

GlueXSensitiveDetectorDIRC &GlueXSensitiveDetectorDIRC::operator=(const
								  GlueXSensitiveDetectorDIRC &src)
{
  *(G4VSensitiveDetector*)this = src;
  return *this;
}

GlueXSensitiveDetectorDIRC::~GlueXSensitiveDetectorDIRC() 
{
  --instanceCount;
}

void GlueXSensitiveDetectorDIRC::Initialize(G4HCofThisEvent* hce)
{
}

G4bool GlueXSensitiveDetectorDIRC::ProcessHits(G4Step* step, 
                                               G4TouchableHistory* unused)
{
  
  const G4ThreeVector &pin = step->GetPreStepPoint()->GetMomentum();
  const G4ThreeVector &xin = step->GetPreStepPoint()->GetPosition();
  const G4ThreeVector &xout = step->GetPostStepPoint()->GetPosition();
  double Ein = step->GetPreStepPoint()->GetTotalEnergy();
  double tin = step->GetPreStepPoint()->GetGlobalTime();
  double tout = step->GetPostStepPoint()->GetGlobalTime();
  G4ThreeVector x = (xin + xout) / 2;
  G4ThreeVector dx = xout - xin;
  double t = (tin + tout) / 2;

  const G4VTouchable* touch = step->GetPreStepPoint()->GetTouchable();
  const G4TouchableHistory* touch_hist = (G4TouchableHistory*)touch;
  const G4AffineTransform &local_from_global = touch->GetHistory()->GetTopTransform();
  G4ThreeVector xlocal = local_from_global.TransformPoint(x);
  
  // For particles that range out inside the active volume, the
  // "out" time may sometimes be set to something enormously high.
  // This screws up the hit. Check for this case here by looking
  // at tout and making sure it is less than 1 second. If it's
  // not, then just use tin for "t".

  if (tout > 1.0*s) t = tin;

  // Post the hit to the points list in the
  // order of appearance in the event simulation.

  G4Track *track = step->GetTrack();
  // G4int trackID = track->GetTrackID();
  GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
    track->GetUserInformation();
  int itrack = trackinfo->GetGlueXTrackID();

  // radiator volume
  if (touch->GetVolume()->GetName() == "QZBL") {
    if (trackinfo->GetGlueXHistory() == 0 && itrack > 0 && xin.dot(pin) > 0) {
      int pdgtype = track->GetDynamicParticle()->GetPDGcode();
      int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
      
      GlueXHitDIRCBar barhit;
      barhit.E_GeV = Ein/GeV;
      barhit.t_ns = t/ns;
      barhit.x_cm = x[0]/cm;
      barhit.y_cm = x[1]/cm;
      barhit.z_cm = x[2]/cm;
      barhit.px_GeV = pin[0]/GeV;
      barhit.py_GeV = pin[1]/GeV;
      barhit.pz_GeV = pin[2]/GeV;
      barhit.pdg = g3type;
      barhit.bar = touch_hist->GetReplicaNumber(0)/4; // each bar is glued from 4 pieces
      barhit.track = itrack;
      fHitsBar.push_back(barhit);
    }
    return true;
  }

  // MCP's pixel volume
  if (touch->GetVolume()->GetName() == "PIXV"){
    if ((int)fHitsMcp.size() < MAX_HITS){

      GlueXHitDIRCMcp mcphit;
      mcphit.E_GeV = Ein/GeV;
      mcphit.t_ns = t/ns;
      mcphit.x_cm = x[0]/cm;
      mcphit.y_cm = x[1]/cm;
      mcphit.z_cm = x[2]/cm;
      
      G4double box = touch_hist->GetReplicaNumber(3);   // [1,2]
      G4double mcp = touch_hist->GetReplicaNumber(1)-1; // [0,101]
      G4double pix = touch_hist->GetReplicaNumber(0)-1; // [0,63]

      mcphit.ch = box*mcp*64+pix;      
      fHitsMcp.push_back(mcphit);
    }else {
      G4cerr << "GlueXSensitiveDetectorDIRC::ProcessHits error: "
             << "max hit count " << MAX_HITS << " exceeded, truncating!"
             << G4endl;
    }
  }
  return true;
}

void GlueXSensitiveDetectorDIRC::EndOfEvent(G4HCofThisEvent*)
{
  if (fHitsBar.size() == 0 || fHitsMcp.size() == 0)
    return;
  
  if (verboseLevel > 1) { 
    G4cout << G4endl
	   << "--------> Hits Collection: in this event there are "
	   << fHitsBar.size() << " bar hits:"
	   << G4endl;
    for(unsigned int h=0; h<fHitsBar.size(); h++)
      fHitsBar[h].Print();
    
    G4cout << G4endl
	   << "--------> Hits Collection: in this event there are "
	   << fHitsMcp.size() << " MCP hits: "
	   << G4endl;
    for(unsigned int h=0; h<fHitsMcp.size(); h++)
      fHitsMcp[h].Print();
  }

  // pack hits into ouptut hddm record
 
  G4EventManager* mgr = G4EventManager::GetEventManager();
  G4VUserEventInformation* info = mgr->GetUserInformation();
  hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
  if (record == 0) {
    G4cerr << "GlueXSensitiveDetectorDIRC::EndOfEvent error - "
	   << "hits seen but no output hddm record to save them into, "
	   << "cannot continue!" << G4endl;
    exit(1);
  }

  if (record->getPhysicsEvents().size() == 0) record->addPhysicsEvents();
  if (record->getHitViews().size() == 0) record->getPhysicsEvent().addHitViews();
  hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
  if (hitview.getDIRCs().size() == 0) hitview.addDIRCs();
  hddm_s::DIRC &dirc = hitview.getDIRC();

  // Collect and output the BarHits
  for(unsigned int h=0; h<fHitsBar.size(); h++){
    hddm_s::DircTruthBarHitList bhit = dirc.addDircTruthBarHits(1);
    bhit(0).setE(fHitsBar[h].E_GeV);
    bhit(0).setT(fHitsBar[h].t_ns);
    bhit(0).setX(fHitsBar[h].x_cm);
    bhit(0).setY(fHitsBar[h].y_cm);
    bhit(0).setZ(fHitsBar[h].z_cm);
    bhit(0).setPx(fHitsBar[h].px_GeV);
    bhit(0).setPy(fHitsBar[h].py_GeV);
    bhit(0).setPz(fHitsBar[h].pz_GeV);
    bhit(0).setPdg(fHitsBar[h].pdg);
    bhit(0).setBar(fHitsBar[h].bar);
    bhit(0).setTrack(fHitsBar[h].track);
  }

  // Collect and output the DircTruthPoints
  for(unsigned int h=0; h<fHitsMcp.size(); h++){
    hddm_s::DircTruthMcpHitList mhit = dirc.addDircTruthMcpHits(1);
    mhit(0).setE(fHitsMcp[h].E_GeV);
    mhit(0).setT(fHitsMcp[h].t_ns);
    mhit(0).setX(fHitsMcp[h].x_cm);
    mhit(0).setY(fHitsMcp[h].y_cm);
    mhit(0).setZ(fHitsMcp[h].z_cm);
    mhit(0).setCh(fHitsMcp[h].ch);
    mhit(0).setKey_bar(fHitsMcp[h].key_bar);
  }

  fHitsBar.clear();
  fHitsMcp.clear();
}

int GlueXSensitiveDetectorDIRC::GetIdent(std::string div, 
                                         const G4VTouchable *touch)
{
  const HddsG4Builder* bldr = GlueXDetectorConstruction::GetBuilder();
  std::map<std::string, std::vector<int> >::const_iterator iter;
  std::map<std::string, std::vector<int> > *identifiers;
  int max_depth = touch->GetHistoryDepth();
  for (int depth = 0; depth < max_depth; ++depth) {
    G4VPhysicalVolume *pvol = touch->GetVolume(depth);
    G4LogicalVolume *lvol = pvol->GetLogicalVolume();
    int volId = fVolumeTable[lvol];
    if (volId == 0) {
      volId = bldr->getVolumeId(lvol);
      fVolumeTable[lvol] = volId;
    }
    identifiers = &Refsys::fIdentifierTable[volId];
    if ((iter = identifiers->find(div)) != identifiers->end()) {
      if (dynamic_cast<G4PVPlacement*>(pvol))
	return iter->second[pvol->GetCopyNo() - 1];
      else
	return iter->second[pvol->GetCopyNo()];
    }
  }
  return -1;
}
