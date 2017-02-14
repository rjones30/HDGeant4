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
int GlueXSensitiveDetectorDIRC::MAX_HITS = 100;

// Minimum hit time difference for two hits on the same tube
double GlueXSensitiveDetectorDIRC::TWO_HIT_TIME_RESOL = 50*ns;

// Minimum photoelectron count for a hit
double GlueXSensitiveDetectorDIRC::THRESH_PE = 2.;

int GlueXSensitiveDetectorDIRC::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorDIRC::fMutex = G4MUTEX_INITIALIZER;

std::map<G4LogicalVolume*, int> GlueXSensitiveDetectorDIRC::fVolumeTable;

GlueXSensitiveDetectorDIRC::GlueXSensitiveDetectorDIRC(const G4String& name)
 : G4VSensitiveDetector(name),
   fFlashesMap(0), fPointsMap(0)
{
   collectionName.insert("DIRCFlashesCollection");
   collectionName.insert("DIRCPointsCollection");

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

GlueXSensitiveDetectorDIRC::GlueXSensitiveDetectorDIRC(
                     const GlueXSensitiveDetectorDIRC &src)
 : G4VSensitiveDetector(src),
   fFlashesMap(src.fFlashesMap), fPointsMap(src.fPointsMap)
{
   ++instanceCount;
}

GlueXSensitiveDetectorDIRC &GlueXSensitiveDetectorDIRC::operator=(const
                                         GlueXSensitiveDetectorDIRC &src)
{
   *(G4VSensitiveDetector*)this = src;
   fFlashesMap = src.fFlashesMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorDIRC::~GlueXSensitiveDetectorDIRC() 
{
   --instanceCount;
}

void GlueXSensitiveDetectorDIRC::Initialize(G4HCofThisEvent* hce)
{
   fFlashesMap = new 
              GlueXHitsMapDIRCflash(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapDIRCpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fFlashesMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorDIRC::ProcessHits(G4Step* step, 
                                               G4TouchableHistory* unused)
{
   double dEsum = step->GetTotalEnergyDeposit();
   if (dEsum == 0)
      return false;

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
   const G4AffineTransform &local_from_global = touch->GetHistory()
                                                     ->GetTopTransform();
   G4ThreeVector xlocal = local_from_global.TransformPoint(x);
  
   // For particles that range out inside the active volume, the
   // "out" time may sometimes be set to something enormously high.
   // This screws up the hit. Check for this case here by looking
   // at tout and making sure it is less than 1 second. If it's
   // not, then just use tin for "t".

   if (tout > 1.0*s)
      t = tin;

   // Post the hit to the points list in the
   // order of appearance in the event simulation.

   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (touch->GetVolume()->GetName() == "RDCD") {
      if (trackinfo->GetGlueXHistory() == 0 && itrack > 0 && xin.dot(pin) > 0) {
         G4int key = fPointsMap->entries();
         GlueXHitDIRCpoint* lastPoint = (*fPointsMap)[key - 1];
         if (lastPoint == 0 || lastPoint->track_ != trackID ||
             fabs(lastPoint->t_ns - t/ns) > 0.1 ||
             fabs(lastPoint->x_cm - x[0]/cm) > 2. ||
             fabs(lastPoint->y_cm - x[1]/cm) > 2. ||
             fabs(lastPoint->z_cm - x[2]/cm) > 2.)
         {
            GlueXHitDIRCpoint* newPoint = new GlueXHitDIRCpoint();
            fPointsMap->add(key, newPoint);
            int pdgtype = track->GetDynamicParticle()->GetPDGcode();
            int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
            newPoint->ptype_G3 = g3type;
            newPoint->track_ = trackID;
            newPoint->trackID_ = itrack;
            newPoint->primary_ = (track->GetParentID() == 0);
            newPoint->t_ns = t/ns;
            newPoint->x_cm = x[0]/cm;
            newPoint->y_cm = x[1]/cm;
            newPoint->z_cm = x[2]/cm;
            newPoint->px_GeV = pin[0]/GeV;
            newPoint->py_GeV = pin[1]/GeV;
            newPoint->pz_GeV = pin[2]/GeV;
            newPoint->E_GeV = Ein/GeV;
         }
      }
      // The DIRC is a special type of detector in that the volume used
      // to generate truth points is different from the volume where the
      // DIRC hits are generated. This is not true in the present pseudo
      // hits code, but as soon as phototube array hits are implemented,
      // the following line should be uncommented.
      //return true;
   }

   // Post the hit to the hits map, ordered by time
   // This code makes pseudo-hits called "flashes", and
   // needs to be replaced with actual simulated DIRC hits.

   if (dEsum > 0) {
      // int bar = GetIdent("bar", touch);
      // int key = GlueXHitDIRCflash::GetKey(bar);
      int key = 0;
      GlueXHitDIRCflash *flash = (*fFlashesMap)[key];
      if (flash == 0) {
         flash = new GlueXHitDIRCflash();
         fFlashesMap->add(key, flash);
      }

      // Add the hit to the hits vector, maintaining strict time ordering

      int merge_hit = 0;
      std::vector<GlueXHitDIRCflash::hitinfo_t>::iterator hiter;
      for (hiter = flash->hits.begin(); hiter != flash->hits.end(); ++hiter)
      {
         //if (fabs(hiter->t_ns*ns - t) < TWO_HIT_TIME_RESOL) {
         //   merge_hit = 1;
         //   break;
         //}
         if (hiter->t_ns*ns > t) {
            break;
         }
      }
      if (merge_hit) { // more like overwrite
         if (hiter->t_ns*ns > t) {
            hiter->E_GeV = Ein/GeV;
            hiter->t_ns = t/ns;
            hiter->x_cm = x[0]/cm;
            hiter->y_cm = x[1]/cm;
            hiter->z_cm = x[2]/cm;
         }
      }
      else if ((int)flash->hits.size() < MAX_HITS)	{ // create new hit 
         hiter = flash->hits.insert(hiter, GlueXHitDIRCflash::hitinfo_t());
         hiter->E_GeV = Ein/GeV;
         hiter->t_ns = t/ns;
         hiter->x_cm = x[0]/cm;
         hiter->y_cm = x[1]/cm;
         hiter->z_cm = x[2]/cm;
      }
      else {
         G4cerr << "GlueXSensitiveDetectorDIRC::ProcessHits error: "
             << "max hit count " << MAX_HITS << " exceeded, truncating!"
             << G4endl;
      }
   }
   return true;
}

void GlueXSensitiveDetectorDIRC::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitDIRCflash*> *flashes = fFlashesMap->GetMap();
   std::map<int,GlueXHitDIRCpoint*> *points = fPointsMap->GetMap();
   if (flashes->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitDIRCflash*>::iterator siter;
   std::map<int,GlueXHitDIRCpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << flashes->size() << " flashes registered in the DIRC: "
             << G4endl;
      for (siter = flashes->begin(); siter != flashes->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the DIRC: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
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

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getDIRCs().size() == 0)
      hitview.addDIRCs();
   hddm_s::DIRC &dirc = hitview.getDIRC();

   // Collect and output the dircTruthHits
   for (siter = flashes->begin(); siter != flashes->end(); ++siter) {
      std::vector<GlueXHitDIRCflash::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=1; ih < hits.size(); ++ih) {
         if (fabs(hits[ih].t_ns - hits[ih-1].t_ns) < 0.1 &&
             fabs(hits[ih].x_cm - hits[ih-1].x_cm) < 0.1 &&
             fabs(hits[ih].y_cm - hits[ih-1].y_cm) < 0.1 &&
             fabs(hits[ih].z_cm - hits[ih-1].z_cm) < 0.1)
         {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            hddm_s::DircTruthHitList thit = dirc.addDircTruthHits(1);
            thit(0).setE(hits[ih].E_GeV);
            thit(0).setT(hits[ih].t_ns);
            thit(0).setX(hits[ih].x_cm);
            thit(0).setY(hits[ih].y_cm);
            thit(0).setZ(hits[ih].z_cm);
         }
      }
   }

   // Collect and output the DircTruthPoints
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::DircTruthPointList point = dirc.addDircTruthPoints(1);
      point(0).setPrimary(piter->second->primary_);
      point(0).setPtype(piter->second->ptype_G3);
      point(0).setPx(piter->second->px_GeV);
      point(0).setPy(piter->second->py_GeV);
      point(0).setPz(piter->second->pz_GeV);
      point(0).setE(piter->second->E_GeV);
      point(0).setX(piter->second->x_cm);
      point(0).setY(piter->second->y_cm);
      point(0).setZ(piter->second->z_cm);
      point(0).setT(piter->second->t_ns);
      point(0).setTrack(piter->second->track_);
      hddm_s::TrackIDList tid = point(0).addTrackIDs();
      tid(0).setItrack(piter->second->trackID_);
   }
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
