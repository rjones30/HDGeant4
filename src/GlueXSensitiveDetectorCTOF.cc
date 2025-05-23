//
// GlueXSensitiveDetectorCTOF - class implementation
//
// author: staylor at jlab.org
// version: october 25, 2021

#include "GlueXSensitiveDetectorCTOF.hh"
#include "GlueXDetectorConstruction.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "HddmOutput.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include <JANA/JApplication.h>
#include <JANA/Calibrations/JCalibrationManager.h>


// Cutoff on the number of allowed hits per bar
int GlueXSensitiveDetectorCTOF::MAX_HITS_PER_BAR = 25;

// Cutoff on the maximum time of flight
double GlueXSensitiveDetectorCTOF::MAX_TOF = 1000*ns;

// Light propagation parameters in tof bars
double GlueXSensitiveDetectorCTOF::ATTENUATION_LENGTH = 150.*cm;
double GlueXSensitiveDetectorCTOF::C_EFFECTIVE = 15*cm/ns;

// Geometric parameters of CTOF scintillators
double GlueXSensitiveDetectorCTOF::FULL_BAR_LENGTH = 120.*cm;

// Minimum hit time difference for two hits on the same bar
double GlueXSensitiveDetectorCTOF::TWO_HIT_TIME_RESOL = 25*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorCTOF::THRESH_MEV = 0.;

int GlueXSensitiveDetectorCTOF::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorCTOF::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorCTOF::GlueXSensitiveDetectorCTOF(const G4String& name)
 : G4VSensitiveDetector(name),
   fBarHitsMap(0), fPointsMap(0)
{
   collectionName.insert("CTOFBarHitsCollection");
   collectionName.insert("CTOFPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the CTOF, you must delete all old
   // objects of this class and create new ones.

   /* Copied over from FTOF, but currently no ccdb entries for CTOF
   G4AutoLock barrier(&fMutex);
   if (instanceCount++ == 0) {
      int runno = HddmOutput::getRunNo();
      extern JApplication *japp;
      if (japp == 0) {
         G4cerr << "Error in GlueXSensitiveDetector constructor - "
                << "jana global DApplication object not set, "
                << "cannot continue." << G4endl;
         exit(-1);
      }
      
      JCalibration *jcalib = japp->GetService<JCalibrationManager>()->GetJCalibration(runno);
      std::map<string, float> tof_parms;
      jcalib->Get("TOF/tof_parms", tof_parms);
      ATTENUATION_LENGTH = tof_parms.at("TOF_ATTEN_LENGTH")*cm;
      C_EFFECTIVE = tof_parms.at("TOF_C_EFFECTIVE")*cm/ns;
      FULL_BAR_LENGTH = tof_parms.at("TOF_PADDLE_LENGTH")*cm;
      TWO_HIT_TIME_RESOL = tof_parms.at("TOF_TWO_HIT_RESOL")*ns;
      THRESH_MEV = tof_parms.at("TOF_THRESH_MEV");
      MAX_HITS_PER_BAR = tof_parms.at("TOF_MAX_PAD_HITS");

      G4cout << "CTOF: ALL parameters loaded from ccdb" << G4endl;
   }
   */
}

GlueXSensitiveDetectorCTOF::GlueXSensitiveDetectorCTOF(
                     const GlueXSensitiveDetectorCTOF &src)
 : G4VSensitiveDetector(src),
   fBarHitsMap(src.fBarHitsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorCTOF &GlueXSensitiveDetectorCTOF::operator=(const
                                         GlueXSensitiveDetectorCTOF &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fBarHitsMap = src.fBarHitsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorCTOF::~GlueXSensitiveDetectorCTOF() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorCTOF::Initialize(G4HCofThisEvent* hce)
{
   fBarHitsMap = new
              GlueXHitsMapCTOFbar(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapCTOFpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fBarHitsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorCTOF::ProcessHits(G4Step* step, 
                                               G4TouchableHistory* ROhist)
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

   int barIndex = GetIdent("column", touch);

   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (trackinfo->GetGlueXHistory() == 0 && itrack > 0) {
      G4int key = fPointsMap->entries();
      GlueXHitCTOFpoint* lastPoint = (*fPointsMap)[key - 1];
      if (lastPoint == 0 || lastPoint->track_ != trackID ||
          fabs(lastPoint->t_ns - t/ns) > 0.1 ||
          fabs(lastPoint->x_cm - x[0]/cm) > 2. ||
          fabs(lastPoint->y_cm - x[1]/cm) > 2. ||
          fabs(lastPoint->z_cm - x[2]/cm) > 2.)
      {
         GlueXHitCTOFpoint newPoint;
         newPoint.ptype_G3 = g3type;
         newPoint.track_ = trackID;
         newPoint.trackID_ = itrack;
         newPoint.primary_ = (track->GetParentID() == 0);
         newPoint.t_ns = t/ns;
         newPoint.x_cm = x[0]/cm;
         newPoint.y_cm = x[1]/cm;
         newPoint.z_cm = x[2]/cm;
         newPoint.px_GeV = pin[0]/GeV;
         newPoint.py_GeV = pin[1]/GeV;
         newPoint.pz_GeV = pin[2]/GeV;
         newPoint.E_GeV = Ein/GeV;
         fPointsMap->add(key, newPoint);
      }
   }

   // Post the hit to the hits map, ordered by plane,bar,end index

   if (dEsum > 0) {
      int key = GlueXHitCTOFbar::GetKey(barIndex);
      GlueXHitCTOFbar *counter = (*fBarHitsMap)[key];
      if (counter == 0) {
         GlueXHitCTOFbar newcounter(barIndex);
         fBarHitsMap->add(key, newcounter);
         counter = (*fBarHitsMap)[key];
      }

      double dist = x[1];
      double dxtop = FULL_BAR_LENGTH / 2 - dist;
      double dxbottom = FULL_BAR_LENGTH / 2 + dist;

      // Calculate time at the PMT "normalized" to the center, 
      // so a hit in the center will have time "t" at both PMTs.
      double ttop = t + dxtop / C_EFFECTIVE;
      double tbottom = t + dxbottom / C_EFFECTIVE;
 
      // calculate energy seen by PM for this track step using attenuation factor
      double dEtop = dEsum * exp(-dxtop / ATTENUATION_LENGTH);
      double dEbottom = dEsum * exp(-dxbottom / ATTENUATION_LENGTH);

      // Add the hit to the hits vector, maintaining strict time ordering

      if (dEtop/MeV > 0 && ttop < MAX_TOF) {
         // add the hit on end=0 (top end of the bar)
         int merge_hit = 0;
         std::vector<GlueXHitCTOFbar::hitinfo_t>::iterator hiter;
         for (hiter = counter->hits.begin(); 
              hiter != counter->hits.end(); ++hiter)
         {
            if (hiter->end_ == 0) {
               if (fabs(hiter->t_ns*ns - ttop) < TWO_HIT_TIME_RESOL) {
                  merge_hit = 1;
                  break;
               }
               else if (hiter->t_ns*ns > ttop) {
                  break;
               }
            }
         }

         if (merge_hit) {
            // sum the charge, do energy weighting of the time --disabled for bad behavior --rtj
            //hiter->t_ns = (hiter->dE_GeV * hiter->t_ns + dEtop/GeV * ttop/ns) /
            //(hiter->dE_GeV += dEtop/GeV);
 
            // Use the time from the earlier hit but add the charge
            hiter->dE_GeV += dEtop/GeV;
            if (hiter->t_ns*ns > ttop) {
                hiter->t_ns = ttop/ns;
            }
         }
         else {
            // create new hit 
            hiter = counter->hits.insert(hiter, GlueXHitCTOFbar::hitinfo_t());
            hiter->end_ = 0;
            hiter->t_ns = ttop/ns;
            hiter->dE_GeV = dEtop/GeV;
         }
      }

      if (dEbottom/MeV > 0 && tbottom < MAX_TOF) {
         // add the hit on end=1 (bottom end of the bar)
         int merge_hit = 0;
         std::vector<GlueXHitCTOFbar::hitinfo_t>::iterator hiter;
         for (hiter = counter->hits.begin(); 
              hiter != counter->hits.end(); ++hiter)
         {
            if (hiter->end_ == 1) {
               if (fabs(hiter->t_ns*ns - tbottom) < TWO_HIT_TIME_RESOL) {
                  merge_hit = 1;
                  break;
               }
               else if (hiter->t_ns*ns > tbottom) {
                  break;
               }
            }
         }
         if (merge_hit) {
            // sum the charge, do energy weighting of the time --disabled for bad behavior --rtj
            // hiter->t_ns = (hiter->dE_GeV * hiter->t_ns +
            //               dEbottom/GeV * tbottom/ns) /
            //(hiter->dE_GeV += dEbottom/GeV);

            // Use the time from the earlier hit but add the charge
            hiter->dE_GeV += dEbottom/GeV;
            if (hiter->t_ns*ns > tbottom) {
               hiter->t_ns = tbottom/ns;
            }
         }
         else {
            // create new hit 
            hiter = counter->hits.insert(hiter, GlueXHitCTOFbar::hitinfo_t());
            hiter->end_ = 1;
            hiter->t_ns = tbottom/ns;
            hiter->dE_GeV = dEbottom/GeV;
         }
      }
   }
   return true;
}

void GlueXSensitiveDetectorCTOF::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitCTOFbar*> *bars = fBarHitsMap->GetMap();
   std::map<int,GlueXHitCTOFpoint*> *points = fPointsMap->GetMap();
   if (bars->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitCTOFbar*>::iterator siter;
   std::map<int,GlueXHitCTOFpoint*>::iterator piter;

   if (verboseLevel > 1) 
     { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << bars->size() << " bars with hits in the CTOF: "
             << G4endl;
      for (siter = bars->begin(); siter != bars->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the CTOF: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorCTOF::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getCppTOFs().size() == 0)
      hitview.addCppTOFs();
   hddm_s::CppTOF &cppTof = hitview.getCppTOF();

   // Collect and output the ctofTruthHits
   for (siter = bars->begin(); siter != bars->end(); ++siter) {
      std::vector<GlueXHitCTOFbar::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].dE_GeV*1e3 <= THRESH_MEV) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::CtofCounterList counter = cppTof.addCtofCounters(1);
         counter(0).setBar(siter->second->bar_);
         // first the end=0 hits
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            if (hits[ih].end_ == 0) {
               hddm_s::CtofTruthHitList thit = counter(0).addCtofTruthHits(1);
               thit(0).setEnd(hits[ih].end_);
               thit(0).setDE(hits[ih].dE_GeV);
               thit(0).setT(hits[ih].t_ns);
            }
         }
         // followed by the end=1 hits
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            if (hits[ih].end_ == 1) {
               hddm_s::CtofTruthHitList thit = counter(0).addCtofTruthHits(1);
               thit(0).setEnd(hits[ih].end_);
               thit(0).setDE(hits[ih].dE_GeV);
               thit(0).setT(hits[ih].t_ns);
            }
         }
         int hitscount = counter(0).getCtofTruthHits().size();
         if (hitscount > 2 * MAX_HITS_PER_BAR) {
            G4cerr << "GlueXSensitiveDetectorCTOF::EndOfEvent warning: "
                   << "max hit count " << 2 * MAX_HITS_PER_BAR << " exceeded, "
                   << hitscount - 2 * MAX_HITS_PER_BAR << " hits discarded."
                   << G4endl;
            counter(0).deleteCtofTruthHits(-1, 2 * MAX_HITS_PER_BAR);
         }
      }
   }

   // Collect and output the ctofTruthPoints
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::CtofTruthPointList point = cppTof.addCtofTruthPoints(1);
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

int GlueXSensitiveDetectorCTOF::GetIdent(std::string div, 
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
         int copyNum = touch->GetCopyNumber(depth);
         copyNum += (dynamic_cast<G4PVPlacement*>(pvol))? -1 : 0;
         return iter->second[copyNum];
      }
   }
   return -1;
}
