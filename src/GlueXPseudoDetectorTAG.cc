//
// class implementation for GlueXPseudoDetectorTAG
//
// author: richard.t.jones at uconn.edu
// version: october 20, 2016

#include "GlueXPseudoDetectorTAG.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserOptions.hh"
#include "G4EventManager.hh"

#include <JANA/JApplication.h>

double GlueXPseudoDetectorTAG::TRIGGER_TIME_SIGMA = 10.*ns;
double GlueXPseudoDetectorTAG::TIME_REFERENCE_PLANE_Z = 65.*cm;
double GlueXPseudoDetectorTAG::BEAM_BUCKET_PERIOD = 4.*ns;
double GlueXPseudoDetectorTAG::BEAM_VELOCITY = 2.99792458e8*m/s;

int GlueXPseudoDetectorTAG::HODO_MAX_HITS = 5000;
int GlueXPseudoDetectorTAG::MICRO_MAX_HITS = 5000;
double GlueXPseudoDetectorTAG::HODO_TWO_HIT_TIME_RESOL = 25.*ns;
double GlueXPseudoDetectorTAG::MICRO_TWO_HIT_TIME_RESOL = 25.*ns;
double GlueXPseudoDetectorTAG::MICRO_HIT_DE = 3.5*MeV;
double GlueXPseudoDetectorTAG::HODO_HIT_DE = 0.55*MeV;

double GlueXPseudoDetectorTAG::MICRO_LIMITS_ERANGE[2];
int GlueXPseudoDetectorTAG::MICRO_CHANNEL_NUMBER[GlueXPseudoDetectorTAG::MICRO_NCHANNELS];
double GlueXPseudoDetectorTAG::MICRO_CHANNEL_EMIN[GlueXPseudoDetectorTAG::MICRO_NCHANNELS];
double GlueXPseudoDetectorTAG::MICRO_CHANNEL_EMAX[GlueXPseudoDetectorTAG::MICRO_NCHANNELS];

double GlueXPseudoDetectorTAG::HODO_LIMITS_ERANGE[2];
int GlueXPseudoDetectorTAG::HODO_CHANNEL_NUMBER[GlueXPseudoDetectorTAG::HODO_NCHANNELS];
double GlueXPseudoDetectorTAG::HODO_CHANNEL_EMIN[GlueXPseudoDetectorTAG::HODO_NCHANNELS];
double GlueXPseudoDetectorTAG::HODO_CHANNEL_EMAX[GlueXPseudoDetectorTAG::HODO_NCHANNELS];

double GlueXPseudoDetectorTAG::TAGGER_TMIN_NS = -200.*ns;
double GlueXPseudoDetectorTAG::TAGGER_TMAX_NS = +200.*ns;

int GlueXPseudoDetectorTAG::instanceCount = 0;
G4Mutex GlueXPseudoDetectorTAG::fMutex = G4MUTEX_INITIALIZER;

GlueXPseudoDetectorTAG::GlueXPseudoDetectorTAG(int run_number)
 : fRunNo(0)
{
   if (run_number != 0)
      setRunNo(run_number);
   instanceCount++;

   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXPseudoDetectorTAG constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   std::map<int,double> trefpars;
   if (user_opts->Find("TREFSIGMA", trefpars)) {
      TRIGGER_TIME_SIGMA = trefpars[1]*ns;
   }
}

GlueXPseudoDetectorTAG::~GlueXPseudoDetectorTAG()
{
   --instanceCount;
}

GlueXPseudoDetectorTAG::GlueXPseudoDetectorTAG(GlueXPseudoDetectorTAG &src)
{}

GlueXPseudoDetectorTAG& GlueXPseudoDetectorTAG::operator=(GlueXPseudoDetectorTAG &src)
{
   return *this;
}

inline void GlueXPseudoDetectorTAG::setRunNo(int runno)
{
   fRunNo = runno;
   G4AutoLock barrier(&fMutex);
   extern int run_number;
   extern jana::JApplication *japp;
   if (japp == 0) {
      G4cerr << "Error in GlueXPseudoDetectorTAG constructor - "
             << "jana global DApplication object not set, "
             << "cannot continue." << G4endl;
      exit(-1);
   }
   run_number = runno;
   jana::JCalibration *jcalib = japp->GetJCalibration(run_number);
   std::map<string, float> rf_parms;
   jcalib->Get("PHOTON_BEAM/RF/beam_period", rf_parms);
   BEAM_BUCKET_PERIOD = rf_parms.at("beam_period")*ns;
   std::map<string, float> beam_parms;
   jcalib->Get("PHOTON_BEAM/endpoint_energy", beam_parms);
   double endpoint_energy = beam_parms.at("PHOTON_BEAM_ENDPOINT_ENERGY")*GeV;
   std::vector<std::map<string, float> > micro_parms;
   jcalib->Get("PHOTON_BEAM/microscope/scaled_energy_range", micro_parms);
   double Emin = 1e9;
   double Emax = -1e9;
   for (unsigned int i=0; i < micro_parms.size(); ++i) {
      MICRO_CHANNEL_NUMBER[i] = micro_parms[i]["column"];
      MICRO_CHANNEL_EMIN[i] = micro_parms[i]["xlow"] * endpoint_energy;
      MICRO_CHANNEL_EMAX[i] = micro_parms[i]["xhigh"] * endpoint_energy;
      Emin = (MICRO_CHANNEL_EMIN[i] < Emin)? MICRO_CHANNEL_EMIN[i] : Emin;
      Emax = (MICRO_CHANNEL_EMAX[i] > Emax)? MICRO_CHANNEL_EMAX[i] : Emax;
   }
   MICRO_LIMITS_ERANGE[0] = Emin;
   MICRO_LIMITS_ERANGE[1] = Emax;
   std::vector<std::map<string, float> > hodo_parms;
   jcalib->Get("PHOTON_BEAM/hodoscope/scaled_energy_range", hodo_parms);
   for (unsigned int i=0; i < hodo_parms.size(); ++i) {
      HODO_CHANNEL_NUMBER[i] = hodo_parms[i]["column"];
      HODO_CHANNEL_EMIN[i] = hodo_parms[i]["xlow"] * endpoint_energy;
      HODO_CHANNEL_EMAX[i] = hodo_parms[i]["xhigh"] * endpoint_energy;
      Emin = (HODO_CHANNEL_EMIN[i] < Emin)? HODO_CHANNEL_EMIN[i] : Emin;
      Emax = (HODO_CHANNEL_EMAX[i] > Emax)? HODO_CHANNEL_EMAX[i] : Emax;
   }
   HODO_LIMITS_ERANGE[0] = Emin;
   HODO_LIMITS_ERANGE[1] = Emax;

   G4cout << "TAGGER: all parameters loaded from ccdb" << G4endl;
}

int GlueXPseudoDetectorTAG::addTaggerHit(G4ThreeVector &vertex,
                                         double energy,
                                         double time, 
                                         int bg)
{
   // look up which tagger channel is hit, if any

   int micro_channel = -1;
   int hodo_channel = -1;
   double E = energy;
   if (E > MICRO_LIMITS_ERANGE[0] && E < MICRO_LIMITS_ERANGE[1]) {
      int i = MICRO_NCHANNELS * (E - MICRO_LIMITS_ERANGE[0]) /
              (MICRO_LIMITS_ERANGE[1] - MICRO_LIMITS_ERANGE[0]);
      while (E < MICRO_CHANNEL_EMIN[i])
         --i;
      while (E > MICRO_CHANNEL_EMAX[i])
         ++i;
      if (E >= MICRO_CHANNEL_EMIN[i] && E <= MICRO_CHANNEL_EMAX[i]) {
         E = (MICRO_CHANNEL_EMIN[i] + MICRO_CHANNEL_EMAX[i]) / 2;
         micro_channel = MICRO_CHANNEL_NUMBER[i];
      }
   }
   else if (E > HODO_LIMITS_ERANGE[0] && E < HODO_LIMITS_ERANGE[1]) {
      int i = HODO_NCHANNELS * (E - HODO_LIMITS_ERANGE[0]) /
              (HODO_LIMITS_ERANGE[1] - HODO_LIMITS_ERANGE[0]);
      while (E < HODO_CHANNEL_EMIN[i])
         --i;
      while (E > HODO_CHANNEL_EMAX[i])
         ++i;
      if (E >= HODO_CHANNEL_EMIN[i] && E <= HODO_CHANNEL_EMAX[i]) {
         E = (HODO_CHANNEL_EMIN[i] + HODO_CHANNEL_EMAX[i]) / 2;
         hodo_channel = HODO_CHANNEL_NUMBER[i];
      }
   }

   double t = time - (vertex[2] - TIME_REFERENCE_PLANE_Z) / BEAM_VELOCITY;
   t = floor(t / BEAM_BUCKET_PERIOD + 0.5) * BEAM_BUCKET_PERIOD;

   // pack hit into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXPseudoDetectorTAG::addTaggerHit error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getTaggers().size() == 0)
      hitview.addTaggers();
   hddm_s::Tagger &tagger = hitview.getTagger();

   // post microscope hit to the microscope hits tree

   if (micro_channel > -1) {
      hddm_s::MicroChannelList columns = tagger.getMicroChannels();
      hddm_s::MicroChannelList::iterator citer;
      for (citer = columns.begin(); citer != columns.end(); ++citer) {
         if (citer->getColumn() == micro_channel)
            break;
      }
      if (citer == columns.end()) {
         columns = tagger.addMicroChannels();
         citer = columns.begin();
         citer->setColumn(micro_channel);
         citer->setRow(0);
         citer->setE(E/GeV);
      }
      hddm_s::TaggerTruthHitList hits = citer->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator hiter;
      int hit = 0;
      for (hiter = hits.begin(); hiter != hits.end(); ++hiter, ++hit) {
         if (fabs(t - hiter->getT()*ns) < MICRO_TWO_HIT_TIME_RESOL)
            break;
         else if (t > hiter->getT()*ns) {
            hits = citer->addTaggerTruthHits(1, ++hit);
            hits(0).setE(energy/GeV);
            hits(0).setT(1e99);
            hits(0).setDE(0);
            hits(0).setBg(bg);
            break;
         }
      }
      if (hiter != hits.end()) {              // merge with a prior hit
         if (t < hiter->getT()*ns) {
            hiter->setT(t/ns);
            hiter->setE(energy/GeV);
            hiter->setBg(bg);
         }
         hiter->setDE(hiter->getDE() + MICRO_HIT_DE/GeV);
      }
      else {                                 // make a new hit
         hiter->setT(t/ns);
         hiter->setE(energy/GeV);
         hiter->setDE(MICRO_HIT_DE/GeV);
         hiter->setBg(bg);
      }
   }
 
   // post hodoscope hit to the hodoscope hits tree

   if (hodo_channel > -1) {
      hddm_s::HodoChannelList counters = tagger.getHodoChannels();
      hddm_s::HodoChannelList::iterator citer;
      for (citer = counters.begin(); citer != counters.end(); ++citer) {
         if (citer->getCounterId() == hodo_channel)
            break;
      }
      if (citer == counters.end()) {
         counters = tagger.addHodoChannels();
         citer = counters.begin();
         citer->setCounterId(hodo_channel);
         citer->setE(E/GeV);
      }
      hddm_s::TaggerTruthHitList hits = citer->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator hiter;
      int hit = 0;
      for (hiter = hits.begin(); hiter != hits.end(); ++hiter, ++hit) {
         if (fabs(t - hiter->getT()*ns) < HODO_TWO_HIT_TIME_RESOL)
            break;
         else if (t > hiter->getT()*ns) {
            hits = citer->addTaggerTruthHits(1, ++hit);
            hits(0).setE(energy/GeV);
            hits(0).setT(1e99);
            hits(0).setDE(0);
            hits(0).setBg(bg);
            break;
         }
      }
      if (hiter != hits.end()) {              // merge with a prior hit
         if (t < hiter->getT()*ns) {
            hiter->setT(t/ns);
            hiter->setE(energy/GeV);
            hiter->setBg(bg);
         }
         hiter->setDE(hiter->getDE() + HODO_HIT_DE/GeV);
      }
      else {                                 // make a new hit
         hiter->setT(t/ns);
         hiter->setE(energy/GeV);
         hiter->setDE(HODO_HIT_DE/GeV);
         hiter->setBg(bg);
      }
   }
   return (micro_channel > -1 || hodo_channel > -1);
}

int GlueXPseudoDetectorTAG::addRFsync(double tsync)
{
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXPseudoDetectorTAG::addRFsync error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getRFtimes().size() == 0)
      hitview.addRFtimes();
   hddm_s::RFtime &rftime = hitview.getRFtime();
   rftime.setTsync(tsync);
   return 1;
}
