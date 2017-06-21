//
// GlueXUserEventInformation - class implementation
//
// author: richard.t.jones at uconn.edu
// version: september 24, 2016

#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "HddmOutput.hh"
#include "HddsG4Builder.hh"
#include "Randomize.hh"

int GlueXUserEventInformation::fWriteNoHitEvents = 0;

GlueXUserEventInformation::GlueXUserEventInformation(hddm_s::HDDM *hddmevent)
 : fKeepEvent(true),
   fNprimaries(0),
   fNvertices(0)
{
   if (hddmevent == 0) {
      fOutputRecord = new hddm_s::HDDM();
      fOutputRecord->addPhysicsEvents();
   }
   else {
      fOutputRecord = hddmevent;
      fNprimaries = fOutputRecord->getProducts().size();
      fNvertices = fOutputRecord->getVertices().size();
   }
   SetRandomSeeds();
}

GlueXUserEventInformation::~GlueXUserEventInformation()
{
   static bool geometry_record_written = false;
   if (fOutputRecord != 0) {
      if (! geometry_record_written) {
         hddm_s::GeometryList geom = fOutputRecord->addGeometrys();
         geom(0).setMd5simulation(last_md5_checksum);
         geometry_record_written = true;
      }
      if (fKeepEvent) {
         hddm_s::PhysicsEventList pev = fOutputRecord->getPhysicsEvents();
         int runno = HddmOutput::getRunNo();
         if (runno > 0)
            pev(0).setRunNo(runno);
         if (pev(0).getEventNo() == 0) {
            pev(0).setEventNo(HddmOutput::incrementEventNo());
         }
         if (fWriteNoHitEvents || pev(0).getHitViews().size() > 0) {
            HddmOutput::WriteOutputHDDM(*fOutputRecord);
         }
      }
      delete fOutputRecord;
   }
}

void GlueXUserEventInformation::AddBeamParticle(int geanttype, double t0,
                                                const G4ThreeVector &pos, 
                                                const G4ThreeVector &mom)
{
   hddm_s::PhysicsEventList pev = fOutputRecord->getPhysicsEvents();
   hddm_s::ReactionList rea = pev(0).getReactions();
   if (rea.size() == 0) {
      rea = pev(0).addReactions();
   }
   hddm_s::BeamList beam = rea(0).addBeams();
   beam(0).setType((Particle_t)geanttype);
   hddm_s::MomentumList pmo = beam(0).addMomenta();
   pmo(0).setPx(mom[0]/GeV);
   pmo(0).setPy(mom[1]/GeV);
   pmo(0).setPz(mom[2]/GeV);
   double mass = GlueXPrimaryGeneratorAction::GetMass(geanttype);
   double E = sqrt(mass*mass + mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
   pmo(0).setE(E/GeV);
}

void GlueXUserEventInformation::AddBeamParticle(int geanttype, double t0,
                                                const G4ThreeVector &pos, 
                                                const G4ThreeVector &mom,
                                                const G4ThreeVector &pol)
{
   AddBeamParticle(geanttype, t0, pos, mom);
   hddm_s::BeamList beam = fOutputRecord->getBeams();
   hddm_s::PolarizationList polar = beam(0).addPolarizations();
   polar(0).setPx(pol[0]);
   polar(0).setPy(pol[1]);
   polar(0).setPz(pol[2]);
}

void GlueXUserEventInformation::AddTargetParticle(int geanttype, double t0,
                                                  const G4ThreeVector &pos,
                                                  const G4ThreeVector &mom)
{
   hddm_s::PhysicsEventList pev = fOutputRecord->getPhysicsEvents();
   hddm_s::ReactionList rea = pev(0).getReactions();
   if (rea.size() == 0) {
      rea = pev(0).addReactions();
   }
   hddm_s::TargetList beam = rea(0).addTargets();
   beam(0).setType((Particle_t)geanttype);
   hddm_s::MomentumList pmo = beam(0).addMomenta();
   pmo(0).setPx(mom[0]/GeV);
   pmo(0).setPy(mom[1]/GeV);
   pmo(0).setPz(mom[2]/GeV);
   double mass = GlueXPrimaryGeneratorAction::GetMass(geanttype);
   double E = sqrt(mass*mass + mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
   pmo(0).setE(E/GeV);
}

void GlueXUserEventInformation::AddTargetParticle(int geanttype, double t0,
                                                  const G4ThreeVector &pos, 
                                                  const G4ThreeVector &mom,
                                                  const G4ThreeVector &pol)
{
   AddTargetParticle(geanttype, t0, pos, mom);
   hddm_s::TargetList beam = fOutputRecord->getTargets();
   hddm_s::PolarizationList polar = beam(0).addPolarizations();
   polar(0).setPx(pol[0]);
   polar(0).setPy(pol[1]);
   polar(0).setPz(pol[2]);
}

void GlueXUserEventInformation::AddPrimaryVertex(const G4PrimaryVertex &vertex)
{
   hddm_s::PhysicsEventList pev = fOutputRecord->getPhysicsEvents();
   hddm_s::ReactionList rea = pev(0).getReactions();
   if (rea.size() == 0) {
      rea = pev(0).addReactions();
   }
   hddm_s::VertexList ver = rea(0).addVertices();
   hddm_s::OriginList ori = ver(0).addOrigins();
   ori(0).setVx(vertex.GetX0()/cm);
   ori(0).setVy(vertex.GetY0()/cm);
   ori(0).setVz(vertex.GetZ0()/cm);
   ori(0).setT(vertex.GetT0()/ns);
   int np = vertex.GetNumberOfParticle();
   hddm_s::ProductList pro = ver(0).addProducts(np);
   for (int ip=0; ip < np; ++ip) {
      G4PrimaryParticle *part = vertex.GetPrimary(ip);
      int gluexID = AssignNextGlueXTrackID();
      int pdgtype = part->GetPDGcode();
      int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
      pro(ip).setType((Particle_t)g3type);
      pro(ip).setPdgtype(pdgtype);
      pro(ip).setId(gluexID);
      pro(ip).setParentid(0);
      pro(ip).setMech(0);
      hddm_s::MomentumList pmo = pro(ip).addMomenta();
      double px = part->GetPx()/GeV;
      double py = part->GetPy()/GeV;
      double pz = part->GetPz()/GeV;
      double mass = part->GetMass()/GeV;
      double E = sqrt(mass*mass + px*px + py*py + pz*pz);
      pmo(0).setPx(px);
      pmo(0).setPy(py);
      pmo(0).setPz(pz);
      pmo(0).setE(E);
      double polx = part->GetPolX();
      double poly = part->GetPolY();
      double polz = part->GetPolZ();
      if (polx != 0 || poly != 0 || polz != 0) {
         hddm_s::PolarizationList polar = pro(ip).addPolarizations();
         polar(0).setPx(polx);
         polar(0).setPy(poly);
         polar(0).setPz(polz);
      }
      ++fNprimaries;
   }
   ++fNvertices;
}

void GlueXUserEventInformation::AddSecondaryVertex(
                                const G4TrackVector &secondaries,
                                int parentID, int mech)
{
   if (secondaries.size() == 0)
      return;

   hddm_s::PhysicsEventList pev = fOutputRecord->getPhysicsEvents();
   hddm_s::ReactionList rea = pev(0).getReactions();
   if (rea.size() == 0) {
      rea = pev(0).addReactions();
   }
   hddm_s::VertexList ver = rea(0).addVertices();
   hddm_s::OriginList ori = ver(0).addOrigins();
   G4ThreeVector vertex(secondaries[0]->GetPosition());
   ori(0).setVx(secondaries[0]->GetPosition()[0]/cm);
   ori(0).setVy(secondaries[0]->GetPosition()[1]/cm);
   ori(0).setVz(secondaries[0]->GetPosition()[2]/cm);
   ori(0).setT(secondaries[0]->GetGlobalTime()/ns);
   int np = secondaries.size();
   hddm_s::ProductList pro = ver(0).addProducts(np);
   for (int ip=0; ip < np; ++ip) {
       GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                  secondaries[ip]->GetUserInformation();
       if (trackinfo == 0) {
          G4cerr << "GlueXUserEventInformation::AddSecondaryVertex error - "
                 << "track found without any UserTrackInformation attached "
                 << "in secondaries list, cannot continue, aborting!"
                 << G4endl;
       }
       int gluexID = trackinfo->GetGlueXTrackID();
       int pdgtype = secondaries[ip]->GetDefinition()->GetPDGEncoding();
       int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
       pro(ip).setType((Particle_t)g3type);
       pro(ip).setPdgtype(pdgtype);
       pro(ip).setId(gluexID);
       pro(ip).setParentid(parentID);
       pro(ip).setDecayVertex(fNvertices);
       pro(ip).setMech(mech);
       hddm_s::MomentumList pmo = pro(ip).addMomenta();
       double px = secondaries[ip]->GetMomentum()[0]/GeV;
       double py = secondaries[ip]->GetMomentum()[1]/GeV;
       double pz = secondaries[ip]->GetMomentum()[2]/GeV;
       double mass = secondaries[ip]->GetDefinition()->GetPDGMass()/GeV;
       double E = sqrt(mass*mass + px*px + py*py + pz*pz);
       pmo(0).setPx(px);
       pmo(0).setPy(py);
       pmo(0).setPz(pz);
       pmo(0).setE(E);
       double polx = secondaries[ip]->GetPolarization()[0];
       double poly = secondaries[ip]->GetPolarization()[1];
       double polz = secondaries[ip]->GetPolarization()[2];
       if (polx != 0 || poly != 0 || polz != 0) {
          hddm_s::PolarizationList polar = pro(ip).addPolarizations();
          polar(0).setPx(polx);
          polar(0).setPy(poly);
          polar(0).setPz(polz);
       }
      ++fNprimaries;
   }
   ++fNvertices;
}

void GlueXUserEventInformation::AddMCtrajectoryPoint(const G4Step &step,
                                                     int save_option)
{
   // The meaning of the save_option argument is as follows.
   //  save_option = 0  don't store trajectory info
   //  save_option = 1  store birth and death points of primary tracks
   //  save_option = 2  store birth and death points of all particles
   //  save_option = 3  store full trajectory of primary tracks
   //  save_option = 4  store full trajectory of primary tracks and
   //                   birth/death points of secondaries
   //  save_option = 5  store full trajectory for all particles

   hddm_s::McTrajectoryList traj = fOutputRecord->getMcTrajectorys();
   if (traj.size() == 0) {
      hddm_s::PhysicsEventList pev = fOutputRecord->getPhysicsEvents();
      if (pev.size() == 0) 
         pev = fOutputRecord->addPhysicsEvents();
      hddm_s::HitViewList view = pev(0).getHitViews();
      if (view.size() == 0) 
         view = pev(0).addHitViews();
      traj = view(0).addMcTrajectorys();
   }

   const G4Track *track = step.GetTrack();
   GlueXUserTrackInformation *trackinfo;
   trackinfo = (GlueXUserTrackInformation*)track->GetUserInformation();
   int isBorn = (track->GetCurrentStepNumber() == 1);
   int isDead = (track->GetTrackStatus() == fStopAndKill);
   int isPrimary = (trackinfo->GetGlueXTrackID() > 0);
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
   int itrack = trackinfo->GetGlueXTrackID();
   G4StepPoint *xin = step.GetPreStepPoint();
   G4StepPoint *xout = step.GetPostStepPoint();
   double radlen = xin->GetMaterial()->GetRadlen();
   if ((save_option == 1 && isPrimary && isBorn) ||
       (save_option == 2 && isBorn) ||
       (save_option == 3 && isPrimary) ||
       (save_option == 4 && (isPrimary || isBorn)) ||
       (save_option == 5 && isBorn))
   {
      G4ThreeVector mom(xin->GetMomentum());
      G4ThreeVector pos(xin->GetPosition());
      hddm_s::McTrajectoryPointList trajpt = traj(0).addMcTrajectoryPoints();
      trajpt(0).setStep(0);
      trajpt(0).setDE(0);
      trajpt(0).setE(xin->GetKineticEnergy()/GeV);
      trajpt(0).setMech(0);
      trajpt(0).setPart(g3type);
      trajpt(0).setPrimary_track(abs(itrack));
      trajpt(0).setPx(mom[0]/GeV);
      trajpt(0).setPy(mom[1]/GeV);
      trajpt(0).setPz(mom[2]/GeV);
      trajpt(0).setRadlen(radlen/cm);
      trajpt(0).setT(xin->GetGlobalTime()/ns*1e-9);
      trajpt(0).setTrack(track->GetTrackID());
      trajpt(0).setX(pos[0]/cm);
      trajpt(0).setY(pos[1]/cm);
      trajpt(0).setZ(pos[2]/cm);
   }
   if ((save_option == 1 && isPrimary && isDead) ||
       (save_option == 2 && isDead) ||
       (save_option == 3 && isPrimary) ||
       (save_option == 4 && (isPrimary || isDead)) ||
       (save_option == 5))
   {
      G4ProcessType mechtype(xout->GetProcessDefinedStep()->GetProcessType());
      int mech4c[2];
      switch (mechtype) {
         case fTransportation:
            snprintf((char*)mech4c, 5, "GEOM");
            break;
         case fElectromagnetic:
            snprintf((char*)mech4c, 5, "ELEC");
            break;
         case fOptical:
            snprintf((char*)mech4c, 5, "OPTI");
            break;
         case fHadronic:
            snprintf((char*)mech4c, 5, "HADR");
            break;
         case fPhotolepton_hadron:
            snprintf((char*)mech4c, 5, "PHOT");
            break;
         case fDecay:
            snprintf((char*)mech4c, 5, "DCAY");
            break;
         case fGeneral:
            snprintf((char*)mech4c, 5, "GENE");
            break;
         case fParameterisation:
            snprintf((char*)mech4c, 5, "PARA");
            break;
         case fUserDefined:
            snprintf((char*)mech4c, 5, "USER");
            break;
         default:
            snprintf((char*)mech4c, 5, "STOP");
            break;
      }
      G4ThreeVector mom(xout->GetMomentum());
      G4ThreeVector pos(xout->GetPosition());
      double steplen;
      if (save_option == 1 || save_option == 2 ||
          (save_option == 4 && !isPrimary))
      {
         steplen = track->GetTrackLength();
      }
      else {
         steplen = step.GetStepLength();
      }
      hddm_s::McTrajectoryPointList trajpt = traj(0).addMcTrajectoryPoints();
      trajpt(0).setStep(steplen/cm);
      trajpt(0).setDE(step.GetTotalEnergyDeposit()/GeV);
      trajpt(0).setE(xout->GetKineticEnergy()/GeV);
      trajpt(0).setMech(mech4c[0]);
      trajpt(0).setPart(g3type);
      trajpt(0).setPrimary_track(abs(itrack));
      trajpt(0).setPx(mom[0]/GeV);
      trajpt(0).setPy(mom[1]/GeV);
      trajpt(0).setPz(mom[2]/GeV);
      trajpt(0).setRadlen(radlen/cm);
      trajpt(0).setT(xout->GetGlobalTime()/ns*1e-9);
      trajpt(0).setTrack(track->GetTrackID());
      trajpt(0).setX(pos[0]/cm);
      trajpt(0).setY(pos[1]/cm);
      trajpt(0).setZ(pos[2]/cm);
   }
}

void GlueXUserEventInformation::SetRandomSeeds()
{
   hddm_s::PhysicsEventList pev = fOutputRecord->getPhysicsEvents();
   hddm_s::ReactionList rea = pev(0).getReactions();
   if (rea.size() == 0) {
      rea = pev(0).addReactions();
   }
   hddm_s::RandomList rnd = rea(0).getRandoms();
   if (rnd.size() > 0) {
      long int seed[2];
      seed[0] = rnd(0).getSeed1();
      seed[1] = rnd(0).getSeed2();
      G4Random::setTheSeeds(seed);
   }
   else {
      const long int *seed = G4Random::getTheSeeds();
      rnd = rea(0).addRandoms();
      rnd(0).setSeed1(seed[0]);
      rnd(0).setSeed2(seed[1]);
      rnd(0).setSeed3(709975946 + pev(0).getEventNo());
      rnd(0).setSeed4(912931182 + pev(0).getEventNo());
#if VERBOSE_RANDOMS
      G4cout << "New event with starting seeds " 
             << seed[0] << ", " << seed[1] << G4endl;
#endif
   }
}

double GlueXUserEventInformation::GetBeamPhotonEnergy()
{
   hddm_s::BeamList beam = fOutputRecord->getBeams();
   if (beam.size() > 0)
      return beam(0).getMomentum().getE()*GeV;
   return 0;
}

int GlueXUserEventInformation::AssignNextGlueXTrackID(const G4Track* track)
{
   int lastG4 = 0;
   int lastGX = 0;
   std::map<int,int>::iterator iter;
   for (iter = fGlueXTrackID.begin(); iter != fGlueXTrackID.end(); ++iter) {
      lastG4 = (iter->first > lastG4)? iter->first : lastG4;
      lastGX = (iter->second > lastGX)? iter->second : lastGX;
   }
   if (track == 0) {
      fGlueXTrackID[++lastG4] = ++lastGX;
   }
   else {
      fGlueXTrackID[track->GetTrackID()] = ++lastGX;
   }
   return lastGX;
}

int GlueXUserEventInformation::GetGlueXTrackID(int g4ID)
{
   if (fGlueXTrackID.find(g4ID) != fGlueXTrackID.end())
      return fGlueXTrackID[g4ID];
   return 0;
}

int GlueXUserEventInformation::GetGlueXTrackID(const G4Track *track)
{
   int g4Id = track->GetTrackID();
   if (fGlueXTrackID.find(g4Id) != fGlueXTrackID.end())
      return fGlueXTrackID[g4Id];
   int parId = track->GetParentID();
   if (fGlueXTrackID.find(parId) != fGlueXTrackID.end())
      return fGlueXTrackID[parId] * 1000000 + g4Id;
   return 1000000 + g4Id;
}

void GlueXUserEventInformation::SetGlueXTrackID(int g4ID, int gluexID)
{
   fGlueXTrackID[g4ID] = gluexID;
}

void GlueXUserEventInformation::Print() const
{
   G4cout << "GlueXUserEventInformation: hddm_s=" << fOutputRecord
          << G4endl;
}

int GlueXUserEventInformation::AssignBCALincidentID(const G4Track* track)
{
   BCALincidentParticle ipart;
   ipart.pos = track->GetPosition();
   ipart.mom = track->GetMomentum();
   ipart.E = track->GetTotalEnergy();
   ipart.ptype = track->GetDynamicParticle()->GetPDGcode();
   ipart.trackID = track->GetTrackID();
   fBCALincidentParticle.push_back(ipart);
   return fBCALincidentParticle.size();
}

const BCALincidentParticle *GlueXUserEventInformation::
                            GetBCALincidentParticle(int incidentID)
{
   if (incidentID < 0 || incidentID >= (int)fBCALincidentParticle.size())
      return 0;
   return &fBCALincidentParticle[incidentID];
}
