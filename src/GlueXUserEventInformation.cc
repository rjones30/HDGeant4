//
// GlueXUserEventInformation - class implementation
//
// author: richard.t.jones at uconn.edu
// version: september 24, 2016

#include "GlueXUserEventInformation.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserOptions.hh"
#include "HddmOutput.hh"
#include "HddsG4Builder.hh"

GlueXUserEventInformation::GlueXUserEventInformation(hddm_s::HDDM *hddmevent)
 : fKeepEvent(true)
{
   if (hddmevent == 0) {
      fOutputRecord = new hddm_s::HDDM();
      fOutputRecord->addPhysicsEvents();
      fNprimaries = 0;
   }
   else {
      fOutputRecord = hddmevent;
      hddm_s::ProductList products = fOutputRecord->getProducts();
      fNprimaries = products.size();
   }
   hddm_s::PhysicsEventList pev = fOutputRecord->getPhysicsEvents();
   pev(0).setRunNo(HddmOutput::getRunNo());
   if (pev(0).getEventNo() == 0) {
      pev(0).setEventNo(HddmOutput::incrementEventNo());
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
         HddmOutput::WriteOutputHDDM(*fOutputRecord);
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
       int pdgtype = part->GetPDGcode();
       int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
       pro(ip).setType((Particle_t)g3type);
       pro(ip).setPdgtype(pdgtype);
       pro(ip).setId(ip);
       pro(ip).setParentid(0);
       pro(ip).setMech(0);
       hddm_s::MomentumList pmo = pro(ip).addMomenta();
       double px = part->GetPx()/GeV;
       double py = part->GetPy()/GeV;
       double pz = part->GetPz()/GeV;
       double mass = part->GetMass()/GeV;
       double E = sqrt(mass*mass + px*px + py*py + pz*pz);
       pmo(ip).setPx(px);
       pmo(ip).setPy(py);
       pmo(ip).setPz(pz);
       pmo(ip).setE(E);
       double polx = part->GetPolX();
       double poly = part->GetPolY();
       double polz = part->GetPolZ();
       if (polx != 0 || poly != 0 || polz != 0) {
          hddm_s::PolarizationList polar = pro(ip).addPolarizations();
          polar(0).setPx(px);
          polar(0).setPy(py);
          polar(0).setPz(pz);
       }
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

void GlueXUserEventInformation::Print() const
{
   G4cout << "GlueXUserEventInformation: hddm_s=" << fOutputRecord
          << G4endl;
}
