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

GlueXUserEventInformation::GlueXUserEventInformation(int geanttype, double t0,
                                                     const G4ThreeVector &pos, 
                                                     const G4ThreeVector &mom)
 : fKeepEvent(true)
{
   fOutputRecord = new hddm_s::HDDM();
   hddm_s::PhysicsEventList pev = fOutputRecord->addPhysicsEvents();
   pev(0).setRunNo(HddmOutput::getRunNo());
   pev(0).setEventNo(HddmOutput::incrementEventNo());
   hddm_s::ReactionList rea = pev(0).addReactions();
   hddm_s::VertexList ver = rea(0).addVertices();
   hddm_s::OriginList ori = ver(0).addOrigins();
   ori(0).setVx(pos[0]/cm);
   ori(0).setVy(pos[1]/cm);
   ori(0).setVz(pos[2]/cm);
   ori(0).setT(t0/ns);
   hddm_s::ProductList pro = ver(0).addProducts();
   pro(0).setType((Particle_t)geanttype);
   pro(0).setPdgtype(GlueXPrimaryGeneratorAction::ConvertGeant3ToPdg(geanttype));
   pro(0).setId(1);
   pro(0).setParentid(0);
   pro(0).setMech(0);
   hddm_s::MomentumList pmo = pro(0).addMomenta();
   pmo(0).setPx(mom[0]/GeV);
   pmo(0).setPy(mom[1]/GeV);
   pmo(0).setPz(mom[2]/GeV);
   double mass = GlueXPrimaryGeneratorAction::GetMass(geanttype);
   double E = sqrt(mass*mass + mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
   pmo(0).setE(E/GeV);
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

void GlueXUserEventInformation::SetRandomSeeds()
{
   hddm_s::PhysicsEventList pev = fOutputRecord->getPhysicsEvents();
   hddm_s::ReactionList rea = pev(0).getReactions();
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
