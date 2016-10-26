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
   SetRandomSeeds();
}

GlueXUserEventInformation::GlueXUserEventInformation(int geanttype,
                                                     G4ThreeVector &pos, 
                                                     G4ThreeVector &mom)
 : fKeepEvent(true)
{
   fOutputRecord = new hddm_s::HDDM();
   hddm_s::PhysicsEventList pev = fOutputRecord->addPhysicsEvents();
   hddm_s::ReactionList rea = pev(0).addReactions();
   hddm_s::VertexList ver = rea(0).addVertices();
   hddm_s::OriginList ori = ver(0).addOrigins();
   ori(0).setVx(pos[0]);
   ori(0).setVy(pos[1]);
   ori(0).setVz(pos[2]);
   hddm_s::ProductList pro = ver(0).addProducts();
   pro(0).setType((Particle_t)geanttype);
   pro(0).setPdgtype(GlueXPrimaryGeneratorAction::ConvertGeant3ToPdg(geanttype));
   pro(0).setId(1);
   pro(0).setParentid(0);
   pro(0).setMech(0);
   hddm_s::MomentumList pmo = pro(0).addMomenta();
   pmo(0).setPx(mom[0]);
   pmo(0).setPy(mom[1]);
   pmo(0).setPz(mom[2]);
   double mass = GlueXPrimaryGeneratorAction::GetMass(geanttype) / GeV;
   pmo(0).setE(sqrt(mass*mass + mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]));
   SetRandomSeeds();
}

GlueXUserEventInformation::~GlueXUserEventInformation()
{
   static bool geometry_record_written = false;
   if (fOutputRecord != 0 && ! geometry_record_written) {
      hddm_s::GeometryList geom = fOutputRecord->addGeometrys();
      geom(0).setMd5simulation(last_md5_checksum);
      geometry_record_written = true;
   }

   if (fOutputRecord != 0 && fKeepEvent) {
      hddm_s::PhysicsEvent &event = fOutputRecord->getPhysicsEvent(0);
      event.setRunNo(HddmOutput::getRunNo());
      if (event.getEventNo() == 0) {
         int evno = HddmOutput::getEventNo();
         event.setEventNo(++evno);
         HddmOutput::setEventNo(evno);
      }
      HddmOutput::WriteOutputHDDM(*fOutputRecord);
      delete fOutputRecord;
   }
}

void GlueXUserEventInformation::SetRandomSeeds()
{
   hddm_s::ReactionList rea = fOutputRecord->getReactions();
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
      rnd(0).setSeed3(709975946 + HddmOutput::getEventNo());
      rnd(0).setSeed4(912931182 + HddmOutput::getEventNo());
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
