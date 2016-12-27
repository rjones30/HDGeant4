//
// GlueXHitPSCpoint - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016

#include "GlueXHitPSCpoint.hh"

G4ThreadLocal G4Allocator<GlueXHitPSCpoint>* GlueXHitPSCpointAllocator = 0;

int GlueXHitPSCpoint::operator==(const GlueXHitPSCpoint &right) const
{
   if (arm_        == right.arm_        &&
       module_     == right.module_     &&
       E_GeV       == right.E_GeV       &&
       dEdx_GeV_cm == right.dEdx_GeV_cm &&
       primary_    == right.primary_    &&
       ptype_G3    == right.ptype_G3    &&
       px_GeV      == right.px_GeV      &&
       py_GeV      == right.py_GeV      &&
       pz_GeV      == right.pz_GeV      &&
       x_cm        == right.x_cm        &&
       y_cm        == right.y_cm        &&
       z_cm        == right.z_cm        &&
       t_ns        == right.t_ns        &&
       track_      == right.track_      &&
       trackID_    == right.trackID_    )
   {
      return 1;
   }
   return 0;
}

GlueXHitPSCpoint &GlueXHitPSCpoint::operator+=(const GlueXHitPSCpoint &right)
{
   G4cerr << "Error in GlueXHitPSCpoint::operator+= - "
          << "illegal attempt to merge two TruthPoint objects in the psc!"
          << G4endl;
   return *this;
}

void GlueXHitPSCpoint::Draw() const
{
   // not yet implemented
}

void GlueXHitPSCpoint::Print() const
{
   G4cout << "GlueXHitPSCpoint:" << G4endl
          << "   arm = " << arm_ << G4endl
          << "   module = " << module_ << G4endl
          << "   track = " << track_ << G4endl
          << "   trackID = " << trackID_ << G4endl
          << "   E = " << E_GeV << " GeV" << G4endl
          << "   dEdx = " << dEdx_GeV_cm << " GeV/cm" << G4endl
          << "   primary = " << primary_ << G4endl
          << "   ptype = " << ptype_G3 << G4endl
          << "   px = " << px_GeV << " GeV/c" << G4endl
          << "   py = " << py_GeV << " GeV/c" << G4endl
          << "   pz = " << pz_GeV << " GeV/c" << G4endl
          << "   x = " << x_cm << " cm" << G4endl
          << "   y = " << y_cm << " cm" << G4endl
          << "   z = " << z_cm << " cm" << G4endl
          << "   t = " << t_ns << " ns" << G4endl
          << G4endl;
}

void printallhits(GlueXHitsMapPSCpoint *hitsmap)
{
   std::map<int, GlueXHitPSCpoint*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitPSCpoint*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}
