//
// GlueXHitCDCstraw - class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 28, 2015

#include "GlueXHitCDCstraw.hh"

G4ThreadLocal G4Allocator<GlueXHitCDCstraw>* GlueXHitCDCstrawAllocator = 0;

GlueXHitCDCstraw::GlueXHitCDCstraw(G4int ring, G4int sector)
 : G4VHit(),
   ring_(ring),
   sector_(sector)
{}

int GlueXHitCDCstraw::operator==(const GlueXHitCDCstraw &right) const
{
   if (ring_ != right.ring_  || sector_ !=  right.sector_)
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].q_fC     == right.hits[ih].q_fC    &&
          hits[ih].t_ns     == right.hits[ih].t_ns    &&
          hits[ih].d_cm     == right.hits[ih].d_cm    &&
          hits[ih].itrack_  == right.hits[ih].itrack_ &&
          hits[ih].t0_ns    == right.hits[ih].t0_ns   &&
          hits[ih].z_cm     == right.hits[ih].z_cm    &&
          hits[ih].ptype_G3 == right.hits[ih].ptype_G3)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitCDCstraw &GlueXHitCDCstraw::operator+=(const GlueXHitCDCstraw &right)
{
   if (ring_ != right.ring_  || sector_ !=  right.sector_) {
      G4cerr << "Error in GlueXHitCDCstraw::operator+=() - "
             << "illegal attempt to merge hits from two different straws!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitCDCstraw::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitCDCstraw::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; hiter != hits.end(); ++hiter) {
         if (hiter->t_ns > hitsrc->t_ns)
            break;
      }
      hiter = hits.insert(hiter, GlueXHitCDCstraw::hitinfo_t());
      hiter->q_fC = hitsrc->q_fC;
      hiter->t_ns = hitsrc->t_ns;
      hiter->d_cm = hitsrc->d_cm;
      hiter->itrack_ = hitsrc->itrack_;
      hiter->ptype_G3 = hitsrc->ptype_G3;
      hiter->t0_ns = hitsrc->t0_ns;
      hiter->z_cm = hitsrc->z_cm;
   }
   return *this;
}

void GlueXHitCDCstraw::Draw() const
{
   // not yet implemented
}

void GlueXHitCDCstraw::Print() const
{
   G4cout << "GlueXHitCDCstraw: "
          << "   ring = " << ring_
          << "   sector = " << sector_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   q = " << hiter->q_fC << " fC"
             << "   t = " << hiter->t_ns << " ns"
             << "   d = " << hiter->d_cm << " cm"
             << "   itrack = " << hiter->itrack_
             << "   ptype = " << hiter->ptype_G3
             << "   t0 = " << hiter->t0_ns << " ns"
             << "   z = " << hiter->z_cm << " cm"
             << G4endl;
   }
}

void printallhits(GlueXHitsMapCDCstraw *hitsmap)
{
   std::map<int, GlueXHitCDCstraw*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitCDCstraw*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}
