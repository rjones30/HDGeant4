//
// GlueXHitFTOFbar - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016

#include "GlueXHitFTOFbar.hh"

G4ThreadLocal G4Allocator<GlueXHitFTOFbar>* GlueXHitFTOFbarAllocator = 0;

GlueXHitFTOFbar::GlueXHitFTOFbar(G4int plane, G4int bar)
 : G4VHit(),
   plane_(plane),
   bar_(bar)
{}

int GlueXHitFTOFbar::operator==(const GlueXHitFTOFbar &right) const
{
   if (plane_ !=  right.plane_ || bar_ != right.bar_ )
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].end_     == right.hits[ih].end_    &&
          hits[ih].dE_GeV   == right.hits[ih].dE_GeV  &&
          hits[ih].t_ns     == right.hits[ih].t_ns    &&
          hits[ih].itrack_  == right.hits[ih].itrack_ &&
          hits[ih].ptype_G3 == right.hits[ih].ptype_G3 &&
          hits[ih].px_GeV   == right.hits[ih].px_GeV  &&
          hits[ih].py_GeV   == right.hits[ih].py_GeV  &&
          hits[ih].pz_GeV   == right.hits[ih].pz_GeV  &&
          hits[ih].E_GeV    == right.hits[ih].E_GeV   &&
          hits[ih].x_cm     == right.hits[ih].x_cm    &&
          hits[ih].y_cm     == right.hits[ih].y_cm    &&
          hits[ih].z_cm     == right.hits[ih].z_cm    &&
          hits[ih].dist_cm  == right.hits[ih].dist_cm )
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitFTOFbar &GlueXHitFTOFbar::operator+=(const GlueXHitFTOFbar &right)
{
   if (plane_ !=  right.plane_ || bar_ != right.bar_) {
      G4cerr << "Error in GlueXHitFTOFbar::operator+=() - "
             << "illegal attempt to merge hits from two different bars!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitFTOFbar::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitFTOFbar::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; hiter != hits.end(); ++hiter) {
         if (hiter->t_ns > hitsrc->t_ns)
            break;
      }
      hiter = hits.insert(hiter, GlueXHitFTOFbar::hitinfo_t());
      hiter->end_   = hitsrc->end_;
      hiter->dE_GeV = hitsrc->dE_GeV;
      hiter->t_ns = hitsrc->t_ns;
      hiter->itrack_ = hitsrc->itrack_;
      hiter->ptype_G3 = hitsrc->ptype_G3;
      hiter->px_GeV = hitsrc->px_GeV;
      hiter->py_GeV = hitsrc->py_GeV;
      hiter->pz_GeV = hitsrc->pz_GeV;
      hiter->E_GeV = hitsrc->E_GeV;
      hiter->x_cm = hitsrc->x_cm;
      hiter->y_cm = hitsrc->y_cm;
      hiter->z_cm = hitsrc->z_cm;
      hiter->dist_cm = hitsrc->dist_cm;
   }
   return *this;
}

void GlueXHitFTOFbar::Draw() const
{
   // not yet implemented
}

void GlueXHitFTOFbar::Print() const
{
   G4cout << "GlueXHitFTOFbar: "
          << "   plane = " << plane_ 
          << ",  bar = " << bar << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   dE = " << hiter->dE_GeV << " GeV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   itrack = " << hiter->itrack_ << G4endl
             << "   ptype = " << hiter->ptype_G3 << G4endl
             << "   px = " << hiter->px_GeV << " GeV/c" << G4endl
             << "   py = " << hiter->py_GeV << " GeV/c" << G4endl
             << "   pz = " << hiter->pz_GeV << " GeV/c" << G4endl
             << "   E = " << hiter->E_GeV << " GeV" << G4endl
             << "   x = " << hiter->x_cm << " cm" << G4endl
             << "   y = " << hiter->y_cm << " cm" << G4endl
             << "   z = " << hiter->z_cm << " cm" << G4endl
             << "   dist = " << hiter->dist_cm << " cm" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapFTOFbar *hitsmap)
{
   std::map<int, GlueXHitFTOFbar*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitFTOFbar*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}
