//
// GlueXHitDIRCWob - class implementation
//
// created on: 05.04.2017
// original author: r.dzhygadlo at gsi.de

#include "GlueXHitDIRCWob.hh"

G4ThreadLocal G4Allocator<GlueXHitDIRCWob>* GlueXHitDIRCWobAllocator = 0;

GlueXHitDIRCWob::GlueXHitDIRCWob():G4VHit()
{
}

void GlueXHitDIRCWob::Draw() const
{
   // not yet implemented
}

void GlueXHitDIRCWob::Print() const
{
   G4cout << "GlueXHitDIRCWob:" << G4endl
          << "   track = " << track << G4endl
          << "   normalId = " << normalId << " " << G4endl
          << G4endl;
}
