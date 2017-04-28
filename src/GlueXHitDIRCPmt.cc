//
// GlueXHitDIRCPmt.cc
//
// created on: 05.04.2017
// original author: r.dzhygadlo at gsi.de
// contributors:

#include "GlueXHitDIRCPmt.hh"

G4ThreadLocal G4Allocator<GlueXHitDIRCPmt>* GlueXHitDIRCPmtAllocator = 0;

GlueXHitDIRCPmt::GlueXHitDIRCPmt()
  : G4VHit()
{}

void GlueXHitDIRCPmt::Draw() const
{
  // not yet implemented
}

void GlueXHitDIRCPmt::Print() const
{
  G4cout << "   E = " << E_GeV << " GeV" << G4endl
	 << "   t = " << t_ns << " ns" << G4endl
	 << "   x = " << x_cm << " cm" << G4endl
	 << "   y = " << y_cm << " cm" << G4endl
	 << "   z = " << z_cm << " cm" << G4endl
	 << "   ch = " << ch << " " << G4endl
	 << "   key_bar = " << key_bar << " " << G4endl
	 << G4endl;
}
