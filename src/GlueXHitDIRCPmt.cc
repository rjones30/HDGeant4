//
// GlueXHitDIRCPmt.cc
//
// created on: 05.04.2017
// original author: r.dzhygadlo at gsi.de
// contributors:

#include "GlueXHitDIRCPmt.hh"

G4ThreadLocal G4Allocator<GlueXHitDIRCPmt>* GlueXHitDIRCPmtAllocator = 0;

GlueXHitDIRCPmt::GlueXHitDIRCPmt(const GlueXHitDIRCPmt &src)
{
  E_GeV = src.E_GeV;
  t_ns = src.t_ns;
  t_fixed_ns = src.t_fixed_ns;
  x_cm = src.x_cm;
  y_cm = src.y_cm;
  z_cm = src.z_cm;
  ch = src.ch;
  key_bar = src.key_bar;
  path = src.path;
  refl = src.refl;
  bbrefl = src.bbrefl;
}

void GlueXHitDIRCPmt::Draw() const
{
  // not yet implemented
}

void GlueXHitDIRCPmt::Print() const
{
  G4cout << "   E = " << E_GeV << " GeV" << G4endl
	 << "   t = " << t_ns << " ns" << G4endl
	 << "   t_fixed = " << t_fixed_ns << " ns" << G4endl
	 << "   x = " << x_cm << " cm" << G4endl
	 << "   y = " << y_cm << " cm" << G4endl
	 << "   z = " << z_cm << " cm" << G4endl
	 << "   path = " << path << " " << G4endl
	 << "   refl = " << refl << " " << G4endl
	 << "   bbrefl = " << bbrefl << " " << G4endl
	 << "   ch = " << ch << " " << G4endl
	 << "   key_bar = " << key_bar << " " << G4endl
	 << G4endl;
}
