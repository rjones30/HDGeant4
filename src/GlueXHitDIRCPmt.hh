//
// GlueXHitDIRCPmt.h
//
// created on: 05.04.2017
// original author: r.dzhygadlo at gsi.de

#ifndef GlueXHitDIRCPmt_h
#define GlueXHitDIRCPmt_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitDIRCPmt : public G4VHit
{
public:
  GlueXHitDIRCPmt() {}
  GlueXHitDIRCPmt(const GlueXHitDIRCPmt &src);
  
  void *operator new(size_t);
  void operator delete(void *aHit);

  void Draw() const;
  void Print() const;

  G4double E_GeV;      // photon energy [GeV]
  G4double t_ns;       // detection time [ns]
  G4double t_fixed_ns; // fixed pathlength time [ns]
  G4double x_cm;       // x coordinate where hit was created
  G4double y_cm;       // y coordinate where hit was created
  G4double z_cm;       // z coordinate where hit was created
  int64_t  path;       // photon's path id in the optical box
  G4int    refl;       // number of reflections in the oprical box
  G4bool   bbrefl;     // reflected off far end mirror of bar box
  G4int    ch;         // PMT channel of the hit
  G4int    key_bar;    // key of the corresponding bar hit
  G4int    track;      // index of the MC track
};

typedef G4THitsMap<GlueXHitDIRCPmt> GlueXHitsMapDIRCPmt;

extern G4ThreadLocal G4Allocator<GlueXHitDIRCPmt>* GlueXHitDIRCPmtAllocator;

inline void* GlueXHitDIRCPmt::operator new(size_t)
{
  if (!GlueXHitDIRCPmtAllocator)
    GlueXHitDIRCPmtAllocator = new G4Allocator<GlueXHitDIRCPmt>;
  return (void *) GlueXHitDIRCPmtAllocator->MallocSingle();
}

inline void GlueXHitDIRCPmt::operator delete(void *aHit)
{
  GlueXHitDIRCPmtAllocator->FreeSingle((GlueXHitDIRCPmt*) aHit);
}

#endif
