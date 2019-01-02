//
// GlueXHitDIRCBar - class implementation
//
// created on: 05.04.2017
// original author: r.dzhygadlo at gsi.de

#ifndef GlueXHitDIRCBar_h
#define GlueXHitDIRCBar_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitDIRCBar : public G4VHit
{
public:
  GlueXHitDIRCBar() {}
  GlueXHitDIRCBar(const GlueXHitDIRCBar &src);

  void *operator new(size_t);
  void operator delete(void *aHit);

  void Draw() const;
  void Print() const;

  G4double E_GeV;      // track particle total energy (GeV)
  G4double t_ns;       // pulse leading-edge time (ns)
  G4double x_cm;       // x coordinate where ch. track hits the radiator
  G4double y_cm;        
  G4double z_cm;        
  G4double px_GeV;     // px component of the track momentum
  G4double py_GeV;     
  G4double pz_GeV;
  G4int pdg;           // PDG of the particle
  G4int bar;           // index of the bar
  G4int track;         // index of the MC track 
};

typedef G4THitsMap<GlueXHitDIRCBar> GlueXHitsMapDIRCBar;

extern G4ThreadLocal G4Allocator<GlueXHitDIRCBar>* GlueXHitDIRCBarAllocator;

inline void* GlueXHitDIRCBar::operator new(size_t)
{
  if (!GlueXHitDIRCBarAllocator)
    GlueXHitDIRCBarAllocator = new G4Allocator<GlueXHitDIRCBar>;
  return (void *) GlueXHitDIRCBarAllocator->MallocSingle();
}

inline void GlueXHitDIRCBar::operator delete(void *aHit)
{
  GlueXHitDIRCBarAllocator->FreeSingle((GlueXHitDIRCBar*) aHit);
}

#endif
