//
// GlueXHitDIRCWob - class implementation
//
// created on: 05.04.2017
// original author: r.dzhygadlo at gsi.de

#ifndef GlueXHitDIRCWob_h
#define GlueXHitDIRCWob_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitDIRCWob : public G4VHit
{
public:
  GlueXHitDIRCWob();

  void *operator new(size_t);
  void operator delete(void *aHit);

  void Draw() const;
  void Print() const;

  G4double normalId;   // uniq identifier of the normal
  G4int track;         // index of the MC track 
};

typedef G4THitsMap<GlueXHitDIRCWob> GlueXHitsMapDIRCWob;

extern G4ThreadLocal G4Allocator<GlueXHitDIRCWob>* GlueXHitDIRCWobAllocator;

inline void* GlueXHitDIRCWob::operator new(size_t)
{
  if (!GlueXHitDIRCWobAllocator)
    GlueXHitDIRCWobAllocator = new G4Allocator<GlueXHitDIRCWob>;
  return (void *) GlueXHitDIRCWobAllocator->MallocSingle();
}

inline void GlueXHitDIRCWob::operator delete(void *aHit)
{
  GlueXHitDIRCWobAllocator->FreeSingle((GlueXHitDIRCWob*) aHit);
}

#endif
