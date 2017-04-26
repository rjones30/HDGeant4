//
// GlueXHitDIRCMcp.h
//
// created on: 05.04.2017
// original author: r.dzhygadlo at gsi.de

#ifndef GlueXHitDIRCMcp_h
#define GlueXHitDIRCMcp_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitDIRCMcp : public G4VHit
{
public:
  GlueXHitDIRCMcp();
  
  void *operator new(size_t);
  void operator delete(void *aHit);

  void Draw() const;
  void Print() const;

  G4double E_GeV;      // photon energy [GeV]
  G4double t_ns;       // detection time [ns]
  G4double x_cm;       // x coordinate where hit was created
  G4double y_cm;       // y coordinate where hit was created
  G4double z_cm;       // z coordinate where hit was created
  G4int    ch;         // MCP channel of the hit
  G4int    key_bar;    // key of the corresponding bar hit
    
};

typedef G4THitsMap<GlueXHitDIRCMcp> GlueXHitsMapDIRCMcp;

extern G4ThreadLocal G4Allocator<GlueXHitDIRCMcp>* GlueXHitDIRCMcpAllocator;

inline void* GlueXHitDIRCMcp::operator new(size_t)
{
  if (!GlueXHitDIRCMcpAllocator)
    GlueXHitDIRCMcpAllocator = new G4Allocator<GlueXHitDIRCMcp>;
  return (void *) GlueXHitDIRCMcpAllocator->MallocSingle();
}

inline void GlueXHitDIRCMcp::operator delete(void *aHit)
{
  GlueXHitDIRCMcpAllocator->FreeSingle((GlueXHitDIRCMcp*) aHit);
}

#endif
