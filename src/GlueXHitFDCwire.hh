//
// GlueXHitFDCwire - class header
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitFDCwire_h
#define GlueXHitFDCwire_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"
#include "G4SystemOfUnits.hh"

class GlueXHitFDCwire : public G4VHit
{
 public:
   GlueXHitFDCwire(G4int chamber, G4int wire);
   int operator==(const GlueXHitFDCwire &right) const;
   GlueXHitFDCwire &operator+=(const GlueXHitFDCwire &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int chamber_;            // FDC chamber, count from upstream starting at 1
   G4int wire_;               // wire number, ordered by u starting at 1

   struct hitinfo_t {
      G4double dE_arb;     // pulse integral (arb. units)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double d_cm;       // distance (cm) of closest cluster to the wire
      G4double itrack_;    // track index of first particle making this hit
      G4double ptype_G3;   // G3 type of first particle making this hit
      G4double t0_ns;      // time of passage of the track making this hit
      G4double v_cm;       // v coordinate of the hit in local refsys
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(chamber_, wire_); }
   static G4int GetKey(G4int chamber, G4int wire) {
      return (chamber << 20) + (2 << 10) + wire;
   }
};

typedef G4THitsMap<GlueXHitFDCwire> GlueXHitsMapFDCwire;

extern G4ThreadLocal G4Allocator<GlueXHitFDCwire>* GlueXHitFDCwireAllocator;

inline void* GlueXHitFDCwire::operator new(size_t)
{
   if (!GlueXHitFDCwireAllocator)
      GlueXHitFDCwireAllocator = new G4Allocator<GlueXHitFDCwire>;
   return (void *) GlueXHitFDCwireAllocator->MallocSingle();
}

inline void GlueXHitFDCwire::operator delete(void *aHit)
{
   GlueXHitFDCwireAllocator->FreeSingle((GlueXHitFDCwire*) aHit);
}

#endif
