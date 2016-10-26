//
// GlueXHitFTOFbar - class header
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitFTOFbar_h
#define GlueXHitFTOFbar_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"
#include "G4SystemOfUnits.hh"

class GlueXHitFTOFbar : public G4VHit
{
 public:
   GlueXHitFTOFbar(G4int plane, G4int bar, G4int end);
   int operator==(const GlueXHitFTOFbar &right) const;
   GlueXHitFTOFbar &operator+=(const GlueXHitFTOFbar &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int plane_;           // =0 (downstream, vertical) or =1 (up, horizontal)
   G4int bar_;             // bar number, bottom-top, south-north (see hdds)
   G4int end_;             // =0 (top, north) or =1 (bottom, south)

   struct hitinfo_t {
      G4double dE_GeV;     // energy deposition (GeV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double itrack_;    // track index of first particle making this hit
      G4double ptype_G3;   // G3 type of first particle making this hit
      G4double px_GeV;     // hit track momentum x component (GeV)
      G4double py_GeV;     // hit track momentum y component (GeV)
      G4double pz_GeV;     // hit track momentum z component (GeV)
      G4double E_GeV;      // hit track total energy (GeV)
      G4double x_cm;       // x coordinate of the hit in global refsys (cm)
      G4double y_cm;       // x coordinate of the hit in global refsys (cm)
      G4double z_cm;       // z coordinate of the hit in global refsys (cm)
      G4double dist_cm;    // distance of hit (from end?) in cm
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(plane_, bar_, end_); }
   static G4int GetKey(G4int plane, G4int bar, G4int end) {
      return (plane << 20) + (bar << 10) + end;
   }
};

typedef G4THitsMap<GlueXHitFTOFbar> GlueXHitsMapFTOFbar;

extern G4ThreadLocal G4Allocator<GlueXHitFTOFbar>* GlueXHitFTOFbarAllocator;

inline void* GlueXHitFTOFbar::operator new(size_t)
{
   if (!GlueXHitFTOFbarAllocator)
      GlueXHitFTOFbarAllocator = new G4Allocator<GlueXHitFTOFbar>;
   return (void *) GlueXHitFTOFbarAllocator->MallocSingle();
}

inline void GlueXHitFTOFbar::operator delete(void *aHit)
{
   GlueXHitFTOFbarAllocator->FreeSingle((GlueXHitFTOFbar*) aHit);
}

#endif
