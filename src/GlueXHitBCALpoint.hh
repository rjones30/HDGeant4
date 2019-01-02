//
// GlueXHitBCALpoint - class header
//
// author: richard.t.jones at uconn.edu
// version: october 24, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitBCALpoint_h
#define GlueXHitBCALpoint_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitBCALpoint : public G4VHit
{
 public:
   GlueXHitBCALpoint() {}
   GlueXHitBCALpoint(const GlueXHitBCALpoint &src);
   int operator==(const GlueXHitBCALpoint &right) const;
   GlueXHitBCALpoint &operator+=(const GlueXHitBCALpoint &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data
 
   G4double E_GeV;       // total energy (GeV) of this shower particle at this point
   G4double phi_rad;     // phi angle (rad) of show particle at this point
   G4bool primary_;      // true if track belongs to from a primary particle
   G4int ptype_G3;       // G3 type of particle making this shower
   G4double px_GeV;      // momentum (GeV/c) of shower particle at point, x component
   G4double py_GeV;      // momentum (GeV/c) of shower particle at point, y component
   G4double pz_GeV;      // momentum (GeV/c) of shower particle at point, z component
   G4double r_cm;        // global sqrt(x*x + y*y) of shower particle at point (cm)
   G4double z_cm;        // global z coordinate of shower particle at point (cm)
   G4double t_ns;        // time of shower particle crossing at point (ns)
   G4int track_;         // Geant4 track ID of particle making this shower
   G4int trackID_;       // GlueX-assigned track ID of particle making this shower

   G4int GetKey() const { return (track_ << 20) + int(t_ns * 100); }
};

typedef G4THitsMap<GlueXHitBCALpoint> GlueXHitsMapBCALpoint;

extern G4ThreadLocal G4Allocator<GlueXHitBCALpoint>* GlueXHitBCALpointAllocator;

inline void* GlueXHitBCALpoint::operator new(size_t)
{
   if (!GlueXHitBCALpointAllocator)
      GlueXHitBCALpointAllocator = new G4Allocator<GlueXHitBCALpoint>;
   return (void *) GlueXHitBCALpointAllocator->MallocSingle();
}

inline void GlueXHitBCALpoint::operator delete(void *aHit)
{
   GlueXHitBCALpointAllocator->FreeSingle((GlueXHitBCALpoint*) aHit);
}

#endif
