//
// GlueXHitGCALpoint - class header
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitGCALpoint_h
#define GlueXHitGCALpoint_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitGCALpoint : public G4VHit
{
 public:
   GlueXHitGCALpoint() {}
   GlueXHitGCALpoint(const GlueXHitGCALpoint &src);
   int operator==(const GlueXHitGCALpoint &right) const;
   GlueXHitGCALpoint &operator+=(const GlueXHitGCALpoint &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data
 
   G4double E_GeV;       // total energy (GeV) of this shower particle at this point
   G4bool primary_;      // true if track belongs to from a primary particle
   G4int ptype_G3;       // G3 type of particle making this shower
   G4double px_GeV;      // momentum (GeV/c) of shower particle at point, x component
   G4double py_GeV;      // momentum (GeV/c) of shower particle at point, y component
   G4double pz_GeV;      // momentum (GeV/c) of shower particle at point, z component
   G4double r_cm;        // global sqrt(x*x + y*y) of shower particle at point (cm)
   G4double phi_rad;     // global y coordinate of shower particle at point (cm)
   G4double z_cm;        // global z coordinate of shower particle at point (cm)
   G4double t_ns;        // time of shower particle crossing at point (ns)
   G4int track_;         // Geant4 track ID of particle making this shower
   G4int trackID_;       // GlueX-assigned track ID of particle making this shower

   G4int GetKey() const { return (track_ << 20) + int(t_ns * 100); }
};

typedef G4THitsMap<GlueXHitGCALpoint> GlueXHitsMapGCALpoint;

extern G4ThreadLocal G4Allocator<GlueXHitGCALpoint>* GlueXHitGCALpointAllocator;

inline void* GlueXHitGCALpoint::operator new(size_t)
{
   if (!GlueXHitGCALpointAllocator)
      GlueXHitGCALpointAllocator = new G4Allocator<GlueXHitGCALpoint>;
   return (void *) GlueXHitGCALpointAllocator->MallocSingle();
}

inline void GlueXHitGCALpoint::operator delete(void *aHit)
{
   GlueXHitGCALpointAllocator->FreeSingle((GlueXHitGCALpoint*) aHit);
}

#endif
