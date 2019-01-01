//
// GlueXHitSTCpoint - class header
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitSTCpoint_h
#define GlueXHitSTCpoint_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitSTCpoint : public G4VHit
{
 public:
   GlueXHitSTCpoint() {}
   GlueXHitSTCpoint(const GlueXHitSTCpoint &src);
   int operator==(const GlueXHitSTCpoint &right) const;
   GlueXHitSTCpoint &operator+=(const GlueXHitSTCpoint &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data
 
   G4double E_GeV;       // total energy (GeV) of this track at this point
   G4double dEdx_GeV_cm; // dE/dx (GeV/cm) of this track inside straw
   G4double phi_rad;     // phi angle (rad) of track at this point
   G4bool primary_;      // true if track belongs to from a primary particle
   G4int ptype_G3;       // G3 type of particle making this track
   G4double px_GeV;      // momentum (GeV/c) of track at point, x component
   G4double py_GeV;      // momentum (GeV/c) of track at point, y component
   G4double pz_GeV;      // momentum (GeV/c) of track at point, z component
   G4double r_cm;        // global sqrt(x*x + y*y) of track at point (cm)
   G4double z_cm;        // global z coordinate of track at point (cm)
   G4double t_ns;        // time of track crossing at point (ns)
   G4int sector_;        // azimuthal sector number of hit counter
   G4int track_;         // Geant4 track ID of particle making this track
   G4int trackID_;       // GlueX-assigned track ID of particle making this track

   G4int GetKey() const { return (track_ << 20) + int(t_ns * 100); }
};

typedef G4THitsMap<GlueXHitSTCpoint> GlueXHitsMapSTCpoint;

extern G4ThreadLocal G4Allocator<GlueXHitSTCpoint>* GlueXHitSTCpointAllocator;

inline void* GlueXHitSTCpoint::operator new(size_t)
{
   if (!GlueXHitSTCpointAllocator)
      GlueXHitSTCpointAllocator = new G4Allocator<GlueXHitSTCpoint>;
   return (void *) GlueXHitSTCpointAllocator->MallocSingle();
}

inline void GlueXHitSTCpoint::operator delete(void *aHit)
{
   GlueXHitSTCpointAllocator->FreeSingle((GlueXHitSTCpoint*) aHit);
}

#endif
