//
// GlueXHitCDCpoint - class header
//
// author: richard.t.jones at uconn.edu
// version: august 28, 2015
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitCDCpoint_h
#define GlueXHitCDCpoint_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitCDCpoint : public G4VHit
{
 public:
   GlueXHitCDCpoint() {}
   GlueXHitCDCpoint(const GlueXHitCDCpoint &src);
   int operator==(const GlueXHitCDCpoint &right) const;
   GlueXHitCDCpoint &operator+=(const GlueXHitCDCpoint &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data
 
   G4double dEdx_GeV_cm; // dE/dx (GeV/cm) of this track inside straw
   G4double dradius_cm;  // track distance (cm) of closest approach to wire
   G4double phi_rad;     // phi angle (rad) of track closest approach to wire
   G4bool primary_;      // true if track belongs to from a primary particle
   G4int ptype_G3;       // G3 type of particle making this track
   G4double px_GeV;      // momentum (GeV/c) of track at point, x component
   G4double py_GeV;      // momentum (GeV/c) of track at point, y component
   G4double pz_GeV;      // momentum (GeV/c) of track at point, z component
   G4double r_cm;        // global sqrt(x*x + y*y) of track at point (cm)
   G4double z_cm;        // global z coordinate of track at point (cm)
   G4double t_ns;        // time of track crossing at point (ns)
   G4int track_;         // Geant4 track ID of particle making this track
   G4int trackID_;       // GlueX-assigned track ID of particle making this track
   G4int sector_;        // CDC straw number within ring where point lies
   G4int ring_;          // CDC ring number where point lies

   G4int GetKey() const { return (track_ << 20) + int(t_ns * 100); }
};

typedef G4THitsMap<GlueXHitCDCpoint> GlueXHitsMapCDCpoint;

extern G4ThreadLocal G4Allocator<GlueXHitCDCpoint>* GlueXHitCDCpointAllocator;

inline void* GlueXHitCDCpoint::operator new(size_t)
{
   if (!GlueXHitCDCpointAllocator)
      GlueXHitCDCpointAllocator = new G4Allocator<GlueXHitCDCpoint>;
   return (void *) GlueXHitCDCpointAllocator->MallocSingle();
}

inline void GlueXHitCDCpoint::operator delete(void *aHit)
{
   GlueXHitCDCpointAllocator->FreeSingle((GlueXHitCDCpoint*) aHit);
}

#endif
