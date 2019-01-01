//
// GlueXHitPSpoint - class header
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitPSpoint_h
#define GlueXHitPSpoint_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitPSpoint : public G4VHit
{
 public:
   GlueXHitPSpoint() {}
   GlueXHitPSpoint(const GlueXHitPSpoint &src);
   int operator==(const GlueXHitPSpoint &right) const;
   GlueXHitPSpoint &operator+=(const GlueXHitPSpoint &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data
 
   G4int arm_;           // arm=0 for right/south, arm=1 for left/north
   G4int column_;        // tile number, starts at 1, low x to high
   G4double E_GeV;       // total energy (GeV) of this track at this point
   G4double dEdx_GeV_cm; // dE/dx (GeV/cm) of this track inside straw
   G4bool primary_;      // true if track belongs to from a primary particle
   G4int ptype_G3;       // G3 type of particle making this track
   G4double px_GeV;      // momentum (GeV/c) of track at point, x component
   G4double py_GeV;      // momentum (GeV/c) of track at point, y component
   G4double pz_GeV;      // momentum (GeV/c) of track at point, z component
   G4double x_cm;        // global x coordinate of track at point (cm)
   G4double y_cm;        // global y coordinate of track at point (cm)
   G4double z_cm;        // global z coordinate of track at point (cm)
   G4double t_ns;        // time of track crossing at point (ns)
   G4int track_;         // Geant4 track ID of particle making this track
   G4int trackID_;       // GlueX-assigned track ID of particle making this track

   G4int GetKey() const { return (track_ << 20) + int(t_ns * 100); }
};

typedef G4THitsMap<GlueXHitPSpoint> GlueXHitsMapPSpoint;

extern G4ThreadLocal G4Allocator<GlueXHitPSpoint>* GlueXHitPSpointAllocator;

inline void* GlueXHitPSpoint::operator new(size_t)
{
   if (!GlueXHitPSpointAllocator)
      GlueXHitPSpointAllocator = new G4Allocator<GlueXHitPSpoint>;
   return (void *) GlueXHitPSpointAllocator->MallocSingle();
}

inline void GlueXHitPSpoint::operator delete(void *aHit)
{
   GlueXHitPSpointAllocator->FreeSingle((GlueXHitPSpoint*) aHit);
}

#endif
