//
// GlueXHitDIRCpoint - class header
//
// author: richard.t.jones at uconn.edu
// version: november 26, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitDIRCpoint_h
#define GlueXHitDIRCpoint_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitDIRCpoint : public G4VHit
{
 public:
   GlueXHitDIRCpoint() {}
   GlueXHitDIRCpoint(const GlueXHitDIRCpoint &src);
   int operator==(const GlueXHitDIRCpoint &right) const;
   GlueXHitDIRCpoint &operator+=(const GlueXHitDIRCpoint &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data
 
   G4double E_GeV;       // total energy (GeV) of this track at this point
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

typedef G4THitsMap<GlueXHitDIRCpoint> GlueXHitsMapDIRCpoint;

extern G4ThreadLocal G4Allocator<GlueXHitDIRCpoint>* GlueXHitDIRCpointAllocator;

inline void* GlueXHitDIRCpoint::operator new(size_t)
{
   if (!GlueXHitDIRCpointAllocator)
      GlueXHitDIRCpointAllocator = new G4Allocator<GlueXHitDIRCpoint>;
   return (void *) GlueXHitDIRCpointAllocator->MallocSingle();
}

inline void GlueXHitDIRCpoint::operator delete(void *aHit)
{
   GlueXHitDIRCpointAllocator->FreeSingle((GlueXHitDIRCpoint*) aHit);
}

#endif
