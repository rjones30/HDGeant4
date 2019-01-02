//
// GlueXHitFDCpoint - class header
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitFDCpoint_h
#define GlueXHitFDCpoint_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitFDCpoint : public G4VHit
{
 public:
   GlueXHitFDCpoint() {}
   GlueXHitFDCpoint(int chamber);
   GlueXHitFDCpoint(const GlueXHitFDCpoint &src);
   int operator==(const GlueXHitFDCpoint &right) const;
   GlueXHitFDCpoint &operator+=(const GlueXHitFDCpoint &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   G4int chamber_;        // FDC chamber, count from upstream starting at 1

   G4double E_GeV;       // energy of primary track at point
   G4double dEdx_GeV_cm; // dE/dx (GeV/cm) of this track inside straw
   G4double dradius_cm;  // track distance (cm) of closest approach to wire
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

typedef G4THitsMap<GlueXHitFDCpoint> GlueXHitsMapFDCpoint;

extern G4ThreadLocal G4Allocator<GlueXHitFDCpoint>* GlueXHitFDCpointAllocator;

inline void* GlueXHitFDCpoint::operator new(size_t)
{
   if (!GlueXHitFDCpointAllocator)
      GlueXHitFDCpointAllocator = new G4Allocator<GlueXHitFDCpoint>;
   return (void *) GlueXHitFDCpointAllocator->MallocSingle();
}

inline void GlueXHitFDCpoint::operator delete(void *aHit)
{
   GlueXHitFDCpointAllocator->FreeSingle((GlueXHitFDCpoint*) aHit);
}

#endif
