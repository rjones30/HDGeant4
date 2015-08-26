//
// GlueXParticleGun class header
//
// author: richard.t.jones at uconn.edu
// version: july 28, 2015
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state. It is
// invoked from the GlueXPrimaryGeneratorAction class, which is
// responsible for managing the interlocks to ensure thread safety.

#ifndef GlueXParticleGun_H
#define GlueXParticleGun_H

#include <G4ParticleGun.hh>
#include <G4ParticleDefinition.hh>

class GlueXParticleGun: public G4ParticleGun
{
 public:
   GlueXParticleGun() : G4ParticleGun() {}
   GlueXParticleGun(int numberofparticles)
    : G4ParticleGun(numberofparticles) {}
   GlueXParticleGun(G4ParticleDefinition *particleDef, int numberofparticles=1)
    : G4ParticleGun(particleDef, numberofparticles) {}

   void Reset() {
      particle_energy = 0.0;
   }
};

#endif
