
// 4/9/2017  D. Lawrence
// When trying to run on OS X 10.11 I got a run time link error
// for G4MTHepRandom::getTheTableSeeds(long*, int). The routine
// CLHEP::HepRandom::getTheTableSeeds(long*, int) was defined in
// the G4 libraries. Unfortunately, it seems the G4MTHepRandom
// only has a pointer to a HepRandomEngine which does not have
// this method.
//
// Looking at G4MTHepRandom in the G4 source, it oddly declares
// this in the header file, but does not define it in the .cc
// file along with all of the other static methods (???).
//
// This just defines the static method modeled after the other
// static methods in G4MTHepRandom

#include <stdlib.h>
#include <iostream>
using namespace std;

#include <G4MTHepRandom.hh>
#include <CLHEP/Random/Random.h>

#ifdef __APPLE__

void G4MTHepRandom::getTheTableSeeds (G4long* seeds, G4int index)
{
	static CLHEP::HepRandom *hr = NULL;
	if(!hr) hr = new CLHEP::HepRandom(getTheEngine());
	hr->getTheTableSeeds(seeds, index);
}

#endif // __APPLE__
