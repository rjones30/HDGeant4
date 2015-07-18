//
// GlueXPhysicsList class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#ifndef GlueXPhysicsList_h
#define GlueXPhysicsList_h 1

#include <GlueXDetectorConstruction.hh>

#include <G4VUserPhysicsList.hh>
#include <QGSP_FTFP_BERT.hh>
#include "globals.hh"

class GlueXPhysicsList: public QGSP_FTFP_BERT
{
 public:
   GlueXPhysicsList(GlueXDetectorConstruction *geometry);
};

#endif
