//
// GlueXUserEventInformation class header
//
// author: richard.t.jones at uconn.edu
// version: aug 1, 2015
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.

#ifndef _GLUEXUSEREVENTINFORMATION_
#define _GLUEXUSEREVENTINFORMATION_

#include <iostream>
#include <vector>

#include "G4UImanager.hh"
#include "G4VUserEventInformation.hh"

#include <HDDM/hddm_s.hpp>

class GlueXUserEventInformation: public G4VUserEventInformation
{
 public:
   GlueXUserEventInformation(hddm_s::HDDM *hddmevent=NULL)
    : hddmevent(hddmevent),
      keep_event(true),
      Nprimaries(0)
   {
      if (hddmevent) {
         // Count number of primaries already in event
         const hddm_s::ProductList &products = hddmevent->getProducts();
         Nprimaries = products.size();
      }
   }

   GlueXUserEventInformation(int geanttype, G4ThreeVector &pos, 
                                            G4ThreeVector &mom)
    : hddmevent(NULL),
      Nprimaries(1)
   {
      vertex_t v(geanttype, pos, mom);
      vertices.push_back(v);
   }

   ~GlueXUserEventInformation() {
      if (hddmevent != NULL)
         delete hddmevent;
   }

   void Print() const {
      G4cout << "GlueXUserEventInformation: hddm_s=" << hddmevent 
             << G4endl;
   }

   hddm_s::HDDM *hddmevent;
   bool keep_event;

   // These hold info. on particles whose parameters are generated within GlueXsim
   // This complicated structure matches what is already in HDDM and will allow
   // more complicated things like pi0 decay by GEANT4 or Lambda decay to be
   // captured here (eventually).

   class partdef_t{
    public:
      partdef_t(int geanttype, G4ThreeVector &mom,
                int id=1, int parentid=0, int mech=0)
       : geanttype(geanttype),
         mom(mom),
         id(id),
         parentid(parentid),
         mech(mech)
      {}

      int geanttype;
      G4ThreeVector mom;
      int id;
      int parentid;
      int mech;
   };

   class vertex_t{
    public:
       vertex_t(int geanttype, G4ThreeVector &pos, G4ThreeVector &mom)
        : pos(pos)
       {
          partdefs.push_back(partdef_t(geanttype, mom));
       }

       G4ThreeVector pos;
       std::vector<partdef_t> partdefs;
   };

   std::vector<vertex_t> vertices;
   int Nprimaries;
};

#endif // _GLUEXUSEREVENTINFORMATION_
