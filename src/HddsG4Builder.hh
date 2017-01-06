//
// HddsG4Builder - class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state.

#ifndef _HDDSG4BUILDER_
#define _HDDSG4BUILDER_

#include "XString.hpp"
#include "XParsers.hpp"
#include "hddsCommon.hpp"

#include <G4Element.hh>
#include <G4Material.hh>
#include <G4MagneticField.hh>
#include <G4FieldManager.hh>
#include <G4LogicalVolume.hh>
#include <G4RotationMatrix.hh>

#include <GlueXMagneticField.hh>

#include <map>

class HddsG4Builder : public CodeWriter
{
 public:
   HddsG4Builder();
   HddsG4Builder(const HddsG4Builder &src);
   HddsG4Builder &operator=(const HddsG4Builder &src);
   ~HddsG4Builder();

   int createMaterial(DOMElement* el);   // generate code for materials
   int createSolid(DOMElement* el,
                   Refsys& ref);    	 // generate code for solids
   int createRotation(Refsys& ref);      // generate code for rotations
   int createRegion(DOMElement* el,
                    Refsys& ref);        // generate code for regions
   int createVolume(DOMElement* el,
                    Refsys& ref);   	 // generate code for placement
   int createDivision(XString& divStr,
                      Refsys& ref); 	 // generate code for divisions
   void createSetFunctions(DOMElement* el,
                  const XString& ident); // generate code for properties
   void createGetFunctions(DOMElement* el,
                  const XString& ident); // generate code for identifiers
   void createMapFunctions(DOMElement* el,
                  const XString& ident); // generate code for field maps

   G4LogicalVolume* getWorldVolume(int parallel=0) const;
                                         // return ptr to world volume
   G4LogicalVolume* getVolume(const G4String volname) const;
                                         // look up logical volume by name
   int getVolumeId(G4LogicalVolume* vol) const;
                                         // reverse-find in fLogicalVolumes
   const std::map<int,G4LogicalVolume*> getSensitiveVolumes() const;
                                         // read-only access to fSensitiveVolumes

   void translate(DOMElement* topel);	 // invokes the main translator

 private:
   typedef std::pair<int,int> vpair_t;

   std::map<vpair_t,G4LogicalVolume*>::iterator
   addNewLayer(int volume_id, int layer); // generate code for geometry layers
   void addReflections(int volume_id);    // propagate volume to other layers

   int fWorldVolume;

   // many-to-one maps from volume id to attribute pointer
   std::map<int,G4Element*> fElements;
   std::map<int,G4Material*> fMaterials;
   std::map<int,G4MagneticField*> fMagneticRegions;
   std::map<int,G4FieldManager*> fFieldManagers;
   std::map<int,G4RotationMatrix*> fRotations;

   // one-to-one map from (volume id, layer index) to logical volume,
   // keeps track of each named volume as it is initially populated
   // (layer=0), and then any subsequent reflections, indexed by layer
   std::map<vpair_t,G4LogicalVolume*> fLogicalVolumes;

   // one-to-one map from (volume id, copy number) to physical volume,
   // keeps track of each named volume as it is initially placed,
   // whatever the layer; placement of reflections is not stored here
   std::map<vpair_t,G4VPhysicalVolume*> fPhysicalVolumes;

   // one-to-one map from volume id to logical volume for sensitive volumes
   std::map<int,G4LogicalVolume*> fSensitiveVolumes;

   // iterator arrays indexed by volume id into the above tables of
   // logical and physical volumes, about the current working state
   typedef std::map<int,std::map<vpair_t,G4LogicalVolume*>::iterator> 
           LogicalVolume_map_iter;
   typedef std::map<int,std::map<vpair_t,G4VPhysicalVolume*>::iterator> 
           PhysicalVolume_map_iter;
   LogicalVolume_map_iter fCurrentMother;
   PhysicalVolume_map_iter fCurrentPlacement;
   PhysicalVolume_map_iter fCurrentDivision;
   std::map<int,double> fCurrentPhiCenter;
};

#endif
