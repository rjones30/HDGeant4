//
// HddsG4Builder - class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//

#ifndef _HDDSG4BUILDER_
#define _HDDSG4BUILDER_

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLStringTokenizer.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>

using namespace xercesc;

#include "XString.hpp"
#include "XParsers.hpp"
#include "hddsCommon.hpp"

#include <G4Element.hh>
#include <G4Material.hh>
#include <G4MagneticField.hh>
#include <G4LogicalVolume.hh>
#include <G4RotationMatrix.hh>
#include <G4PVPlacement.hh>
#include <G4PVDivision.hh>
#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4Tubs.hh>
#include <G4Trap.hh>
#include <G4Sphere.hh>
#include <G4EllipticalTube.hh>
#include <G4Polycone.hh>
#include <G4Polyhedra.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4VisAttributes.hh>

#include <GlueXMagneticField.hh>

#include <map>
#include <vector>

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
                    Refsys& ref);  	 // generate code for placement
   int createDivision(XString& divStr,
                      Refsys& ref);	 // generate code for divisions
   void createSetFunctions(DOMElement* el,
                  const XString& ident); // generate code for properties
   void createGetFunctions(DOMElement* el,
                  const XString& ident); // generate code for identifiers
   void createMapFunctions(DOMElement* el,
                  const XString& ident); // generate code for field maps

   G4LogicalVolume* getWorldVolume(int parallel=0);
                                         // return ptr to world volume
   void translate(DOMElement* topel);	 // invokes the translator

 private:
   typedef std::pair<int,int> vpair_t;

   std::map<vpair_t,G4LogicalVolume*>::iterator
   addNewLayer(int volume_id, int layer); // generate code for geometry layers

   int fWorldVolume;
   std::map<int,G4Element*> fElements;
   std::map<int,G4Material*> fMaterials;
   std::map<int,G4MagneticField*> fMagneticRegions;
   std::map<vpair_t,G4LogicalVolume*> fLogicalVolumes;
   std::map<vpair_t,G4VPhysicalVolume*> fPhysicalVolumes;
   std::map<int,G4RotationMatrix*> fRotations;
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
