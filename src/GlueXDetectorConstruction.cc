//
// GlueXDetectorConstruction class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
 
#include "GlueXDetectorConstruction.hh"
#include "GlueXDetectorMessenger.hh"
#include "GlueXMagneticField.hh"

#include "G4Box.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

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

#include <string.h>
#include <libgen.h>
#include <errno.h>

#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

#define APP_NAME "HDGeant4"

GlueXDetectorConstruction::GlueXDetectorConstruction(G4String hddsFile)
: fMaxStep(0),
  fUniformField(0),
  fpMagneticField(0)
{
   fpDetectorMessenger = new GlueXDetectorMessenger(this);

   // Read the geometry description from a HDDS document
   // (Hall D Detector Specification) and build a Geant4
   // model of it in memory.

   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException& toCatch)
   {
      XString message(toCatch.getMessage());
      std::cerr
           << APP_NAME << " - error during xercesc initialization!"
           << std::endl << S(message) << std::endl;
      return;
   }

   XString xmlFile;
   if (hddsFile.size() > 0)
   {
      xmlFile = hddsFile.c_str();
   }
   else if (getenv("JANA_GEOMETRY_URL"))
   {
      hddsFile = getenv("JANA_GEOMETRY_URL");
      if (hddsFile.index("xmlfile://") == 0)
      {
         hddsFile.remove(0,10);
         xmlFile = hddsFile.c_str();
      }
      else
      {
         std::cerr
           << APP_NAME << " - unsupported protocol for JANA_GEOMETRY_URL"
           << std::endl << getenv("JANA_GEOMETRY_URL") << std::endl;
         return;
      }
   }
   else {
      std::cerr
           << APP_NAME << " - no hdds geometry file specified!"
           << std::endl;
      return;
   }
   
   int size=100;
   char *saved_cwd = new char[size];
   while (getcwd(saved_cwd, size) == 0)
   {
      delete [] saved_cwd;
      saved_cwd = new char[size *= 2];
   }
   char *dirpath = new char[hddsFile.size() + 2];
   chdir(dirname(strcpy(dirpath, hddsFile.c_str())));
   DOMDocument* document = buildDOMDocument(xmlFile,false);
   if (document == 0)
   {
      std::cerr
           << APP_NAME << " - error parsing HDDS document, "
           << "cannot continue" << std::endl;
      return;
   }
   chdir(saved_cwd);
   delete [] saved_cwd;
   delete [] dirpath;

   DOMNode* docEl;
   try {
      docEl = document->getDocumentElement();
   }
   catch (DOMException& e) {
      std::cerr << APP_NAME << " - Woops " << e.msg << std::endl;
      return;
   }

   DOMElement* rootEl = document->getElementById(X("everything"));
   if (rootEl == 0)
   {
      std::cerr
           << APP_NAME << " - error scanning HDDS document, " << std::endl
           << "  no element named \"everything\" found" << std::endl;
      return;
   }

   fHddsBuilder.translate(rootEl);

   XMLPlatformUtils::Terminate();
}

GlueXDetectorConstruction::
GlueXDetectorConstruction(const GlueXDetectorConstruction &src)
 : fHddsBuilder(src.fHddsBuilder)
{
     fMaxStep = src.fMaxStep;
     fUniformField = src.fUniformField;
     fpMagneticField = src.fpMagneticField; // shallow copy, sharing magfield
     fpDetectorMessenger = src.fpDetectorMessenger;  // shallow copy, sharing
}

GlueXDetectorConstruction::~GlueXDetectorConstruction()
{
  delete fpDetectorMessenger;             
}

void GlueXDetectorConstruction::SetUniformField(G4double field_T)
{
  // This method embeds the entire Hall D (spectrometer + beam line)
  // in an external uniform magnetic field. Used only for testing,
  // this is not how to simulate the GlueX spectrometer field!

  fUniformField = field_T;
  std::cout << "setting up a new field manager with a uniform field "
            << " of " << field_T << " Tesla."
            << std::endl;
  G4ThreeVector B(0,0,fUniformField);
  G4AffineTransform xform;
  fpMagneticField = new GlueXUniformMagField(B,tesla,xform);
}

void GlueXDetectorConstruction::SetMaxStep(G4double step_mm)
{
  fMaxStep = (step_mm > 0.)? step_mm*mm : 0;
}

G4VPhysicalVolume* GlueXDetectorConstruction::Construct()
{
   G4LogicalVolume *worldvol = fHddsBuilder.getWorldVolume();
   std::cout << "Root geometry volume " << worldvol->GetName();
   worldvol->SetName("World");
   std::cout << " configured as " << worldvol->GetName() << std::endl;
   worldvol->SetVisAttributes(new G4VisAttributes(false));
   return new G4PVPlacement(0, G4ThreeVector(), worldvol, "World", 0, 0, 0);
}

int GlueXDetectorConstruction::GetParallelWorldCount()
{
   for (int worlds=0; ; ++worlds) {
      if (fHddsBuilder.getWorldVolume(worlds) == 0) {
         return --worlds;
      }
   }
   return 0;
}

G4String GlueXDetectorConstruction::GetParallelWorldName(int paraIndex)
{
   std::stringstream worldStr;
   worldStr << "Parallel World " << paraIndex;
   return worldStr.str();
}

G4LogicalVolume*
 GlueXDetectorConstruction::GetParallelWorldVolume(int paraIndex)
{
   return fHddsBuilder.getWorldVolume(paraIndex);
}

GlueXParallelWorld::GlueXParallelWorld(const GlueXParallelWorld &src)
 : G4VUserParallelWorld(src.fWorldName), fTopVolume(src.fTopVolume)
{ }

GlueXParallelWorld::~GlueXParallelWorld() { }

void GlueXParallelWorld::Construct()
{
   G4VPhysicalVolume* ghostWorld = GetWorld();
   G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
   for (int child = fTopVolume->GetNoDaughters() - 1; child >= 0; --child) {
      G4VPhysicalVolume* pvol = fTopVolume->GetDaughter(child);
      fTopVolume->RemoveDaughter(pvol);
      worldLogical->AddDaughter(pvol);
   }
   if (worldLogical->GetNoDaughters()) {
      std::cout << "Additional geometry layer configured as "
                << ghostWorld->GetName() << std::endl;
   }
   else {
      std::cout << "Additional geometry layer configured as "
                << ghostWorld->GetName() << " is EMPTY!" << std::endl
                << "This is probably due to a geometry bug -- "
                << "PLEASE REPORT THIS TO THE AUTHORS" << std::endl;
   }
}
