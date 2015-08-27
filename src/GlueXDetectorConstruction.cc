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
#include "G4MTRunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4HelixMixedStepper.hh"
#include "G4ClassicalRK4.hh"

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

pthread_mutex_t *GlueXDetectorConstruction::fMutex = 0;
std::list<GlueXDetectorConstruction*> GlueXDetectorConstruction::fInstance;

GlueXDetectorConstruction::GlueXDetectorConstruction(G4String hddsFile)
: fMaxStep(0),
  fUniformField(0),
  fpMagneticField(0)
{
   pthread_mutex_t *mutex = new pthread_mutex_t;
   if (fMutex == 0) {
      fMutex = mutex;
      pthread_mutex_init(mutex, 0);
   }
   if (fMutex != mutex) {
      pthread_mutex_destroy(mutex);
      delete mutex;
   }
   pthread_mutex_lock(fMutex);
   fInstance.push_back(this);
   pthread_mutex_unlock(fMutex);

   // Initialize the class that implements custom interactive
   // commands under the command prefix /hdgeant4/.
 
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
   pthread_mutex_lock(fMutex);
   fInstance.push_back(this);
   pthread_mutex_unlock(fMutex);
   fMaxStep = src.fMaxStep;
   fUniformField = src.fUniformField;
   fpMagneticField = src.fpMagneticField; // shallow copy, sharing magfield
   fpDetectorMessenger = src.fpDetectorMessenger;  // shallow copy, sharing
}

GlueXDetectorConstruction::~GlueXDetectorConstruction()
{
   pthread_mutex_lock(fMutex);
   fInstance.remove(this);
   int remaining = fInstance.size();
   pthread_mutex_unlock(fMutex);
   if (remaining == 0) {
      pthread_mutex_destroy(fMutex);
      delete fMutex;
      fMutex = 0;
   }
   delete fpDetectorMessenger;             
}

const GlueXDetectorConstruction *GlueXDetectorConstruction::GetInstance()
{
   // Generally one only needs a single instance of this object
   // per process, and this static method lets any component in the
   // application obtain the primary instance, if any. If none has
   // yet been constructed, it returns zero.

   if (fInstance.size() > 0)
      return *fInstance.begin();
   else
      return 0;
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

void GlueXDetectorConstruction::ConstructSDandField()
{
   G4RunManager::RMType rmtype = G4RunManager::GetRunManager()->
                                               GetRunManagerType();
   if (rmtype == G4RunManager::sequentialRM)
      return;
   else
      CloneF();
}

void GlueXDetectorConstruction::CloneF()
{
   typedef std::map<G4FieldManager*,G4FieldManager*> FMtoFMmap;
   typedef std::pair<G4FieldManager*,G4FieldManager*> FMpair;
   FMtoFMmap masterToWorker;
   G4LogicalVolumeStore* const logVolStore = G4LogicalVolumeStore::GetInstance();
   assert(logVolStore != NULL);
   for (G4LogicalVolumeStore::const_iterator iter = logVolStore->begin();
        iter != logVolStore->end();
        ++iter)
   {
      G4LogicalVolume *lvol = *iter;
      G4FieldManager* masterFM = lvol->GetFieldManager();
      G4FieldManager* clonedFM = 0;
      if (masterFM) {
         FMtoFMmap::iterator fmFound = masterToWorker.find(masterFM);
         if (fmFound == masterToWorker.end()) {

            // First time we see this FM, let's clone and remember...

            G4ChordFinder *cfinder = masterFM->GetChordFinder();
            G4MagInt_Driver *midriver = cfinder->GetIntegrationDriver();
            double stepMinimum = midriver->GetHmin();
            G4MagIntegratorStepper *stepper = midriver->GetStepper();
            G4EquationOfMotion *eqn = stepper->GetEquationOfMotion();
            const G4Field *field = masterFM->GetDetectorField();

            G4MagneticField *field_copy;
            if (dynamic_cast<const GlueXUniformMagField*>(field)) {
               GlueXUniformMagField &orig = *(GlueXUniformMagField*)field;
               field_copy = new GlueXUniformMagField(orig);
            }
            else if (dynamic_cast<const GlueXMappedMagField*>(field)) {
               GlueXMappedMagField &orig = *(GlueXMappedMagField*)field;
               field_copy = new GlueXMappedMagField(orig);
            }
            else if (dynamic_cast<const GlueXComputedMagField*>(field)) {
               GlueXComputedMagField &orig = *(GlueXComputedMagField*)field;
               field_copy = new GlueXComputedMagField(orig);
            }
            else if (dynamic_cast<const G4UniformMagField*>(field)) {
               G4UniformMagField &orig = *(G4UniformMagField*)field;
               field_copy = new G4UniformMagField(orig);
            }
            else {
               std::cerr << "GlueXDetectorConstruction::CloneF error - "
                         << "unknown G4MagneticField class found "
                         << "attached to detector volume " << lvol->GetName()
                         << ", cannot continue!" << std::endl;
               exit(1);
            }
            G4Mag_UsualEqRhs *eqn_copy = new G4Mag_UsualEqRhs(field_copy);
            G4MagIntegratorStepper *stepper_copy;
            if (dynamic_cast<const G4ExactHelixStepper*>(stepper)) {
               G4ExactHelixStepper &orig = *(G4ExactHelixStepper*)stepper;
               stepper_copy = new G4ExactHelixStepper(eqn_copy);
            }
            else if (dynamic_cast<const G4ClassicalRK4*>(stepper)) {
               G4ClassicalRK4 &orig = *(G4ClassicalRK4*)stepper;
               stepper_copy = new G4ClassicalRK4(eqn_copy);
            }
            else if (dynamic_cast<const G4HelixMixedStepper*>(stepper)) {
               G4HelixMixedStepper &orig = *(G4HelixMixedStepper*)stepper;
               stepper_copy = new G4HelixMixedStepper(eqn_copy);
            }
            else {
               std::cerr << "GlueXDetectorConstruction::CloneF error - "
                         << "unknown G4MagIntegratorStepper class found "
                         << "attached to detector volume " << lvol->GetName()
                         << ", cannot continue!" << std::endl;
               exit(1);
            }
            G4ChordFinder *cfinder_copy = new G4ChordFinder(field_copy,
                                                            stepMinimum,
                                                            stepper_copy);
            clonedFM = new G4FieldManager(field_copy, cfinder_copy);
            masterToWorker[masterFM] = clonedFM;
         }
         else {

            // We have already seen this FM attached to a different 
            // LogicalVolume, let's re-use the previous clone
 
            clonedFM = (*fmFound).second;
         }
         lvol->SetFieldManager(clonedFM, false);
      }
   }
}

int GlueXDetectorConstruction::GetParallelWorldCount() const
{
   for (int worlds=0; ; ++worlds) {
      if (fHddsBuilder.getWorldVolume(worlds) == 0) {
         return --worlds;
      }
   }
   return 0;
}

G4String GlueXDetectorConstruction::GetParallelWorldName(int paraIndex) const
{
   std::stringstream worldStr;
   worldStr << "Parallel World " << paraIndex;
   return worldStr.str();
}

G4LogicalVolume*
 GlueXDetectorConstruction::GetParallelWorldVolume(int paraIndex) const
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
