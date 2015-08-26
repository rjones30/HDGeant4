//
// GlueXDetectorConstruction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#ifndef GlueXDetectorConstruction_h
#define GlueXDetectorConstruction_h 1

#include <list>
#include <pthread.h>

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include <G4VUserParallelWorld.hh>
#include <GlueXMagneticField.hh>
#include <HddsG4Builder.hh>

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class GlueXDetectorMessenger;

class GlueXDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
     GlueXDetectorConstruction(G4String hddsFile="");
     GlueXDetectorConstruction(const GlueXDetectorConstruction &src);
     GlueXDetectorConstruction &operator=(const GlueXDetectorConstruction &src);
    ~GlueXDetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
     void SetUniformField(G4double field_T);
     void SetMaxStep (G4double step_mm);     

     G4double GetUniformField(G4double unit) const {
        return fUniformField * tesla / unit;
     }
     G4double GetMaxStep(G4double unit) const {
        return fMaxStep * mm / unit;
     }

     int GetParallelWorldCount() const;
     G4String GetParallelWorldName(int paraIndex) const;
     G4LogicalVolume* GetParallelWorldVolume(int paraIndex) const;
     
     static const GlueXDetectorConstruction* GetInstance();

  private:
     G4double fMaxStep;		// maximum step size for tracking
     G4double fUniformField;  	// optional uniform field (for testing)
     G4MagneticField* fpMagneticField; // pointer to the field manager
     GlueXDetectorMessenger* fpDetectorMessenger;  // pointer to the Messenger
     HddsG4Builder fHddsBuilder; // hdds translator object instance

     static pthread_mutex_t *fMutex;
     static std::list<GlueXDetectorConstruction*> fInstance;
};

class GlueXParallelWorld : public G4VUserParallelWorld
{
 public:
   GlueXParallelWorld(G4String worldName, G4LogicalVolume *topVolume)
    : G4VUserParallelWorld(worldName), fTopVolume(topVolume)
   {}
   GlueXParallelWorld(const GlueXParallelWorld &src);
   GlueXParallelWorld &operator=(const GlueXParallelWorld &src);
   ~GlueXParallelWorld();

   virtual void Construct(); 

 private:
   G4LogicalVolume *fTopVolume;
};

#endif
