//
// GlueXMagneticField class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state. The
// getter methods of these classes are thread-safe and provide
// concurrent interlock-free access to the magnetic fields from
// all worker threads.

#ifndef GlueXMagneticField_H
#define GlueXMagneticField_H

#include <vector>
#include <string>

#include <G4UniformMagField.hh>
#include <G4AffineTransform.hh>
#include <HDGEOMETRY/DMagneticFieldMap.h>
#include <HDGEOMETRY/DMagneticFieldMapPS.h>

class GlueXUniformMagField: public G4UniformMagField
{
 public:
   GlueXUniformMagField(G4ThreeVector B, G4double unit,
                        const G4AffineTransform &xform);
   GlueXUniformMagField(const GlueXUniformMagField &src);
   GlueXUniformMagField& operator=(const GlueXUniformMagField &src);
   ~GlueXUniformMagField();  

   virtual void SetMagFieldTransform(const G4AffineTransform &xform);
   virtual const G4AffineTransform& GetMagFieldTransform() const;

   virtual void SetMagField(G4ThreeVector B, G4double unit);
   virtual G4ThreeVector GetMagField(G4double unit) const;

 private:
   G4AffineTransform fXform;   // converts field map coordinates to region
};

class GlueXMappedMagField: public G4MagneticField
{
 public:
   GlueXMappedMagField(G4double Bmax, G4double unit,
                       const G4AffineTransform &xform);
   GlueXMappedMagField(const GlueXMappedMagField &src);
   GlueXMappedMagField& operator=(const GlueXMappedMagField &src);
   ~GlueXMappedMagField();  

   virtual void SetMagFieldTransform(const G4AffineTransform &xform);
   virtual const G4AffineTransform& GetMagFieldTransform() const;
 
   virtual int AddCartesianGrid(const int axsamples[4],
                                const int axorder[4],
                                const int axsense[4],
                                const G4double axlower[4],
                                const G4double axupper[4]);
   virtual int AddCylindricalGrid(const int axsamples[4],
                                  const int axorder[4],
                                  const int axsense[4],
                                  const G4double axlower[4],
                                  const G4double axupper[4]);
   virtual int ReadMapFile(const char *mapS);

   virtual G4ThreeVector GetMagField(const G4double point[4],
                                     G4double unit) const;
   virtual void  GetFieldValue(const G4double point[4],
                                     G4double *Bfield ) const;
   double GetBmax(G4double unit) const { return fBmax * fUnit/unit; }

 private:
   G4double fUnit;             // converts stored map into G4 field units
   G4double fBmax;             // max value of mapped field (may be useful)
   G4AffineTransform fXform;   // converts field map coordinates to region
   G4AffineTransform fXfinv;   // converts region coordinates to map coords
   int fGridtype;              // none (0), cartesian (1), or cylindrical (2)
   int fDim[3];                // sample grid dimensions (x,y,z) or (r,phi,z)
   int fOrder[3];              // specifies how the 3D grid is strung out into
                               // a 1D array, eg. cartesian (2,3,1) indicates
                               // index = fDim[0] * (fDim[2] * iy + iz) + ix

   struct field_map_grid_t {
      std::vector<int> sense;       // field map component sign reversal flags
      std::vector<G4double> lower;  // grid starting coordinate, by dimension
      std::vector<G4double> upper;  // grid ending coordinate, by dimension
                                    // in cm for length, radians for angle
   };
   std::vector<struct field_map_grid_t> fGrid;

   union field_map_entry_t {
      struct cartesian_map_entry_t {
         G4double x, y, z;
      } cart;
      struct cylindrical_map_entry_t {
         G4double r, phi, z;
      } cyl;
   };
   std::vector<union field_map_entry_t> fFieldMap;

   int lookup_field(int i1, int i2, int i3, G4double mapvalue[3]) const;
};

class GlueXComputedMagField: public G4MagneticField
{
 public:
   GlueXComputedMagField(G4double Bmax, G4double unit,
                         const G4AffineTransform &xform);
   GlueXComputedMagField(const GlueXComputedMagField &src);
   GlueXComputedMagField& operator=(const GlueXComputedMagField &src);
   ~GlueXComputedMagField();  

   void SetMagFieldTransform(const G4AffineTransform &xform);
   const G4AffineTransform& GetMagFieldTransform() const;
 
   void SetFunction(std::string function);
   std::string GetFunction() const;

   virtual G4ThreeVector GetMagField(const G4double point[4],
                                     G4double unit) const;
   virtual void  GetFieldValue(const G4double point[4],
                                     G4double *Bfield ) const;
   double GetBmax(G4double unit) const { return fBmax * fUnit/unit; }

 private:
   G4double fUnit;                 // converts stored map into G4 field units
   G4double fBmax;                 // max value of mapped field (may be useful)
   G4AffineTransform fXform;       // converts field map coordinates to region
   G4AffineTransform fXfinv;       // converts region coordinates to map coords
   std::string fFunction;          // keyword string identifying lookup field
   DMagneticFieldMap *fJanaFieldMap; // object from JANA framework that actually
                                   // contains the solenoid field map, normally
                                   // read from the ccdb database
   DMagneticFieldMapPS *fJanaFieldMapPS; // object from JANA framework that 
                                   // contains the pair spectrometer field map,
                                   // normally read from the ccdb database
};

#endif
