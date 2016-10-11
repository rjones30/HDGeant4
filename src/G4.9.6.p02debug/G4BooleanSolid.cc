//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// Implementation for the abstract base class for solids created by boolean 
// operations between other solids
//
// History:
//
// 10.09.98 V.Grichine, created
//
// --------------------------------------------------------------------

#include "G4BooleanSolid.hh"
#include "G4VSolid.hh"
#include "G4Polyhedron.hh"
#include "HepPolyhedronProcessor.h"
#include "Randomize.hh"

//////////////////////////////////////////////////////////////////
//
// Constructor

G4BooleanSolid::G4BooleanSolid( const G4String& pName,
                                G4VSolid* pSolidA ,
                                G4VSolid* pSolidB   ) :
  G4VSolid(pName), fAreaRatio(0.), fStatistics(1000000), fCubVolEpsilon(0.001),
  fAreaAccuracy(-1.), fCubicVolume(0.), fSurfaceArea(0.),
  fpPolyhedron(0), createdDisplacedSolid(false)
{
  fPtrSolidA = pSolidA ;
  fPtrSolidB = pSolidB ;
}

//////////////////////////////////////////////////////////////////
//
// Constructor

G4BooleanSolid::G4BooleanSolid( const G4String& pName,
                                      G4VSolid* pSolidA ,
                                      G4VSolid* pSolidB ,
                                      G4RotationMatrix* rotMatrix,
                                const G4ThreeVector& transVector    ) :
  G4VSolid(pName), fAreaRatio(0.), fStatistics(1000000), fCubVolEpsilon(0.001),
  fAreaAccuracy(-1.), fCubicVolume(0.), fSurfaceArea(0.),
  fpPolyhedron(0), createdDisplacedSolid(true)
{
  fPtrSolidA = pSolidA ;
  fPtrSolidB = new G4DisplacedSolid("placedB",pSolidB,rotMatrix,transVector) ;
}

//////////////////////////////////////////////////////////////////
//
// Constructor

G4BooleanSolid::G4BooleanSolid( const G4String& pName,
                                      G4VSolid* pSolidA ,
                                      G4VSolid* pSolidB ,
                                const G4Transform3D& transform    ) :
  G4VSolid(pName), fAreaRatio(0.), fStatistics(1000000), fCubVolEpsilon(0.001),
  fAreaAccuracy(-1.), fCubicVolume(0.), fSurfaceArea(0.),
  fpPolyhedron(0), createdDisplacedSolid(true)
{
  fPtrSolidA = pSolidA ;
  fPtrSolidB = new G4DisplacedSolid("placedB",pSolidB,transform) ;
}

///////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4BooleanSolid::G4BooleanSolid( __void__& a )
  : G4VSolid(a), fPtrSolidA(0), fPtrSolidB(0), fAreaRatio(0.),
    fStatistics(1000000), fCubVolEpsilon(0.001), 
    fAreaAccuracy(-1.), fCubicVolume(0.), fSurfaceArea(0.),
    fpPolyhedron(0), createdDisplacedSolid(false)
{
}

///////////////////////////////////////////////////////////////
//
// Destructor deletes transformation contents of the created displaced solid

G4BooleanSolid::~G4BooleanSolid() 
{
  if(createdDisplacedSolid)
  {
    ((G4DisplacedSolid*)fPtrSolidB)->CleanTransformations();
  }
  delete fpPolyhedron;
}

///////////////////////////////////////////////////////////////
//
// Copy constructor

G4BooleanSolid::G4BooleanSolid(const G4BooleanSolid& rhs)
  : G4VSolid (rhs), fPtrSolidA(rhs.fPtrSolidA), fPtrSolidB(rhs.fPtrSolidB),
    fAreaRatio(rhs.fAreaRatio),
    fStatistics(rhs.fStatistics), fCubVolEpsilon(rhs.fCubVolEpsilon),
    fAreaAccuracy(rhs.fAreaAccuracy), fCubicVolume(rhs.fCubicVolume),
    fSurfaceArea(rhs.fSurfaceArea), fpPolyhedron(0),
    createdDisplacedSolid(rhs.createdDisplacedSolid)
{
}

///////////////////////////////////////////////////////////////
//
// Assignment operator

G4BooleanSolid& G4BooleanSolid::operator = (const G4BooleanSolid& rhs) 
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4VSolid::operator=(rhs);

  // Copy data
  //
  fPtrSolidA= rhs.fPtrSolidA; fPtrSolidB= rhs.fPtrSolidB;
  fAreaRatio= rhs.fAreaRatio;
  fStatistics= rhs.fStatistics; fCubVolEpsilon= rhs.fCubVolEpsilon;
  fAreaAccuracy= rhs.fAreaAccuracy; fCubicVolume= rhs.fCubicVolume;
  fSurfaceArea= rhs.fSurfaceArea; fpPolyhedron= 0;
  createdDisplacedSolid= rhs.createdDisplacedSolid;

  return *this;
}  

///////////////////////////////////////////////////////////////
//
// If Solid is made up from a Boolean operation of two solids,
//   return the corresponding solid (for no=0 and 1)
// If the solid is not a "Boolean", return 0

const G4VSolid* G4BooleanSolid::GetConstituentSolid(G4int no) const
{
  const G4VSolid*  subSolid=0;
  if( no == 0 )  
    subSolid = fPtrSolidA;
  else if( no == 1 ) 
    subSolid = fPtrSolidB;
  else
  {
    DumpInfo();
    G4Exception("G4BooleanSolid::GetConstituentSolid()",
                "GeomSolids0002", FatalException, "Invalid solid index.");
  }

  return subSolid;
}

///////////////////////////////////////////////////////////////
//
// If Solid is made up from a Boolean operation of two solids,
//   return the corresponding solid (for no=0 and 1)
// If the solid is not a "Boolean", return 0

G4VSolid* G4BooleanSolid::GetConstituentSolid(G4int no)
{
  G4VSolid*  subSolid=0;
  if( no == 0 )  
    subSolid = fPtrSolidA;
  else if( no == 1 ) 
    subSolid = fPtrSolidB;
  else
  {
    DumpInfo();
    G4Exception("G4BooleanSolid::GetConstituentSolid()",
                "GeomSolids0002", FatalException, "Invalid solid index.");
  }

  return subSolid;
}

//////////////////////////////////////////////////////////////////////////
//
// Returns entity type

G4GeometryType G4BooleanSolid::GetEntityType() const 
{
  return G4String("G4BooleanSolid");
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4BooleanSolid::StreamInfo(std::ostream& os) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for Boolean solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: " << GetEntityType() << "\n"
     << " Parameters of constituent solids: \n"
     << "===========================================================\n";
  fPtrSolidA->StreamInfo(os);
  fPtrSolidB->StreamInfo(os);
  os << "===========================================================\n";

  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// Returns a point (G4ThreeVector) randomly and uniformly selected
// on the solid surface
//

G4ThreeVector G4BooleanSolid::GetPointOnSurface() const
{
  G4double rand;
  G4ThreeVector p;

  do
  {
    rand = G4UniformRand();

    if (rand < GetAreaRatio()) { p = fPtrSolidA->GetPointOnSurface(); }
    else                       { p = fPtrSolidB->GetPointOnSurface(); }
  } while (Inside(p) != kSurface);

  return p;
}

//////////////////////////////////////////////////////////////////////////
//
// Returns polyhedron for visualization

G4Polyhedron* G4BooleanSolid::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
    }
  return fpPolyhedron;
}

//////////////////////////////////////////////////////////////////////////
//
// Stacks polyhedra for processing. Returns top polyhedron.

G4Polyhedron*
G4BooleanSolid::StackPolyhedron(HepPolyhedronProcessor& processor,
                                const G4VSolid* solid) const
{
  HepPolyhedronProcessor::Operation operation;
  const G4String& type = solid->GetEntityType();
  if (type == "G4UnionSolid")
    { operation = HepPolyhedronProcessor::UNION; }
  else if (type == "G4IntersectionSolid")
    { operation = HepPolyhedronProcessor::INTERSECTION; }
  else if (type == "G4SubtractionSolid")
    { operation = HepPolyhedronProcessor::SUBTRACTION; }
  else
  {
    std::ostringstream message;
    message << "Solid - " << solid->GetName()
            << " - Unrecognised composite solid" << G4endl
            << " Returning NULL !";
    G4Exception("StackPolyhedron()", "GeomSolids1001", JustWarning, message);
    return 0;
  }

  G4Polyhedron* top = 0;
  const G4VSolid* solidA = solid->GetConstituentSolid(0);
  const G4VSolid* solidB = solid->GetConstituentSolid(1);

  if (solidA->GetConstituentSolid(0))
  {
    top = StackPolyhedron(processor, solidA);
  }
  else
  {
    top = solidA->GetPolyhedron();
  }

  G4Polyhedron* operand = solidB->GetPolyhedron();

#ifdef DEBUG_SECTIONPLANE
G4cout << "Here we are analyzing a " << type << std::endl
       << "  between volume " << solidA->GetName()
       << "  with " << top->GetNoFacets() << " facets:" << std::endl;
#ifdef DEBUG_VERTEX_LISTING
for (int v=1; v <= top->GetNoVertices(); ++v) {
   G4Point3D vtx = top->GetVertex(v);
   G4cout << std::setw(12) << vtx[0] << std::setw(12) << vtx[1]
          << std::setw(12) << vtx[2] << std::endl;
}
#endif
G4cout << "  and volume " << solidB->GetName()
       << "  with " << operand->GetNoFacets() << " facets:" << std::endl;
#ifdef DEBUG_VERTEX_LISTING
for (int v=1; v <= operand->GetNoVertices(); ++v) {
   G4Point3D vtx = operand->GetVertex(v);
   G4cout << std::setw(12) << vtx[0] << std::setw(12) << vtx[1]
          << std::setw(12) << vtx[2] << std::endl;
}
#endif
#ifdef DEBUG_FACETS_DETAIL
const G4VSolid *sold = solidB;
const G4Polyhedron *solp = operand;
G4cout << "Facets detail for volume " << sold->GetName() << " :" << std::endl;
for (int f=1; f <= solp->GetNoFacets(); ++f) {
   int Nedge;
   int ivertex[4];
   int edgeFlag[4];
   int opface[4];
   G4Point3D vertex[4];
   G4Normal3D normal[4];
   solp->GetFacet(f,Nedge,ivertex,edgeFlag,opface);
   solp->GetFacet(f,Nedge,vertex,edgeFlag,normal);
   G4cout << "  facet " << f << " has " << Nedge << " edges:" << std::endl;
   for (int e=1; e <= Nedge; ++e) {
      G4cout << "    " << e << ": vertex " << vertex[e-1]
             << ", adjoins facet " << opface[e-1]
             << ", edge flag " << edgeFlag[e-1]
             << ", has normal:" << std::endl
             << std::setw(12) << normal[e-1][0]
             << std::setw(12) << normal[e-1][1]
             << std::setw(12) << normal[e-1][2] << std::endl;
   }
}
#endif
#endif

  processor.push_back (operation, *operand);

  return top;
}
