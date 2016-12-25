//
// PairConversionGenerator class header
//
// author: richard.t.jones at uconn.edu
// version: december 17, 2016
//
// notes:
//
// This class computes differential pair production rates with full
// polarization degrees of freedom taken into account for both the
// photon and the leptons. Both incoherent pair production with a
// recoil electron and nuclear/coherent atomic scattering without a
// free recoil electron are included.
//
// dependencies:
//
// 1. This class depends on the Dirac++ package by the same author.
// 2. The Dirac++ package uses extended precision floating point arithmetic
//    called "long double" in c++. On Intel cpu architectures, this is 
//    implemented using 80-bit floating-point registers known as "x87".
//    Extended precision is needed to keep rounding errors in the sum over
//    cancelling Feynman amplitudes from producing excessive rounding
//    errors in the computation of the differential cross section.
//
// units:
// Any length is in m; energy,momentum,mass in GeV (c=1); angles in
// radians; time in seconds; cross section in barns.

#ifdef USING_DIRACXX
#ifndef PairConversionGenerator_h
#define PairConversionGenerator_h 1

#include <string>
#include <vector>

#include "Complex.h"
#include "TPhoton.h"
#include "TLepton.h"
#include "TCrossSection.h"
#include "TLorentzBoost.h"

#if BOOST_PYTHON_WRAPPING
#include <boost/python.hpp>
#endif

class PairConversionGenerator {
 public:
   PairConversionGenerator();
   ~PairConversionGenerator();

   double FFatomic(double qRecoil);
   double DiffXS_pair(const TPhoton &gIn, 
                      const TLepton &pOut, const TLepton &eOut);
   double DiffXS_triplet(const TPhoton &gIn, const TLepton &pOut,
                         const TLepton &eOut2, const TLepton &eOut3);

   const TThreeVectorReal &GetPolarization();
   void SetPolarization(double polx, double poly, double polz);
   void SetPolarization(const TThreeVectorReal &pol);
   void SetPolarization(double pol[3]);
   unsigned int GetConverterZ();
   void SetConverterZ(unsigned int Z);

 protected:
   TThreeVectorReal fPolar;   // incident photon polarization (Stokes parameterization)
   unsigned int fConverterZ;  // atomic number of converter material

 private:
   PairConversionGenerator(const PairConversionGenerator &src);
   PairConversionGenerator &operator=(const PairConversionGenerator &src);
};

inline const TThreeVectorReal &PairConversionGenerator::GetPolarization() {
   return fPolar;
}

inline void PairConversionGenerator::SetPolarization(double polx,
                                                     double poly,
                                                     double polz)
{
   fPolar[0] = polx;
   fPolar[1] = poly;
   fPolar[2] = polz;
}

inline void PairConversionGenerator::SetPolarization(const TThreeVectorReal &pol)
{
   fPolar = pol;
}

inline void PairConversionGenerator::SetPolarization(double pol[3])
{
   fPolar[0] = pol[0];
   fPolar[1] = pol[1];
   fPolar[2] = pol[2];
}

inline unsigned int PairConversionGenerator::GetConverterZ()
{
   return fConverterZ;
}

inline void PairConversionGenerator::SetConverterZ(unsigned int Z)
{
   fConverterZ = Z;
}

#endif
#endif
