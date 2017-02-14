//
// GlueXBremsstrahlungGenerator class implementation
//
// author: richard.t.jones at uconn.edu
// version: february 4, 2017
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state. It is
// invoked from the GlueXGeneratorAction class, which is
// responsible for managing the interlocks to ensure thread safety.
//
// Notes:
// This is an auxilliary class that is not used in the main production
// version of HDGeant4. It was written to provide an alternate means to
// calculate statistical properties of the bremsstrahlung beam from the
// underlying QED processes. The primary photon beam generator class is
// GlueXPhotonBeamGenerator, which uses analytic distribution formulas
// to generate distributions in photon energy, angle, and polarization.
// By contrast, this GlueXBremsstrahlungGenerator directly computes the
// QED amplitude for a specific choise for the event kinematics, and
// then employs importance sampling to achieve a true reproduction of
// the parent distribution without needing to know any formulas that
// describe them. This latter method has the advantage that one can
// use it to sample any property of the reaction (eg. the polarization
// of the recoil electron, correlation of the photon polarization with
// emission angle, etc.) simply by histogramming it over a sample.

#include <GlueXBremsstrahlungGenerator.hh>

#include <iostream>
#include <math.h>

#define sqr(x) ((x)*(x))

#if USING_DIRACXX

#include <TPhoton.h>
#include <TLepton.h>
#include <TLorentzBoost.h>
#include <TCrossSection.h>

#include <TFile.h>
#include <TTree.h>

ImportanceSampler GlueXBremsstrahlungGenerator::fDummyPDFx;

TTree *fTree = 0;
TFile *treefile = 0;

GlueXBremsstrahlungGenerator::GlueXBremsstrahlungGenerator()
 : fBeamEnergy(12.), fMinEnergy(3.)
{
   if (treefile == 0) {
      treefile = new TFile("beamtree.root", "update");
      if (fTree) {
         delete fTree;
         fTree = 0;
      }
   }
   if (fTree == 0) {
      fTree = new TTree("beam", "polarized bremstrahlung generator");
      fTree->Branch("Ebeam", &Ebeam, "Ebeam/D");
      fTree->Branch("qT", &qT, "qT/D");
      fTree->Branch("qTphi", &qTphi, "qTphi/D");
      fTree->Branch("qL", &qL, "qL/D");
      fTree->Branch("kstar", &kstar, "kstar/D");
      fTree->Branch("M2", &M2, "M2/D");
      fTree->Branch("phi", &phi, "phi/D");
      fTree->Branch("u", &u[0], "u[5]/D");
      fTree->Branch("k", &k[0], "k[4]/D");
      fTree->Branch("pR", &pR[0], "pR[4]/D");
      fTree->Branch("polar_0_90", &polar_0_90, "polar_0_90/D");
      fTree->Branch("polar_45_135", &polar_45_135, "polar_45_135/D");
      fTree->Branch("diffXS", &diffXS, "diffXS/D");
      fTree->Branch("weight", &weight, "weight/D");
   }

#if defined DO_IMPORTANCE_SAMPLE
   if (fTDummyPDFx.density.size() == 0) {
      std::cout << "GlueXBremsstrahlungGenerator constructor: "
                << "Setting up importance sampling tables, please wait... "
                << std::flush;
      prepareImportanceSamplingPDFs();
      std::cout << "finished." << std::endl;
   }
#endif
}

GlueXBremsstrahlungGenerator::~GlueXBremsstrahlungGenerator()
{
   if (fTree && treefile) {
      fTree->Write();
      delete fTree;
      fTree = 0;
   }
   if (treefile) {
      delete treefile;
      treefile = 0;
   }
}

void GlueXBremsstrahlungGenerator::GenerateBeamPhotons(int nevents)
{
   // This method implements the generator as a Monte Carlo loop
   // that runs until nevents bremsstrahlung events have been
   // generated between photon energy fMinEnergy and fBeamEnergy.
   // The results are stored in a ROOT tree called "beam".
   //
   // The key input for this generator is the differential cross section
   // for atomic bremsstrahlung that is computed using the TCrossSection
   // class provided by the Dirac++ package. The Bremsstrahlung method
   // of this class takes as its arguments Dirac++ objects describing
   // the incident electron, the recoil electron and the emitted photon,
   // and returns the differential cross section in microbarns/GeV^4/r.
   // The differential measure is d(sigma)/(dk dphi d^3 q) where k is
   // the energy of the bremsstrahlung photon and phi is the azimuthal angle
   // of the photon.  The polar angle of the photon is fixed by kinematics.
   // It is assumed that eIn.Mom()[0] = eOut.Mom()[0]+gOut.Mom()[0], that
   // the energy carried away by the recoil is zero in the laboratory frame,
   // but it is not checked. This result must be multiplied by the square
   // of the form factor of the atomic target and then integrated over
   // recoil kinematics to generate the photon spectrum.
   //
   // The way this is done below is as follows. First the recoil momentum
   // q is split into two transverse components collectively called qT
   // and longitudinal component qL, then qL is transformed into a Lorentz
   // invariant M2 = (P0 + q)^2 where P0 is the incident electron momentum
   // 4-vector and q is the four-vector (0, qTx, qTy, qL) of the momentum
   // exchange with the target atom: qL = -(M2 -mElectron^2 + qT^2)/(2 E0)
   // for a ultra-relativistic incident electron energy E0. The physical
   // meaning of M2 is the invariant mass-squared of the final electron-
   // photon system. Boosting the photon into this system, its lab energy
   // gets translated into an emission angle theta_cm with respect to the
   // boost vector. In combination with phi_cm = phi, generating a random
   // direction in the M2 rest frame becomes equivalent to generating a
   // random energy for the photon in the lab. The advantage of these new
   // coordinates is that their kinematic bounds are now independent from
   // one another. In terms of these new coordinates, the measure is now
   //
   //  dk dphi d^3 q = [(M2 - mElectron^2)/(4 M2)] d^2 qT dM2 dOmega_cm
   //
   // By this approach, the final photon and recoil electron kinematics
   // are generated first in a reaction rest frame, and then boosted into
   // the lab as a final step.

   int good_event;
   for (int event=0; event < nevents; event += good_event) {
      good_event = 0; 
      fRandom.flatArray(5, u);
      Ebeam = fBeamEnergy;
      weight = 1;

      // generate u = q0^2 / (q0^2 + qT^2) ~ Uniform[0,1]
      double q0 = 30.e-6; // 30 keV typical atomic scale
      qT = q0 * sqrt(1/u[0] - 1);
      qTphi = u[1] * 2*M_PI;
      weight *= sqr((sqr(q0) + sqr(qT)) / q0);

      // generate u = m0^2 / (m0^2 + M2 -mElectron^2) ~ Uniform[0,1]
      double m0sqr = sqr(0.3e-3); // 300 keV
      M2 = sqr(mElectron) + m0sqr * (1/u[2] - 1);
      qL = (sqr(mElectron) - M2 - sqr(qT)) / (2 * Ebeam);
      weight *= sqr(m0sqr + M2 - sqr(mElectron)) / m0sqr;

      // generate cm angles theta,phi uniform on the unit sphere
      double costheta = 2 * u[3] - 1;
      double sintheta = sqrt(1 - sqr(costheta));
      phi = u[4] * 2*M_PI;
      weight *= 4*M_PI;

      // compute the weight on the chosen measure
      weight *= (M2 - sqr(mElectron)) / (4 * M2);

      // compute the reaction kinematics based on what was generated above
      TFourVectorReal pin(Ebeam, 0, 0, sqrt(sqr(Ebeam) - sqr(mElectron)));
      TFourVectorReal q(0, qT * cos(qTphi), qT * sin(qTphi), qL);
      TThreeVectorReal pcm(sintheta * cos(phi), sintheta * sin(phi), costheta);
      TLorentzBoost toLab(pin + q);
      toLab.Invert();
      kstar = (M2 - sqr(mElectron)) / (2 * sqrt(M2));
      TFourVectorReal kout(kstar, kstar * pcm);
      kout.Boost(toLab);
      TFourVectorReal pout(sqrt(sqr(mElectron) + sqr(kstar)), -kstar * pcm);
      pout.Boost(toLab);

      if (kout[0] < fMinEnergy)
         continue;

      TLepton ein(pin, mElectron);
      TLepton eout(pout, mElectron);
      TPhoton gout(kout);
      TThreeVectorReal zeroVector(0,0,0);
      TThreeVectorReal posXhat(1,0,-kout[1]/kout[3]);
      TThreeVectorReal posUhat(1,1,-(kout[1] + kout[2])/kout[3]);
      posUhat /= sqrt(2.);
      ein.SetPol(zeroVector);
      eout.AllPol();
      gout.AllPol();
      double FF2 = sqr(AtomicFormFactor(q.LengthSqr()));
      diffXS = FF2 * TCrossSection::Bremsstrahlung(ein, eout, gout);
      gout.SetPol(posXhat);
      double polar0 = FF2 * TCrossSection::Bremsstrahlung(ein, eout, gout);
      gout.SetPol(posUhat);
      double polar45 = FF2 * TCrossSection::Bremsstrahlung(ein, eout, gout);

      for (int i=0; i < 4; ++i) {
         k[i] = kout[i];
         pR[i] = pout[i];
      }
      polar_0_90 = polar0 / diffXS;
      polar_45_135 = polar45 / diffXS;
      fTree->Fill();
      good_event = 1;
      if (event % 1000 == 0) {
         std::cout << fTree->GetEntries() << " events written\r"
                   << std::flush;
      }
   }
}

double GlueXBremsstrahlungGenerator::AtomicFormFactor(double q2)
{
   // This is a generic dipole model of a light atom, atomic number Z. 
   // The argument q2 should be in GeV^2, and the result is a
   // dimensionless coefficient that multiplies the charge of
   // the electron. Both electron cloud and nuclear contributions
   // to the atomic electric form factor are included.

   double Z = 13; // aluminum
   double beta = 111 * pow(Z, -1/3.) / mElectron;	// ff cutoff in /GeV
   double Fff = 1 / (1 + sqr(beta) * q2);
   return Z * (1 - Fff);
}

void GlueXBremsstrahlungGenerator::prepareImportanceSamplingPDFs()
{}

#endif
