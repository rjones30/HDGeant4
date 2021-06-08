//
// class implementation for GlueXBeamConversionProcess
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//

#include "GlueXBeamConversionProcess.hh"
#include "GlueXPhotonBeamGenerator.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserOptions.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4PairProductionRelModel.hh"

#define VERBOSE_PAIRS_SPLITTING 1
#define DO_PAIRCOH_IMPORTANCE_SAMPLE 1
#define USE_ADAPTIVE_SAMPLER 1

// If you set this flag to 1 then all beam photons that reach
// the TPOL converter target will convert to e+e- pairs inside,
// otherwise the standard pair conversion probabilities apply.
#define FORCED_PTAR_PAIR_CONVERSION 0

// If you set this flag to 1 then all beam photons that reach
// the LIH2 target will undergo Bethe-Heitler conversion to e+e-
// pairs inside, otherwise the standard pair conversion
// probabilities apply.
#define FORCED_LIH2_PAIR_CONVERSION 0

#include <G4SystemOfUnits.hh>
#include "G4PhysicalConstants.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4RunManager.hh"
#include "G4TrackVector.hh"

#include <stdio.h>
#include <iomanip>

#ifdef USING_DIRACXX
#include <TLorentzBoost.h>
#include <TPhoton.h>
#include <TLepton.h>
#include <TCrossSection.h>

#endif

void unif01(int n, double *u) {
   G4Random::getTheEngine()->flatArray(n,u);
}

G4double GlueXBeamConversionProcess::fBHpair_mass_min = 0;

std::vector<AdaptiveSampler*> GlueXBeamConversionProcess::fAdaptiveSamplerRegistry;

G4Mutex GlueXBeamConversionProcess::fMutex = G4MUTEX_INITIALIZER;
int GlueXBeamConversionProcess::fStopBeamBeforeConverter = 0;
int GlueXBeamConversionProcess::fStopBeamAfterConverter = 0;
int GlueXBeamConversionProcess::fStopBeamAfterTarget = 0;
int GlueXBeamConversionProcess::fConfigured = 0;

GlueXBeamConversionProcess::GlueXBeamConversionProcess(const G4String &name, 
                                                       G4ProcessType aType)
 : G4VEmProcess(name, aType),
#  ifdef USING_DIRACXX
   fPairsGeneration(0),
#  endif
   fPaircohPDF(0),
   fTripletPDF(0),
   fAdaptiveSampler(0),
   isInitialised(false)
{
   verboseLevel = 0;

   G4AutoLock barrier(&fMutex);
   if (fConfigured)
      return;

   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXBeamConversionProcess constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   fStopBeamBeforeConverter = 0;
   fStopBeamAfterConverter = 0;
   fStopBeamAfterTarget = 0;

   std::map<int,std::string> infile;
   std::map<int,double> beampars;
   std::map<int,std::string> genbeampars;
   if (!user_opts->Find("INFI", infile) &&
       user_opts->Find("BEAM", beampars) &&
       user_opts->Find("GENBEAM", genbeampars))
   {
      if (genbeampars.find(1) != genbeampars.end() &&
         (genbeampars[1] == "POSTCOL" ||
          genbeampars[1] == "postcol" ||
          genbeampars[1] == "Postcol" ||
          genbeampars[1] == "PostCol" ))
      {
         fStopBeamBeforeConverter = 1;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "POSTCONV" ||
               genbeampars[1] == "postconv" ||
               genbeampars[1] == "Postconv" ||
               genbeampars[1] == "PostConv" ))
      {
         fStopBeamAfterConverter = 1;
      }
      else if (genbeampars.find(1) != genbeampars.end() &&
              (genbeampars[1] == "BHgen" ||
               genbeampars[1] == "BHGEN" ||
               genbeampars[1] == "BHGen" ||
               genbeampars[1] == "bhgen" ))
      {
         fStopBeamAfterTarget = 1;
         if (genbeampars.find(2) != genbeampars.end()) {
            fBHpair_mass_min = std::atof(genbeampars[2].c_str());
         }
      }
   }

   fConfigured = 1;

   if (verboseLevel > 0) {
       G4cout << GetProcessName() << " is created " << G4endl
	      << "    Stop beam before converter? "
	      << (fStopBeamBeforeConverter? "yes" : "no") << G4endl
	      << "    Stop beam after converter? "
	      << (fStopBeamAfterConverter? "yes" : "no") << G4endl
	      << "    Stop beam after target? "
	      << (fStopBeamAfterTarget? "yes" : "no") << G4endl;
   }
}

GlueXBeamConversionProcess::~GlueXBeamConversionProcess()
{
   if (fPaircohPDF)
      delete fPaircohPDF;
   if (fTripletPDF)
      delete fTripletPDF;
#ifdef USING_DIRACXX
   if (fPairsGeneration)
      delete fPairsGeneration;
#endif
   G4AutoLock barrier(&fMutex);
   int nsamplers = fAdaptiveSamplerRegistry.size();
   if (nsamplers > 0) {
      fAdaptiveSampler = fAdaptiveSamplerRegistry[0];
      for (int i=nsamplers-1; i > 0; --i) {
         std::stringstream astatefile;
         astatefile << "BHgen_thread_" << i << ".astate";
         fAdaptiveSamplerRegistry[i]->saveState(astatefile.str());
         fAdaptiveSampler->mergeState(astatefile.str());
         delete fAdaptiveSamplerRegistry[i];
         fAdaptiveSamplerRegistry.erase(fAdaptiveSamplerRegistry.begin() + i);
      }
      fAdaptiveSamplerRegistry[0]->saveState("BHgen_stats.astate");
      delete fAdaptiveSamplerRegistry[0];
      fAdaptiveSamplerRegistry.clear();
   }
}

G4bool GlueXBeamConversionProcess::IsApplicable(const G4ParticleDefinition& p)
{
   return (&p == G4Gamma::Gamma());
}

void GlueXBeamConversionProcess::InitialiseProcess(const G4ParticleDefinition*)
{
   if (!isInitialised) {
      isInitialised = true;
      G4EmParameters* param = G4EmParameters::Instance();
      G4double emin = std::max(param->MinKinEnergy(), 2*electron_mass_c2);
      SetMinKinEnergy(emin);

      if (!EmModel(0)) {
         SetEmModel(new G4PairProductionRelModel());
      }
      EmModel(0)->SetLowEnergyLimit(emin);
      AddEmModel(1, EmModel(0));

#ifdef USING_DIRACXX
      std::vector<double> Z;
      std::vector<double> A;
      std::vector<double> w;
      Z.push_back(4);    // assumes 4Be converter target for TPOL
      A.push_back(9.012);
      w.push_back(1);
      fPairsGeneration = new PairConversionGeneration(Z, A, w);

#if defined DO_TRIPLET_IMPORTANCE_SAMPLE || defined DO_PAIRCOH_IMPORTANCE_SAMPLE
      if (fTripletPDF == 0 || fPaircohPDF == 0) {
         G4cout << "GlueXBeamConversionProcess::InitialiseProcess:"
                << G4endl
                << "   Setting up cross section tables, please wait... "
                << std::flush;
         prepareImportanceSamplingPDFs();
         G4cout << "finished." << G4endl;
      }
#endif

#if USE_ADAPTIVE_SAMPLER
      if (fAdaptiveSampler == 0) {
         prepareAdaptiveSampler();
      }
#endif

#endif
   }
}

G4double GlueXBeamConversionProcess::MinPrimaryEnergy(const G4ParticleDefinition*,
                                                      const G4Material*)
{
  return 2 * electron_mass_c2;
}

void GlueXBeamConversionProcess::PrintInfo()
{}         


void GlueXBeamConversionProcess::ProcessDescription(std::ostream& out) const
{
  out << "  Dirac++ Gamma conversion (forced)";
  G4VEmProcess::ProcessDescription(out);
}


G4double GlueXBeamConversionProcess::GetMeanFreePath(
                                     const G4Track &track, 
                                     G4double previousStepSize,
                                     G4ForceCondition *condition)
{
   return 100*cm;
}

G4double GlueXBeamConversionProcess::PostStepGetPhysicalInteractionLength(
                                     const G4Track &track,
                                     G4double previousStepSize,
                                     G4ForceCondition *condition)
{
   const G4Step *step = G4ParallelWorldProcess::GetHyperStep();
   G4VPhysicalVolume *pvol = step->GetPostStepPoint()->GetPhysicalVolume();
   if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "PTAR" &&
       (FORCED_PTAR_PAIR_CONVERSION || 
        fStopBeamBeforeConverter ||
        fStopBeamAfterConverter ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      *condition = Forced;
      return 100*cm;
   }
   else if (track.GetTrackID() == 1 && pvol && pvol->GetName() == "LIH2" &&
       (FORCED_LIH2_PAIR_CONVERSION || 
        fStopBeamAfterTarget ))
   {
      fPIL = G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             track, previousStepSize, condition);
      *condition = Forced;
      return 100*cm;
   }
   *condition = NotForced;
   return 1e99;
}

G4VParticleChange *GlueXBeamConversionProcess::PostStepDoIt(
                                               const G4Track &track, 
                                               const G4Step &step)
{
   pParticleChange->Initialize(track);
   GlueXUserEventInformation *eventinfo;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   eventinfo = (GlueXUserEventInformation*)event->GetUserInformation();
   if (fStopBeamBeforeConverter) {
      double tvtx = step.GetPreStepPoint()->GetGlobalTime();
      G4ThreeVector vtx = step.GetPreStepPoint()->GetPosition();
      G4ThreeVector mom = step.GetPreStepPoint()->GetMomentum();
      G4ThreeVector pol = step.GetPreStepPoint()->GetPolarization();
      G4PrimaryVertex* vertex = new G4PrimaryVertex(vtx, tvtx);
      G4PrimaryParticle* photon = new G4PrimaryParticle(G4Gamma::Definition(),
                                                        mom[0], mom[1], mom[2]);
      photon->SetPolarization(pol);
      vertex->SetPrimary(photon);
      eventinfo->AddPrimaryVertex(*vertex);
      eventinfo->AddBeamParticle(1, tvtx, vtx, mom, pol);

      if (verboseLevel > 0) {
         G4cout << "GlueXBeamConversionProcess: beam particle stopped"
                << " before converter, stored in ouptut primary vertex."
                << G4endl;
      }
   }
   else if (fStopBeamAfterConverter) {
      double tvtx = step.GetPreStepPoint()->GetGlobalTime();
      G4ThreeVector vtx = step.GetPreStepPoint()->GetPosition();
      G4ThreeVector mom = step.GetPreStepPoint()->GetMomentum();
      G4ThreeVector pol = step.GetPreStepPoint()->GetPolarization();
      eventinfo->AddBeamParticle(1, tvtx, vtx, mom, pol);
      GenerateBeamPairConversion(step);

      if (verboseLevel > 0) {
         G4cout << "GlueXBeamConversionProcess: beam particle stopped"
                << " at converter exit, pair conversion forced."
                << G4endl;
      }
   }
   else if (fStopBeamAfterTarget) {
      double tvtx = step.GetPreStepPoint()->GetGlobalTime();
      G4ThreeVector vtx = step.GetPreStepPoint()->GetPosition();
      G4ThreeVector mom = step.GetPreStepPoint()->GetMomentum();
      G4ThreeVector pol = step.GetPreStepPoint()->GetPolarization();
      eventinfo->AddBeamParticle(1, tvtx, vtx, mom, pol);
      double BetheHeitler_fraction = 0.5;
      double maxTripletMass2 = 2 * 0.511*MeV *
                               step.GetPreStepPoint()->GetKineticEnergy();
      if (fBHpair_mass_min*GeV > sqrt(maxTripletMass2))
         BetheHeitler_fraction = 1;
      double weight_factor;
      if (G4UniformRand() < BetheHeitler_fraction) {
         if (verboseLevel > 0) {
            G4cout << "GlueXBeamConversionProcess: beam particle stopped"
                   << " in the LiH2 target, Bethe-Heitler conversion forced."
                   << G4endl;
         }
         GenerateBetheHeitlerProcess(step);
         weight_factor = 1 / BetheHeitler_fraction;
      }
      else {
         if (verboseLevel > 0) {
            G4cout << "GlueXBeamConversionProcess: beam particle stopped"
                   << " in the LiH2 target, triplet conversion forced."
                   << G4endl;
         }
         GenerateTripletProcess(step);
         weight_factor = 1 / (1. - BetheHeitler_fraction);
      }
      GlueXUserEventInformation *event_info;
      const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
      event_info = (GlueXUserEventInformation*)event->GetUserInformation();
      hddm_s::ReactionList rea = event_info->getOutputRecord()->getReactions();
      rea(0).setWeight(rea(0).getWeight() * weight_factor);
   }
   else {
      GenerateBeamPairConversion(step);

      if (verboseLevel > 0) {
         G4cout << "GlueXBeamConversionProcess: beam particle stopped"
                << " unexpectedly, pair conversion forced. But WHY?"
                << G4endl;
      }
   }
   pParticleChange->ProposeTrackStatus(fStopAndKill);
   eventinfo->SetKeepEvent(1);
   return pParticleChange;
}

void GlueXBeamConversionProcess::prepareImportanceSamplingPDFs()
{
   // Construct lookup tables representing the PDFs used for
   // importance-sampling the gamma pair conversion kinematics.

   // Compute 2D histogram containing rate as a function of u0,u1
   // where u0 is the (originally uniform [0,1]) random number used
   // to generate Mpair and u1 generates qR. The algorithm assumes
   // that the mapping u0->Mpair and u1->qR used here is the same
   // as the mmpaing used in the process generation methods, currently
   // GenerateBeamPairConversion and GenerateTargetConversionProcess.

   if (fPaircohPDF == 0 || fTripletPDF == 0) {
      fPaircohPDF = new ImportanceSampler;
      fTripletPDF = new ImportanceSampler;
   }

   fPaircohPDF->Pcut = 60;
   fTripletPDF->Pcut = 15;

#if USING_DIRACXX

   LDouble_t kin = 9.; // GeV
   TPhoton g0;
   TLepton p1(mElectron);
   TLepton e2(mElectron);
   TLepton e3(mElectron);
   TThreeVectorReal p;
   g0.SetMom(p.SetPolar(kin,0,0));
   g0.SetPol(TThreeVector(0,0,0));
   p1.AllPol();
   e2.AllPol();
   e3.AllPol();

   LDouble_t Epos = kin / 2;
   LDouble_t Mmin = 2 * mElectron;
   LDouble_t Mcut = 5e-3; // 5 MeV cutoff parameter
   LDouble_t um0 = 1 + sqr(Mcut / Mmin);
   LDouble_t qRcut = 1e-3; // 1 MeV/c cutoff parameter

   int Nbins = 50;
   fTripletPDF->Psum = 0;
   fPaircohPDF->Psum = 0;
   for (int i0=0; i0 < Nbins; ++i0) {
      LDouble_t u0 = (i0 + 0.5) / Nbins;
      LDouble_t um = pow(um0, u0);
      LDouble_t Mpair = Mcut / sqrt(um - 1);
      LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mElectron);
      LDouble_t qRmin = sqr(Mpair) / (2 * kin);
      LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
      LDouble_t weight0 = sqr(Mpair) * (sqr(Mcut) + sqr(Mpair));
      for (int i1=0; i1 < Nbins; ++i1) {
         LDouble_t u1 = (i1 + 0.5) / Nbins;
         LDouble_t uq = pow(uq0, u1);
         LDouble_t qR = 2 * qRcut * uq / (1 - sqr(uq));
         LDouble_t weight = weight0 * sqr(qR) * sqrt(sqr(qRcut) + sqr(qR));
         LDouble_t E3 = sqrt(sqr(qR) + sqr(mElectron));
         LDouble_t E12 = kin + mElectron - E3;
         if (k12star2 < 0 || E12 < Mpair) {
            fPaircohPDF->density.push_back(0);
            fTripletPDF->density.push_back(0);
            fPaircohPDF->integral.push_back(fPaircohPDF->Psum);
            fTripletPDF->integral.push_back(fTripletPDF->Psum);
            continue;
         }
         LDouble_t k12star = sqrt(k12star2);
         LDouble_t q12mag = sqrt(sqr(E12) - sqr(Mpair));
         LDouble_t costhetastar = (Epos - E12 / 2) * 
                                  Mpair / (k12star * q12mag);
         LDouble_t costhetaR = (sqr(Mpair) / 2 + (kin + mElectron) *
                                (E3 - mElectron)) / (kin * qR);
         if (fabs(costhetastar) > 1 || fabs(costhetaR) > 1) {
            fPaircohPDF->density.push_back(0);
            fTripletPDF->density.push_back(0);
            fPaircohPDF->integral.push_back(fPaircohPDF->Psum);
            fTripletPDF->integral.push_back(fTripletPDF->Psum);
            continue;
         }
         LDouble_t qRlong = qR * costhetaR;
         LDouble_t qRperp = sqrt(sqr(qR) - sqr(qRlong));
         TThreeVectorReal q3(0, qRperp, qRlong);
         LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
         TThreeVectorReal k12(k12star * sinthetastar, 0, 
                              k12star * costhetastar);
         TFourVectorReal q1(Mpair / 2, k12);
         TFourVectorReal q2(Mpair / 2, -k12);
         TLorentzBoost tolab(q3[1] / E12, q3[2] / E12, (q3[3] - kin) / E12);
         q1.Boost(tolab);
         q2.Boost(tolab);
         p1.SetMom(q1);
         e2.SetMom(q2);
         e3.SetMom(q3);
         LDouble_t tripXS = fPairsGeneration->DiffXS_triplet(g0,p1,e2,e3);
         LDouble_t pairXS = fPairsGeneration->DiffXS_pair(g0,p1,e2);
         fTripletPDF->Psum += fTripletPDF->Pmax = tripXS * weight;
         fPaircohPDF->Psum += fPaircohPDF->Pmax = pairXS * weight;
         fTripletPDF->density.push_back(fTripletPDF->Pmax);
         fPaircohPDF->density.push_back(fPaircohPDF->Pmax);
         fPaircohPDF->integral.push_back(fPaircohPDF->Psum);
         fTripletPDF->integral.push_back(fTripletPDF->Psum);
      }
   }

   LDouble_t du2 = 1. / sqr(Nbins);
   for (int i0=0, index=0; i0 < Nbins; ++i0) {
      for (int i1=0; i1 < Nbins; ++i1, ++index) {
         LDouble_t randvar = i0 + (i1 + 0.5) / Nbins;
         fTripletPDF->randvar.push_back(randvar);
         fPaircohPDF->randvar.push_back(randvar);
         fTripletPDF->density[index] /= fTripletPDF->Psum * du2;
         fPaircohPDF->density[index] /= fPaircohPDF->Psum * du2;
         fTripletPDF->integral[index] /= fTripletPDF->Psum;
         fPaircohPDF->integral[index] /= fPaircohPDF->Psum;
      }
   }
#endif

   fTripletPDF->Pmax = 0;
   fTripletPDF->Psum = 0;
   fPaircohPDF->Pmax = 0;
   fPaircohPDF->Psum = 0;
}

void GlueXBeamConversionProcess::GenerateBeamPairConversion(const G4Step &step)
{
   // Unlike the other GenerateXXX methods in this class, this method should
   // be invoked after tracking of an event is already underway. Its purpose
   // is to generate pair conversion / triplet production by gamma rays in a
   // converter target with the full polarization-dependent QED differential
   // cross section. Both incoherent (triplet) and nuclear + atomic coherent
   // scattering is simulated by this method. It should be called from your
   // G4SteppingAction::UserSteppingAction method to force pair conversion
   // in special simulations dedicated to the study of the analyzing power
   // of this reaction. The incident gamma ray is stopped and the pair (plus
   // recoil electron, if any) are added to the event stack as new primary
   // particles.

#ifndef USING_DIRACXX

   G4cerr << "GlueXBeamConversionProcess::GenerateBeamPairConversion error:"
          << G4endl
          << "  You have enabled pair/triplet conversion in the PTAR target,"
          << G4endl
          << "  but your have built HDGeant4 without support for the Dirac++"
		  << G4endl
          << "  library. Either rebuild with Dirac++ support or else jack up"
		  << G4endl
          << "  the value of constant FORCED_PTAR_PAIR_CONVERSION_THRESHOLD"
		  << G4endl
          << "  in GlueXSteppingAction.cc and try again. Aborting this run..."
		  << G4endl;
   exit(1);

#else

   TPhoton gIn;
   TLepton p1(mElectron);
   TLepton e2(mElectron);
   TLepton e3(mElectron);
   p1.AllPol();
   e2.AllPol();
   e3.AllPol();
   const G4Track *track = step.GetTrack();
   LDouble_t kin = track->GetKineticEnergy()/GeV;

   // If we are below pair production threshold, do nothing
   const double mTarget = step.GetPreStepPoint()->GetMaterial()->GetA() * 0.932;
   if (kin < 2 * mElectron * (1 + mElectron / mTarget))
      return;

   G4ThreeVector mom(track->GetMomentum());
   TThreeVectorReal mom0(mom[0]/GeV, mom[1]/GeV, mom[2]/GeV);
   gIn.SetMom(mom0);
   G4ThreeVector pol(track->GetPolarization());
   TThreeVectorReal pol0(pol[0], pol[1], pol[2]);
   gIn.SetPlanePolarization(pol0, pol0.Length());

   // Define an angle and axis that rotates zhat into the direction
   // of the incident gamma, so that the generated kinematics is
   // defined with the incident gamma aligned with zhat, and then
   // rotated at the end into the final spatial direction.
   TThreeVectorReal rockaxis(mom0);
   rockaxis.Cross(TThreeVectorReal(0,0,1));
   double rockangle = asin(rockaxis.Length() / mom0.Length());
   rockaxis /= rockaxis.Length();

   while (true) {
      LDouble_t weight = 1;

      // Generate uniform in E+, phi12, phiR
      LDouble_t Epos = kin * G4UniformRand();
      while (Epos < mElectron) {
         Epos = kin * G4UniformRand();
      }
      weight *= kin - mElectron;
      LDouble_t phi12 = 2*M_PI * G4UniformRand();
      weight *= 2*M_PI;
      LDouble_t phiR = 2*M_PI * G4UniformRand();
      weight *= 2*M_PI;

      LDouble_t u0 = G4UniformRand();
      LDouble_t u1 = G4UniformRand();

#if DO_TRIPLET_IMPORTANCE_SAMPLE

      int i = fTripletPDF->search(u1);
      LDouble_t fi = fTripletPDF->density[i];
      LDouble_t ui = fTripletPDF->integral[i];
      LDouble_t ri = fTripletPDF->randvar[i];
      LDouble_t xi = ri - floor(ri);
      LDouble_t dx = (xi > 0.5)? ri - fTripletPDF->randvar[i-1]:
                              fTripletPDF->randvar[i+1] - ri;
      u1 = xi + dx / 2 - (ui - u1) / (fi * dx);
      u0 = (u0 + floor(ri)) * dx;
      weight /= fi;

#elif DO_PAIRCOH_IMPORTANCE_SAMPLE

      int i = fPaircohPDF->search(u1);
      LDouble_t fi = fPaircohPDF->density[i];
      LDouble_t ui = fPaircohPDF->integral[i];
      LDouble_t ri = fPaircohPDF->randvar[i];
      LDouble_t xi = ri - floor(ri);
      LDouble_t dx = (xi > 0.5)? ri - fPaircohPDF->randvar[i-1]:
                              fPaircohPDF->randvar[i+1] - ri;
      u1 = xi + dx / 2 - (ui - u1) / (fi * dx);
      u0 = (u0 + floor(ri)) * dx;
      weight /= fi;

#endif
   
      // Generate Mpair as 1 / (M [M^2 + Mcut^2])
      LDouble_t Mmin = 2 * mElectron;
      LDouble_t Mcut = 0.005;  // GeV
      LDouble_t um0 = 1 + sqr(Mcut / Mmin);
      LDouble_t um = pow(um0, u0);
      LDouble_t Mpair = Mcut / sqrt(um - 1 + 1e-99);
      weight *= Mpair * (sqr(Mcut) + sqr(Mpair)) * log(um0) / (2 * sqr(Mcut));
   
      // Generate qR^2 with weight 1 / [qR^2 sqrt(qRcut^2 + qR^2)]
      LDouble_t qRmin = sqr(Mpair) /(2 * kin);
      LDouble_t qRcut = 1e-3; // GeV
      LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
      LDouble_t uq = pow(uq0, u1);
      LDouble_t qR = 2 * qRcut * uq / (1 - sqr(uq));
      LDouble_t qR2 = qR * qR;
      weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut)) * (-2 * log(uq0));
   
      // Include overall measure Jacobian factor
      weight *= Mpair / (2 * kin);
   
      // Generate with importance sampling
      LDouble_t Striplet = fTripletPDF->Npassed * 
                        (fTripletPDF->Ntested / (fTripletPDF->Psum + 1e-99));
      LDouble_t Spaircoh = fPaircohPDF->Npassed *
                        (fPaircohPDF->Ntested / (fPaircohPDF->Psum + 1e-99));
      if (Striplet < Spaircoh) {                     // try incoherent generation
         ++fTripletPDF->Ntested;
   
         // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
         LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mElectron);
         if (k12star2 < 0) {
            // try again (this should never happen!)
            continue;
         }
         LDouble_t k12star = sqrt(k12star2);
         LDouble_t E3 = sqrt(qR2 + sqr(mElectron));
         LDouble_t E12 = kin + mElectron - E3;
         if (E12 < Mpair) {
            // no kinematic solution because E12 < Mpair, try again
            continue;
         }
         LDouble_t q12mag = sqrt(sqr(E12) - sqr(Mpair));
         LDouble_t costhetastar = (Epos - E12 / 2) * Mpair / (k12star * q12mag);
         if (Epos > E12 - mElectron) {
            // no kinematic solution because Epos > E12 - mElectron, try again
            continue;
         }
         else if (fabs(costhetastar) > 1) {
            // no kinematic solution because |costhetastar| > 1, try again
            continue;
         }
   
         // Solve for the recoil electron kinematics
         LDouble_t costhetaR = (sqr(Mpair) / 2 + (kin + mElectron) *
                             (E3 - mElectron)) / (kin * qR);
         if (fabs(costhetaR) > 1) {
            // no kinematic solution because |costhetaR| > 1, try again
            continue;
         }
         LDouble_t sinthetaR = sqrt(1 - sqr(costhetaR));
         TFourVectorReal q3(E3, qR * sinthetaR * cos(phiR),
                                qR * sinthetaR * sin(phiR),
                                qR * costhetaR);
   
         // Boost the pair momenta into the lab
         LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
         TThreeVectorReal k12(k12star * sinthetastar * cos(phi12),
                              k12star * sinthetastar * sin(phi12),
                              k12star * costhetastar);
         TFourVectorReal q1(Mpair / 2, k12);
         TFourVectorReal q2(Mpair / 2, -k12);
         TLorentzBoost toLab(q3[1] / E12, q3[2] / E12, (q3[3] - kin) / E12);
         q1.Boost(toLab);
         q2.Boost(toLab);
   
         // To avoid double-counting, return zero if recoil electron
         // momentum is greater than the momentum of the pair electron.
         if (q2.Length() < qR) {
            // recoil/pair electrons switched, try again
            continue;
         }

         // Rotate final-state momenta into the lab frame
         p1.SetMom(q1.Rotate(rockaxis, rockangle));
         e2.SetMom(q2.Rotate(rockaxis, rockangle));
         e3.SetMom(q3.Rotate(rockaxis, rockangle));

         // Check 4-momentum conservation
         TFourVectorReal pIn, pFi;
         TFourVectorReal::SetResolution(1e-10);
         pIn = gIn.Mom() + TFourVectorReal(mElectron,0,0,0);
         pFi = p1.Mom() + e2.Mom() + e3.Mom();
         if (pIn != pFi) {
            std::cout << "Warning in GenerateBeamPairConversion - "
                      << "momentum conservation violated." << std::endl
                      << "   pIn = "; 
            pIn.Print();
            std::cout << "   pFi = ";
            pFi.Print();
            std::cout << "   pIn - pFi = ";
            (pIn-pFi).Print();
         }

         // Compute the differential cross section (barns/GeV^4)
         // returned as d(sigma)/(dE+ dphi+ d^3qR)
         LDouble_t diffXS = fPairsGeneration->DiffXS_triplet(gIn, p1, e2, e3);

         // Use keep/discard algorithm
         LDouble_t Pfactor = diffXS * weight;
         if (Pfactor > fTripletPDF->Pmax)
            fTripletPDF->Pmax = Pfactor;
         if (Pfactor > fTripletPDF->Pcut) {
            G4cout << "Warning in GenerateBeamPairConversion - Pfactor " 
                   << Pfactor << " exceeds fTripletPDF->Pcut = " 
                   << fTripletPDF->Pcut << G4endl
                   << "  present qR = " << qR << G4endl
                   << "  present Mpair = " << Mpair << G4endl
                   << "  present maximum Pfactor = "
                   << fTripletPDF->Pmax << G4endl
                   << "  current generator efficiency = "
                   << fTripletPDF->Npassed /
                      (fTripletPDF->Ntested + 1e-99)
                   << G4endl;
         }
         fTripletPDF->Psum += Pfactor;
         if (G4UniformRand() * fTripletPDF->Pcut > Pfactor) {
            continue;
         }
         ++fTripletPDF->Npassed;
         break;
      }

      else {                          // try coherent generation
         ++fPaircohPDF->Ntested;
   
         // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
         // assuming that the atomic target absorbs zero energy
         LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mElectron);
         if (k12star2 < 0) {
            // try again (this should never happen!)
            continue;
         }
         LDouble_t k12star = sqrt(k12star2);
         LDouble_t Eele = kin - Epos;
         if (kin < Mpair) {
            // no kinematic solution because kin < Mpair, try again
            continue;
         }
         else if (Eele < mElectron) {
            // no kinematic solution because Eele < mElectron, try again
            continue;
         }
         LDouble_t q12mag = sqrt(sqr(kin) - sqr(Mpair));
         LDouble_t costhetastar = (Epos - kin / 2) * Mpair / (k12star * q12mag);
         if (fabs(costhetastar) > 1) {
            // no kinematic solution because |costhetastar| > 1, try again
            continue;
         }
   
         // Solve for the recoil kinematics kinematics
         LDouble_t costhetaR = (sqr(Mpair) + qR2) / (2 * kin * qR);
         if (fabs(costhetaR) > 1) {
            // no kinematic solution because |costhetaR| > 1, try again
            continue;
         }
         LDouble_t sinthetaR = sqrt(1 - sqr(costhetaR));
         TThreeVectorReal q3(qR * sinthetaR * cos(phiR),
                             qR * sinthetaR * sin(phiR),
                             qR * costhetaR);
   
         // Boost the pair momenta into the lab
         LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
         TThreeVectorReal k12(k12star * sinthetastar * cos(phi12),
                              k12star * sinthetastar * sin(phi12),
                              k12star * costhetastar);
         TFourVectorReal q1(Mpair / 2, k12);
         TFourVectorReal q2(Mpair / 2, -k12);
         TLorentzBoost toLab(q3[1] / kin, q3[2] / kin, (q3[3] - kin) / kin);
         q1.Boost(toLab);
         q2.Boost(toLab);
   
         // Compute the differential cross section (barnes/GeV^4)
         // returned as d(sigma)/(dE+ dphi+ d^3qR)
         p1.SetMom(q1.Rotate(rockaxis, rockangle));
         e2.SetMom(q2.Rotate(rockaxis, rockangle));
         e3.SetMom(TThreeVectorReal(0,0,0));
         LDouble_t diffXS = fPairsGeneration->DiffXS_pair(gIn, p1, e2);
   
         // Use keep/discard algorithm
         LDouble_t Pfactor = diffXS * weight;
         if (Pfactor > fPaircohPDF->Pmax)
            fPaircohPDF->Pmax = Pfactor;
         if (Pfactor > fPaircohPDF->Pcut) {
            G4cout << "Warning in GenerateBeamPairConversion - Pfactor " 
                   << Pfactor << " exceeds fPaircohPDF->Pcut = " 
                   << fPaircohPDF->Pcut << G4endl
                   << "  present qR = " << qR << G4endl
                   << "  present Mpair = " << Mpair << G4endl
                   << "  present maximum Pfactor = "
                   << fPaircohPDF->Pmax << G4endl
                   << "  current generator efficiency = "
                   << fPaircohPDF->Npassed /
                      (fPaircohPDF->Ntested + 1e-99)
                   << G4endl;
         }
         fPaircohPDF->Psum += Pfactor;
         if (G4UniformRand() * fPaircohPDF->Pcut > Pfactor) {
            continue;
         }
         ++fPaircohPDF->Npassed;
         break;
      }
   }

   // Generate new primaries for the pair
   double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
   double steplength = pParticleChange->GetTrueStepLength();
   G4ThreeVector direction(track->GetMomentumDirection());
   G4ThreeVector x0(track->GetPosition());
   double t0 = track->GetGlobalTime();
   double uvtx = G4UniformRand();
   double lvtx = steplength + fPIL * log(1 - uvtx * (1 - exp(-steplength / fPIL)));
   x0 -= lvtx * direction;
   t0 -= lvtx / beamVelocity;
   G4PrimaryVertex vertex(x0, t0);
   G4ParticleDefinition *positron = G4Positron::Definition();
   G4ParticleDefinition *electron = G4Electron::Definition();
   G4ThreeVector psec1(p1.Mom()[1]*GeV, p1.Mom()[2]*GeV, p1.Mom()[3]*GeV);
   G4ThreeVector psec2(e2.Mom()[1]*GeV, e2.Mom()[2]*GeV, e2.Mom()[3]*GeV);
   G4ThreeVector psec3(e3.Mom()[1]*GeV, e3.Mom()[2]*GeV, e3.Mom()[3]*GeV);
   G4DynamicParticle *sec1 = new G4DynamicParticle(positron, psec1);
   G4DynamicParticle *sec2 = new G4DynamicParticle(electron, psec2);
   G4DynamicParticle *sec3 = new G4DynamicParticle(electron, psec3);
   G4TrackVector secondaries;
   secondaries.push_back(new G4Track(sec1, t0, x0));
   secondaries.push_back(new G4Track(sec2, t0, x0));
   if (e3.Mom().Length() > 0) {
      secondaries.push_back(new G4Track(sec3, t0, x0));
   }

   GlueXUserEventInformation *event_info;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   event_info = (GlueXUserEventInformation*)event->GetUserInformation();
   G4TrackVector::iterator iter;
   for (iter = secondaries.begin(); iter != secondaries.end(); ++iter) {
      GlueXUserTrackInformation *trackinfo = new GlueXUserTrackInformation();
      if (event_info) {
         trackinfo->SetGlueXTrackID(event_info->AssignNextGlueXTrackID());
      }
      (*iter)->SetUserInformation(trackinfo);
      if (fStopBeamAfterConverter == 0) {
         pParticleChange->AddSecondary(*iter);
      }
   }

   // append secondary vertex to MC record
   if (event_info) {
      int mech[2];
      char *cmech = (char*)mech;
      snprintf(cmech, 5, "%c%c%c%c", 'B', 'E', 'A', 'M');
      event_info->AddSecondaryVertex(secondaries, 1, mech[0]);
   }

#if VERBOSE_PAIRS_SPLITTING
   if ((fTripletPDF->Npassed + fPaircohPDF->Npassed) % 500 == 0) {
      G4cout << std::setprecision(5)
             << "triplet cross section is " << fTripletPDF->Pcut *
             fTripletPDF->Npassed / (fTripletPDF->Ntested + 1e-99) 
             << " +/- " << fTripletPDF->Pcut *
             sqrt(fTripletPDF->Npassed) / (fTripletPDF->Ntested + 1e-99) 
             << " barns, efficiency is " 
             << fTripletPDF->Npassed / (fTripletPDF->Ntested + 1e-99)
             << G4endl
             << "pair cross section is " << fPaircohPDF->Pcut *
             fPaircohPDF->Npassed / (fPaircohPDF->Ntested + 1e-99) 
             << " +/- " << fPaircohPDF->Pcut *
             sqrt(fPaircohPDF->Npassed) / (fPaircohPDF->Ntested + 1e-99) 
             << " barns, efficiency is " 
             << fPaircohPDF->Npassed / (fPaircohPDF->Ntested + 1e-99) 
             << G4endl
             << "counts are "
             << fTripletPDF->Npassed << " / " << fPaircohPDF->Npassed
             << " = "
             << fTripletPDF->Npassed / (fPaircohPDF->Npassed + 1e-99)
             << G4endl;
   }
#endif

#endif
}

void GlueXBeamConversionProcess::GenerateBetheHeitlerProcess(const G4Step &step)
{
   // Unlike the other GenerateXXX methods in this class, this method should
   // be invoked after tracking of an event is already underway. Its purpose
   // is to generate Bethe-Heitler pair conversion of the beam photon in the
   // primary GlueX target with the full polarization-dependent QED differential
   // cross section. A minimum invariant mass of the pair may be specified in
   // the member variable fBHpair_mass_min. It should be called from your
   // G4SteppingAction::UserSteppingAction method to force pair conversion
   // in special simulations dedicated to the study of the Bethe-Heitler
   // process. The incident gamma ray is stopped and the e+/e- pair vertex,
   // plus the recoil proton, is added to the output event record.

#ifndef USING_DIRACXX

   G4cerr << "GlueXBeamConversionProcess::GenerateBetheHeitlerProcess error:"
          << G4endl
          << "  You have enabled Bethe-Heitler conversion in the LIH2 target,"
          << G4endl
          << "  but your have built HDGeant4 without support for the Dirac++"
		  << G4endl
          << "  library. Either rebuild with Dirac++ support or else turn off"
		  << G4endl
          << "  this process. Aborting this run..."
		  << G4endl;
   exit(1);

#else

   TPhoton gIn;
   TLepton eOut(mElectron);
   TLepton pOut(mElectron);
   TLepton nIn(mProton);
   TLepton nOut(mProton);
   LDouble_t diffXS(0);
   eOut.AllPol();
   pOut.AllPol();
   nOut.AllPol();
   const TThreeVectorReal zero3vector(0,0,0);
   nIn.SetPol(TThreeVectorReal(zero3vector));
   nIn.SetMom(TThreeVectorReal(zero3vector));
   const G4Track *track = step.GetTrack();
   LDouble_t kin = track->GetKineticEnergy()/GeV;

   // If we are below pair production threshold, do nothing
   if (kin < 2 * mElectron * (1 + mElectron / mProton))
      return;

   G4ThreeVector mom(track->GetMomentum());
   TThreeVectorReal mom0(mom[0]/GeV, mom[1]/GeV, mom[2]/GeV);
   gIn.SetMom(mom0);
   G4ThreeVector pol(track->GetPolarization());
   TThreeVectorReal pol0(pol[0], pol[1], pol[2]);
   gIn.SetPlanePolarization(pol0, pol0.Length());

   // Define an angle and axis that rotates zhat into the direction
   // of the incident gamma, so that the generated kinematics is
   // defined with the incident gamma aligned with zhat, and then
   // rotated at the end into the final spatial direction.
   TThreeVectorReal rockaxis(mom0);
   rockaxis.Cross(TThreeVectorReal(0,0,1));
   double rockangle = asin(rockaxis.Length() / mom0.Length());
   rockaxis /= rockaxis.Length();

   // The beam photon energy is already generated by the cobrems
   // generator. Supply a normalized beam energy as the leading
   // elemnt of the u-vector in the call to sample(), and let it
   // return the remaining 5 random numbers needed to fix the
   // kinematics of the pair-production reaction.
   double E0_GeV = GlueXPrimaryGeneratorAction::getBeamEndpointEnergy()/GeV;
   double u[6];
   u[0] = kin / E0_GeV;
#if USE_ADAPTIVE_SAMPLER
   double weight = fAdaptiveSampler->sample(u);
#else
   unif01(5, u+1);
   double weight=1;
#endif

   // Generate uniform in E+, phi12, phiR
   LDouble_t Epos = kin * u[1];
   weight *= kin;
   LDouble_t phi12 = 2*M_PI * u[2];
   weight *= 2*M_PI;
   LDouble_t phiR = 2*M_PI * u[3];
   weight *= 2*M_PI;

   // Generate Mpair as 1 / (M [M^2 + Mcut^2])
   LDouble_t Mthresh = 2 * mElectron;
   LDouble_t Mcut = 0.005;  // GeV
   LDouble_t um0 = 1 + sqr(Mcut / Mthresh);
   LDouble_t umax = 1;
   if (fBHpair_mass_min > Mthresh)
      umax = log(1 + sqr(Mcut / fBHpair_mass_min)) / log(um0);
   LDouble_t um = pow(um0, u[4] * umax);
   LDouble_t Mpair = Mcut / sqrt(um - 1 + 1e-99);
   weight *= Mpair * (sqr(Mcut) + sqr(Mpair)) * log(um0) / (2 * sqr(Mcut));
   weight *= umax;

   // Generate qR^2 with weight 1 / [qR^2 sqrt(qRcut^2 + qR^2)]
   LDouble_t qRmin = sqr(Mpair) /(2 * kin);
   LDouble_t qRcut = 1e-3; // GeV
   LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
   LDouble_t uq = pow(uq0, u[5]);
   LDouble_t qR = 2 * qRcut * uq / (1 - sqr(uq));
   LDouble_t qR2 = qR * qR;
   weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut)) * (-2 * log(uq0));

   // Include overall measure Jacobian factor
   weight *= Mpair / (2 * kin);

   try {
      if (Epos < mElectron) {
         throw std::runtime_error("positron energy less than its rest mass.");
      }

      // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
      LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mElectron);
      if (k12star2 < 0) {
         throw std::runtime_error("no kinematic solution because k12star2 < 0");
      }
      LDouble_t k12star = sqrt(k12star2);
      LDouble_t Erec = sqrt(qR2 + sqr(mProton));
      LDouble_t E12 = kin + mProton - Erec;
      if (E12 < Mpair) {
         throw std::runtime_error("no kinematic solution because E12 < Mpair");
      }
      LDouble_t q12mag = sqrt(sqr(E12) - sqr(Mpair));
      LDouble_t costhetastar = (Epos - E12 / 2) * Mpair / (k12star * q12mag);
      if (Epos > E12 - mElectron) {
         throw std::runtime_error("no kinematic solution because Epos > E12 - mElectron");
      }
      else if (fabs(costhetastar) > 1) {
         throw std::runtime_error("no kinematic solution because |costhetastar| > 1");
      }

      // Solve for the recoil nucleon kinematics
      LDouble_t costhetaR = (sqr(Mpair) / 2 + (kin + mProton) *
                            (Erec - mProton)) / (kin * qR);
      if (fabs(costhetaR) > 1) {
         throw std::runtime_error("no kinematic solution because |costhetaR| > 1");
      }
      LDouble_t sinthetaR = sqrt(1 - sqr(costhetaR));
      TFourVectorReal q3(Erec, qR * sinthetaR * cos(phiR),
                               qR * sinthetaR * sin(phiR),
                               qR * costhetaR);

      // Boost the pair momenta into the lab
      TLorentzBoost toLab(q3[1] / E12, q3[2] / E12, (q3[3] - kin) / E12);
      LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
      TThreeVectorReal k12(k12star * sinthetastar * cos(phi12),
                           k12star * sinthetastar * sin(phi12),
                           k12star * costhetastar);
      TFourVectorReal q1(Mpair / 2, k12);
      TFourVectorReal q2(Mpair / 2, -k12);
      q1.Boost(toLab);
      q2.Boost(toLab);

      // Rotate final state momenta to the lab frame
      pOut.SetMom(q1.Rotate(rockaxis, rockangle));
      eOut.SetMom(q2.Rotate(rockaxis, rockangle));
      nOut.SetMom(q3.Rotate(rockaxis, rockangle));

      // Check 4-momentum conservation
      TFourVectorReal pIn(gIn.Mom() + nIn.Mom());
      TFourVectorReal pFi(pOut.Mom() + eOut.Mom() + nOut.Mom());
      TFourVectorReal::SetResolution(1e-10);
      if (pIn != pFi) {
         std::cout << "Warning in GenerateBetheHeitlerConversion - "
                   << "momentum conservation violated." << std::endl
                   << "   pIn = "; 
         pIn.Print();
         std::cout << "       = ";
         gIn.Mom().Print();
         std::cout << "       + ";
         nIn.Mom().Print();
         std::cout << "   pFi = ";
         pFi.Print();
         std::cout << "       = ";
         pOut.Mom().Print();
         std::cout << "       + ";
         eOut.Mom().Print();
         std::cout << "       + ";
         nOut.Mom().Print();
         std::cout << "   pIn - pFi = ";
         (pIn-pFi).Print();
      }

      // Compute the polarized differential cross section (barnes/GeV^4)
      // returned as d(sigma)/(dE+ dphi+ d^3qR)
      LDouble_t t = sqr(Erec - mProton) - qR2;
      LDouble_t tau = -t / sqr(2 * mProton);
      LDouble_t proton_magnetic_moment = 2.793;
      LDouble_t F1_timelike = 1;
      LDouble_t F2_timelike = 0;
      LDouble_t F1_spacelike = (1 / sqr(1 - t/0.71)) / (1 + tau) *
                               (1 + proton_magnetic_moment * tau);
      LDouble_t F2_spacelike = (1 / sqr(1 - t/0.71)) / (1 + tau) * 
                               (proton_magnetic_moment - 1);
      diffXS = TCrossSection::BetheHeitlerNucleon(gIn, nIn,
                                                  pOut, eOut, nOut,
                                                  F1_spacelike,
                                                  F2_spacelike,
                                                  F1_timelike,
                                                  F2_timelike);
#if USE_ADAPTIVE_SAMPLER
      fAdaptiveSampler->feedback(u, weight * diffXS);
#endif
   }
   catch (const std::exception &e) {

      // These events have no cross section, but do not discard
      // them because they are needed to get the right MC integral.
 
      nOut.SetMom(TThreeVectorReal(0,0,1e-12));
      eOut.SetMom(TThreeVectorReal(0,0,1e-12));
      pOut.SetMom(TThreeVectorReal(0,0,1e-12));
#if USE_ADAPTIVE_SAMPLER
      fAdaptiveSampler->feedback(u, 0);
#endif
      diffXS = 0;
   }

   // Generate the new vertex
   double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
   double steplength = pParticleChange->GetTrueStepLength();
   G4ThreeVector direction(track->GetMomentumDirection());
   G4ThreeVector x0(track->GetPosition());
   double t0 = track->GetGlobalTime();
   double uvtx = G4UniformRand();
   x0 -= uvtx * steplength * direction;
   t0 -= uvtx * steplength / beamVelocity;
   G4PrimaryVertex vertex(x0, t0);
   G4ParticleDefinition *positron = G4Positron::Definition();
   G4ParticleDefinition *electron = G4Electron::Definition();
   G4ParticleDefinition *proton = G4Proton::Definition();
   G4ThreeVector psec1(pOut.Mom()[1]*GeV, pOut.Mom()[2]*GeV, pOut.Mom()[3]*GeV);
   G4ThreeVector psec2(eOut.Mom()[1]*GeV, eOut.Mom()[2]*GeV, eOut.Mom()[3]*GeV);
   G4ThreeVector psec3(nOut.Mom()[1]*GeV, nOut.Mom()[2]*GeV, nOut.Mom()[3]*GeV);
   G4DynamicParticle *sec1 = new G4DynamicParticle(positron, psec1);
   G4DynamicParticle *sec2 = new G4DynamicParticle(electron, psec2);
   G4DynamicParticle *sec3 = new G4DynamicParticle(proton, psec3);
   G4TrackVector secondaries;
   secondaries.push_back(new G4Track(sec1, t0, x0));
   secondaries.push_back(new G4Track(sec2, t0, x0));
   if (nOut.Mom().Length() > 0) {
      secondaries.push_back(new G4Track(sec3, t0, x0));
   }

   GlueXUserEventInformation *event_info;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   event_info = (GlueXUserEventInformation*)event->GetUserInformation();
   G4TrackVector::iterator iter;
   for (iter = secondaries.begin(); iter != secondaries.end(); ++iter) {
      GlueXUserTrackInformation *trackinfo = new GlueXUserTrackInformation();
      if (event_info) {
         trackinfo->SetGlueXTrackID(event_info->AssignNextGlueXTrackID());
      }
      (*iter)->SetUserInformation(trackinfo);
      if (fStopBeamAfterTarget == 0) {
         pParticleChange->AddSecondary(*iter);
      }
   }

   // append secondary vertex to MC record
   if (event_info) {
      int mech[2];
      char *cmech = (char*)mech;
      snprintf(cmech, 5, "%c%c%c%c", 'C', 'O', 'N', 'V');
      event_info->AddSecondaryVertex(secondaries, 1, mech[0]);
      hddm_s::ReactionList rea = event_info->getOutputRecord()->getReactions();
      if (rea.size() > 0) {
         rea(0).setType(221); // Bethe-Heitler process
         rea(0).setWeight(weight * diffXS);
      }
   }

#if USE_ADAPTIVE_SAMPLER
   static int passes = 0;
   ++passes; // thread safety not important here
   int tens = 1000;
   while (passes >= tens*10)
      tens *= 10;
   if (passes % tens == 0) {
      std::cout << "AdaptiveSampler reports efficiency " 
                << fAdaptiveSampler->getEfficiency()
                << std::endl;
      std::stringstream astatefile;
      astatefile << "BHgen_thread_" << G4Threading::G4GetThreadId() << ".astate";
      fAdaptiveSampler->saveState(astatefile.str());
   }
#endif

#endif
}

void GlueXBeamConversionProcess::GenerateTripletProcess(const G4Step &step)
{
   // This method is completemary to GenerateBetheHeitlerProcess, for the case
   // of pair production off a target electron rather than coherently from the
   // atom as a whole. It generates the triplet process in the GlueX liquid
   // hydrogen target with the full polarization-dependent QED differential
   // cross section. A minimum invariant mass of the pair may be specified in
   // the member variable fBHpair_mass_min. It should be called from your
   // G4SteppingAction::UserSteppingAction method to force pair conversion
   // in special simulations dedicated to the study of the triplet proces.
   // The incident gamma ray is stopped and the e+/e- pair vertex plus the
   // recoil electron is added to the output event record.

#ifndef USING_DIRACXX

   G4cerr << "GlueXBeamConversionProcess::GenerateTripletProcess error:"
          << G4endl
          << "  You have enabled triplet conversion in the LIH2 target,"
          << G4endl
          << "  but your have built HDGeant4 without support for the Dirac++"
		  << G4endl
          << "  library. Either rebuild with Dirac++ support or else turn off"
		  << G4endl
          << "  this process. Aborting this run..."
		  << G4endl;
   exit(1);

#else

   TPhoton gIn;
   TLepton eOut(mElectron);
   TLepton pOut(mElectron);
   TLepton eIn(mElectron);
   TLepton eOut3(mElectron);
   LDouble_t diffXS(0);
   eOut.AllPol();
   pOut.AllPol();
   eOut3.AllPol();
   const TThreeVectorReal zero3vector(0,0,0);
   eIn.SetPol(TThreeVectorReal(zero3vector));
   eIn.SetMom(TThreeVectorReal(zero3vector));
   const G4Track *track = step.GetTrack();
   LDouble_t kin = track->GetKineticEnergy()/GeV;

   // If we are below pair production threshold, do nothing
   if (kin < 4 * mElectron)
      return;

   G4ThreeVector mom(track->GetMomentum());
   TThreeVectorReal mom0(mom[0]/GeV, mom[1]/GeV, mom[2]/GeV);
   gIn.SetMom(mom0);
   G4ThreeVector pol(track->GetPolarization());
   TThreeVectorReal pol0(pol[0], pol[1], pol[2]);
   gIn.SetPlanePolarization(pol0, pol0.Length());

   // Define an angle and axis that rotates zhat into the direction
   // of the incident gamma, so that the generated kinematics is
   // defined with the incident gamma aligned with zhat, and then
   // rotated at the end into the final spatial direction.
   TThreeVectorReal rockaxis(mom0);
   rockaxis.Cross(TThreeVectorReal(0,0,1));
   double rockangle = asin(rockaxis.Length() / mom0.Length());
   rockaxis /= rockaxis.Length();

   // The beam photon energy is already generated by the cobrems
   // generator. Supply a normalized beam energy as the leading
   // elemnt of the u-vector in the call to sample(), and let it
   // return the remaining 5 random numbers needed to fix the
   // kinematics of the pair-production reaction.
   double E0_GeV = GlueXPrimaryGeneratorAction::getBeamEndpointEnergy()/GeV;
   double u[6];
   u[0] = kin / E0_GeV;
#if USE_ADAPTIVE_SAMPLER
   double weight = fAdaptiveSampler->sample(u);
#else
   unif01(5, u+1);
   double weight=1;
#endif

   // Generate uniform in E+, phi12, phiR
   LDouble_t Epos = kin * u[1];
   weight *= kin;
   LDouble_t phi12 = 2*M_PI * u[2];
   weight *= 2*M_PI;
   LDouble_t phiR = 2*M_PI * u[3];
   weight *= 2*M_PI;

   // Generate Mpair as 1 / (M [M^2 + Mcut^2])
   LDouble_t Mthresh = 2 * mElectron;
   LDouble_t Mcut = 0.005;  // GeV
   LDouble_t um0 = 1 + sqr(Mcut / Mthresh);
   LDouble_t umax = 1;
   if (fBHpair_mass_min > Mthresh)
      umax = log(1 + sqr(Mcut / fBHpair_mass_min)) / log(um0);
   LDouble_t um = pow(um0, u[4] * umax);
   LDouble_t Mpair = Mcut / sqrt(um - 1 + 1e-99);
   weight *= Mpair * (sqr(Mcut) + sqr(Mpair)) * log(um0) / (2 * sqr(Mcut));
   weight *= umax;

   // Generate qR^2 with weight 1 / [qR^2 sqrt(qRcut^2 + qR^2)]
   LDouble_t qRmin = sqr(Mpair) /(2 * kin);
   LDouble_t qRcut = 1e-3; // GeV
   LDouble_t uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
   LDouble_t uq = pow(uq0, u[5]);
   LDouble_t qR = 2 * qRcut * uq / (1 - sqr(uq));
   LDouble_t qR2 = qR * qR;
   weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut)) * (-2 * log(uq0));

   // Include overall measure Jacobian factor
   weight *= Mpair / (2 * kin);

   try {
      if (Epos < mElectron) {
         throw std::runtime_error("positron energy less than its rest mass.");
      }

      // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
      LDouble_t k12star2 = sqr(Mpair / 2) - sqr(mElectron);
      if (k12star2 < 0) {
         throw std::runtime_error("no kinematic solution because k12star2 < 0");
      }
      LDouble_t k12star = sqrt(k12star2);
      LDouble_t Erec = sqrt(qR2 + sqr(mElectron));
      LDouble_t E12 = kin + mElectron - Erec;
      if (E12 < Mpair) {
         throw std::runtime_error("no kinematic solution because E12 < Mpair");
      }
      LDouble_t q12mag = sqrt(sqr(E12) - sqr(Mpair));
      LDouble_t costhetastar = (Epos - E12 / 2) * Mpair / (k12star * q12mag);
      if (Epos > E12 - mElectron) {
         throw std::runtime_error("no kinematic solution because Epos > E12 - mElectron");
      }
      else if (fabs(costhetastar) > 1) {
         throw std::runtime_error("no kinematic solution because |costhetastar| > 1");
      }

      // Solve for the recoil nucleon kinematics
      LDouble_t costhetaR = (sqr(Mpair) / 2 + (kin + mElectron) *
                            (Erec - mElectron)) / (kin * qR);
      if (fabs(costhetaR) > 1) {
         throw std::runtime_error("no kinematic solution because |costhetaR| > 1");
      }
      LDouble_t sinthetaR = sqrt(1 - sqr(costhetaR));
      TFourVectorReal q3(Erec, qR * sinthetaR * cos(phiR),
                               qR * sinthetaR * sin(phiR),
                               qR * costhetaR);

      // Boost the pair momenta into the lab
      TLorentzBoost toLab(q3[1] / E12, q3[2] / E12, (q3[3] - kin) / E12);
      LDouble_t sinthetastar = sqrt(1 - sqr(costhetastar));
      TThreeVectorReal k12(k12star * sinthetastar * cos(phi12),
                           k12star * sinthetastar * sin(phi12),
                           k12star * costhetastar);
      TFourVectorReal q1(Mpair / 2, k12);
      TFourVectorReal q2(Mpair / 2, -k12);
      q1.Boost(toLab);
      q2.Boost(toLab);

      // Rotate final state momenta to the lab frame
      pOut.SetMom(q1.Rotate(rockaxis, rockangle));
      eOut.SetMom(q2.Rotate(rockaxis, rockangle));
      eOut3.SetMom(q3.Rotate(rockaxis, rockangle));

      // Check 4-momentum conservation
      TFourVectorReal pIn(gIn.Mom() + eIn.Mom());
      TFourVectorReal pFi(pOut.Mom() + eOut.Mom() + eOut3.Mom());
      TFourVectorReal::SetResolution(1e-10);
      if (pIn != pFi) {
         std::cout << "Warning in GenerateTripletConversion - "
                   << "momentum conservation violated." << std::endl
                   << "   pIn = "; 
         pIn.Print();
         std::cout << "       = ";
         gIn.Mom().Print();
         std::cout << "       + ";
         eIn.Mom().Print();
         std::cout << "   pFi = ";
         pFi.Print();
         std::cout << "       = ";
         pOut.Mom().Print();
         std::cout << "       + ";
         eOut.Mom().Print();
         std::cout << "       + ";
         eOut3.Mom().Print();
         std::cout << "   pIn - pFi = ";
         (pIn-pFi).Print();
      }

      // Compute the polarized differential cross section (barnes/GeV^4)
      // returned as d(sigma)/(dE+ dphi+ d^3qR)
      diffXS = TCrossSection::TripletProduction(gIn, eIn, pOut, eOut, eOut3);
#if USE_ADAPTIVE_SAMPLER
      fAdaptiveSampler->feedback(u, weight * diffXS);
#endif
   }
   catch (const std::exception &e) {

      // These events have no cross section, but do not discard
      // them because they are needed to get the right MC integral.
 
      eOut3.SetMom(TThreeVectorReal(0,0,1e-12));
      eOut.SetMom(TThreeVectorReal(0,0,1e-12));
      pOut.SetMom(TThreeVectorReal(0,0,1e-12));
#if USE_ADAPTIVE_SAMPLER
      fAdaptiveSampler->feedback(u, 0);
#endif
      diffXS = 0;
   }

   // Generate the new vertex
   double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
   double steplength = pParticleChange->GetTrueStepLength();
   G4ThreeVector direction(track->GetMomentumDirection());
   G4ThreeVector x0(track->GetPosition());
   double t0 = track->GetGlobalTime();
   double uvtx = G4UniformRand();
   x0 -= uvtx * steplength * direction;
   t0 -= uvtx * steplength / beamVelocity;
   G4PrimaryVertex vertex(x0, t0);
   G4ParticleDefinition *positron = G4Positron::Definition();
   G4ParticleDefinition *electron = G4Electron::Definition();
   G4ThreeVector psec1(pOut.Mom()[1]*GeV, pOut.Mom()[2]*GeV, pOut.Mom()[3]*GeV);
   G4ThreeVector psec2(eOut.Mom()[1]*GeV, eOut.Mom()[2]*GeV, eOut.Mom()[3]*GeV);
   G4ThreeVector psec3(eOut3.Mom()[1]*GeV, eOut3.Mom()[2]*GeV, eOut3.Mom()[3]*GeV);
   G4DynamicParticle *sec1 = new G4DynamicParticle(positron, psec1);
   G4DynamicParticle *sec2 = new G4DynamicParticle(electron, psec2);
   G4DynamicParticle *sec3 = new G4DynamicParticle(electron, psec3);
   G4TrackVector secondaries;
   secondaries.push_back(new G4Track(sec1, t0, x0));
   secondaries.push_back(new G4Track(sec2, t0, x0));
   if (eOut3.Mom().Length() > 0) {
      secondaries.push_back(new G4Track(sec3, t0, x0));
   }

   GlueXUserEventInformation *event_info;
   const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
   event_info = (GlueXUserEventInformation*)event->GetUserInformation();
   G4TrackVector::iterator iter;
   for (iter = secondaries.begin(); iter != secondaries.end(); ++iter) {
      GlueXUserTrackInformation *trackinfo = new GlueXUserTrackInformation();
      if (event_info) {
         trackinfo->SetGlueXTrackID(event_info->AssignNextGlueXTrackID());
      }
      (*iter)->SetUserInformation(trackinfo);
      if (fStopBeamAfterTarget == 0) {
         pParticleChange->AddSecondary(*iter);
      }
   }

   // append secondary vertex to MC record
   if (event_info) {
      int mech[2];
      char *cmech = (char*)mech;
      snprintf(cmech, 5, "%c%c%c%c", 'C', 'O', 'N', 'V');
      event_info->AddSecondaryVertex(secondaries, 1, mech[0]);
      hddm_s::ReactionList rea = event_info->getOutputRecord()->getReactions();
      if (rea.size() > 0) {
         rea(0).setType(231); // triplet process
         rea(0).setWeight(weight * diffXS);
      }
   }

#if USE_ADAPTIVE_SAMPLER
   static int passes = 0;
   ++passes; // thread safety not important here
   int tens = 1000;
   while (passes >= tens*10)
      tens *= 10;
   if (passes % tens == 0) {
      std::cout << "AdaptiveSampler reports efficiency " 
                << fAdaptiveSampler->getEfficiency()
                << std::endl;
      std::stringstream astatefile;
      astatefile << "BHgen_thread_" << G4Threading::G4GetThreadId() << ".astate";
      fAdaptiveSampler->saveState(astatefile.str());
   }
#endif

#endif
}

void GlueXBeamConversionProcess::prepareAdaptiveSampler()
{
   fAdaptiveSampler = new AdaptiveSampler(6, &unif01, 1);
   fAdaptiveSampler->restoreState("BHgen.astate");
   fAdaptiveSampler->setVerbosity(2);
   fAdaptiveSampler->reset_stats();

   G4AutoLock barrier(&fMutex);
   fAdaptiveSamplerRegistry.push_back(fAdaptiveSampler);
}
