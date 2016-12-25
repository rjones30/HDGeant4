//
// class implementation for GlueXPrimaryGeneratorAction
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserOptions.hh"
#include "HddmOutput.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "Randomize.hh"

#include <JANA/jerror.h>
#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>

typedef GlueXPrimaryGeneratorAction::source_type_t source_type_t;
typedef GlueXPrimaryGeneratorAction::single_particle_gun_t particle_gun_t;
typedef GlueXPrimaryGeneratorAction::ImportanceSampler ImportanceSampler;

int GlueXPrimaryGeneratorAction::instanceCount = 0;
source_type_t GlueXPrimaryGeneratorAction::fSourceType = SOURCE_TYPE_NONE;

std::ifstream *GlueXPrimaryGeneratorAction::fHDDMinfile = 0;
hddm_s::istream *GlueXPrimaryGeneratorAction::fHDDMistream = 0;
CobremsGenerator *GlueXPrimaryGeneratorAction::fCobremsGenerator = 0;
#ifdef USING_DIRACXX
PairConversionGenerator *GlueXPrimaryGeneratorAction::fPairsGenerator = 0;
#endif
G4ParticleTable *GlueXPrimaryGeneratorAction::fParticleTable = 0;
GlueXParticleGun *GlueXPrimaryGeneratorAction::fParticleGun = 0;
GlueXPseudoDetectorTAG *GlueXPrimaryGeneratorAction::fTagger = 0;
particle_gun_t GlueXPrimaryGeneratorAction::fGunParticle;

double GlueXPrimaryGeneratorAction::fBeamBackgroundRate = 0;
double GlueXPrimaryGeneratorAction::fBeamBackgroundGateStart = 0;
double GlueXPrimaryGeneratorAction::fBeamBackgroundGateStop = 0;
double GlueXPrimaryGeneratorAction::fBeamBucketPeriod = 4. * ns;
double GlueXPrimaryGeneratorAction::fL1triggerTimeSigma = 10 * ns;
double GlueXPrimaryGeneratorAction::fBeamStartZ = -24 * m;
double GlueXPrimaryGeneratorAction::fTargetCenterZ = 65 * cm;
double GlueXPrimaryGeneratorAction::fTargetLength = 29.9746 * cm;
double GlueXPrimaryGeneratorAction::fBeamDiameter = 0.5 * cm;
double GlueXPrimaryGeneratorAction::fBeamVelocity = 2.99792e8 * m/s;

int GlueXPrimaryGeneratorAction::fEventCount = 0;

ImportanceSampler GlueXPrimaryGeneratorAction::fCoherentPDFx; 
ImportanceSampler GlueXPrimaryGeneratorAction::fIncoherentPDFlogx;
ImportanceSampler GlueXPrimaryGeneratorAction::fIncoherentPDFy;
double GlueXPrimaryGeneratorAction::fIncoherentPDFtheta02;
ImportanceSampler GlueXPrimaryGeneratorAction::fPaircohPDF;
ImportanceSampler GlueXPrimaryGeneratorAction::fTripletPDF;

G4Mutex GlueXPrimaryGeneratorAction::fMutex = G4MUTEX_INITIALIZER;
std::list<GlueXPrimaryGeneratorAction*> GlueXPrimaryGeneratorAction::fInstance;

//--------------------------------------------
// GlueXPrimaryGeneratorAction (constructor)
//--------------------------------------------

GlueXPrimaryGeneratorAction::GlueXPrimaryGeneratorAction()
{
   G4AutoLock barrier(&fMutex);
   fInstance.push_back(this);
   ++instanceCount;

   // Initializaton is driven by the control.in file, which
   // gets read and parsed only once, by the first constructor.

   if (fSourceType != SOURCE_TYPE_NONE) {
      return;
   }

   fParticleGun = new GlueXParticleGun();
   fParticleTable = G4ParticleTable::GetParticleTable();
   
   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXPrimaryGeneratorAction constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   std::map<int,std::string> infile;
   std::map<int,double> beampars;
   std::map<int,double> kinepars;

   // Three event source options are supported:
   // 1) external generator, hddm input stream source
   // 2) internal coherent bremsstrahlung beam generator
   // 3) internal particle gun generator
 
   if (user_opts->Find("INFILE", infile) ||
       user_opts->Find("INFI", infile))
   {
      fHDDMinfile = new std::ifstream(infile[1].c_str());
      if (!fHDDMinfile->is_open()) {
         G4cerr << "GlueXPrimaryGeneratorAction error: "
                << "Unable to open HDDM input file: " << infile[1]
                << G4endl;
         exit(-1);
      }
      fHDDMistream = new hddm_s::istream(*fHDDMinfile);
      G4cout << "Opened input file: " << infile[1] << G4endl;
      fSourceType = SOURCE_TYPE_HDDM;
   }

   else if (user_opts->Find("BEAM", beampars))
   {
      double beamE0 = beampars[1];
      double beamEpeak = beampars[2];
      double beamEmin = (beampars[3] > 0)? beampars[3] : 0.120;
      double radColDist = (beampars[4] > 0)? beampars[4] : 76.;
      double colDiam = (beampars[5] > 0)? beampars[5] : 0.0034;
      double beamEmit = (beampars[6] > 0)? beampars[6] : 2.5e-9;
      double radThick = (beampars[7] > 0)? beampars[7] : 20e-6;

      if (beamE0 == 0 || beamEpeak == 0) {
         G4cerr << "GlueXPrimaryGeneratorAction error: "
                << "BEAM card specified in control.in but required values "
                << "Emax and/or Epeak are missing, cannot continue."
                << G4endl;
         exit(-1);
      }

      // CobremsGenerator has its own standard units that it uses:
      //  length : m
      //  angles : radians
      //  energy : GeV
      //  time   : s
      //  current: microAmps
 
      fCobremsGenerator = new CobremsGenerator(beamE0, beamEpeak);
      fCobremsGenerator->setPhotonEnergyMin(beamEmin);
      fCobremsGenerator->setCollimatorDistance(radColDist);
      fCobremsGenerator->setCollimatorDiameter(colDiam);
      fCobremsGenerator->setBeamEmittance(beamEmit);
      fCobremsGenerator->setTargetThickness(radThick);
      fPairsGenerator = new PairConversionGenerator();

      // These cutoffs should be set empirically, as low as possible
      // for good efficiency, but not too low so as to avoid excessive
      // warnings about Pcut violations.
      fCoherentPDFx.Pcut = .003;
      fIncoherentPDFlogx.Pcut = .003;
      fPaircohPDF.Pcut = 10;
      fTripletPDF.Pcut = 2.5;

      std::map<int, double> bgratepars;
      std::map<int, double> bggatepars;
      if (user_opts->Find("BGRATE", bgratepars) &&
          user_opts->Find("BGGATE", bggatepars))
      {
         fBeamBackgroundRate = bgratepars[1] * 1/s;
         fBeamBackgroundGateStart = bgratepars[1] * ns;
         fBeamBackgroundGateStop = bgratepars[2] * ns;
         if (fBeamBackgroundRate > 0 &&
             fBeamBackgroundGateStart >= fBeamBackgroundGateStop)
         {
            G4cerr << "GlueXPrimaryGeneratorAction error: "
                   << "BGRATE is non-zero, but the time window specified "
                   << "in BGGATE is invalid."
                   << G4endl;
            exit(-1);
         }
      }
      fSourceType = SOURCE_TYPE_COBREMS_GEN;
   }

   else if (user_opts->Find("KINE", kinepars))
   {
      if (kinepars[1] == 1000) {
         fGunParticle.geantType = 0;
         fGunParticle.pdgType = 999999;
         fGunParticle.partDef = fParticleTable->FindParticle("geantino");
      }
      else if (kinepars[1] == 1001) {
         fGunParticle.geantType = 0;
         fGunParticle.pdgType = 999999;
         fGunParticle.partDef = fParticleTable->FindParticle("chargedgeantino");
      }
      else {
         if (kinepars[1] > 100)
            fGunParticle.geantType = kinepars[1] - 100;
         else
            fGunParticle.geantType = kinepars[1];
         fGunParticle.pdgType = ConvertGeant3ToPdg(fGunParticle.geantType);
         fGunParticle.partDef = fParticleTable->FindParticle(fGunParticle.pdgType);
      }
      if (fGunParticle.partDef == 0) {   
         G4cerr << "GlueXPrimaryGeneratorAction constructor error - "
                << "Unknown GEANT particle type: " << kinepars[1] 
                << " was specified in the control.in file." << G4endl;
         exit(-1);
      }
      fParticleGun->SetParticleDefinition(fGunParticle.partDef);

      double x(0), y(0), z(65 * cm);
      std::map<int,double> scappars;
      if (user_opts->Find("SCAP", scappars)) {
         x = scappars[1] * cm;
         y = scappars[2] * cm;
         z = scappars[3] * cm;
      }
      fGunParticle.pos.set(x,y,z);
      std::map<int,double> tgtwidthpars;
      if (user_opts->Find("tgtwidth", tgtwidthpars)) {
         fGunParticle.deltaR = tgtwidthpars[1] * cm;
         fGunParticle.deltaZ = tgtwidthpars[2] * cm;
      }
      else {
         fGunParticle.deltaR = 0;
         fGunParticle.deltaZ = 0;
      }

      fGunParticle.mom = kinepars[2] * GeV;
      if (kinepars[1] > 100) {
         fGunParticle.theta = kinepars[3] * degree;
         fGunParticle.phi = kinepars[4] * degree;
         fGunParticle.deltaMom = kinepars[5] * GeV;
         fGunParticle.deltaTheta = kinepars[6] * degree;
         fGunParticle.deltaPhi = kinepars[7] * degree;
      }
      else {
         fGunParticle.deltaMom = 0;
         fGunParticle.theta = 90 * degree;
         fGunParticle.deltaTheta = 180 * degree;
         fGunParticle.phi = 0;
         fGunParticle.deltaPhi = 360 * degree;
      }
      fSourceType = SOURCE_TYPE_PARTICLE_GUN;
   }

   std::map<int,double> trefsigma;
   if (user_opts->Find("trefsigma", trefsigma)) {
      fL1triggerTimeSigma = trefsigma[1] * ns;
   }
   else {
      fL1triggerTimeSigma = 10 * ns;
   }
}

GlueXPrimaryGeneratorAction::GlueXPrimaryGeneratorAction(const
                             GlueXPrimaryGeneratorAction &src)
 : G4VUserPrimaryGeneratorAction(src)
{
   G4AutoLock barrier(&fMutex);
   fInstance.push_back(this);
   ++instanceCount;
}

GlueXPrimaryGeneratorAction &GlueXPrimaryGeneratorAction::operator=(const
                             GlueXPrimaryGeneratorAction &src)
{
   *(G4VUserPrimaryGeneratorAction*)this = src;
   return *this;
}

//--------------------------------------------
// ~GlueXPrimaryGeneratorAction (destructor)
//--------------------------------------------

GlueXPrimaryGeneratorAction::~GlueXPrimaryGeneratorAction()
{
   G4AutoLock barrier(&fMutex);
   fInstance.remove(this);
   if (--instanceCount == 0) {
      if (fHDDMistream)
         delete fHDDMistream;
      if (fHDDMinfile)
         delete fHDDMinfile;
      if (fCobremsGenerator)
         delete fCobremsGenerator;
      if (fPairsGenerator)
         delete fPairsGenerator;
      delete fParticleGun;
   }
}

const GlueXPrimaryGeneratorAction *GlueXPrimaryGeneratorAction::GetInstance()
{
   // Generally one only needs a single instance of this object
   // per process, and this static method lets any component in the
   // application obtain the primary instance, if any. If none has
   // yet been constructed, it returns zero.

   if (fInstance.size() > 0)
      return *fInstance.begin();
   return 0;
}

void GlueXPrimaryGeneratorAction::prepareCobremsImportanceSamplingPDFs()
{
   // Construct lookup tables representing the PDFs used for
   // importance-sampling the coherent bremsstrahlung kinematics.

   const int Ndim = 500;
   double Emin = fCobremsGenerator->getPhotonEnergyMin() * GeV;
   double Emax = fCobremsGenerator->getBeamEnergy() * GeV;
   double sum;

   // Compute approximate PDF for dNc/dx
   double xmin = Emin / Emax;
   double dx = (1 - xmin) / Ndim;
   double xarr[Ndim], yarr[Ndim];
   for (int i=0; i < Ndim; ++i) {
      xarr[i] = xmin + (i + 0.5) * dx;
      yarr[i] = fCobremsGenerator->Rate_dNcdxdp(xarr[i], pi/4);
      yarr[i] = (yarr[i] > 0)? yarr[i] : 0;
   }
   fCobremsGenerator->applyBeamCrystalConvolution(Ndim, xarr, yarr);
   sum = 0;
   for (int i=0; i < Ndim; ++i) {
      sum += yarr[i];
      fCoherentPDFx.randvar.push_back(xarr[i]);
      fCoherentPDFx.density.push_back(yarr[i]);
      fCoherentPDFx.integral.push_back(sum);
   }
   for (int i=0; i < Ndim; ++i) {
      fCoherentPDFx.density[i] /= sum * dx;
      fCoherentPDFx.integral[i] /= sum;
   }

   // Compute approximate PDF for dNi/dlogx
   double logxmin = log(xmin);
   double dlogx = -logxmin / Ndim;
   double dNidlogx;
   sum = 0;
   for (int i=0; i < Ndim; ++i) {
      double logx = logxmin + (i + 0.5) * dlogx;
      double x = exp(logx);
      dNidlogx = fCobremsGenerator->Rate_dNidxdt2(x, 0) * x;
      dNidlogx = (dNidlogx > 0)? dNidlogx : 0;
      sum += dNidlogx;
      fIncoherentPDFlogx.randvar.push_back(logx);
      fIncoherentPDFlogx.density.push_back(dNidlogx);
      fIncoherentPDFlogx.integral.push_back(sum);
   }
   for (int i=0; i < Ndim; ++i) {
      fIncoherentPDFlogx.density[i] /= sum * dlogx;
      fIncoherentPDFlogx.integral[i] /= sum;
   }
 
   // Compute approximate PDF for dNi/dy
   fIncoherentPDFtheta02 = 1.8;
   double ymin = 1e-3;
   double dy = (1 - ymin) / Ndim;
   double dNidxdy;
   sum = 0;
   for (int i=0; i < Ndim; ++i) {
      double y = ymin + (i + 0.5) * dy;
      double theta2 = fIncoherentPDFtheta02 * (1 / y - 1);
      dNidxdy = fCobremsGenerator->Rate_dNidxdt2(0.5, theta2) *
                fIncoherentPDFtheta02 / (y*y);
      dNidxdy = (dNidxdy > 0)? dNidxdy : 0;
      sum += dNidxdy;
      fIncoherentPDFy.randvar.push_back(y);
      fIncoherentPDFy.density.push_back(dNidxdy);
      fIncoherentPDFy.integral.push_back(sum);
   }
   for (int i=0; i < Ndim; ++i) {
      fIncoherentPDFy.density[i] /= sum * dy;
      fIncoherentPDFy.integral[i] /= sum;
   }
   fCoherentPDFx.Pmax = 0;
   fCoherentPDFx.Psum = 0;
   fIncoherentPDFlogx.Pmax = 0;
   fIncoherentPDFlogx.Psum = 0;
   fIncoherentPDFy.Pmax = 0;
   fIncoherentPDFy.Psum = 0;
}

void GlueXPrimaryGeneratorAction::preparePairsImportanceSamplingPDFs()
{
   // Construct lookup tables representing the PDFs used for
   // importance-sampling the gamma pair conversion kinematics.

   // Compute 2D histogram containing rate as a function of u0,u1
   // where u0 is the (originally uniform [0,1]) random number used
   // to generate Mpair and u1 generates qR. The algorithm succeeds
   // because the mapping u0->Mpair and u1->qR used here is the
   // same as is used in GenerateBeamPairConversion.

   double kin = 9.; // GeV
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

   double Epos = kin / 2;
   double Mmin = 2 * mElectron;
   double Mcut = 5e-3; // 5 MeV cutoff parameter
   double um0 = 1 + sqr(Mcut / Mmin);
   double qRcut = 1e-3; // 1 MeV/c cutoff parameter

   int Nbins = 50;
   fTripletPDF.Psum = 0;
   fPaircohPDF.Psum = 0;
   for (int i0=0; i0 < Nbins; ++i0) {
      double u0 = (i0 + 0.5) / Nbins;
      double um = pow(um0, u0);
      double Mpair = Mcut / sqrt(um - 1);
      double k12star2 = sqr(Mpair / 2) - sqr(mElectron);
      double qRmin = sqr(Mpair) / (2 * kin);
      double uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
      double weight0 = sqr(Mpair) * (sqr(Mcut) + sqr(Mpair));
      for (int i1=0; i1 < Nbins; ++i1) {
         double u1 = (i1 + 0.5) / Nbins;
         double uq = pow(uq0, u1);
         double qR = 2 * qRcut * uq / (1 - sqr(uq));
         double weight = weight0 * sqr(qR) * sqrt(sqr(qRcut) + sqr(qR));
         double E3 = sqrt(sqr(qR) + sqr(mElectron));
         double E12 = kin + mElectron - E3;
         if (k12star2 < 0 || E12 < Mpair) {
            fPaircohPDF.density.push_back(0);
            fTripletPDF.density.push_back(0);
            fPaircohPDF.integral.push_back(fPaircohPDF.Psum);
            fTripletPDF.integral.push_back(fTripletPDF.Psum);
            continue;
         }
         double k12star = sqrt(k12star2);
         double q12mag = sqrt(sqr(E12) - sqr(Mpair));
         double costhetastar = (Epos - E12 / 2) * 
                               Mpair / (k12star * q12mag);
         double costhetaR = (sqr(Mpair) / 2 + (kin + mElectron) *
                             (E3 - mElectron)) / (kin * qR);
         if (fabs(costhetastar) > 1 || fabs(costhetaR) > 1) {
            fPaircohPDF.density.push_back(0);
            fTripletPDF.density.push_back(0);
            fPaircohPDF.integral.push_back(fPaircohPDF.Psum);
            fTripletPDF.integral.push_back(fTripletPDF.Psum);
            continue;
         }
         double qRlong = qR * costhetaR;
         double qRperp = sqrt(sqr(qR) - sqr(qRlong));
         TThreeVectorReal q3(0, qRperp, qRlong);
         double sinthetastar = sqrt(1 - sqr(costhetastar));
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
         double tripXS = fPairsGenerator->DiffXS_triplet(g0,p1,e2,e3);
         double pairXS = fPairsGenerator->DiffXS_pair(g0,p1,e2);
         fTripletPDF.Psum += fTripletPDF.Pmax = tripXS * weight;
         fPaircohPDF.Psum += fPaircohPDF.Pmax = pairXS * weight;
         fTripletPDF.density.push_back(fTripletPDF.Pmax);
         fPaircohPDF.density.push_back(fPaircohPDF.Pmax);
         fPaircohPDF.integral.push_back(fPaircohPDF.Psum);
         fTripletPDF.integral.push_back(fTripletPDF.Psum);
      }
   }

   double du2 = 1. / sqr(Nbins);
   for (int i0=0, index=0; i0 < Nbins; ++i0) {
      for (int i1=0; i1 < Nbins; ++i1, ++index) {
         double randvar = i0 + (i1 + 0.5) / Nbins;
         fTripletPDF.randvar.push_back(randvar);
         fPaircohPDF.randvar.push_back(randvar);
         fTripletPDF.density[index] /= fTripletPDF.Psum * du2;
         fPaircohPDF.density[index] /= fPaircohPDF.Psum * du2;
         fTripletPDF.integral[index] /= fTripletPDF.Psum;
         fPaircohPDF.integral[index] /= fPaircohPDF.Psum;
      }
   }
   fTripletPDF.Pmax = 0;
   fTripletPDF.Psum = 0;
   fPaircohPDF.Pmax = 0;
   fPaircohPDF.Psum = 0;
}

unsigned int ImportanceSampler::search(double u, const std::vector<double> &list)
{
   // Perform a fast search through non-decreasing list
   // for the index of the first element not less than u.

   int imin = -1;
   int imax = list.size() - 1;
   while (imax > imin + 1) {
      int i = (imin + imax) / 2;
      if (list[i] >= u)
         imax = i;
      else
         imin = i;
   }
   return imax;
}

unsigned int ImportanceSampler::search(double u) const
{
   return search(u, integral);
}

//--------------------------------------------
// GeneratePrimaries
//--------------------------------------------

void GlueXPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   G4AutoLock barrier(&fMutex);

   switch(fSourceType){
      case SOURCE_TYPE_HDDM:
         GeneratePrimariesHDDM(anEvent);
         break;
      case SOURCE_TYPE_COBREMS_GEN:
         GeneratePrimariesCobrems(anEvent);
         break;
      case SOURCE_TYPE_PARTICLE_GUN:
         GeneratePrimariesParticleGun(anEvent);
         break;
      default:
         G4cout << "No event source selected, cannot continue!" << G4endl;
         exit(-1);
   }
}   

//--------------------------------------------
// GeneratePrimariesParticleGun
//--------------------------------------------

void GlueXPrimaryGeneratorAction::GeneratePrimariesParticleGun(G4Event* anEvent)
{
   // Unbelievably, GEANT4's G4ParticleGun class insists on printing
   // a message whenever the momentum or energy is changed, unless
   // the other is 0. Here, we reset the particle gun energy using 
   // our own derived class. (Sheesh!!)
   fParticleGun->Reset();

   // place and smear the particle gun origin
   G4ThreeVector pos(fGunParticle.pos);
   if (fGunParticle.deltaR > 0) {
      double dx, dy;
      while (true) {
         double rx = G4UniformRand() - 0.5;
         double ry = G4UniformRand() - 0.5;
         if (rx*rx + ry*ry <= 0.25) {
            dx = rx * 2 * fGunParticle.deltaR;
            dy = ry * 2 * fGunParticle.deltaR;
            break;
         }
      }
      pos += G4ThreeVector(dx, dy, 0);
   }
   if (fGunParticle.deltaZ > 0) {
      double dz = (G4UniformRand() - 0.5) * fGunParticle.deltaZ;
      pos += G4ThreeVector(0, 0, dz);
   }
   fParticleGun->SetParticlePosition(pos);

   // Assign and optionally smear the particle momentum
   double p = fGunParticle.mom;
   double thetap = fGunParticle.theta;
   double phip = fGunParticle.phi;
   if (fGunParticle.deltaMom > 0)
      p += (G4UniformRand() - 0.5) * fGunParticle.deltaMom;
   if (fGunParticle.deltaTheta > 0)
      thetap += (G4UniformRand() - 0.5) * fGunParticle.deltaTheta;
   if (fGunParticle.deltaPhi > 0)
      phip += (G4UniformRand() - 0.5) * fGunParticle.deltaPhi;
   G4ThreeVector mom(p * sin(thetap) * cos(phip),
                     p * sin(thetap) * sin(phip),
                     p * cos(thetap));
   fParticleGun->SetParticleMomentum(mom);

   // Sync the particle gun generator to the beam bunch clock
   double tvtx = (pos[2] - fTargetCenterZ) / fBeamVelocity;
   tvtx -= GenerateTriggerTime();
   fParticleGun->SetParticleTime(tvtx);

   // Set the event number and fire the gun
   anEvent->SetEventID(++fEventCount);
   fParticleGun->GeneratePrimaryVertex(anEvent);

   // Store generated particle info so it can be written to output file
   GlueXUserEventInformation *event_info = new GlueXUserEventInformation();
   event_info->AddPrimaryVertex(*anEvent->GetPrimaryVertex());
   anEvent->SetUserInformation(event_info);
}

//--------------------------------------------
// GeneratePrimariesHDDM
//--------------------------------------------

void GlueXPrimaryGeneratorAction::GeneratePrimariesHDDM(G4Event* anEvent)
{
   if (! fHDDMinfile->good()) {
      anEvent->SetEventAborted();
      return;
   }

   hddm_s::HDDM *hddmevent = new hddm_s::HDDM;
   try {
      *fHDDMistream >> *hddmevent;
   }
   catch(std::exception e) {
      G4cout << e.what() << G4endl;
      anEvent->SetEventAborted();
      return;
   }

   // Store generated event info so it can be written to output file
   ++fEventCount;
   GlueXUserEventInformation *event_info;
   event_info = new GlueXUserEventInformation(hddmevent);
   anEvent->SetUserInformation(event_info);

   // Unpack generated event and prepare initial state for simulation
   int Nprimaries = 0;
   hddm_s::VertexList vertices = hddmevent->getVertices();
   if (vertices.size() == 0) {
      G4cout << "No vertices in input HDDM event!" << G4endl;
      anEvent->SetEventAborted();
      return;
   }
   hddm_s::VertexList::iterator it_vertex;
   for (it_vertex = vertices.begin();
        it_vertex != vertices.end(); ++it_vertex)
   {
      anEvent->SetEventID(it_vertex->getEventNo());
      hddm_s::Origin &origin = it_vertex->getOrigin();
      double x = origin.getVx() * cm;
      double y = origin.getVy() * cm;
      double z = origin.getVz() * cm;
      double t = origin.getT() * ns;
      if (x == 0 && y == 0 && z == 0) {
         while (true) {
            x = G4UniformRand() - 0.5;
            y = G4UniformRand() - 0.5;
            if (x*x + y*y <= 0.25) {
               x *= fBeamDiameter;
               y *= fBeamDiameter;
            }
         }
         z = fTargetCenterZ + (G4UniformRand() - 0.5) * fTargetLength;
         origin.setVx(x/cm);
         origin.setVy(y/cm);
         origin.setVz(z/cm);
      }
      if (t == 0) {
         t = (z - fTargetCenterZ) / fBeamVelocity;
         t -= GenerateTriggerTime();
         origin.setT(t/ns);
      }
      G4ThreeVector pos(x, y, z);
      G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, t);
      hddm_s::ProductList &products = it_vertex->getProducts();
      hddm_s::ProductList::iterator it_product;
      for (it_product = products.begin();
           it_product != products.end(); ++it_product)
      {
         // ignore intermediaries in the MC record
         if (it_product->getType() <= 0)
           continue;

         int g3type = it_product->getType();
         int pdgtype = it_product->getPdgtype();
         G4ParticleDefinition *part;
         if (pdgtype > 0 && pdgtype < 999999) {
            part = fParticleTable->FindParticle(pdgtype);
         }
         else if (g3type > 0) {
            pdgtype = ConvertGeant3ToPdg(g3type);
            part = fParticleTable->FindParticle(pdgtype);
#if FORCE_PARTICLE_TYPE_CHARGED_GEANTINO
            part = fParticleTable->FindParticle("chargedgeantino");
#endif
         }
         else {
            G4cerr << "Unknown particle found in input MC record, "
                   << "geant3 type " << g3type 
                   << ", PDG type " << pdgtype
                   << ", failing over to geantino!"
                   << G4endl;
            part = fParticleTable->FindParticle("geantino");
         }
         hddm_s::Momentum &momentum = it_product->getMomentum();
         double px = momentum.getPx() * GeV;
         double py = momentum.getPy() * GeV;
         double pz = momentum.getPz() * GeV;
         double Etot = momentum.getE() * GeV;
         vertex->SetPrimary(new G4PrimaryParticle(part, px, py, pz, Etot));
         ++Nprimaries;
      }
      anEvent->AddPrimaryVertex(vertex);
   }
   
   if (Nprimaries == 0) {
      G4cerr << "Number of primaries in event is zero!!" << G4endl;
      anEvent->SetEventAborted();
   }

   // Superimpose any request background minimum-bias beam interactions

   if (fBeamBackgroundRate > 0) {
      double t = fBeamBackgroundGateStart;
      while (true) {
         t += -log(G4UniformRand()) / fBeamBackgroundRate;
         if (t > fBeamBackgroundGateStop)
            break;
         GenerateBeamPhoton(anEvent, t);
      }
   }
}

void GlueXPrimaryGeneratorAction::GeneratePrimariesCobrems(G4Event* anEvent)
{
   GenerateBeamPhoton(anEvent, 0);
   ++fEventCount;
}

void GlueXPrimaryGeneratorAction::GenerateBeamPhoton(G4Event* anEvent,
                                                     double t0)
{
   // Generates a single beam photon according to the coherent bremsstrahlung
   // model defined by class CobremsGenerator.  The photon begins its lifetime
   // just upstream of the primary collimator (WARNING: position is hard-wired
   // in the code below) and is tracked by the simulation from there forward.
   // Its time t0 should identify its beam bucket, ie. the time the photon
   // would reach the midplane of the target. To enable beam motion spreading,
   // define the beam box size below.

   // The algorithm below generates coherent bremsstrahlung photons using a
   // importance-sampling technique. This algorithm requires that we prepare
   // an approximate probability density function for the generated photons.
   // The function is not in general equal to the true physical PDF, which
   // varies from event to event depending on the direction of the incident
   // electron, and also the exact angle of the crystal planes at the point
   // of scattering which moves because of the mosaic spread of the crystal.
   // The important thing is that the approximate PDF be reasonably close to
   // the average over all beam particles and the crystal mosaic, and that
   // deviations from event to event are sufficiently small that rejection
   // sampling can be used to take them into account with high efficiency.
   //
   // The kinematics of bremsstrahlung are described by three independent
   // variables (x, theta, phi) where x is the photon energy in units of the
   // incident electron energy, and theta,phi are the polar,azimuthal angles
   // of the photon in a lab frame tilted so that the incident electron comes
   // in along the z axis. Polar angle theta is represented by dimensionless
   // variable y = theta0^2 / (theta^2 + theta0^2) where contant theta0 is
   // chosen to optimize the uniformity of the PDF in y. On each event,
   // a new random tuple (x, phi, y) is generated on the interval x:[0,1],
   // phi:[0,2pi], y:[0,1] using a split-and-recombine strategy. One side 
   // of the split covers the coherent process and the other side covers the
   // incoherent one.
   //
   //  1) coherent process - the PDF here is continuous in x,phi according
   //     the dNc/(dx dphi), and the dependence on y is a sequence of delta 
   //     functions corresponding to the different planes that contribute to
   //     the scattering at the given value of x. Here we take advantage of
   //     the fact that the marginal distribution dNc/dx is proportional to
   //     dNc/(dx dphi) at phi=pi/4. This allows us to decompose the generation
   //     into two stages, first generating x from dNc/dx and then generating
   //     phi from dNc/(dx dphi) at fixed x. The x generation step is performed
   //     using importance sampling based on the average PDF stored in table
   //     fCoherentPDF, followed by rejection sampling based on the value of
   //     dNc/(dx dphi) computed for the particular kinematics of each event.
   //     The y value is obtained by sampling the weighted list of q2 values
   //     that contributed the to q-sum in the calculation of dNc/(dx dphi).
   //
   //  2) incoherent process - the PDF here is continuous in x,phi,y but it
   //     is uniform in phi, so it is effectively a 2D distribution. Here we
   //     take advantage of the fact that x and y are independent variables
   //     to a good approximation, which allows us to generate x using
   //     importance sampling from an approximation to dNi/(dx dtheta^2) at
   //     theta=0 and y ~ uniform [0,1], then employ rejection sampling based
   //     on the exact PDF dNi/(dx dtheta2) to get a true sample.
   //
   // Recombination after the split is very simple. First we integrate over
   // phi in both cases to obtain values dNc/dx and dNi/(dx dy). It turns
   // out that in spite of the fact that the y-dependence is discrete in the
   // coherent case and continuous in the incoherent case, the sum over the
   // probabilities for all values of y in dNc/dx is always normalized to 1
   // independently for all values of x. Hence we can treat y as a psuedo
   // coordinate y' ~ Unif[0,1] and form a 2D PDF dNc/(dx dy') which is
   // numerically equal to dNc/dx, do the rejection sampling in combination
   // with that applied to dNi/(dx dy) and then replace the fake variable y'
   // with the true y that was sampled as described above.

   if (fCoherentPDFx.density.size() == 0) {
      prepareCobremsImportanceSamplingPDFs();
   }

   double phiMosaic = twopi * G4UniformRand();
   double rhoMosaic = sqrt(-2 * log(G4UniformRand()));
   rhoMosaic *= fCobremsGenerator->getTargetCrystalMosaicSpread() * radian;
   double thxMosaic = rhoMosaic * cos(phiMosaic);
   double thyMosaic = rhoMosaic * sin(phiMosaic);

   double xemittance = fCobremsGenerator->getBeamEmittance() * m*radian;
   double yemittance = xemittance / 2.5; // nominal, should be checked
   double xspotsize = fCobremsGenerator->getCollimatorSpotrms() * m;
   double yspotsize = xspotsize; // nominal, should be checked
   double phiBeam = twopi * G4UniformRand();
   double rhoBeam = sqrt(-2 * log(G4UniformRand()));
   double thxBeam = (xemittance / xspotsize) * rhoBeam * cos(phiBeam);
   double thyBeam = (yemittance / yspotsize) * rhoBeam * sin(phiBeam);

   double raddz = fCobremsGenerator->getTargetThickness() * m;
   double varMS = fCobremsGenerator->Sigma2MS(raddz/m * G4UniformRand());
   double rhoMS = sqrt(-2 * varMS * log(G4UniformRand()));
   double phiMS = twopi * G4UniformRand();
   double thxMS = rhoMS * cos(phiMS);
   double thyMS = rhoMS * sin(phiMS);

   double targetThetax = fCobremsGenerator->getTargetThetax() * radian;
   double targetThetay = fCobremsGenerator->getTargetThetay() * radian;
   double targetThetaz = fCobremsGenerator->getTargetThetaz() * radian;
   double thetax = thxBeam + thxMS - targetThetax - thxMosaic;
   double thetay = thyBeam + thyMS - targetThetay - thyMosaic;
   double thetaz = -targetThetaz;
   fCobremsGenerator->setTargetOrientation(thetax, thetay, thetaz);

   // Generate with importance sampling
   double x = 0;
   double phi = 0;
   double theta2 = 0;
   double polarization = 0;
   double Scoherent = fCoherentPDFx.Npassed * 
                     (fCoherentPDFx.Ntested / (fCoherentPDFx.Psum + 1e-99));
   double Sincoherent = fIncoherentPDFlogx.Npassed *
                       (fIncoherentPDFlogx.Ntested /
                       (fIncoherentPDFlogx.Psum + 1e-99));
   if (Scoherent < Sincoherent) {
      while (true) {                             // try coherent generation
         ++fCoherentPDFx.Ntested;

         double u = G4UniformRand();
         int i = fCoherentPDFx.search(u);
         double fi = fCoherentPDFx.density[i];
         double ui = fCoherentPDFx.integral[i];
         double xi = fCoherentPDFx.randvar[i];
         double dx = (i > 0)? xi - fCoherentPDFx.randvar[i-1]:
                              fCoherentPDFx.randvar[i+1] - xi;
         x = xi + dx / 2 - (ui - u) / fi;
         double dNcdxPDF = fi;
         double dNcdx = twopi * fCobremsGenerator->Rate_dNcdxdp(x, pi / 4);
         double Pfactor = dNcdx / dNcdxPDF;
         if (Pfactor > fCoherentPDFx.Pmax)
            fCoherentPDFx.Pmax = Pfactor;
         if (Pfactor > fCoherentPDFx.Pcut) {
            G4cout << "Warning in GenerateBeamPhoton - Pfactor " << Pfactor
                   << " exceeds fCoherentPDFx.Pcut = " << fCoherentPDFx.Pcut
                   << G4endl
                   << "  present x = " << x << G4endl
                   << "  present maximum Pfactor = "
                   << fCoherentPDFx.Pmax << G4endl
                   << "  current generator efficiency = "
                   << fCoherentPDFx.Npassed /
                      (fCoherentPDFx.Ntested + 1e-99)
                   << G4endl;
         }
         fCoherentPDFx.Psum += Pfactor;
         if (G4UniformRand() * fCoherentPDFx.Pcut > Pfactor) {
            continue;
         }
         ++fCoherentPDFx.Npassed;

         double freq;
         double fmax = dNcdx / pi;
         while (true) {
            phi = twopi * G4UniformRand();
            freq = fCobremsGenerator->Rate_dNcdxdp(x, phi);
            if (G4UniformRand() * fmax < freq)
               break;
         }
         double uq = freq * G4UniformRand();
         int j = ImportanceSampler::search(uq, fCobremsGenerator->fQ2weight);
         theta2 = fCobremsGenerator->fQ2theta2[j];
         polarization = fCobremsGenerator->Polarization(x, theta2);
         break;
      }
   }
   else {
      while (true) {                           // try incoherent generation
         ++fIncoherentPDFlogx.Ntested;

         double ux = G4UniformRand();
         int i = fIncoherentPDFlogx.search(ux);
         double fi = fIncoherentPDFlogx.density[i];
         double ui = fIncoherentPDFlogx.integral[i];
         double logxi = fIncoherentPDFlogx.randvar[i];
         double dlogx = (i > 0)? logxi - fIncoherentPDFlogx.randvar[i-1]:
                                 fIncoherentPDFlogx.randvar[i+1] - logxi;
         double logx = logxi + dlogx / 2 - (ui - ux) / fi;
         x = exp(logx);
         double dNidxdyPDF = fi / x;
         double uy = G4UniformRand();
         int j = fIncoherentPDFy.search(uy);
         double fj = fIncoherentPDFy.density[j];
         double uj = fIncoherentPDFy.integral[j];
         double yj = fIncoherentPDFy.randvar[j];
         double dy = (j > 0)? yj - fIncoherentPDFy.randvar[j-1]:
                             fIncoherentPDFy.randvar[j+1] - yj;
         double y = yj + dy / 2 - (uj - uy) / fj;
         dNidxdyPDF *= fj;
         theta2 = fIncoherentPDFtheta02 * (1 / (y + 1e-99) - 1);
         double dNidxdy = fCobremsGenerator->Rate_dNidxdt2(x, theta2) *
                          fIncoherentPDFtheta02 / (y*y + 1e-99);
         double Pfactor = dNidxdy / dNidxdyPDF;
         if (Pfactor > fIncoherentPDFlogx.Pmax)
            fIncoherentPDFlogx.Pmax = Pfactor;
         if (Pfactor > fIncoherentPDFlogx.Pcut) {
            G4cout << "Warning in GenerateBeamPhoton - Pfactor " << Pfactor
                   << " exceeds fIncoherentPDFlogx.Pcut = " 
                   << fIncoherentPDFlogx.Pcut
                   << G4endl
                   << "  present x = " << x << G4endl
                   << "  present y = " << y << G4endl
                   << "  present maximum Pfactor = "
                   << fIncoherentPDFlogx.Pmax << G4endl
                   << "  current generator efficiency = "
                   << fIncoherentPDFlogx.Npassed /
                      (fIncoherentPDFlogx.Ntested + 1e-99)
                   << G4endl;
         }
         fIncoherentPDFlogx.Psum += Pfactor;
         if (G4UniformRand() * fIncoherentPDFlogx.Pcut > Pfactor) {
            continue;
         }
         ++fIncoherentPDFlogx.Npassed;

         phi = twopi * G4UniformRand();
         polarization = 0;
         break;
      }
   }

   // Put the radiator back the way your found it
   fCobremsGenerator->setTargetOrientation(targetThetax,
                                           targetThetay,
                                           targetThetaz);

   // Define the particle kinematics and polarization in lab coordinates
   G4ParticleDefinition *part = fParticleTable->FindParticle("gamma");
   double Emax = fCobremsGenerator->getBeamEnergy() * GeV;
   double Erms = fCobremsGenerator->getBeamErms() * GeV;
   double Ebeam = Emax + Erms * G4RandGauss::shoot();
   double theta = sqrt(theta2) * electron_mass_c2 / Emax;
   double alphax = thxBeam + thxMS + theta * cos(phi);
   double alphay = thyBeam + thyMS + theta * sin(phi);
   double pabs = Ebeam * x;
   double px = pabs * alphax;
   double py = pabs * alphay;
   double pz = sqrt(pabs*pabs - px*px - py*py);
   double colphi = twopi * G4UniformRand();
   double vspotrms = fCobremsGenerator->getCollimatorSpotrms() * m;
   double colrho = vspotrms * sqrt(-2 * log(G4UniformRand()));
   double colDist = fCobremsGenerator->getCollimatorDistance() * m;
   double radx = colrho * cos(colphi) - colDist * thxBeam;
   double rady = colrho * sin(colphi) - colDist * thyBeam;
   double colx = radx + colDist * alphax;
   double coly = rady + colDist * alphay;
#if defined BEAM_BOX_SIZE
   colx += BEAM_BOX_SIZE * (G4UniformRand() - 0.5);
   coly += BEAM_BOX_SIZE * (G4UniformRand() - 0.5);
#endif
   G4ThreeVector vtx(colx, coly, fBeamStartZ);
   G4ThreeVector pol(0, polarization, -polarization * py / pz);
   G4ThreeVector mom(px, py, pz);

   // If beam photon is primary particle, use it to initialize event info
   GlueXUserEventInformation *event_info;
   int bg = 1;
   double tvtx;
   if (t0 == 0) {
      event_info = new GlueXUserEventInformation();
      anEvent->SetUserInformation(event_info);
      tvtx = (vtx[2] - fTargetCenterZ) / fBeamVelocity;
      tvtx -= GenerateTriggerTime();
      event_info->AddBeamParticle(1, tvtx, vtx, mom, pol);
      bg = 0;
   }
   else {
      event_info = (GlueXUserEventInformation*)anEvent->GetUserInformation();
      assert (event_info != 0);
      tvtx = fBeamBucketPeriod * floor(t0 / fBeamBucketPeriod + 0.5);
      tvtx += (vtx[2] - fTargetCenterZ) / fBeamVelocity;
   }

   // Generate new primary for the beam photon
   G4PrimaryVertex* vertex = new G4PrimaryVertex(vtx, tvtx);
   G4PrimaryParticle* photon = new G4PrimaryParticle(part, px, py, pz);
   photon->SetPolarization(pol);
   vertex->SetPrimary(photon);
   anEvent->AddPrimaryVertex(vertex);

   // If bg beam particle, append to MC record
   if (bg) {
      event_info->AddPrimaryVertex(*vertex);
   }

   // Register a tagger hit for each beam photon
   int runNo = HddmOutput::getRunNo();
   if (fTagger == 0) {
      fTagger = new GlueXPseudoDetectorTAG(runNo);
   }
   else if (fTagger->getRunNo() != runNo) {
      delete fTagger;
      fTagger = new GlueXPseudoDetectorTAG(runNo);
   }
   double ttag = tvtx + (fTargetCenterZ - vtx[2]) / fBeamVelocity;
   fTagger->addTaggerPhoton(anEvent, vtx, pabs, ttag, bg);

#if VERBOSE_COBREMS_SPLITTING
   if (fIncoherentPDFlogx.Npassed / 100 * 100 == fIncoherentPDFlogx.Npassed) {
      G4cout << "coherent rate is "
             << fCoherentPDFx.Psum / (fCoherentPDFx.Ntested + 1e-99)
             << ", efficiency is "
             << fCoherentPDFx.Npassed / (fCoherentPDFx.Ntested + 1e-99)
             << G4endl
             << "incoherent rate is "
             << fIncoherentPDFlogx.Psum / (fIncoherentPDFlogx.Ntested + 1e-99)
             << ", efficiency is "
             << fIncoherentPDFlogx.Npassed / (fIncoherentPDFlogx.Ntested + 1e-99)
             << G4endl
             << "counts are "
             << fCoherentPDFx.Npassed << " / " << fIncoherentPDFlogx.Npassed
             << " = "
             << fCoherentPDFx.Npassed / (fIncoherentPDFlogx.Npassed + 1e-99)
             << G4endl;
   }
#endif
}

void GlueXPrimaryGeneratorAction::GenerateBeamPairConversion(const G4Step* step)
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

   G4AutoLock barrier(&fMutex);

#ifndef USING_DIRACXX

   G4cerr << "GlueXPrimaryGeneratorAction::GeneratorBeamPairConversion error:"
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

#if defined DO_TRIPLET_IMPORTANCE_SAMPLE || defined DO_PAIRCOH_IMPORTANCE_SAMPLE
   if (fTripletPDF.density.size() == 0) {
      G4cout << "GlueXPrimaryGeneratorAction::GenerateBeamPairConversion:"
             << G4endl
             << "   Setting up cross section tables, please wait... "
             << std::flush;
      preparePairsImportanceSamplingPDFs();
      G4cout << "finished." << G4endl;
   }
#endif

   TPhoton gIn;
   TLepton p1(mElectron);
   TLepton e2(mElectron);
   TLepton e3(mElectron);
   p1.AllPol();
   e2.AllPol();
   e3.AllPol();
   const G4Track *track = step->GetTrack();
   double kin = track->GetKineticEnergy()/GeV;
   G4ThreeVector mom(track->GetMomentum());
   TThreeVectorReal mom0(mom[0]/GeV, mom[1]/GeV, mom[2]/GeV);
   gIn.SetMom(mom0);
   G4ThreeVector pol(track->GetPolarization());
   TThreeVectorReal pol0(pol[0], pol[1], pol[2]);
   gIn.SetPol(pol0);

   // Define an angle and axis that rotates zhat into the direction
   // of the incidentn gamma, so that the generated kinematics is
   // defined with the incident gamma aligned with zhat, and then
   // rotated at the end into the final spatial direction.
   TThreeVectorReal rockaxis(mom0);
   rockaxis.Cross(TThreeVectorReal(0,0,1));
   double rockangle = rockaxis.Length() / kin;
   rockaxis /= rockaxis.Length();

   while (true) {
      double weight = 1;

      // Generate uniform in E+, phi12, phiR
      double Epos = kin * G4UniformRand();
      while (Epos < mElectron) {
         Epos = kin * G4UniformRand();
      }
      weight *= kin - mElectron;
      double phi12 = 2*M_PI * G4UniformRand();
      weight *= 2*M_PI;
      double phiR = 2*M_PI * G4UniformRand();
      weight *= 2*M_PI;

      double u0 = G4UniformRand();
      double u1 = G4UniformRand();

#if DO_TRIPLET_IMPORTANCE_SAMPLE

      int i = fTripletPDF.search(u1);
      double fi = fTripletPDF.density[i];
      double ui = fTripletPDF.integral[i];
      double ri = fTripletPDF.randvar[i];
      double xi = ri - floor(ri);
      double dx = (xi > 0.5)? ri - fTripletPDF.randvar[i-1]:
                              fTripletPDF.randvar[i+1] - ri;
      u1 = xi + dx / 2 - (ui - u1) / (fi * dx);
      u0 = (u0 + floor(ri)) * dx;
      weight /= fi;

#elif DO_PAIRCOH_IMPORTANCE_SAMPLE

      int i = fPaircohPDF.search(u1);
      double fi = fPaircohPDF.density[i];
      double ui = fPaircohPDF.integral[i];
      double ri = fPaircohPDF.randvar[i];
      double xi = ri - floor(ri);
      double dx = (xi > 0.5)? ri - fPaircohPDF.randvar[i-1]:
                              fPaircohPDF.randvar[i+1] - ri;
      u1 = xi + dx / 2 - (ui - u1) / (fi * dx);
      u0 = (u0 + floor(ri)) * dx;
      weight /= fi;

#endif
   
      // Generate Mpair as 1 / (M [M^2 + Mcut^2])
      double Mmin = 2 * mElectron;
      double Mcut = 0.005;  // GeV
      double um0 = 1 + sqr(Mcut / Mmin);
      double um = pow(um0, u0);
      double Mpair = Mcut / sqrt(um - 1 + 1e-99);
      weight *= Mpair * (sqr(Mcut) + sqr(Mpair)) * log(um0) / (2 * sqr(Mcut));
   
      // Generate qR^2 with weight 1 / [qR^2 sqrt(qRcut^2 + qR^2)]
      double qRmin = sqr(Mpair) /(2 * kin);
      double qRcut = 1e-3; // GeV
      double uq0 = qRmin / (qRcut + sqrt(sqr(qRcut) + sqr(qRmin)));
      double uq = pow(uq0, u1);
      double qR = 2 * qRcut * uq / (1 - sqr(uq));
      double qR2 = qR * qR;
      weight *= qR2 * sqrt(1 + qR2 / sqr(qRcut)) * (-2 * log(uq0));
   
      // Include overall measure Jacobian factor
      weight *= Mpair / (2 * kin);
   
      // Generate with importance sampling
      double Striplet = fTripletPDF.Npassed * 
                        (fTripletPDF.Ntested / (fTripletPDF.Psum + 1e-99));
      double Spaircoh = fPaircohPDF.Npassed *
                        (fPaircohPDF.Ntested / (fPaircohPDF.Psum + 1e-99));
      if (Striplet < Spaircoh) {                     // try incoherent generation
         ++fTripletPDF.Ntested;
   
         // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
         double k12star2 = sqr(Mpair / 2) - sqr(mElectron);
         if (k12star2 < 0) {
            // try again (this should never happen!)
            continue;
         }
         double k12star = sqrt(k12star2);
         double E3 = sqrt(qR2 + sqr(mElectron));
         double E12 = kin + mElectron - E3;
         if (E12 < Mpair) {
            // no kinematic solution because E12 < Mpair, try again
            continue;
         }
         double q12mag = sqrt(sqr(E12) - sqr(Mpair));
         double costhetastar = (Epos - E12 / 2) * Mpair / (k12star * q12mag);
         if (Epos > E12 - mElectron) {
            // no kinematic solution because Epos > E12 - mElectron, try again
            continue;
         }
         else if (fabs(costhetastar) > 1) {
            // no kinematic solution because |costhetastar| > 1, try again
            continue;
         }
   
         // Solve for the recoil electron kinematics
         double costhetaR = (sqr(Mpair) / 2 + (kin + mElectron) *
                             (E3 - mElectron)) / (kin * qR);
         if (fabs(costhetaR) > 1) {
            // no kinematic solution because |costhetaR| > 1, try again
            continue;
         }
         double sinthetaR = sqrt(1 - sqr(costhetaR));
         TFourVectorReal q3(E3, qR * sinthetaR * cos(phiR),
                                qR * sinthetaR * sin(phiR),
                                qR * costhetaR);
   
         // Boost the pair momenta into the lab
         double sinthetastar = sqrt(1 - sqr(costhetastar));
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

         // Compute the differential cross section (barnes/GeV^4)
         // returned as d(sigma)/(dE+ dphi+ d^3qR)
         p1.SetMom(q1.Rotate(rockaxis, rockangle));
         e2.SetMom(q2.Rotate(rockaxis, rockangle));
         e3.SetMom(q3.Rotate(rockaxis, rockangle));
         double diffXS = fPairsGenerator->DiffXS_triplet(gIn, p1, e2, e3);
   
         // Use keep/discard algorithm
         double Pfactor = diffXS * weight;
         if (Pfactor > fTripletPDF.Pmax)
            fTripletPDF.Pmax = Pfactor;
         if (Pfactor > fTripletPDF.Pcut) {
            G4cout << "Warning in GenerateBeamPairConversion - Pfactor " 
                   << Pfactor << " exceeds fTripletPDF.Pcut = " 
                   << fTripletPDF.Pcut << G4endl
                   << "  present qR = " << qR << G4endl
                   << "  present Mpair = " << Mpair << G4endl
                   << "  present maximum Pfactor = "
                   << fTripletPDF.Pmax << G4endl
                   << "  current generator efficiency = "
                   << fTripletPDF.Npassed /
                      (fTripletPDF.Ntested + 1e-99)
                   << G4endl;
         }
         fTripletPDF.Psum += Pfactor;
         if (G4UniformRand() * fTripletPDF.Pcut > Pfactor) {
            if (fTripletPDF.Ntested / 100 * 100 == fTripletPDF.Ntested) {
               G4cout << "rate is " << fTripletPDF.Psum / fTripletPDF.Ntested
                      << G4endl;
            }
            continue;
         }
         ++fTripletPDF.Npassed;
         break;
      }

      else {                          // try coherent generation
         ++fPaircohPDF.Ntested;
   
         // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
         // assuming that the atomic target absorbs zero energy
         double k12star2 = sqr(Mpair / 2) - sqr(mElectron);
         if (k12star2 < 0) {
            // try again (this should never happen!)
            continue;
         }
         double k12star = sqrt(k12star2);
         double Eele = kin - Epos;
         if (kin < Mpair) {
            // no kinematic solution because kin < Mpair, try again
            continue;
         }
         else if (Eele < mElectron) {
            // no kinematic solution because Eele < mElectron, try again
            continue;
         }
         double q12mag = sqrt(sqr(kin) - sqr(Mpair));
         double costhetastar = (Epos - kin / 2) * Mpair / (k12star * q12mag);
         if (fabs(costhetastar) > 1) {
            // no kinematic solution because |costhetastar| > 1, try again
            continue;
         }
   
         // Solve for the recoil kinematics kinematics
         double costhetaR = (sqr(Mpair) + qR2) / (2 * kin * qR);
         if (fabs(costhetaR) > 1) {
            // no kinematic solution because |costhetaR| > 1, try again
            continue;
         }
         double sinthetaR = sqrt(1 - sqr(costhetaR));
         TThreeVectorReal q3(qR * sinthetaR * cos(phiR),
                             qR * sinthetaR * sin(phiR),
                             qR * costhetaR);
   
         // Boost the pair momenta into the lab
         double sinthetastar = sqrt(1 - sqr(costhetastar));
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
         double diffXS = fPairsGenerator->DiffXS_pair(gIn, p1, e2);
   
         // Use keep/discard algorithm
         double Pfactor = diffXS * weight;
         if (Pfactor > fPaircohPDF.Pmax)
            fPaircohPDF.Pmax = Pfactor;
         if (Pfactor > fPaircohPDF.Pcut) {
            G4cout << "Warning in GenerateBeamPairConversion - Pfactor " 
                   << Pfactor << " exceeds fPaircohPDF.Pcut = " 
                   << fPaircohPDF.Pcut << G4endl
                   << "  present qR = " << qR << G4endl
                   << "  present Mpair = " << Mpair << G4endl
                   << "  present maximum Pfactor = "
                   << fPaircohPDF.Pmax << G4endl
                   << "  current generator efficiency = "
                   << fPaircohPDF.Npassed /
                      (fPaircohPDF.Ntested + 1e-99)
                   << G4endl;
         }
         fPaircohPDF.Psum += Pfactor;
         if (G4UniformRand() * fPaircohPDF.Pcut > Pfactor) {
            if (fPaircohPDF.Ntested / 100 * 100 == fPaircohPDF.Ntested) {
               G4cout << "rate is " << fPaircohPDF.Psum / fPaircohPDF.Ntested
                      << G4endl;
            }
            continue;
         }
         ++fPaircohPDF.Npassed;
         break;
      }
   }

   // Generate new primaries for the pair
   G4PrimaryVertex* vertex = new G4PrimaryVertex(track->GetPosition(),
                                                 track->GetGlobalTime());
   vertex->SetPrimary(new G4PrimaryParticle(G4Positron::Positron(),
                                            p1.Mom()[1]*GeV,
                                            p1.Mom()[2]*GeV,
                                            p1.Mom()[3]*GeV));
   vertex->SetPrimary(new G4PrimaryParticle(G4Electron::Electron(),
                                            e2.Mom()[1]*GeV,
                                            e2.Mom()[2]*GeV,
                                            e2.Mom()[3]*GeV));
   if (e3.Mom().Length() > 0) {
      vertex->SetPrimary(new G4PrimaryParticle(G4Electron::Electron(),
                                               e3.Mom()[1]*GeV,
                                               e3.Mom()[2]*GeV,
                                               e3.Mom()[3]*GeV));
   }
#if 0
   const G4Event *anEvent = G4RunManager::GetRunManager()->GetCurrentEvent();
   ((G4Event*)anEvent)->AddPrimaryVertex(vertex);

   // If bg beam particle, append to MC record
   GlueXUserEventInformation *event_info;
   event_info = (GlueXUserEventInformation*)anEvent->GetUserInformation();
   if (event_info) {
      event_info->AddPrimaryVertex(*vertex);
   }
#endif

#if VERBOSE_PAIRS_SPLITTING
   if (fTripletPDF.Npassed / 100 * 100 == fTripletPDF.Npassed) {
      G4cout << "triplet rate is "
             << fTripletPDF.Psum / (fTripletPDF.Ntested + 1e-99) 
             << ", efficiency is " 
             << fTripletPDF.Npassed / (fTripletPDF.Ntested + 1e-99)
             << G4endl
             << "pair rate is "
             << fPaircohPDF.Psum / (fPaircohPDF.Ntested + 1e-99) 
             << ", efficiency is " 
             << fPaircohPDF.Npassed / (fPaircohPDF.Ntested + 1e-99) 
             << G4endl
             << "counts are "
             << fTripletPDF.Npassed << " / " << fPaircohPDF.Npassed
             << " = "
             << fTripletPDF.Npassed / (fPaircohPDF.Npassed + 1e-99)
             << G4endl;
   }
#endif

#endif
}

// Convert particle types from Geant3 types to PDG scheme

int GlueXPrimaryGeneratorAction::ConvertGeant3ToPdg(int Geant3type)
{
   // This method was imported from ROOT source file TDatabasePDG.cc

   switch (Geant3type) {
      case 1   : return 22;       // photon
      case 2   : return -11;      // e+
      case 3   : return 11;       // e-
      case 4   : return 12;       // e-neutrino (NB: flavour undefined by Geant)
      case 5   : return -13;      // mu+
      case 6   : return 13;       // mu-
      case 7   : return 111;      // pi0
      case 8   : return 211;      // pi+
      case 9   : return -211;     // pi-
      case 10  : return 130;      // K long
      case 11  : return 321;      // K+
      case 12  : return -321;     // K-
      case 13  : return 2112;     // n
      case 14  : return 2212;     // p
      case 15  : return -2212;    // anti-proton
      case 16  : return 310;      // K short
      case 17  : return 221;      // eta
      case 18  : return 3122;     // Lambda
      case 19  : return 3222;     // Sigma+
      case 20  : return 3212;     // Sigma0
      case 21  : return 3112;     // Sigma-
      case 22  : return 3322;     // Xi0
      case 23  : return 3312;     // Xi-
      case 24  : return 3334;     // Omega- (PB)
      case 25  : return -2112;    // anti-neutron
      case 26  : return -3122;    // anti-Lambda
      case 27  : return -3222;    // Sigma-
      case 28  : return -3212;    // Sigma0
      case 29  : return -3112;    // Sigma+ (PB)*/
      case 30  : return -3322;    // Xi0
      case 31  : return -3312;    // Xi+
      case 32  : return -3334;    // Omega+ (PB)
      case 33  : return -15;      // tau+
      case 34  : return 15;       // tau-
      case 35  : return 411;      // D+
      case 36  : return -411;     // D-
      case 37  : return 421;      // D0
      case 38  : return -421;     // D0
      case 39  : return 431;      // Ds+
      case 40  : return -431;     // anti Ds-
      case 41  : return 4122;     // Lamba_c+
      case 42  : return 24;       // W+
      case 43  : return -24;      // W-
      case 44  : return 23;       // Z
      case 45  : return 1000010020; // deuteron
      case 46  : return 1000010030; // triton
      case 47  : return 1000020040; // alpha
      case 48  : return 0;        // geantino (no PDG type)
      case 49  : return 1000020030; // He3 ion
      case 50  : return 0;        // Cerenkov photon (no PDG type)

      default  :
         G4cout << "Warning in GlueXPrimaryGeneratorAction::"
                   "ConvertGeant3ToPdg - lookup performed on unknown"
                   " Geant3 particle type " << Geant3type << ","
                   " returning 0 for the PDG particle code." << G4endl;
   }
   return 0;
}

// Convert particle types from PDG scheme to Geant3 types

int GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(int PDGtype)
{
   // Invert the table contained in ConvertGeant3ToPdg

   switch (PDGtype) {
      case         22 : return 1;     // photon
      case        -11 : return 2;     // e+
      case         11 : return 3;     // e-
      case         12 : return 4;     // e-neutrino (NB: flavour undefined by Geant)
      case        -13 : return 5;     // mu+
      case         13 : return 6;     // mu-
      case        111 : return 7;     // pi0
      case        211 : return 8;     // pi+
      case       -211 : return 9;     // pi-
      case        130 : return 10;    // K long
      case        321 : return 11;    // K+
      case       -321 : return 12;    // K-
      case       2112 : return 13;    // n
      case       2212 : return 14;    // p
      case      -2212 : return 15;    // anti-proton
      case        310 : return 16;    // K short
      case        221 : return 17;    // eta
      case       3122 : return 18;    // Lambda
      case       3222 : return 19;    // Sigma+
      case       3212 : return 20;    // Sigma0
      case       3112 : return 21;    // Sigma-
      case       3322 : return 22;    // Xi0
      case       3312 : return 23;    // Xi-
      case       3334 : return 24;    // Omega- (PB)
      case      -2112 : return 25;    // anti-neutron
      case      -3122 : return 26;    // anti-Lambda
      case      -3222 : return 27;    // Sigma-
      case      -3212 : return 28;    // Sigma0
      case      -3112 : return 29;    // Sigma+ (PB)*/
      case      -3322 : return 30;    // Xi0
      case      -3312 : return 31;    // Xi+
      case      -3334 : return 32;    // Omega+ (PB)
      case        -15 : return 33;    // tau+
      case         15 : return 34;    // tau-
      case        411 : return 35;    // D+
      case       -411 : return 36;    // D-
      case        421 : return 37;    // D0
      case       -421 : return 38;    // D0
      case        431 : return 39;    // Ds+
      case       -431 : return 40;    // anti Ds-
      case       4122 : return 41;    // Lamba_c+
      case         24 : return 42;    // W+
      case        -24 : return 43;    // W-
      case         23 : return 44;    // Z
      case 1000010020 : return 45;    // deuteron
      case 1000010030 : return 46;    // triton
      case 1000020040 : return 47;    // alpha
      case 1000020030 : return 49;    // He3 ion

      default  :
         G4cout << "Warning in GlueXPrimaryGeneratorAction::"
                   "ConvertPdgToGeant3 - lookup performed on unknown"
                   " PDG particle type " << PDGtype << ","
                   " returning 0 for the Geant3 particle code." << G4endl;
   }
   return 0;
}

double GlueXPrimaryGeneratorAction::getBeamBucketPeriod(int runno)
{
   // Look up the beam bucket period for this run in ccdb
   // unless the user has already set the value by hand.

   if (runno > 0) {
      jana::JCalibration *jcalib = japp->GetJCalibration(runno);
      std::map<std::string, double> result;
      std::string map_key("/PHOTON_BEAM/RF/rf_period");
      if (jcalib->Get(map_key, result)) {
         G4cerr << "Error in GeneratePrimariesHDDM - "
                << "error fetching " << map_key << " from ccdb, "
                << "cannot continue." << G4endl;
         exit(-1);
      }
      else if (result.find("rf_period") != result.end()) {
         fBeamBucketPeriod = result["rf_period"] * ns;
      }
      else {
         G4cerr << "Error in GeneratePrimariesHDDM - "
                << "error finding value for " << map_key
                << " in ccdb, cannot continue." << G4endl;
         exit(-1);
      }
   }
   return fBeamBucketPeriod;
}

double GlueXPrimaryGeneratorAction::GetMassPDG(int PDGtype)
{
   return fParticleTable->FindParticle(PDGtype)->GetPDGMass();
}

double GlueXPrimaryGeneratorAction::GetMass(int Geant3Type)
{
   return GetMassPDG(ConvertGeant3ToPdg(Geant3Type));
}

double GlueXPrimaryGeneratorAction::GenerateTriggerTime()
{
   // The primary interaction vertex time is referenced to a clock
   // whose t=0 is synchronized to the crossing of a beam bunch
   // through the target midplane. This beam bunch may not contain
   // the beam particle whose interaction generated the vertex,
   // but it represents best-guess based on the arrival time of
   // the L1 trigger signal. The spread in the L1 relative to the
   // interacting bunch time is parameterized as a Gaussian.

   extern int run_number;
   static int last_run_number = 0;
   if (run_number != last_run_number) {
      getBeamBucketPeriod(run_number);
      last_run_number = run_number;
   }
   double t0 = fL1triggerTimeSigma * G4RandGauss::shoot();
   return fBeamBucketPeriod * floor(t0 / fBeamBucketPeriod + 0.5);
}
