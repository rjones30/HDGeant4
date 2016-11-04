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
G4ParticleTable *GlueXPrimaryGeneratorAction::fParticleTable = 0;
GlueXParticleGun *GlueXPrimaryGeneratorAction::fParticleGun = 0;
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

G4Mutex GlueXPrimaryGeneratorAction::fMutex = G4MUTEX_INITIALIZER;

//--------------------------------------------
// GlueXPrimaryGeneratorAction (constructor)
//--------------------------------------------

GlueXPrimaryGeneratorAction::GlueXPrimaryGeneratorAction()
{
   G4AutoLock barrier(&fMutex);
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
      prepareCobremsImportanceSamplingPDFs();

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
   if (--instanceCount == 0) {
      if (fHDDMistream)
         delete fHDDMistream;
      if (fHDDMinfile)
         delete fHDDMinfile;
      if (fCobremsGenerator)
         delete fCobremsGenerator;
      delete fParticleGun;
   }
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
   double xarr[Ndim + 1], yarr[Ndim + 1];
   for (int i=0; i <= Ndim; ++i) {
      xarr[i] = xmin + i * dx;
      yarr[i] = twopi * fCobremsGenerator->Rate_dNcdxdp(xarr[i], pi/2);
   }
   fCobremsGenerator->applyBeamCrystalConvolution(Ndim + 1, xarr, yarr);
   sum = 0;
   for (int i=0; i <= Ndim; ++i) {
      sum += (i > 0)? (yarr[i] + yarr[i - 1]) / 2 : 0;
      fCoherentPDFx.randvar.push_back(xarr[i]);
      fCoherentPDFx.density.push_back(yarr[i]);
      fCoherentPDFx.integral.push_back(sum);
   }
   for (int i=0; i <= Ndim; ++i) {
      fCoherentPDFx.density[i] /= sum * dx;
      fCoherentPDFx.integral[i] /= sum;
   }

   // Compute approximate PDF for dNi/dx
   double logxmin = log(xmin);
   double dlogx = -logxmin / Ndim;
   sum = 0;
   for (int i=0; i <= Ndim; ++i) {
      double logx = logxmin + i * dlogx;
      double x = exp(logx);
      double dNidx = fCobremsGenerator->Rate_dNidxdt2(x, 0);
      double dNidlogx = dNidx * x;
      fIncoherentPDFlogx.randvar.push_back(logx);
      fIncoherentPDFlogx.density.push_back(dNidlogx);
      fIncoherentPDFlogx.integral.push_back(sum);
      sum += (i < Ndim)? dNidlogx : 0;
   }
   for (int i=0; i <= Ndim; ++i) {
      fIncoherentPDFlogx.density[i] /= sum * dlogx;
      fIncoherentPDFlogx.integral[i] /= sum;
   }
 
   // Compute approximate PDF for dNi/dy
   fIncoherentPDFtheta02 = 1.8;
   double ymin = 1e-3;
   double dy = (1 - ymin) / Ndim;
   sum = 0;
   for (int i=0; i <= Ndim; ++i) {
      double y = ymin + i * dy;
      double theta2 = fIncoherentPDFtheta02 * (1 / y - 1);
      double dNidxdt2 = fCobremsGenerator->Rate_dNidxdt2(0.5, theta2);
      fIncoherentPDFy.randvar.push_back(y);
      fIncoherentPDFy.density.push_back(dNidxdt2);
      fIncoherentPDFy.integral.push_back(sum);
      sum += (i < Ndim)? dNidxdt2 : 0;
   }
   for (int i=0; i <= Ndim; ++i) {
      fIncoherentPDFy.density[i] /= sum * dy;
      fIncoherentPDFy.integral[i] /= sum;
   }

   // These cutoffs should be set empirically, as low as possible
   // for good efficiency, but not too low so as to avoid excessive
   // warnings about Pcut violations.
   fCoherentPDFx.Pcut = .001;
   fIncoherentPDFlogx.Pcut = 100.;
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

   // Store generated particle info so it can be written to output file
   int type = fGunParticle.geantType;
   GlueXUserEventInformation *event_info;
   event_info = new GlueXUserEventInformation(type, tvtx, pos, mom);
   anEvent->SetUserInformation(event_info);

   // Set the event number and fire the gun
   anEvent->SetEventID(++fEventCount);
   fParticleGun->GeneratePrimaryVertex(anEvent);
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


   double phiMosaic = twopi * G4UniformRand();
   double rhoMosaic = sqrt(-2 * log(G4UniformRand()));
   rhoMosaic *= fCobremsGenerator->getTargetCrystalMosaicSpread() * m*radian;
   double thxMosaic = rhoMosaic * cos(phiMosaic);
   double thyMosaic = rhoMosaic * sin(phiMosaic);

   double xemittance = fCobremsGenerator->getBeamEmittance() * m;
   double yemittance = xemittance / 2.5; // nominal, should be checked
   double xspotsize = fCobremsGenerator->getCollimatorSpotrms() * m;
   double yspotsize = xspotsize; // nominal, should be checked
   double thxBeam = (xemittance / xspotsize) * sqrt(-2 * log(G4UniformRand()));
   double thyBeam = (yemittance / yspotsize) * sqrt(-2 * log(G4UniformRand()));

   double raddz = fCobremsGenerator->getTargetThickness() * m;
   double varMS = fCobremsGenerator->Sigma2MS(raddz * G4UniformRand());
   double thxMS = sqrt(-2 * varMS * log(G4UniformRand()));
   double thyMS = sqrt(-2 * varMS * log(G4UniformRand()));

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
   double Scoherent = fCoherentPDFx.Npassed * fCoherentPDFx.Npassed / 
                      (fCoherentPDFx.Psum + 1e-99);
   double Sincoherent = fIncoherentPDFlogx.Npassed * fIncoherentPDFlogx.Npassed /
                        (fIncoherentPDFlogx.Psum + 1e-99);
   if (Scoherent < Sincoherent) {
      while (true) {                             // try coherent generation
         double dNcdxPDF;
         double u = G4UniformRand();
         for (unsigned int i=1; i < fCoherentPDFx.randvar.size(); ++i) {
            if (u <= fCoherentPDFx.integral[i]) {
               double x0 = fCoherentPDFx.randvar[i - 1];
               double x1 = fCoherentPDFx.randvar[i];
               double f0 = fCoherentPDFx.density[i - 1];
               double f1 = fCoherentPDFx.density[i];
               double u0 = fCoherentPDFx.integral[i - 1];
               double u1 = fCoherentPDFx.integral[i];
               x = (x0 * (u1 - u) + x1 * (u - u0)) / (u1 - u0);
               dNcdxPDF = (f0 * (u1 - u) + f1 * (u - u0)) / (u1 - u0);
               break;
            }
         }
         double dNcdx = twopi * fCobremsGenerator->Rate_dNcdxdp(x, pi / 4);
         double Pfactor = dNcdx / dNcdxPDF;
         if (Pfactor > fCoherentPDFx.Pmax)
            fCoherentPDFx.Pmax = Pfactor;
         if (Pfactor > fCoherentPDFx.Pcut) {
            G4cout << "Warning in GenerateBeamPhoton - Pfactor " << Pfactor
                   << " exceeds fCoherentPDFx.Pcut = " << fCoherentPDFx.Pcut
                   << ", please increase." << G4endl;
         }
         if (G4UniformRand() * fCoherentPDFx.Pcut > Pfactor) {
            ++fCoherentPDFx.Nfailed;
            continue;
         }
         fCoherentPDFx.Psum += Pfactor;
         ++fCoherentPDFx.Npassed;

         double fmax = dNcdx / pi;
         while (true) {
            phi = twopi * G4UniformRand();
            double f = fCobremsGenerator->Rate_dNcdxdp(x, phi);
            if (G4UniformRand() * fmax < f)
               break;
         }
         double uq = G4UniformRand();
         for (unsigned int i=0; i < fCobremsGenerator->fQ2theta2.size(); ++i) {
            if (uq <= fCobremsGenerator->fQ2weight[i]) {
               theta2 = fCobremsGenerator->fQ2theta2[i];
               break;
            }
         }
         polarization = fCobremsGenerator->Polarization(x, theta2);
         break;
      }
   }
   else {
      while (true) {                           // try incoherent generation
         double dNidxdyPDF;
         double u = G4UniformRand();
         for (unsigned int i=1; i < fIncoherentPDFlogx.randvar.size(); ++i) {
            if (u <= fIncoherentPDFlogx.integral[i]) {
               double logx0 = fIncoherentPDFlogx.randvar[i - 1];
               double logx1 = fIncoherentPDFlogx.randvar[i];
               double f0 = fIncoherentPDFlogx.density[i - 1];
               double f1 = fIncoherentPDFlogx.density[i];
               double u0 = fIncoherentPDFlogx.integral[i - 1];
               double u1 = fIncoherentPDFlogx.integral[i];
               double logx = (logx0 * (u1 - u) + logx1 * (u - u0)) / (u1 - u0);
               dNidxdyPDF = (f0 * (u1 - u) + f1 * (u - u0)) / (u1 - u0);
               x = exp(logx);
               break;
            }
         }
         double y;
         double uy = G4UniformRand();
         for (unsigned int i=1; i < fIncoherentPDFy.randvar.size(); ++i) {
            if (uy <= fIncoherentPDFy.integral[i]) {
               double y0 = fIncoherentPDFy.randvar[i - 1];
               double y1 = fIncoherentPDFy.randvar[i];
               double f0 = fIncoherentPDFy.density[i - 1];
               double f1 = fIncoherentPDFy.density[i];
               double u0 = fIncoherentPDFy.integral[i - 1];
               double u1 = fIncoherentPDFy.integral[i];
               y = (y0 * (u1 - uy) + y1 * (uy - u0)) / (u1 - u0);
               dNidxdyPDF *= (f0 * (u1 - uy) + f1 * (uy - u0)) / (u1 - u0);
               break;
            }
         }
         theta2 = fIncoherentPDFtheta02 * (1 / (y + 1e-99) - 1);
         double dNidxdy = fCobremsGenerator->Rate_dNidxdt2(x, theta2) *
                          fIncoherentPDFtheta02 / (y*y + 1e-99);
         double Pfactor = dNidxdy / dNidxdyPDF;
         if (Pfactor > fIncoherentPDFlogx.Pmax)
            fIncoherentPDFlogx.Pmax = Pfactor;
         if (Pfactor > fIncoherentPDFlogx.Pcut) {
            G4cout << "Warning in GenerateBeamPhoton - Pfactor " << Pfactor
                   << " exceeds fIncoherentPDFlogx.Pcut = " 
                   << fIncoherentPDFlogx.Pcut << ", please increase."
                   << G4endl;
         }
         if (G4UniformRand() * fIncoherentPDFlogx.Pcut > Pfactor) {
            ++fIncoherentPDFlogx.Nfailed;
            continue;
         }
         ++fIncoherentPDFlogx.Psum += Pfactor;
         ++fIncoherentPDFlogx.Npassed;

         phi = twopi * G4UniformRand();
         polarization = 0;
         break;
      }
   }

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

   // If beam photon is primary particle, use it to initialize event info
   int bg = 1;
   double tvtx;
   if (t0 == 0) {
      tvtx = (vtx[2] - fTargetCenterZ) / fBeamVelocity;
      tvtx -= GenerateTriggerTime();
      G4ThreeVector mom(px, py, pz);
      GlueXUserEventInformation *event_info;
      event_info = new GlueXUserEventInformation(1, tvtx, vtx, mom);
      anEvent->SetUserInformation(event_info);
      bg = 0;
   }
   else {
      tvtx = fBeamBucketPeriod * floor(t0 / fBeamBucketPeriod + 0.5);
      tvtx += (vtx[2] - fTargetCenterZ) / fBeamVelocity;
   }

   // Register a tagger hit for each beam photon
   double ttag = tvtx + (fTargetCenterZ - vtx[2]) / fBeamVelocity;
   HddmOutput::getTagger().addTaggerPhoton(anEvent, vtx, pabs, ttag, bg);
   //
   // Generate a new primary for the beam photon
   G4PrimaryVertex* vertex = new G4PrimaryVertex(vtx, tvtx);
   G4PrimaryParticle* photon = new G4PrimaryParticle(part, px, py, pz);
   photon->SetPolarization(pol);
   vertex->SetPrimary(photon);
   anEvent->AddPrimaryVertex(vertex);
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
