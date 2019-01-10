//
// class implementation for GlueXPrimaryGeneratorAction
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXPrimaryGenerator.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserOptions.hh"
#include "G4OpticalPhoton.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

typedef GlueXPrimaryGeneratorAction::source_type_t source_type_t;
typedef GlueXPrimaryGeneratorAction::single_particle_gun_t particle_gun_t;

int GlueXPrimaryGeneratorAction::instanceCount = 0;
source_type_t GlueXPrimaryGeneratorAction::fSourceType = SOURCE_TYPE_NONE;

std::ifstream *GlueXPrimaryGeneratorAction::fHDDMinfile = 0;
hddm_s::istream *GlueXPrimaryGeneratorAction::fHDDMistream = 0;

CobremsGeneration *GlueXPrimaryGeneratorAction::fCobremsGeneration = 0;
GlueXPhotonBeamGenerator *GlueXPrimaryGeneratorAction::fPhotonBeamGenerator = 0;
G4ParticleTable *GlueXPrimaryGeneratorAction::fParticleTable = 0;
particle_gun_t GlueXPrimaryGeneratorAction::fGunParticle;

GlueXParticleGun *GlueXPrimaryGeneratorAction::fParticleGun = 0;
GlueXPrimaryGenerator *GlueXPrimaryGeneratorAction::fPrimaryGenerator = 0;

double GlueXPrimaryGeneratorAction::fBeamBackgroundRate = 0;
double GlueXPrimaryGeneratorAction::fBeamBackgroundGateStart = 0;
double GlueXPrimaryGeneratorAction::fBeamBackgroundGateStop = 0;
double GlueXPrimaryGeneratorAction::fL1triggerTimeSigma = 10 * ns;
double GlueXPrimaryGeneratorAction::fRFreferencePlaneZ = 65 * cm;
double GlueXPrimaryGeneratorAction::fTargetCenterZ = 65 * cm;
double GlueXPrimaryGeneratorAction::fTargetLength = 29.9746 * cm;

int GlueXPrimaryGeneratorAction::fEventCount = 0;

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

   // get positions for LUT from XML geometry
   std::map<int, int> dirclutpars;
   if (instanceCount == 1) {
      
      if (user_opts->Find("DIRCLUT", dirclutpars)) {
         
         extern int run_number;
         extern jana::JApplication *japp;
         if (japp == 0) {
            G4cerr << "Error in GlueXPrimaryGeneratorAction constructor - "
              << "jana global DApplication object not set, "
              << "cannot continue." << G4endl;
            exit(-1);
         }
         jana::JGeometry *jgeom = japp->GetJGeometry(run_number);
         if (japp == 0) {   // dummy
            jgeom = 0;
            G4cout << "DIRC: ALL parameters loaded from ccdb" << G4endl;
         }
         
         vector<double>DIRC;
         vector<double>DRCC;
         vector<double>DCML00_XYZ;
         vector<double>DCML01_XYZ;
         vector<double>DCML10_XYZ;
         vector<double>DCML11_XYZ;
         vector<double>WNGL00_XYZ;
         vector<double>WNGL01_XYZ;
         vector<double>WNGL10_XYZ;
         vector<double>WNGL11_XYZ;
         vector<double>OWDG_XYZ;
         jgeom->Get("//section/composition/posXYZ[@volume='DIRC']/@X_Y_Z", DIRC);
         jgeom->Get("//composition[@name='DIRC']/posXYZ[@volume='DRCC']/@X_Y_Z", DRCC);
         jgeom->Get("//composition[@name='DRCC']/posXYZ[@volume='DCML00']/@X_Y_Z", DCML00_XYZ);
         jgeom->Get("//composition[@name='DRCC']/posXYZ[@volume='DCML01']/@X_Y_Z", DCML01_XYZ);
         jgeom->Get("//composition[@name='DRCC']/posXYZ[@volume='DCML10']/@X_Y_Z", DCML10_XYZ);
         jgeom->Get("//composition[@name='DRCC']/posXYZ[@volume='DCML11']/@X_Y_Z", DCML11_XYZ);
         jgeom->Get("//composition[@name='DCML00']/posXYZ[@volume='WNGL']/@X_Y_Z", WNGL00_XYZ);
         jgeom->Get("//composition[@name='DCML01']/posXYZ[@volume='WNGL']/@X_Y_Z", WNGL01_XYZ);
         jgeom->Get("//composition[@name='DCML10']/posXYZ[@volume='WNGL']/@X_Y_Z", WNGL10_XYZ);
         jgeom->Get("//composition[@name='DCML11']/posXYZ[@volume='WNGL']/@X_Y_Z", WNGL11_XYZ);
         jgeom->Get("//composition[@name='DCML11']/posXYZ[@volume='WNGL']/@X_Y_Z", WNGL11_XYZ);
         jgeom->Get("//trd[@name='OWDG']/@Xmp_Ymp_Z", OWDG_XYZ);
         DIRC_LUT_Z = (DIRC[2] + DRCC[2] + DCML01_XYZ[2] + 0.8625) * cm;
         DIRC_QZBL_DY = 3.5 * cm;   // nominal width to generate LUT
         DIRC_QZBL_DZ = 1.725 * cm; // nominal thickness to generate LUT
         DIRC_OWDG_DZ = OWDG_XYZ[4];

         // set array of bar positions
         for (int i=0; i<48; i++) {
            vector<double>DCBR_XYZ;
            if (i<12) {
               std::stringstream geomDCML10;
               geomDCML10 << "//composition[@name='DCML10']/posXYZ[@volume='DCBR" 
                     << std::setfill('0') << std::setw(2) << i << "']/@X_Y_Z"; 
               jgeom->Get(geomDCML10.str(), DCBR_XYZ);
               DIRC_BAR_Y[i] = (DCML10_XYZ[1] - DCBR_XYZ[1]) * cm;
               DIRC_LUT_X[i] = (DIRC[0] + DRCC[0] + DCML10_XYZ[0] - WNGL10_XYZ[0] + DIRC_OWDG_DZ) * cm;
            }
            else if (i<24) {
               std::stringstream geomDCML11;
               geomDCML11 << "//composition[@name='DCML11']/posXYZ[@volume='DCBR" 
                     << std::setfill('0') << std::setw(2) << i << "']/@X_Y_Z"; 
               jgeom->Get(geomDCML11.str(), DCBR_XYZ);
               DIRC_BAR_Y[i] = (DCML11_XYZ[1] - DCBR_XYZ[1]) * cm;
               DIRC_LUT_X[i] = (DIRC[0] + DRCC[0] + DCML11_XYZ[0] - WNGL11_XYZ[0] + DIRC_OWDG_DZ) * cm;
            }
            else if (i<36) {
               std::stringstream geomDCML01;
               geomDCML01 << "//composition[@name='DCML01']/posXYZ[@volume='DCBR" 
                     << std::setfill('0') << std::setw(2) << i << "']/@X_Y_Z"; 
               jgeom->Get(geomDCML01.str(), DCBR_XYZ);
               DIRC_BAR_Y[i] = (DCML01_XYZ[1] + DCBR_XYZ[1]) * cm;
               DIRC_LUT_X[i] = (DIRC[0] + DRCC[0] + DCML01_XYZ[0] + WNGL01_XYZ[0] - DIRC_OWDG_DZ) * cm;
            }
            else if (i<48) {
               std::stringstream geomDCML00;
               geomDCML00 << "//composition[@name='DCML00']/posXYZ[@volume='DCBR" 
                     << std::setfill('0') << std::setw(2) << i << "']/@X_Y_Z"; 
               jgeom->Get(geomDCML00.str(), DCBR_XYZ);
               DIRC_BAR_Y[i] = (DCML00_XYZ[1] + DCBR_XYZ[1]) * cm;
               DIRC_LUT_X[i] = (DIRC[0] + DRCC[0] + DCML00_XYZ[0] + WNGL00_XYZ[0] - DIRC_OWDG_DZ) * cm;
            }            
         }
      }
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
      fPrimaryGenerator = new GlueXPrimaryGenerator(fHDDMistream);
      fSourceType = SOURCE_TYPE_HDDM;
   }

   else if (user_opts->Find("BEAM", beampars))
   {
      fSourceType = SOURCE_TYPE_COBREMS_GEN;
   }

   else if (user_opts->Find("DIRCLUT", dirclutpars))
   {
      fGunParticle.geantType = 0;
      fGunParticle.pdgType = 999999;
      fGunParticle.partDef = fParticleTable->FindParticle("opticalphoton");
      fGunParticle.deltaR = 0;
      fGunParticle.deltaZ = 0;
      fGunParticle.mom = 3.18 * eV;

      fGunParticle.deltaMom = 0;
      fGunParticle.deltaTheta = 0;
      fGunParticle.deltaPhi = 0;
      fParticleGun->SetParticleDefinition(fGunParticle.partDef);
 
      fSourceType = SOURCE_TYPE_PARTICLE_GUN;
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
      else if (kinepars[1] == 1050) {
         fGunParticle.geantType = 0;
         fGunParticle.pdgType = 999999;
         fGunParticle.partDef = fParticleTable->FindParticle("opticalphoton");
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
      fGunParticle.plogOption = 0;
      std::map<int,int> plogpars;
      if (user_opts->Find("PLOG", plogpars)) {
         fGunParticle.plogOption = plogpars[1];
      }
      fGunParticle.tlogOption = 0;
      std::map<int,int> tlogpars;
      if (user_opts->Find("PLOG", tlogpars)) {
         fGunParticle.tlogOption = tlogpars[1];
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

   if (user_opts->Find("BEAM", beampars)) {
      double beamE0 = beampars[1] * GeV;
      double beamEpeak = beampars[2] * GeV;
      double beamEmin = ((beampars[3] > 0)? beampars[3] : 0.120) * GeV;
      double radColDist = ((beampars[4] > 0)? beampars[4] : 76.) * m;
      double colDiam = ((beampars[5] > 0)? beampars[5] : 0.0034) * m;
      double beamEmit = ((beampars[6] > 0)? beampars[6] : 2.5e-9) * m;
      double radThick = ((beampars[7] > 0)? beampars[7] : 20e-6) * m;

      if (beamE0 == 0) {
         G4cerr << "GlueXPrimaryGeneratorAction error: "
                << "BEAM card specified in control.in but required values "
                << "Emax and/or Epeak are missing, cannot continue."
                << G4endl;
         exit(-1);
      }

      // CobremsGeneration has its own standard units that it uses:
      //  length : m
      //  angles : radians
      //  energy : GeV
      //  time   : s
      //  current: microAmps
 
      fCobremsGeneration = new CobremsGeneration(beamE0/GeV, beamEpeak/GeV);
      fCobremsGeneration->setPhotonEnergyMin(beamEmin/GeV);
      fCobremsGeneration->setCollimatorDistance(radColDist/m);
      fCobremsGeneration->setCollimatorDiameter(colDiam/m);
      fCobremsGeneration->setBeamEmittance(beamEmit/(m*radian));
      fCobremsGeneration->setTargetThickness(radThick/m);
      fPhotonBeamGenerator = new GlueXPhotonBeamGenerator(fCobremsGeneration);

      std::map<int, double> bgratepars;
      std::map<int, double> bggatepars;
      if (user_opts->Find("BGRATE", bgratepars) &&
          user_opts->Find("BGGATE", bggatepars))
      {
         fBeamBackgroundRate = bgratepars[1] * 1/ns;
         fBeamBackgroundGateStart = bggatepars[1] * ns;
         fBeamBackgroundGateStop = bggatepars[2] * ns;
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
      if (fPrimaryGenerator)
         delete fPrimaryGenerator;
      if (fCobremsGeneration)
         delete fCobremsGeneration;
      if (fPhotonBeamGenerator)
         delete fPhotonBeamGenerator;
      delete fParticleGun;
      if (fHDDMistream)
         delete fHDDMistream;
      if (fHDDMinfile)
         delete fHDDMinfile;
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

//--------------------------------------------
// GeneratePrimaries
//--------------------------------------------

void GlueXPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   G4AutoLock barrier(&fMutex);

   switch (fSourceType) {
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

   //   std::cout<<"GlueXPrimaryGeneratorAction:: GENERATE PRIMARIES PARTICLE GUN"<<std::endl;

   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   std::map<int,int> dirclutpars; 

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
   if (fGunParticle.deltaMom > 0) {
      if (fGunParticle.plogOption) {
         double pmin = p - fGunParticle.deltaMom / 2;
         double pmax = p + fGunParticle.deltaMom / 2;
         pmin = (pmin > 0)? pmin : 1e-6*GeV;
         p = pmin * pow(pmax/pmin, G4UniformRand());
      }
      else {
         p += (G4UniformRand() - 0.5) * fGunParticle.deltaMom;
      }
   }
   if (fGunParticle.deltaTheta > 0) {
      if (fGunParticle.plogOption) {
         double thetamin = thetap - fGunParticle.deltaTheta / 2;
         double thetamax = thetap + fGunParticle.deltaTheta / 2;
         thetamin = (thetamin > 0)? thetamin : 1e-6*degree;
         thetap = thetamin * pow(thetamax/thetamin, G4UniformRand());
      }
      else {
         thetap += (G4UniformRand() - 0.5) * fGunParticle.deltaTheta;
      }
   }
   if (fGunParticle.deltaPhi > 0)
      phip += (G4UniformRand() - 0.5) * fGunParticle.deltaPhi;

   // Special case of Cherenkov photon gun for DIRC Look Up Tables (LUT)
   if (user_opts->Find("DIRCLUT", dirclutpars)) {

      // array of bar y-positions for LUT from JGeometry
      double y = 0.; // no shift
      double x = DIRC_LUT_X[dirclutpars[1]];
      double z = DIRC_LUT_Z;

      G4ThreeVector vec(0,0,1);
      double rand1 = G4UniformRand();
      double rand2 = G4UniformRand();
      vec.setTheta(acos(rand1));
      vec.setPhi(2*M_PI*rand2);
      vec.rotateY(M_PI/2.);
      y = DIRC_BAR_Y[dirclutpars[1]];
      if (dirclutpars[1] < 24) {
        vec.rotateY(M_PI);
      }
     
      // spread over end of bar in y and z
      y += DIRC_QZBL_DY/2.0 - DIRC_QZBL_DY*G4UniformRand();
      z += DIRC_QZBL_DZ/2.0 - DIRC_QZBL_DZ*G4UniformRand(); 

      thetap = vec.theta();
      phip = vec.phi();
      fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z));
   }

   G4ThreeVector mom(p * sin(thetap) * cos(phip),
                     p * sin(thetap) * sin(phip),
                     p * cos(thetap));
   fParticleGun->SetParticleMomentum(mom);
   fParticleGun->SetParticleTime(0);

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
   if (fPrimaryGenerator != 0) {
      double beamDiameter = GlueXPhotonBeamGenerator::getBeamDiameter();
      double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
      double x, y;
      while (true) {
         x = G4UniformRand() - 0.5;
         y = G4UniformRand() - 0.5;
         if (x*x + y*y <= 0.25) {
            x *= beamDiameter;
            y *= beamDiameter;
            break;
         }
      }
      double z = fTargetCenterZ + (G4UniformRand() - 0.5) * fTargetLength;
      fPrimaryGenerator->SetParticlePosition(G4ThreeVector(x,y,z));
      double ttag = GlueXPhotonBeamGenerator::GenerateTriggerTime(anEvent);
      double trel = (z - fRFreferencePlaneZ) / beamVelocity;
      fPrimaryGenerator->SetParticleTime(trel + ttag);
      fPrimaryGenerator->GeneratePrimaryVertex(anEvent);
      GlueXUserEventInformation *eventinfo;
      eventinfo = (GlueXUserEventInformation*)anEvent->GetUserInformation();
      if (eventinfo) {
         double E = eventinfo->GetBeamPhotonEnergy();
         GlueXPhotonBeamGenerator::GenerateTaggerHit(anEvent, E, ttag);
      }
      else {
         return;
      }
   }
   ++fEventCount;

   // Superimpose any requested background minimum-bias beam interactions

   if (fBeamBackgroundRate > 0 && fPhotonBeamGenerator != 0) {
      double t = fBeamBackgroundGateStart;
      while (true) {
         t += -log(G4UniformRand()) / fBeamBackgroundRate;
         if (t > fBeamBackgroundGateStop)
            break;
         fPhotonBeamGenerator->GenerateBeamPhoton(anEvent, t);
      }
   }
}

void GlueXPrimaryGeneratorAction::GeneratePrimariesCobrems(G4Event* anEvent)
{
   if (fPhotonBeamGenerator != 0) {
      fPhotonBeamGenerator->GeneratePrimaryVertex(anEvent);
   }
   ++fEventCount;
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
      case 27  : return -3222;    // anti-Sigma-
      case 28  : return -3212;    // anti-Sigma0
      case 29  : return -3112;    // anti-Sigma+
      case 30  : return -3322;    // anti-Xi0
      case 31  : return -3312;    // anti-Xi+
      case 32  : return -3334;    // anti-Omega+ 
      case 45  : return 1000010020; // deuteron
      case 46  : return 1000010030; // triton
      case 47  : return 1000020040; // alpha
      case 48  : return 0;        // geantino (no PDG type)
      case 49  : return 1000020030; // He3 ion
      case 50  : return 0;        // Cerenkov photon (no PDG type)

      // These are "private" Geant3 types that were defined in hdgeant
      case 33  : return 223;      // omega(782)
      case 34  : return 333;      // phi(1020)
      case 35  : return 331;      // etaPrime(958)
      case 36  : return 0;        // unused
      case 37  : return 0;        // unused
      case 38  : return 0;        // unused
      case 39  : return 0;        // unused
      case 40  : return 0;        // unused
      case 41  : return 0;        // unused
      case 42  : return 213;      // rho(770)+
      case 43  : return -213;     // rho(770)-
      case 44  : return 113;      // rho(770)0

      case 82  : return 2224;     // Delta++
      case 83  : return 443;      // Jpsi
      case 84  : return 441;      // Eta_c
      case 85  : return 10441;    // Chi_c0
      case 86  : return 20443;    // Chi_c1
      case 87  : return 445;      // Chi_c2
      case 88  : return 100443;   // Psi2s
      case 89  : return 421;      // D0
      case 90  : return 411;      // D+
      case 91  : return 10421;    // Dstar0
      case 92  : return 10411;    // Dstar+
      case 93  : return 4022;     // Lambda_c+
      case 94  : return -421;     // anti-D0

      case 163  : return 9000111; // a0(980)
      case 164  : return 9010221; // f0(980)
      case 165  : return 313;     // K*(892)0
      case 166  : return 323;     // K*(892)+
      case 167  : return -323;    // K*(892)-
      case 168  : return -313;    // anti-K*(892)0
      case 169  : return 20323;   // K1(1400)+
      case 170  : return -20323;  // K1(1400)-
      case 171  : return 4122;    // b1(1235)+
      case 172  : return 3224;    // Sigma*(1385)+
      case 173  : return 3214;    // Sigma*(1385)0
      case 174  : return 3114;    // Sigma*(1385)-

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
      case          0 : return 50;    // optical photon
      case         22 : return 1;     // photon
      case        -11 : return 2;     // e+
      case         11 : return 3;     // e-
      case        -12 : return 4;     // anti-e-neutrino
      case         12 : return 4;     // e-neutrino
      case        -13 : return 5;     // mu+
      case         13 : return 6;     // mu-
      case        -14 : return 4;     // anti-mu-neutrino
      case         14 : return 4;     // mu-neutrino
      case        -16 : return 4;     // anti-tau-neutrino
      case         16 : return 4;     // tau-neutrino
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
      case       3334 : return 24;    // Omega-
      case      -2112 : return 25;    // anti-neutron
      case      -3122 : return 26;    // anti-Lambda
      case      -3222 : return 27;    // Sigma-
      case      -3212 : return 28;    // Sigma0
      case      -3112 : return 29;    // Sigma+
      case      -3322 : return 30;    // Xi0
      case      -3312 : return 31;    // Xi+
      case      -3334 : return 32;    // Omega+
      case 1000010020 : return 45;    // deuteron
      case 1000010030 : return 46;    // triton
      case 1000020040 : return 47;    // alpha
      case 1000020030 : return 49;    // He3

      // These are "private" Geant3 types that were defined in hdgeant
      case 223        : return 33;    // omega(782)
      case 333        : return 34;    // phi(1020)
      case 331        : return 35;    // etaPrime(958)
      case 213        : return 42;    // rho(770)+
      case -213       : return 43;    // rho(770)-
      case 113        : return 44;    // rho(770)0

      case 2224       : return 82;    // Delta++
      case 443        : return 83;    // Jpsi
      case 441        : return 84;    // Eta_c
      case 10441      : return 85;    // Chi_c0
      case 20443      : return 86;    // Chi_c1
      case 445        : return 87;    // Chi_c2
      case 100443     : return 88;    // Psi2s
      case 421        : return 89;    // D0
      case 411        : return 90;    // D+
      case 10421      : return 91;    // Dstar0
      case 10411      : return 92;    // Dstar+
      case 4022       : return 93;    // Lambda_c+
      case -421       : return 94;    // anti-D0

      case 9000111    : return 163;   // a0(980)
      case 9010221    : return 164;   // f0(980)
      case 313        : return 165;   // K*(892)0
      case 323        : return 166;   // K*(892)+
      case -323       : return 167;   // K*(892)-
      case -313       : return 168;   // anti-K*(892)0
      case 20323      : return 169;   // K1(1400)+
      case -20323     : return 170;   // K1(1400)-
      case 4122       : return 171;   // b1(1235)+
      case 3224       : return 172;   // Sigma*(1385)+
      case 3214       : return 173;   // Sigma*(1385)0
      case 3114       : return 174;   // Sigma*(1385)-

      default  :
         if (PDGtype < 1000000000) {
            G4cout << "Warning in GlueXPrimaryGeneratorAction::"
                      "ConvertPdgToGeant3 - lookup performed on unknown"
                      " PDG particle type " << PDGtype << ","
                      " returning 0 for the Geant3 particle code." << G4endl;
         }
   }
   return 0;
}

double GlueXPrimaryGeneratorAction::GetMassPDG(int PDGtype)
{
   return fParticleTable->FindParticle(PDGtype)->GetPDGMass();
}

double GlueXPrimaryGeneratorAction::GetMass(int Geant3Type)
{
   return GetMassPDG(ConvertGeant3ToPdg(Geant3Type));
}
