//
/// \brief Implementation of the DetectorConstruction class
//

#include <algorithm>
#include <functional>
#include <utility>
#include <typeinfo>

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Para.hh"
#include "G4Trd.hh"
#include "G4Polycone.hh"
#include "G4GenericTrap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "G4Mag_UsualEqRhs.hh"
#include "G4IntegrationDriver.hh"
#include "G4ChordFinder.hh"
#include "G4DormandPrince745.hh"
#include "G4CashKarpRKF45.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4AssemblyVolume.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4UserLimits.hh"
#include "N6SetUp.hh"
#include "LxDetector.hh"
#include "N6MuTracker.hh"
#include "N6SiTracker.hh"
#include "N6TorMagnet.hh"
#include "N6Magnets.hh"
#include "N6BFields.hh"


DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fWorldMaterial(0),fDefaultWorld(true),
 fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),
 fDetectorMessenger(0), fMagnetZPos(0), fMagnetSizeZ(0)
{
  // default parameter values of the calorimeter
//  N6SetUp *lxs = N6SetUp::Instance();
  ComputeCalorParameters();
  
  // materials  
  DefineMaterials();
  SetWorldMaterial   (N6SetUp::Instance()->EnvironmentMaterial);

 
  // create commands for interactive definition of the calorimeter  
  fDetectorMessenger = new DetectorMessenger(this);
}



DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
  for (auto itr : fDetList) {
    if (itr.second) {delete itr.second;} 
  } 
}


void DetectorConstruction::AddSensorSegmentation(const G4String sname, const G4double xsize, 
                                                 const G4double ysize, const G4int nx, const G4int ny)
{
  EventAction *evatmp = new EventAction();
  int res = evatmp->TestSegmentationEncoding(nx, ny);
  delete evatmp;
  if (!res) {
    fSensors[sname] = std::make_tuple (xsize, ysize, nx, ny);
  } else {
    G4String msgstr("Failed to add sensitive detector ");
    msgstr += sname;
    G4Exception("DetectorConstruction::", "AddSensorSegmentation()", FatalException, msgstr.c_str());
  }
}



G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}



void DetectorConstruction::DefineMaterials()
{
  //This function illustrates the possible ways to define materials

  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;

  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;

  //
  // define Elements
  //

  G4Element* H  = new G4Element("Hydrogen",symbol="H",  z= 1, a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
  G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a=  22.99*g/mole);
  G4Element* Ar = new G4Element("Argon",   symbol="Ar", z=18, a=  39.95*g/mole);
  G4Element* I  = new G4Element("Iodine",  symbol="I" , z=53, a= 126.90*g/mole);
  G4Element* Xe = new G4Element("Xenon",   symbol="Xe", z=54, a= 131.29*g/mole);
  G4Element* Be = new G4Element("Beryllium",   symbol="Be", z=4, a= 9.01*g/mole);

  //
  // define simple materials
  //

  new G4Material("H2Liq"    , z= 1, a= 1.01*g/mole, density= 70.8*mg/cm3);
  new G4Material("Beryllium", z= 4, a= 9.01*g/mole, density= 1.848*g/cm3);
  new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);
  new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3);

  G4Material* lAr = new G4Material("liquidArgon", density= 1.390*g/cm3, ncomponents=1);
  lAr->AddElement(Ar, natoms=1);

  new G4Material("Iron",     z=26, a= 55.85*g/mole, density= 7.870*g/cm3);
  new G4Material("Copper",   z=29, a= 63.55*g/mole, density= 8.960*g/cm3);
  new G4Material("Germanium",z=32, a= 72.61*g/mole, density= 5.323*g/cm3);
  new G4Material("Silver",   z=47, a=107.87*g/mole, density= 10.50*g/cm3);
  new G4Material("Gold",     z=79, a=196.97*g/mole, density= 19.32*g/cm3);
  new G4Material("Lead",     z=82, a=207.19*g/mole, density= 11.35*g/cm3);
  
  G4Material* Tungsten = new G4Material("Tungsten", z=74, a=183.85*g/mole, density= 19.30*g/cm3);
  //
  // define a material from elements.   case 1: chemical molecule
  //ConstructCalorimeter

  G4Material* H2O = new G4Material("Water",density= 1.000*g/cm3,ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78*eV);

  G4Material* CH = new G4Material("Plastic",density= 1.04*g/cm3,ncomponents=2);
  CH->AddElement(C, natoms=1);
  CH->AddElement(H, natoms=1);

  G4Material* NaI = new G4Material("NaI", density= 3.67*g/cm3, ncomponents=2);
  NaI->AddElement(Na, natoms=1);
  NaI->AddElement(I , natoms=1);
  NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* Air20 = 
  new G4Material("Air20", density= 1.205*mg/cm3, ncomponents=2,
                   kStateGas, 293.*kelvin, 1.*atmosphere);
  Air20->AddElement(N, fractionmass=0.7);
  Air20->AddElement(O, fractionmass=0.3);

  G4Material* Air20Refractive =
  new G4Material("Air20Refractive", density= 1.205*mg/cm3, ncomponents=2,
                   kStateGas, 293.*kelvin, 1.*atmosphere);
  Air20Refractive->AddElement(N, fractionmass=0.7);
  Air20Refractive->AddElement(O, fractionmass=0.3);

  //Graphite
  //
  G4Material* Graphite = new G4Material("Graphite", density= 1.7*g/cm3, ncomponents=1);
  Graphite->AddElement(C, fractionmass=1.);
  
   //BeO
  // fBeO = fNistMan->FindOrBuildMaterial("G4_BERYLLIUM_OXIDE");;
//  G4Material* BeO;
//  BeO = G4NistManager::Instance()->FindOrBuildMaterial("G4_BERYLLIUM_OXIDE");
   G4Material* BeO = new G4Material("BeO", density= 3.02*g/cm3, ncomponents=2);
   BeO->AddElement(Be, natoms=1);
   BeO->AddElement(O, natoms=1);

  //Havar
  //
  G4Element* Cr = new G4Element("Chrome", "Cr", z=25, a=  51.996*g/mole);
  G4Element* Fe = new G4Element("Iron"  , "Fe", z=26, a=  55.845*g/mole);
  G4Element* Co = new G4Element("Cobalt", "Co", z=27, a=  58.933*g/mole);
  G4Element* Ni = new G4Element("Nickel", "Ni", z=28, a=  58.693*g/mole);
  G4Element* W  = new G4Element("Tungsten","W", z=74, a= 183.850*g/mole);

  G4Material* Havar = new G4Material("Havar", density= 8.3*g/cm3, ncomponents=5);
  Havar->AddElement(Cr, fractionmass=0.1785);
  Havar->AddElement(Fe, fractionmass=0.1822);
  Havar->AddElement(Co, fractionmass=0.4452);
  Havar->AddElement(Ni, fractionmass=0.1310);
  Havar->AddElement(W , fractionmass=0.0631);

  //
  // examples of gas
  //
//  G4Material* ArgonGas = new G4Material("ArgonGas", z=18, a=39.948*g/mole, density= 1.782*mg/cm3,
//                 kStateGas, 293.15*kelvin, 1*atmosphere);

//  G4Material* HeliumGas = new G4Material("HeliumGas", z=2, a=4.003*g/mole, density = 0.166322 *mg/cm3, kStateGas, 293.15*kelvin, 1*atmosphere);

  new G4Material("XenonGas", z=54, a=131.29*g/mole, density= 5.458*mg/cm3,
                 kStateGas, 293.15*kelvin, 1*atmosphere);

  G4Material* CO2 =
    new G4Material("CarbonicGas", density= 1.977*mg/cm3, ncomponents=2);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);

  G4Material* ArCO2 =
    new G4Material("ArgonCO2",   density= 1.8223*mg/cm3, ncomponents=2);
  ArCO2->AddElement (Ar,  fractionmass=0.7844);
  ArCO2->AddMaterial(CO2, fractionmass=0.2156);

  //another way to define mixture of gas per volume
  G4Material* NewArCO2 =
    new G4Material("NewArgonCO2", density= 1.8223*mg/cm3, ncomponents=3);
  NewArCO2->AddElement (Ar, natoms=8);
  NewArCO2->AddElement (C,  natoms=2);
  NewArCO2->AddElement (O,  natoms=4);

  G4Material* ArCH4 = 
    new G4Material("ArgonCH4",    density= 1.709*mg/cm3,  ncomponents=3);
  ArCH4->AddElement (Ar, natoms=93);
  ArCH4->AddElement (C,  natoms=7);
  ArCH4->AddElement (H,  natoms=28);

  G4Material* XeCH = 
    new G4Material("XenonMethanePropane", density= 4.9196*mg/cm3, ncomponents=3,
                   kStateGas, 293.15*kelvin, 1*atmosphere);
  XeCH->AddElement (Xe, natoms=875);
  XeCH->AddElement (C,  natoms=225);
  XeCH->AddElement (H,  natoms=700);

  G4Material* steam = new G4Material("WaterSteam", density= 1.0*mg/cm3, ncomponents=1);
  steam->AddMaterial(H2O, fractionmass=1.);
  steam->GetIonisation()->SetMeanExcitationEnergy(71.6*eV);

  //
  // example of vacuum
  //

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1, a=1.01*g/mole,density, kStateGas,temperature,pressure);

  // Define Stainless steal

  G4Element* Si = new G4Element("Silicon", "Si", z=14, a=28.09 *g/mole);
  G4Element* Mn  = new G4Element("Manganese","Mn", z=25, a=54.938044 *g/mole);

  G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=6);
  StainlessSteel->AddElement(C, fractionmass=0.001);
  StainlessSteel->AddElement(Si, fractionmass=0.007);
  StainlessSteel->AddElement(Cr, fractionmass=0.18);
  StainlessSteel->AddElement(Mn, fractionmass=0.01);
  StainlessSteel->AddElement(Fe, fractionmass=0.712);
  StainlessSteel->AddElement(Ni, fractionmass=0.09);  

  G4Element* Al = 	new G4Element("Aluminum", "Al", 13., 26.9815386 * g/mole);
//  G4Element* Ce = 	new G4Element("Cerium", "Ce", 58., 140.12 * g/mole);
  G4Element* Ca = 	new G4Element("Calcium", "Ca", 20.0, 40.078 * g/mole);

     // Concrete, must  check recipe for concrete

  G4double crdensity = 2.5*g/cm3;
  G4Material* ShieldingConcrete = new G4Material("ShieldingConcrete", crdensity, 6);
  ShieldingConcrete->AddElement(O,  0.52);
  ShieldingConcrete->AddElement(Si, 0.325);
  ShieldingConcrete->AddElement(Ca, 0.06);
  ShieldingConcrete->AddElement(Na, 0.015);
  ShieldingConcrete->AddElement(Fe, 0.04);
  ShieldingConcrete->AddElement(Al, 0.04);

  G4Material* fPstyrene = new G4Material("Polystyrene", density= 1.03*g/cm3, 2);
  fPstyrene->AddElement(C, 8);
  fPstyrene->AddElement(H, 8);
  G4Material* lanex;
  lanex = G4NistManager::Instance()->FindOrBuildMaterial("G4_GADOLINIUM_OXYSULFIDE");
 

 //Honeycomb
  G4Material* Cellulose;
  Cellulose = G4NistManager::Instance()->FindOrBuildMaterial("G4_CELLULOSE_CELLOPHANE");
//  G4Material* Honeycomb = new G4Material("Honeycomb", density= 0.030*g/cm3, ncomponents=2);
  G4Material* Honeycomb = new G4Material("Honeycomb", density= 0.048*g/cm3, ncomponents=2); // Nomex honeycomb
  Honeycomb->AddMaterial(Air, fractionmass=0.47);
  Honeycomb->AddMaterial(Cellulose, fractionmass=0.53);


  density     = 1.250 *mg/cm3;
  temperature = 300.*kelvin;
  pressure    = 1.0*atmosphere;
  G4Material* CO = new G4Material("CarbonMonoxide", density, ncomponents=2, kStateGas, temperature, pressure);
  CO->AddElement(C, natoms=1);
  CO->AddElement(O, natoms=1);

  density     = 8.376E-05 *g/cm3;
  pressure    = 1.0*atmosphere;
  G4Material* H2Gass = new G4Material("HydrogenGass", density, ncomponents=1, kStateGas, temperature, pressure);
  H2Gass->AddElement(H, natoms=2);

  density     = 2.89e-15 *mg/cm3;
  pressure    = 1.0e-6 *pascal;
  G4Material* xfelVacuum = new G4Material("XFELVacuum", density, ncomponents=2, kStateGas, temperature, pressure);
  xfelVacuum->AddMaterial(CO, 0.78);
  xfelVacuum->AddMaterial(H2Gass, 0.22);
}



void DetectorConstruction::ComputeCalorParameters()
{

  if (fDefaultWorld) {
     fWorldSizeZ = 50.0*m;
     fWorldSizeXY= 30.0*m;
  }
}


  
G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{ 
  // Cleanup old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // World
  //
  fSolidWorld = new G4Box("World", fWorldSizeXY/2.0, fWorldSizeXY/2.0, fWorldSizeZ/2.0); 
  fLogicWorld = new G4LogicalVolume(fSolidWorld, fWorldMaterial, "World");
  fPhysiWorld = new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, true);

  ConstructLuxeDetectors();
  ConstructTorMagnetField();
  ConstructIPMagnet();

//   AssignRegions();
// ////////////////////
//   DumpBFieldModel();

  PrintCalorParameters();         
  
  return fPhysiWorld;
}



void DetectorConstruction::ConstructLuxeDetectors()
{
  fDetList["N6MuTracker"] = new N6MuTracker(this);
  fDetList["N6SiTracker"] = new N6SiTracker(this);
  fDetList["N6TorMagnet"] = new N6TorMagnet(this);

  for (auto &nd : fDetList) { nd.second->Construct(); }
}

/*

void DetectorConstruction::ConstructInfrastructure()
{
  N6SetUp *lxs = N6SetUp::Instance();
  G4Material* concreteMaterial = G4NistManager::Instance()->FindOrBuildMaterial("ShieldingConcrete");
  
  G4Box *solidFloor = new G4Box("solidFloor", fWorldSizeXY/2.0, lxs->FloorY/2.0, fWorldSizeZ/2.0);
  G4LogicalVolume *logicFloor = new G4LogicalVolume(solidFloor, concreteMaterial, "logicFloor");
  G4double ypos = lxs->FloorSurfaceYpos - lxs->FloorY/2.0;   // lxs->FloorSurfaceYpos < 0
  new G4PVPlacement(0, G4ThreeVector(0.0, ypos, 0), logicFloor, "Floor", fLogicWorld, false, 0, lxs->OverlapTest);
}

*/

void DetectorConstruction::ConstructTorMagnetField(const G4String magType)
{
  N6SetUp *lxs = N6SetUp::Instance();
  G4ThreeVector trm(0.0, 0.0, lxs->TMagnetZpos);

  G4LogicalVolume *fieldvol = dynamic_cast<N6TorMagnet*>(fDetList["N6TorMagnet"])->GetFieldVolume();

// Add local magnetic field
  G4String magid = "Tor";
  TorBField *torfield = new TorBField(trm);
// Parameter are adjusted to provide the required field shape.
  G4double magR0Max, magR1Max;
  dynamic_cast<N6TorMagnet*>(fDetList["N6TorMagnet"])->GetFieldVolumeParameters(magR0Max, magR1Max);
  torfield->SetParameters(250.0*tesla, 200.0*mm, magR0Max, 200.0*mm, magR1Max, -2000.0*mm, 2000.0*mm, 20.0*mm);
//  torfield->SetParameters(250.0*tesla, 200.0*mm, 1000.0*mm, 200.0*mm, 3000.0*mm, -2000.0*mm, 2000.0*mm, 20.0*mm);

  G4FieldManager* fieldMgr = new G4FieldManager(torfield);
  ConfigureFieldManager(fieldMgr, torfield);
 
  fieldvol->SetFieldManager(fieldMgr, true);
  G4AutoDelete::Register(torfield);
  G4AutoDelete::Register(fieldMgr);

//   torfield->DumpBFieldModel(M_PI/3.0);

}




void DetectorConstruction::ConstructIPMagnet(const G4String magType)
{
  N6SetUp *lxs = N6SetUp::Instance();

  G4String mtype(magType);
//   G4String mtype("MEP48Magnet");
  N6MagnetAssembly *mag = new MEP48MagnetAssembly(mtype);

  G4ThreeVector trm(lxs->DipoleMagnetXpos, lxs->DipoleMagnetZpos, lxs->DipoleMagnetZpos);
  mag->GetAssembly()->MakeImprint(fLogicWorld, trm, 0, 0, lxs->OverlapTest);
//   mag->GetAssembly()->MakeImprint(fLogicWorld, trm, 
//                               new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.5*M_PI), 0, lxs->OverlapTest);
//   mag->CostructSupport(fLogicWorld, trm, "DumpMagnet", true);

}



void DetectorConstruction::AddFieldToLogVolumes(G4FieldManager* fieldMgr, const std::vector<G4String> vname)
{
  G4LogicalVolumeStore *lvstor = G4LogicalVolumeStore::GetInstance();
  for (auto &vitr : vname) {
    G4LogicalVolume  *lv = lvstor->GetVolume (vitr);
    if (lv) {
      lv->SetFieldManager(fieldMgr, true);
      G4cout << "Add magnetic field to " << vitr << G4endl;
    } else G4cout << "!!!!!!!!!!!!!!!! Field was not added to volume " << vitr << " !!!!!!!!!!!!!!!!\n";
  }
}



void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n" << fWorldMaterial    << G4endl;
  G4cout << "\n The  WORLD   is made of "  << G4BestUnit(fWorldSizeZ,"Length")
         << " of " << fWorldMaterial->GetName();
  G4cout << ". The transverse size (XY) of the world is " 
         << G4BestUnit(fWorldSizeXY,"Length") << G4endl;
  G4cout << G4endl;
}


void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fWorldMaterial != pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if(fLogicWorld) fLogicWorld->SetMaterial(fWorldMaterial);
//     if(fLogicMagnet) fLogicMagnet->SetMaterial(fWorldMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}



void DetectorConstruction::SetWorldSizeZ(G4double val)
{
  fWorldSizeZ = val;
  fDefaultWorld = false;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetWorldSizeXY(G4double val)
{
  fWorldSizeXY = val;
  fDefaultWorld = false;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


void DetectorConstruction::ConstructSDandField()
{
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg = new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );
        
    }
}



void DetectorConstruction::ConfigureFieldManager(G4FieldManager *fieldMgr, G4MagneticField *mfield)
{
  G4double minStep = 0.010*mm;
  auto pEquation = new G4Mag_UsualEqRhs(mfield);
//   G4int nvar = pEquation->GetNumberOfVariables();
  G4int nvar = 8;
  auto pStepper = new G4DormandPrince745( pEquation, nvar );
  auto pIntgrationDriver = new G4IntegrationDriver<G4DormandPrince745>(minStep, pStepper, nvar);

  G4ChordFinder *pChordFinder = new G4ChordFinder(pIntgrationDriver);
  fieldMgr->SetChordFinder( pChordFinder );

  fieldMgr->SetMinimumEpsilonStep( 1.0e-5 );
  fieldMgr->SetMaximumEpsilonStep( 1.0e-4 );
  fieldMgr->SetDeltaOneStep( 0.5e-3 * mm );
  fieldMgr->SetDeltaIntersection(0.5e-3 * mm);

//   // magnetic field
//   G4FieldManager fieldMgr = new G4FieldManager();
//   fieldMgr->SetDetectorField(mfield);
//   fieldMgr->CreateChordFinder(mfield);
//   G4bool forceToAllDaughters = true;
//   fieldvol->SetFieldManager(fieldMgr, forceToAllDaughters);

}



void DetectorConstruction::AssignRegions()
{ 
/*  G4Region* region;
  G4String regName;
  G4ProductionCuts* cuts;

  N6SetUp *lxs = N6SetUp::Instance();

  regName = "BDumpRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  if (!region) {
    region = new G4Region(regName);
  }
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(lxs->BeamDumpProductionCut,G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(lxs->BeamDumpProductionCut, G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(lxs->BeamDumpProductionCut, G4ProductionCuts::GetIndex("e+"));
//   cuts->SetProductionCut(lxs->BeamDumpProductionCut); // same cuts for gamma, proton, e- and e+
  region->SetProductionCuts(cuts);
  G4LogicalVolume *bdump = G4LogicalVolumeStore::GetInstance()->GetVolume("logicBeamDump");
  if (bdump) {
    bdump->SetRegion(region);
    region->AddRootLogicalVolume(bdump);
  }

*/
}

/*
void DetectorConstruction::SetAbsorberType(G4String val)
{
  if (val == "foil") {
    N6SetUp::Instance()->GTargetType = N6SetUp::tTargetType::tfoil;
  } else if (val == "wire") {
    N6SetUp::Instance()->GTargetType = N6SetUp::tTargetType::twire;
  } else {
    G4cout << "DetectorConstruction::SetAbsorberType: <" << val << ">"
           << " is not supported!"
           << G4endl;
  }
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
*/


void DetectorConstruction::DumpBFieldModel()
{
  G4cout << "======= Following setting are used to create magnetic fields =======\n"; 
  for (const auto &mitr : fBFieldModelsInfo) {
    const auto &tinfo = mitr.second;
    G4cout << "MagnetID: " << mitr.first << G4endl;
    for (const auto &vitr : tinfo) {
      G4cout << "     " << std::get<0>(vitr) << "  " << std::get<1>(vitr)
             << "  " << std::get<2>(vitr) << "  " << std::get<3>(vitr) << G4endl;
    }
  }
  G4cout << "====================================================================\n"; 
}
