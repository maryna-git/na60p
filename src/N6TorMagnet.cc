//
/// \brief Implementation of the BeamProfiler tracker class
//

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4GenericTrap.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4AssemblyVolume.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4UserLimits.hh"
#include "N6SetUp.hh"
#include "LxDetector.hh"
#include "N6TorMagnet.hh"




void N6TorMagnet::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  CreateMaterial();

  N6SetUp *lxs = N6SetUp::Instance();
  G4Material* environmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* magnetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Material4Toroid");
//   G4Material* CoreMagnetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Epoxy");

///Toroid magnet
  G4double  tMagnetZ = 3350.0 *mm; //3230.0 *mm;

//  G4double  tMagRMax = 2971.0 *mm;
//  G4double  tMagRMin = 1950.0 *mm;
//   G4int  tMagNloops = 12;

  G4double tMagRbend = 250.0* mm; //??
  G4double tMagRinCut = 101.0 * mm;
  G4double tMagCoilThickness = 64.4* mm;
  G4double tMagTrHightWidth = 332.0* mm;

  G4double tMagCoilWidth = 783.0* mm;
  G4double tMagSmallTrbaseH = 1320.0* mm - tMagRinCut - tMagTrHightWidth;

  G4double tMagBigTrbaseH = lxs->tMagRMin - tMagRinCut - tMagTrHightWidth;

  G4double holezin = tMagnetZ - 2.0 * (tMagCoilWidth + tMagRbend);
  G4double sideAngle = atan((tMagBigTrbaseH - tMagSmallTrbaseH)/holezin);
//  G4double tMagContainerRmax = (tMagCoilWidth)*cos(sideAngle); 

//  G4double tMagContainerRmin = tMagRinCut + tMagTrHightWidth + tMagSmallTrbaseH +(tMagCoilWidth + tMagRbend )*tan(M_PI/4.0 - sideAngle/2.0);
//  tMagContainerRmin = sqrt(pow(tMagContainerRmin, 2.0) + pow(tMagCoilThickness/2.0, 2.0));
//  G4double tMagContainerRmax = tMagContainerRmin + tMagnetZ*tan(sideAngle);

  fCoilR0Max = tMagRinCut + tMagTrHightWidth + tMagSmallTrbaseH +(tMagCoilWidth + tMagRbend )*tan(M_PI/4.0 - sideAngle/2.0);
  fCoilR0Max = sqrt(pow(fCoilR0Max, 2.0) + pow(tMagCoilThickness/2.0, 2.0));
  fCoilR1Max = fCoilR0Max + tMagnetZ*tan(sideAngle);

  fCoilSideAngle = sideAngle;
  G4double tMagContainerRmin = fCoilR0Max;
  G4double tMagContainerRmax = fCoilR1Max;
  
  G4Cons *solidTMagnetContainer = new G4Cons("solidTMagnetContainer", tMagRinCut, 1.05*tMagContainerRmin ,
                                                tMagRinCut, 1.05*tMagContainerRmax, tMagnetZ/2.0, 0.0, 2.0*M_PI);

  G4LogicalVolume *logicTMagnetContainer = new G4LogicalVolume(solidTMagnetContainer, environmentMaterial, "logicTMagnetContainer");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->TMagnetZpos), logicTMagnetContainer,
                     "ToroidContainer", fLogicWorld, false, 0, lxs->OverlapTest);
  fFieldVolume = logicTMagnetContainer;

 //Toroid Sector coil
  G4Box *solidTMagSmallTrbase = new G4Box("solidTMagSmallTrbase", tMagCoilThickness/2.0, tMagSmallTrbaseH/2.0,
                                     tMagCoilWidth/2.0);
  G4Box *solidTMagBigTrbase = new G4Box("solidTMagBigTrbase", tMagCoilThickness/2.0, tMagBigTrbaseH/2.0,
                                     tMagCoilWidth/2.0);
  G4Box *solidTMagTrHight = new G4Box("solidTMagTrHight", tMagCoilThickness/2.0, tMagTrHightWidth/2.0,
                                     tMagnetZ/2.0 - tMagTrHightWidth);

  G4double zside = holezin/cos(sideAngle);
  G4Box *solidTMagCoilSide = new G4Box("solidTMagCoilSide", tMagCoilThickness/2.0, tMagCoilWidth/2.0, zside/2.0);

  G4Tubs *solidTMagSmallTrbaseRounding = new G4Tubs("solidTMagSmallTrbaseRounding", tMagRbend, tMagRbend + tMagCoilWidth,
                                                    tMagCoilThickness/2.0, 0.0, 0.5*M_PI - sideAngle);
  G4Tubs *solidTMagBigTrbaseRounding = new G4Tubs("solidTMagBigTrbaseRounding", tMagRbend, tMagRbend + tMagCoilWidth,
                                                       tMagCoilThickness/2.0, 0.0, 0.5*M_PI + sideAngle);
  G4Tubs *solidTMagBottomRounding = new G4Tubs("solidTMagBottomRounding", 0.0, tMagTrHightWidth,
                                                       tMagCoilThickness/2.0, 0.0, 0.5*M_PI);

  std::vector<G4TwoVector> auxvtx{G4TwoVector(-tMagCoilThickness, tMagCoilThickness), G4TwoVector(tMagCoilThickness*2.0, tMagCoilThickness),
	                         G4TwoVector(tMagCoilThickness*2.0, -4.0*tMagCoilThickness), G4TwoVector(-tMagCoilThickness, -tMagCoilThickness),
	                             G4TwoVector(-tMagCoilThickness, tMagCoilThickness), G4TwoVector(tMagCoilThickness*2.0, tMagCoilThickness),
	                         G4TwoVector(tMagCoilThickness*2.0, -4.0*tMagCoilThickness), G4TwoVector(-tMagCoilThickness, -tMagCoilThickness)};
  G4GenericTrap *solidTMagCore = new G4GenericTrap("solidTMagCore", (tMagnetZ-tMagCoilWidth)/2.0, auxvtx);
  
  G4LogicalVolume *logicTMagSmallTrbase = new G4LogicalVolume(solidTMagSmallTrbase, magnetMaterial, "logicTMagSmallTrbase");
  G4LogicalVolume *logicTMagBigTrbase = new G4LogicalVolume(solidTMagBigTrbase, magnetMaterial, "logicTMagBigTrbase");
  G4LogicalVolume *logicTMagTrHight = new G4LogicalVolume(solidTMagTrHight, magnetMaterial, "logicTMagTrHight");
  G4LogicalVolume *logicTMagBottomRounding = new G4LogicalVolume(solidTMagBottomRounding, magnetMaterial, "logicTMagBottomRounding");
  G4LogicalVolume *logicTMagCoilSide = new G4LogicalVolume(solidTMagCoilSide, magnetMaterial, "logicTMagCoilSide");
  G4LogicalVolume *logicTMagSmallTrbaseRounding = new G4LogicalVolume(solidTMagSmallTrbaseRounding, magnetMaterial, "logicTMagSmallTrbaseRounding");
  G4LogicalVolume *logicTMagBigTrbaseRounding = new G4LogicalVolume(solidTMagBigTrbaseRounding, magnetMaterial, "logicTMagBigTrbaseRounding");
  G4LogicalVolume *logicTMagCore = new G4LogicalVolume(solidTMagCore, magnetMaterial, "logicTMagCore");

  G4AssemblyVolume* assemblyMuSector = new G4AssemblyVolume();

  G4int  TMagNsectors = 8; //8;

  G4double ylpos = tMagRinCut + tMagTrHightWidth/2.0;
  G4ThreeVector vtr(0.0, ylpos, 0.0);
  assemblyMuSector->AddPlacedVolume(logicTMagTrHight, vtr, 0);
  vtr.setY(ylpos + tMagTrHightWidth/2.0);
  vtr.setZ(tMagnetZ/2.0 - tMagTrHightWidth);
  assemblyMuSector->AddPlacedVolume(logicTMagBottomRounding, vtr, new G4RotationMatrix(M_PI/2.0, M_PI/2.0, M_PI/2.0));
  vtr.setZ(-vtr.z());
  assemblyMuSector->AddPlacedVolume(logicTMagBottomRounding, vtr, new G4RotationMatrix(M_PI/2.0, -M_PI/2.0, M_PI/2.0));
  
  G4ThreeVector vtraptr(-2.5*tMagCoilThickness, ylpos - tMagTrHightWidth/2.0 + 4.0*tMagCoilThickness, 0.0);
  assemblyMuSector->AddPlacedVolume(logicTMagCore, vtraptr, 0);

  vtr.setY(ylpos + (tMagTrHightWidth + tMagSmallTrbaseH)/2.0);
  vtr.setZ(-0.5*(tMagnetZ - tMagCoilWidth));
  assemblyMuSector->AddPlacedVolume(logicTMagSmallTrbase, vtr, 0);
  vtr.setY(ylpos + (tMagTrHightWidth + tMagBigTrbaseH)/2.0);
  vtr.setZ(-vtr.z());
  assemblyMuSector->AddPlacedVolume(logicTMagBigTrbase, vtr, 0);

  G4double sideypos = ylpos + tMagTrHightWidth/2.0 + tMagSmallTrbaseH
                    + (tMagRbend + 0.5*tMagCoilWidth) * cos(sideAngle) + 0.5*zside*sin(sideAngle);
  G4double sidezpos = -(tMagRbend + tMagCoilWidth/2.0) * sin(sideAngle);
  vtr.setY(sideypos);
  vtr.setZ(sidezpos);
  assemblyMuSector->AddPlacedVolume(logicTMagCoilSide, vtr, new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), -sideAngle));

  vtr.setY(ylpos + tMagTrHightWidth/2.0 + tMagSmallTrbaseH);
  vtr.setZ(tMagCoilWidth + tMagRbend - tMagnetZ/2.0);
  assemblyMuSector->AddPlacedVolume(logicTMagSmallTrbaseRounding, vtr, new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0));
  vtr.setY(ylpos + tMagTrHightWidth/2.0 + tMagBigTrbaseH);
  vtr.setZ(-vtr.z());
  assemblyMuSector->AddPlacedVolume(logicTMagBigTrbaseRounding, vtr, new G4RotationMatrix(G4ThreeVector(0.0, -1.0, 0.0), M_PI/2.0));

  G4ThreeVector sectortr;
  for( G4int i = 0; i < TMagNsectors; i++ )
  {
// Translation of the assembly inside the logicTMagnetContainer

     G4double dphi = i * 2.0*M_PI/static_cast<G4double>(TMagNsectors);
     assemblyMuSector->MakeImprint(logicTMagnetContainer, sectortr, new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), dphi));
  }
//   AddSegmentation();
}



void N6TorMagnet::AddSegmentation()
{
  N6SetUp *lxs = N6SetUp::Instance();
  G4double bsmcalox = lxs->BSMCaloX;
  G4double bsmcaloy = lxs->BSMCaloY;
  G4int ncellx = lxs->BSMCaloNCellX;
  G4int ncelly = lxs->BSMCaloNCellY;

  fDetector->AddSensorSegmentation("BSMCaloLayer", bsmcalox, bsmcaloy, ncellx, ncelly);
}



G4AssemblyVolume* N6TorMagnet::ConstructSupportAssembly()
{
//   N6SetUp *lxs = N6SetUp::Instance();
  G4AssemblyVolume *supportAssembly = new G4AssemblyVolume();
  return supportAssembly;
}



void N6TorMagnet::CreateMaterial()
{
G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  
// 
  G4int ncomponents;
  G4double fractionmass;

  G4NistManager *materials = G4NistManager::Instance();
  
  G4Material* Copper = materials->FindOrBuildMaterial("G4_Cu");
  G4Material* Kapton = materials->FindOrBuildMaterial("G4_KAPTON");
  G4Material* Tungsten = materials->FindOrBuildMaterial("G4_W");
  G4Material* Ni = materials->FindOrBuildMaterial("G4_Ni");
  G4Material* Carbon = materials->FindOrBuildMaterial("G4_C");

//material for toroid
  G4Material* H2Liq = new G4Material("H2Liq"    , z= 1, a= 1.01*g/mole, density= 70.8*mg/cm3);
  G4Material* Aluminium = new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);

  double density4Tor = 0.3* 70.8*mg/cm3 + 0.7*2.700*g/cm3;
  G4Material* Material4Toroid = new G4Material("Material4Toroid", density4Tor, ncomponents=2);
  Material4Toroid->AddMaterial(Aluminium, fractionmass=80.0*perCent);
  Material4Toroid->AddMaterial(H2Liq, fractionmass=20.0*perCent);

  
      // EPOXY - build up from elements
  G4Element* H  =  new G4Element("Hydrogen", symbol="H", z= 1.0, a= 1.01*g/mole);
  G4Element* C  =  new G4Element("Carbon"  , symbol="C", z= 6.0, a= 12.01*g/mole);
  G4Element* O  =  new G4Element("Oxygen"  , symbol="O", z= 8.0, a= 16.00*g/mole);
    
  G4Material* Epoxy = materials->FindOrBuildMaterial("Epoxy");
  if (!Epoxy) {
    Epoxy = new G4Material("Epoxy", density= 1.3*g/cm3, ncomponents=3);
    Epoxy->AddElement(H, fractionmass=0.1310);
    Epoxy->AddElement(C, fractionmass=0.5357);
    Epoxy->AddElement(O, fractionmass=0.3333);
  }
    
  G4double epoxydens =  Epoxy->GetDensity();
  G4double kaptondens = Kapton->GetDensity();
  G4double copperdens = Copper->GetDensity();

  // FrontFanoutMaterial
  G4double Lcal_epoxy_propF    = 1.0;
  G4double Lcal_copper_propF   = 0.5;
  G4double Lcal_epoxy_heightF  = 0.065*mm;
  G4double Lcal_kapton_heightF = 0.050*mm;
  G4double Lcal_copper_heightF = 0.035*mm;
  G4double Lcal_fanoutF_thickness = Lcal_epoxy_heightF + Lcal_kapton_heightF + Lcal_copper_heightF;
  
  G4double epoxydensF  = epoxydens  * Lcal_epoxy_propF;
  G4double copperdensF = copperdens * Lcal_copper_propF;
  G4double epoxyfracF  = Lcal_epoxy_heightF  / Lcal_fanoutF_thickness;
  G4double kaptonfracF = Lcal_kapton_heightF / Lcal_fanoutF_thickness;
  G4double copperfracF = Lcal_copper_heightF / Lcal_fanoutF_thickness;
  G4double frontDensity = (epoxydensF  * epoxyfracF + kaptondens  * kaptonfracF + copperdensF * copperfracF);

  G4Material* FanoutMatF = new G4Material("FrontFanoutMaterial", frontDensity, ncomponents=3);
  FanoutMatF->AddMaterial(Copper, fractionmass=epoxyfracF);
  FanoutMatF->AddMaterial(Kapton, fractionmass=kaptonfracF);
  FanoutMatF->AddMaterial(Epoxy,  fractionmass=copperfracF);
  
  // BackFanoutMaterial
  G4double Lcal_epoxy_propB    = 1.0;
  G4double Lcal_copper_propB   = 1.0;
  G4double Lcal_epoxy_heightB  = 0.060 *mm;
  G4double Lcal_kapton_heightB = 0.065 *mm;
  G4double Lcal_copper_heightB = 0.025 *mm;
  G4double Lcal_fanoutB_thickness = Lcal_epoxy_heightB + Lcal_kapton_heightB + Lcal_copper_heightB;
                                    
  G4double epoxydensB = epoxydens * Lcal_epoxy_propB;
  G4double copperdensB = copperdens * Lcal_copper_propB;
  G4double epoxyfracB  = Lcal_epoxy_heightB  / Lcal_fanoutB_thickness;
  G4double kaptonfracB = Lcal_kapton_heightB / Lcal_fanoutB_thickness;
  G4double copperfracB = Lcal_copper_heightB / Lcal_fanoutB_thickness;
  G4double backDensity = (2.0* epoxydensB  * epoxyfracB + kaptondens  * kaptonfracB + copperdensB * copperfracB);
  G4Material* FanoutMatB = new G4Material("BackFanoutMaterial", backDensity, ncomponents=3);
  FanoutMatB->AddMaterial(Copper, fractionmass=epoxyfracB);
  FanoutMatB->AddMaterial(Kapton, fractionmass=kaptonfracB);
  FanoutMatB->AddMaterial(Epoxy,  fractionmass=copperfracB);
  
  G4Material* C_fiber = new G4Material("ECalCarbonFiber", 1.6*g/cm3, ncomponents=2);
  C_fiber->AddMaterial(Carbon, fractionmass=50.0*perCent);
  C_fiber->AddMaterial(Epoxy, fractionmass=50.0*perCent);
  
  //Absorber MGS
  G4Material* Wabsorber_MGS = new G4Material("Wabsorber_MGS", 17.7*g/cm3, ncomponents=3);
  Wabsorber_MGS->AddMaterial(Tungsten, fractionmass=93.0*perCent);
  Wabsorber_MGS->AddMaterial(Ni, fractionmass=5.25*perCent);
  Wabsorber_MGS->AddMaterial(Copper, fractionmass=1.75*perCent);

  // for the 2014 TB plate
  G4Material* Wabsorber_PL = new G4Material("Wabsorber_PL", 18.0*g/cm3 ,ncomponents=3);
  Wabsorber_PL->AddMaterial(Tungsten, fractionmass=95.0*perCent);
  Wabsorber_PL->AddMaterial(Ni, fractionmass=2.5*perCent);
  Wabsorber_PL->AddMaterial(Copper, fractionmass=2.5*perCent);

  // fiber glass
  G4int natoms;
  G4Element* Si = new G4Element("Silcone"  ,symbol="Si" , z= 14., a= 28.09*g/mole);
  G4Material* fiberglass = new G4Material( "fiberglass",density=2.61*g/cm3,ncomponents=2);
  fiberglass->AddElement(Si, natoms=1);
  fiberglass->AddElement(O, natoms=2);
  // PCBoard material FR4
  G4Material *FR4 = new G4Material("FR4",density=1.93*g/cm3,ncomponents=3);
  FR4->AddMaterial( Epoxy, fractionmass=0.34);
  FR4->AddMaterial( fiberglass, fractionmass=0.56);
  FR4->AddMaterial( Copper, fractionmass=0.1);
        
}

