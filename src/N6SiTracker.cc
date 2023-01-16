//
/// \brief Implementation of the Alpide tracker class
//

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4AssemblyVolume.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4UserLimits.hh"
#include "N6SetUp.hh"
#include "LxDetector.hh"
#include "N6SiTracker.hh"

void N6SiTracker::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  CreateMaterial();

  N6SetUp *lxs = N6SetUp::Instance();
  G4Material* environmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* AlpidFrmConstructMaterial = G4NistManager::Instance()->FindOrBuildMaterial("CarbonFoam");
  G4Material* AlpideConstructMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Silicon");
  G4Material* TargetConstructMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Lead");

  G4Tubs *solidIPContainer = new G4Tubs("solidIPContainer",  0.0,  lxs->DipoleFieldR, lxs->DipoleFieldY/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicIPContainer = new G4LogicalVolume(solidIPContainer, environmentMaterial, "logicIPContainer");
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, 0.0, 0.0), logicIPContainer,
                     "IPContainer", fLogicWorld, false, 0, lxs->OverlapTest);

  ConstructMagnet(logicIPContainer, G4ThreeVector(0.0, lxs->DipoleFieldB, 0.0));

  //Pb targets

  G4double PbTarget_dr_0 = 0.3 *cm;
  G4double PbTarget_dr_1 = 0.1 *cm; 
  G4double PbTarget_dz = 0.15 *cm;
//RCC tgt0       0.0 0.0 -5.3 0.0 0.0 -0.15 0.3 //fluka
  G4Tubs *solidPbTarget_0 = new G4Tubs("solidPbTarget_0",  0.0,  PbTarget_dr_0, PbTarget_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicPbTarget_0 = new G4LogicalVolume(solidPbTarget_0, TargetConstructMaterial, "logicPbTarget_0");
  G4Tubs *solidPbTarget_1 = new G4Tubs("solidPbTarget_1",  0.0,  PbTarget_dr_1, PbTarget_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicPbTarget_1 = new G4LogicalVolume(solidPbTarget_1, TargetConstructMaterial, "logicPbTarget_1");
  G4Tubs *solidPbTarget_2 = new G4Tubs("solidPbTarget_2",  0.0,  PbTarget_dr_1, PbTarget_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicPbTarget_2 = new G4LogicalVolume(solidPbTarget_2, TargetConstructMaterial, "logicPbTarget_2");
  G4Tubs *solidPbTarget_3 = new G4Tubs("solidPbTarget_3",  0.0,  PbTarget_dr_1, PbTarget_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicPbTarget_3 = new G4LogicalVolume(solidPbTarget_3, TargetConstructMaterial, "logicPbTarget_3");
  G4Tubs *solidPbTarget_4 = new G4Tubs("solidPbTarget_4",  0.0,  PbTarget_dr_1, PbTarget_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicPbTarget_4 = new G4LogicalVolume(solidPbTarget_4, TargetConstructMaterial, "logicPbTarget_4");

//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, PbTarget_zpos_0), logicPbTarget_0,
//                      "Target0", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, PbTarget_zpos_1), logicPbTarget_1,
//                      "Target1", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, PbTarget_zpos_2), logicPbTarget_2,
//                      "Target2", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, PbTarget_zpos_3), logicPbTarget_3,
//                      "Target3", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, PbTarget_zpos_4), logicPbTarget_4,
//                      "Target4", fLogicWorld, false, 0, lxs->OverlapTest);

  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, lxs->PbTarget_zpos_0, 0.0), logicPbTarget_0,
                     "Target0", logicIPContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, lxs->PbTarget_zpos_1, 0.0), logicPbTarget_1,
                     "Target1", logicIPContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, lxs->PbTarget_zpos_2, 0.0), logicPbTarget_2,
                     "Target2", logicIPContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, lxs->PbTarget_zpos_3, 0.0), logicPbTarget_3,
                     "Target3", logicIPContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, lxs->PbTarget_zpos_4, 0.0), logicPbTarget_4,
                     "Target4", logicIPContainer, false, 0, lxs->OverlapTest);

// Veretex spectrometer Alpide Trackers

//Frame 
  G4double AlpidFrm_dx = 30.0 *cm;
  G4double AlpidFrm_dy = 30.0 *cm;
  G4double AlpidFrm_dz = 0.50 *cm;
  G4double AlpidFrmHoleR = 0.424 *cm;

  G4double AlpidTrkCut_dx = 12.5 *cm;
  G4double AlpidTrkCut_dy = 13.0 *cm;
  G4double dxycut = 1.0 *cm;

  G4double AlpidTrk_dx = 14.69 *cm;
  G4double AlpidTrk_dy = 14.69 *cm;
  G4double AlpidTrk_dz = 50e-4 *cm;
  G4double AlpiDy = 0.29 *cm;
  G4double AlpiDx = 0.31 *cm;

//   G4double AlpidTrk_zpos_0 = 7.1175 *cm;
//   G4double AlpidTrk_zpos_1 = 15.1175 *cm;
//   G4double AlpidTrk_zpos_2 = 20.1175 *cm;
//   G4double AlpidTrk_zpos_3 = 25.1175 *cm;
//   G4double AlpidTrk_zpos_4 = 38.1175 *cm;

  std::vector<G4double> AlpidTrk_zpos{7.1175 *cm, 15.1175 *cm, 20.1175 *cm, 25.1175 *cm, 38.1175 *cm};
// Pixel supporting frame
  G4Box *solidAlpidFrm = new G4Box("solidAlpidFrm", AlpidFrm_dx/2.0, AlpidFrm_dy/2.0, AlpidFrm_dz/2.0);
  // hole for the beampipe
  G4Tubs *solidAlpidFrm4hole = new G4Tubs("solidAlpidFrm4hole",  0.0,  AlpidFrmHoleR, AlpidFrm_dz, 0.0, 2.0*M_PI);
  G4SubtractionSolid* solidAlpidFrm_hole = new G4SubtractionSolid("solidAlpidFrm_hole", solidAlpidFrm, solidAlpidFrm4hole);
  // holes for alpides
  G4Box *solidAlpidHole = new G4Box("solidAlpidHole", AlpidTrkCut_dx/2.0, AlpidTrkCut_dy/2.0, AlpidFrm_dz);

  std::vector<G4double> shiftx{AlpidTrkCut_dx/2.0+dxycut,  -AlpidTrkCut_dx/2.0, -AlpidTrkCut_dx/2.0-dxycut, AlpidTrkCut_dx/2.0};
  std::vector<G4double> shifty{AlpidTrkCut_dy/2.0, AlpidTrkCut_dy/2.0+dxycut, -AlpidTrkCut_dy/2.0,  -AlpidTrkCut_dy/2.0-dxycut};
  G4SubtractionSolid* solidAlpidFrm_SqrHole3 = solidAlpidFrm_hole;

  for (size_t ii = 0; ii < shiftx.size(); ++ii) {
    G4Transform3D transform0(G4RotationMatrix(), G4ThreeVector(shiftx[ii], shifty[ii], 0.0));
    G4String vname = G4String("solidAlpidFrm_SqrHole") + std::to_string(ii);
    solidAlpidFrm_SqrHole3 = new G4SubtractionSolid(vname, solidAlpidFrm_SqrHole3, solidAlpidHole, transform0);
  }
  G4LogicalVolume *logicAlpidFrm = new G4LogicalVolume(solidAlpidFrm_SqrHole3, AlpidFrmConstructMaterial, "logicAlpidFrm");

// Silicon sensors in the container
  G4Box *solidAlpidTrk = new G4Box("solidAlpidTrk", AlpidTrk_dx/2.0, AlpidTrk_dy/2.0, AlpidTrk_dz/2.0);
  G4LogicalVolume *logicAlpidTrk = new G4LogicalVolume(solidAlpidTrk, AlpideConstructMaterial, "logicAlpidTrk");

  std::vector<G4double> alpdx{AlpidTrk_dx/2.0+AlpiDx,  -AlpidTrk_dx/2.0+AlpiDy, -AlpidTrk_dx/2.0-AlpiDx, AlpidTrk_dx/2.0-AlpiDy};
  std::vector<G4double> alpdy{AlpidTrk_dy/2.0-AlpiDy, AlpidTrk_dy/2.0+AlpiDx, -AlpidTrk_dy/2.0+AlpiDy,  -AlpidTrk_dy/2.0-AlpiDx};

  G4Box *solidAlpidContainer = new G4Box("solidAlpidContainer", AlpidFrm_dx/2.0, AlpidFrm_dy/2.0, AlpidTrk_dz/2.0);
  for (size_t ll = 0; ll < AlpidTrk_zpos.size(); ++ll) {
    G4String stname = G4String("logicPixStnContainer") + std::to_string(ll);
    G4LogicalVolume *logicAlpidContainer = new G4LogicalVolume(solidAlpidContainer, environmentMaterial, stname);
    for (size_t ii = 0; ii < alpdx.size(); ++ii) {
      G4String vname = G4String("PixStn") + std::to_string(ll);
      new G4PVPlacement(0, G4ThreeVector(alpdx[ii], alpdy[ii], 0.0), logicAlpidTrk,
                        vname, logicAlpidContainer, false, ii, lxs->OverlapTest);
    }
    // placement of frames to the IP volume
    G4String stfvname = G4String("PixStnFrame") + std::to_string(ll);
    new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, AlpidTrk_zpos[ll] +  0.5 * (AlpidFrm_dz + AlpidTrk_dz), 0.0),
                      logicAlpidFrm, stfvname, logicIPContainer, false, ll, lxs->OverlapTest);
    // placement of sensor containers to the IP volume
    G4String stpixvname = G4String("PixStnLayer") + std::to_string(ll);
    new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0),
                      G4ThreeVector(0.0, AlpidTrk_zpos[ll], 0.0),
                      logicAlpidContainer, stpixvname, logicIPContainer, false, ll, lxs->OverlapTest);
  }


/* //4*4
  G4int LGidx=0;
  G4SubtractionSolid* solidAlpidFrm_SqrHole;
  G4Box *solidAlpidTrk = new G4Box("solidAlpidTrk", AlpidTrk_dx/2.0, AlpidTrk_dy/2.0, 
                                     AlpidTrk_dz/2.0);
  G4LogicalVolume *logicAlpidTrk = new G4LogicalVolume(solidAlpidTrk, AlpideConstructMaterial, "logicAlpidTrk");
  G4double xy[4] = {(shift - lxs->BPipeR + AlpidTrk_dx/2.0),(shift + lxs->BPipeR  + AlpidTrk_dy/2.0),( -shift - AlpidTrk_dx/2.0 +lxs->BPipeR),(-shift -AlpidTrk_dx/2.0-lxs->BPipeR)};                   
//1
  G4Transform3D transform0(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0), 
                                   G4ThreeVector(xy[0], xy[1], 0.0));
  G4SubtractionSolid* solidAlpidFrm_SqrHole0 = new G4SubtractionSolid("solidAlpidFrm_SqrHole0", solidAlpidFrm_hole, solidAlpidHole, transform0);
 

//2
  G4Transform3D transform1(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0), 
                                   G4ThreeVector(xy[1], xy[2], 0.0));
  G4SubtractionSolid* solidAlpidFrm_SqrHole1 = new G4SubtractionSolid("solidAlpidFrm_SqrHole1", solidAlpidFrm_SqrHole0, solidAlpidHole, transform1);

//3  
  G4Transform3D transform2(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0), 
                                   G4ThreeVector(xy[2], xy[3], 0.0));
  G4SubtractionSolid* solidAlpidFrm_SqrHole2 = new G4SubtractionSolid("solidAlpidFrm_SqrHole2", solidAlpidFrm_SqrHole1, solidAlpidHole, transform2);

//4  
  G4Transform3D transform3(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0), 
                                   G4ThreeVector(xy[3], xy[0], 0.0));
  G4SubtractionSolid* solidAlpidFrm_SqrHole3 = new G4SubtractionSolid("solidAlpidFrm_SqrHole3", solidAlpidFrm_SqrHole2, solidAlpidHole, transform3);

*/
//   G4LogicalVolume *logicAlpidFrm = new G4LogicalVolume(solidAlpidFrm_SqrHole3, AlpidFrmConstructMaterial, "logicAlpidFrm");

//   new G4PVPlacement(transform0, logicAlpidTrk ,
//                      "Alpide0", logicAlpidFrm, false, 0, lxs->OverlapTest);

  
 
  
/*
  for (G4int il = 0; il < 2; ++il) {
    for (G4int jl = 0; jl < 2; ++jl) {
      std::cout << "alpide index: " << il << jl<< std::endl;
 if(!( (il==1&&jl==2)||(il==2&&jl==2)|| (il==2&&jl==1)||(il==1&&jl==1) )){ 
  LGidx=LGidx+1;                                 
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0), 
                    G4ThreeVector(xy[il], xy[jl], 0.0), 
                    logicAlpidTrk, "Alpids", logicAlpidFrm, false, LGidx, lxs->OverlapTest);                   
                    }                               
  			}           
  } */           
 

 
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AlpidTrk_zpos_0), logicAlpidTrk,
//                      "PixStn0", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AlpidTrk_zpos_1), logicAlpidTrk,
//                      "PixStn1", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AlpidTrk_zpos_2), logicAlpidTrk,
//                      "PixStn2", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AlpidTrk_zpos_3), logicAlpidTrk,
//                      "PixStn3", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AlpidTrk_zpos_4), logicAlpidTrk,
//                      "PixStn4", fLogicWorld, false, 0, lxs->OverlapTest);

//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, AlpidTrk_zpos_0, 0.0), logicAlpidFrm,
//                      "PixStn0", logicIPContainer, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, AlpidTrk_zpos_1, 0.0), logicAlpidFrm,
//                      "PixStn1", logicIPContainer, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, AlpidTrk_zpos_2, 0.0), logicAlpidFrm,
//                      "PixStn2", logicIPContainer, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, AlpidTrk_zpos_3, 0.0), logicAlpidFrm,
//                      "PixStn3", logicIPContainer, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, AlpidTrk_zpos_4, 0.0), logicAlpidFrm,
//                      "PixStn4", logicIPContainer, false, 0, lxs->OverlapTest);

//   AddSegmentation();
}


void N6SiTracker::AddSegmentation()
{
  N6SetUp *lxs = N6SetUp::Instance();
  G4double bsmcalox = lxs->BSMCaloX;
  G4double bsmcaloy = lxs->BSMCaloY;
  G4int ncellx = lxs->BSMCaloNCellX;
  G4int ncelly = lxs->BSMCaloNCellY;

  fDetector->AddSensorSegmentation("BSMCaloLayer", bsmcalox, bsmcaloy, ncellx, ncelly);
}



G4AssemblyVolume* N6SiTracker::ConstructSupportAssembly()
{
//   N6SetUp *lxs = N6SetUp::Instance();
  G4AssemblyVolume *supportAssembly = new G4AssemblyVolume();
  return supportAssembly;
}



void N6SiTracker::CreateMaterial()
{
G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  
// 
  G4int ncomponents, nel;
  G4double fractionmass;

  G4NistManager *materials = G4NistManager::Instance();
  
  G4Material* Copper = materials->FindOrBuildMaterial("G4_Cu");
  G4Material* Kapton = materials->FindOrBuildMaterial("G4_KAPTON");
  G4Material* Tungsten = materials->FindOrBuildMaterial("G4_W");
  G4Material* Ni = materials->FindOrBuildMaterial("G4_Ni");
  G4Material* Carbon = materials->FindOrBuildMaterial("G4_C");
//  G4Material* Lead = new G4Material("Lead",     z=82, a=207.19*g/mole, density= 11.35*g/cm3);
  
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
  G4Element* Si = new G4Element("Silicon"  ,symbol="Si" , z= 14., a= 28.09*g/mole);
  G4Material* fiberglass = new G4Material( "fiberglass",density=2.61*g/cm3,ncomponents=2);
  fiberglass->AddElement(Si, natoms=1);
  fiberglass->AddElement(O, natoms=2);
  // PCBoard material FR4
  G4Material *FR4 = new G4Material("FR4",density=1.93*g/cm3,ncomponents=3);
  FR4->AddMaterial( Epoxy, fractionmass=0.34);
  FR4->AddMaterial( fiberglass, fractionmass=0.56);
  FR4->AddMaterial( Copper, fractionmass=0.1);
  
  // Carbon Fiber
  G4Material* CarbonFiber = new G4Material("CarbonFiber",density = 0.145*g/cm3, nel=1);
  CarbonFiber->AddElement(C,1);
   // Carbon Foam Pocofoam® graphite foam
  G4Material* CarbonFoam = new G4Material("CarbonFoam",density = 0.5*g/cm3, nel=1);
  CarbonFoam->AddElement(C,1);
        
}




void N6SiTracker::ConstructMagnet(G4LogicalVolume *mfContainer, const G4ThreeVector mf)
{
//   G4double magnetFieldValue = 14000.0*gauss;
//   G4double magnetFieldValue = 0.0*gauss;

//   G4ThreeVector  fieldVector( 0.0, -magnetFieldValue, 0.0);
  G4ThreeVector  fieldVector(mf);
  G4MagneticField *magField = new G4UniformMagField( fieldVector );
  G4FieldManager  *localFieldMgr = new G4FieldManager (magField);
  localFieldMgr->CreateChordFinder(magField);

  mfContainer->SetFieldManager(localFieldMgr, true);
}







