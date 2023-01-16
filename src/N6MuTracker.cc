//
/// \brief Implementation of the BeamProfiler tracker class
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

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4UserLimits.hh"
#include "N6SetUp.hh"
#include "LxDetector.hh"
#include "N6MuTracker.hh"



void N6MuTracker::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  CreateMaterial();

  N6SetUp *lxs = N6SetUp::Instance();
  G4Material* environmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* chamberConstructMaterial1 = G4NistManager::Instance()->FindOrBuildMaterial("FR4");
  G4Material* chamberConstructMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Honeycomb");  
  G4Material* chamberGasMaterial = G4NistManager::Instance()->FindOrBuildMaterial("ArgonCO2");
  G4Material* DummychamberConstructMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Silicon");

//dummy muon chambers
  G4double MuChamber_drmin = 0.0; //lxs->BPipeR;
  G4double MuChamber_drmax = 350 *cm;
  
  G4double MuChamber_0_drmin = MuChamber_drmin; //6.5329 *cm;
  G4double MuChamber_0_drmax = MuChamber_drmax; //90.0 *cm;
  G4double MuChamber_1_drmin = MuChamber_drmin; //30.000 *cm;
  G4double MuChamber_1_drmax = MuChamber_drmax; //110.0 *cm;
  G4double MuChamber_2_drmin = MuChamber_drmin; //30.0000 *cm;
  G4double MuChamber_2_drmax = MuChamber_drmax; //220.0 *cm;
  G4double MuChamber_3_drmin = MuChamber_drmin; //16.6656 *cm;
  G4double MuChamber_3_drmax = MuChamber_drmax; //240.0 *cm;
  G4double MuChamber_4_drmin = MuChamber_drmin; //26.6649 *cm;
  G4double MuChamber_4_drmax = MuChamber_drmax; //311.0 *cm;  
  G4double MuChamber_5_drmin = MuChamber_drmin; //27.7759 *cm;
  G4double MuChamber_5_drmax = MuChamber_drmax; //311.0 *cm;    
  G4double MuChamber_0_dz = 610e-4 *cm;
  G4double MuChamber_2_dz = 700e-4 *cm;
  G4double MuChamber_4_dz = 1330e-4 *cm;
 
  G4Tubs *solidDummyMuChamber_0 = new G4Tubs("solidDummyMuChamber_0",  MuChamber_0_drmin,  MuChamber_0_drmax, MuChamber_0_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicDummyMuChamber_0 = new G4LogicalVolume(solidDummyMuChamber_0, DummychamberConstructMaterial, "logicDummyMuChamber_0");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ms0_pos), logicDummyMuChamber_0,
                     "MS0", fLogicWorld, false, 0, lxs->OverlapTest);
  G4Tubs *solidDummyMuChamber_1 = new G4Tubs("solidDummyMuChamber_1",  MuChamber_1_drmin,  MuChamber_1_drmax, MuChamber_0_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicDummyMuChamber_1 = new G4LogicalVolume(solidDummyMuChamber_1, DummychamberConstructMaterial, "logicDummyMuChamber_1");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ms1_pos), logicDummyMuChamber_1,
                     "MS1", fLogicWorld, false, 0, lxs->OverlapTest);
  G4Tubs *solidDummyMuChamber_2 = new G4Tubs("solidDummyMuChamber_2",  MuChamber_2_drmin,  MuChamber_2_drmax, MuChamber_2_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicDummyMuChamber_2 = new G4LogicalVolume(solidDummyMuChamber_2, DummychamberConstructMaterial, "logicDummyMuChamber_2");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ms2_pos), logicDummyMuChamber_2,
                     "MS2", fLogicWorld, false, 0, lxs->OverlapTest);
  G4Tubs *solidDummyMuChamber_3 = new G4Tubs("solidDummyMuChamber_3",  MuChamber_3_drmin,  MuChamber_3_drmax, MuChamber_2_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicDummyMuChamber_3 = new G4LogicalVolume(solidDummyMuChamber_3, DummychamberConstructMaterial, "logicDummyMuChamber_3");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ms3_pos), logicDummyMuChamber_3,
                     "MS3", fLogicWorld, false, 0, lxs->OverlapTest); 
  G4Tubs *solidDummyMuChamber_4 = new G4Tubs("solidDummyMuChamber_4",  MuChamber_4_drmin,  MuChamber_4_drmax, MuChamber_4_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicDummyMuChamber_4 = new G4LogicalVolume(solidDummyMuChamber_4, DummychamberConstructMaterial, "logicDummyMuChamber_4");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ms4_pos), logicDummyMuChamber_4,
                     "TrigStn0", fLogicWorld, false, 0, lxs->OverlapTest);
  G4Tubs *solidDummyMuChamber_5 = new G4Tubs("solidDummyMuChamber_5",  MuChamber_5_drmin,  MuChamber_5_drmax, MuChamber_4_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicDummyMuChamber_5 = new G4LogicalVolume(solidDummyMuChamber_5, DummychamberConstructMaterial, "logicDummyMuChamber_5");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ms5_pos), logicDummyMuChamber_5,
                     "TrigStn1", fLogicWorld, false, 0, lxs->OverlapTest);                  
                                                          

//real muon chambers
  G4double inchamber_h = 900.0 *mm;
  G4double inchamber_ldx = 300.0 *mm;
  G4double inchamber_hdx = 700.0 *mm;
//   G4double hcombthickness = 17.0 *mm;
  G4double sidethickness = 7.0 *mm;
  G4double gasthickness = 6.0 *mm;
  G4double rocard_x = 50.0 *mm;

  G4double a = atan2(0.5 * (inchamber_hdx-inchamber_ldx), inchamber_h);
  G4double outchamber_h = inchamber_h + 2.0 * sidethickness;
  G4double outchamber_ldx = inchamber_ldx + 2.0 * (sidethickness/cos(a) - sidethickness * tan(a));
  G4double outchamber_hdx = inchamber_hdx + 2.0 * (sidethickness/cos(a) + sidethickness * tan(a));
  G4double container_ldx = outchamber_ldx + 2.0 * rocard_x/cos(a);
  G4double container_hdx = outchamber_hdx + 2.0 * rocard_x/cos(a);
  G4double outchamber_y = 20.0 *mm;

  G4double r0_pos = 244.0 *mm;

  G4Trd *solidChamberContainer = new G4Trd("solidChamberContainer", container_ldx/2.0, container_hdx/2.0,
                                           outchamber_y/2.0, outchamber_y/2.0, outchamber_h/2.0);
  G4LogicalVolume *logicChamberContainer = new G4LogicalVolume(solidChamberContainer, environmentMaterial,
                                                               "logicChamberContainer");

  G4Trd *solidChamber = new G4Trd("solidChamber", outchamber_ldx/2.0, outchamber_hdx/2.0,
                                  outchamber_y/2.0, outchamber_y/2.0, outchamber_h/2.0);
  G4LogicalVolume *logicChamber = new G4LogicalVolume(solidChamber, chamberConstructMaterial,
                                                               "logicChamber");

  G4Trd *solidChamberGas = new G4Trd("solidChamberGas", inchamber_ldx/2.0, inchamber_hdx/2.0,
                                     gasthickness/2.0, gasthickness/2.0, inchamber_h/2.0);
  G4LogicalVolume *logicChamberGas = new G4LogicalVolume(solidChamberGas, chamberGasMaterial,
                                                               "logicChamberGas");

  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicChamberGas, "ChamberGas", logicChamber, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicChamber, "Chamber", logicChamberContainer, false, 0, lxs->OverlapTest);

  std::vector<G4int> muwheelid{0, 100, 200, 300, 400, 500};
  G4int nchamber0 = 12;
  G4double psi = 11.0*M_PI/180.0;
  G4double zpos_shift = 15.0 *cm;

  for (G4int il = 0; il < nchamber0; ++il) {
    G4double phi = il * 2.0*M_PI/static_cast<G4double>(nchamber0);
    G4double ypos = (r0_pos + outchamber_h/2.0) * cos(phi);
    G4double xpos = (r0_pos + outchamber_h/2.0) * sin(phi);
    G4RotationMatrix *rotlhex = new G4RotationMatrix();
    rotlhex->rotateX(M_PI/2.0);
    rotlhex->rotateY(phi);
    rotlhex->rotateZ(-psi);

    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms0_pos + zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il + muwheelid[0], lxs->OverlapTest);
    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms1_pos - zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il + muwheelid[1], lxs->OverlapTest);
    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms2_pos + zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il + muwheelid[2], lxs->OverlapTest);
    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms3_pos - zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il + muwheelid[3], lxs->OverlapTest);
    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms4_pos + zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il + muwheelid[4], lxs->OverlapTest);
    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms5_pos + zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il + muwheelid[5], lxs->OverlapTest);                   
  }

  G4int nchamber1 = 24;

  for (G4int il = 0; il < nchamber1; ++il) {
    G4double phi = (il-0.5) * 2.0*M_PI/static_cast<G4double>(nchamber1);
    G4double ypos = (r0_pos + 0.5*outchamber_h + inchamber_h) * cos(phi);
    G4double xpos = (r0_pos + 0.5*outchamber_h + inchamber_h) * sin(phi);
    G4RotationMatrix *rotlhex = new G4RotationMatrix();
    rotlhex->rotateX(M_PI/2.0);
    rotlhex->rotateY(phi);
    rotlhex->rotateZ(-psi);

    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms2_pos + zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il+nchamber0 + muwheelid[2], lxs->OverlapTest);
    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms3_pos - zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il+nchamber0 + muwheelid[3], lxs->OverlapTest);
    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms4_pos + zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il+nchamber0 + muwheelid[4], lxs->OverlapTest);
    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms5_pos + zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il+nchamber0 + muwheelid[5], lxs->OverlapTest);                 
  }


  G4int nchamber2 = 48;
  for (G4int il = 0; il < nchamber2; ++il) {
    G4double ddphi = il&1 ? 0.5 : 0.5;
    G4double phi = (il - ddphi) * 2.0*M_PI/static_cast<G4double>(nchamber2);
    G4double ypos = (r0_pos + 0.5*outchamber_h + 2.0*inchamber_h) * cos(phi);
    G4double xpos = (r0_pos + 0.5*outchamber_h + 2.0*inchamber_h) * sin(phi);
    G4RotationMatrix *rotlhex = new G4RotationMatrix();
    rotlhex->rotateX(M_PI/2.0);
    rotlhex->rotateY(phi);
    rotlhex->rotateZ(-psi);
     
 
    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms4_pos + zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il+nchamber0+nchamber1 + muwheelid[4], lxs->OverlapTest);
    new G4PVPlacement (rotlhex, G4ThreeVector(-xpos, ypos, lxs->ms5_pos + zpos_shift), logicChamberContainer, "ChamberContainer",
                     fLogicWorld, false, il+nchamber0+nchamber1 + muwheelid[5], lxs->OverlapTest);                                        
  }

/// Absorbers    

  G4double AbsorberBeO_1_dy = 520.0 *mm;
  G4double AbsorberBeO_2_dy = 1200.0 *mm;
  G4double AbsorberBeO_1_dx = 520.0 *mm; 
  G4double AbsorberBeO_2_dx = 1200.0 *mm; 
  G4double AbsorberBeO_1_dr = 24.7 *mm;
  G4double AbsorberBeO_2_dr = 54.8 *mm;
  G4double AbsorberGr_1_dr = 110.0 *mm;
  G4Material* AbsorberBeOMaterial = G4NistManager::Instance()->FindOrBuildMaterial("BeO");// (lxs->EnvironmentMaterial);//
  G4Material* AbsorberGraphiteMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Graphite"); //(lxs->EnvironmentMaterial);//
  G4Material* AbsorberPlugWMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Tungsten"); //(lxs->EnvironmentMaterial);//

  G4Box *solidAbsorberBeO_1 = new G4Box("solidAbsorberBeO_1", AbsorberBeO_1_dx/2.0, AbsorberBeO_1_dy/2.0, 
                                     lxs->AbsorberBeO_1_dz/2.0);
  G4Box *solidAbsorberBeO_2 = new G4Box("solidAbsorberBeO_2", AbsorberBeO_2_dx/2.0, AbsorberBeO_2_dy/2.0, 
                                     lxs->AbsorberBeO_2_dz/2.0);
  G4Transform3D transformh(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, 0.0));                     
                                     
  G4Tubs *solidAbsorberBeOH_1 = new G4Tubs("solidAbsorberBeOH1", 0.0, AbsorberBeO_1_dr, lxs->AbsorberBeO_1_dz, 0.0, 2.0*M_PI);
  G4Tubs *solidAbsorberBeOH_2 = new G4Tubs("solidAbsorberBeOH2", 0.0, AbsorberBeO_2_dr, lxs->AbsorberBeO_2_dz, 0.0, 2.0*M_PI);
  
  G4SubtractionSolid* solidAbsorberBeO_1S = new G4SubtractionSolid("solidAbsorberBeO_1S", solidAbsorberBeO_1,
                                                               solidAbsorberBeOH_1, transformh); 
  G4SubtractionSolid* solidAbsorberBeO_2S = new G4SubtractionSolid("solidAbsorberBeO_2S", solidAbsorberBeO_2,
                                                               solidAbsorberBeOH_2, transformh);                                                                           
  
  G4LogicalVolume *logicAbsorberBeOH_1 = new G4LogicalVolume(solidAbsorberBeO_1S, AbsorberBeOMaterial, "logicAbsorberBeOH_1");
  G4LogicalVolume *logicAbsorberBeOH_2 = new G4LogicalVolume(solidAbsorberBeO_2S, AbsorberBeOMaterial, "logicAbsorberBeOH_2");
 
 // BeO Absorber Plugs
  G4Tubs *solidAbsorberPlugW_1 = new G4Tubs("solidAbsorberPlugW_1", 0.0, AbsorberBeO_1_dr, lxs->AbsorberBeO_1_dz/2.0, 0.0, 2.0*M_PI);
  G4Tubs *solidAbsorberPlugW_2 = new G4Tubs("solidAbsorberPlugW_2", 0.0, AbsorberBeO_2_dr, lxs->AbsorberBeO_2_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicAbsorberPlugW_1 = new G4LogicalVolume(solidAbsorberPlugW_1, AbsorberPlugWMaterial, "logicAbsorberPlugW_1");
  G4LogicalVolume *logicAbsorberPlugW_2 = new G4LogicalVolume(solidAbsorberPlugW_2, AbsorberPlugWMaterial, "logicAbsorberPlugW_2");
  
 
   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->AbsorberBeOH_1_pos), logicAbsorberBeOH_1,
                     "AbsoBeO_1", fLogicWorld, false, 0, lxs->OverlapTest); 
   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->AbsorberBeOH_2_pos), logicAbsorberBeOH_2,
                     "AbsoBeO_2", fLogicWorld, false, 0, lxs->OverlapTest); 
   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->AbsorberBeOH_1_pos), logicAbsorberPlugW_1,
                     "AbsoPlug1", fLogicWorld, false, 0, lxs->OverlapTest); 
   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->AbsorberBeOH_2_pos), logicAbsorberPlugW_2,
                     "AbsoPlug2", fLogicWorld, false, 0, lxs->OverlapTest);                  

 //Graphite Absorber 1
  G4double AbsorberC_1_dz = 1302.0 *mm; //1300.0 *mm;
  G4double AbsorberC_1_dx = 2600.0 *mm;
  G4double AbsorberC_1_dy = 2600.0 *mm;
 // to make   Graphite Absorber 1 sensitive
 //6.7 X0 so make 7
 G4double AbsorberC_1_dz_1X0 = 186.0 *mm;
 // Graphite Absorber 1 Container
  G4Box *solidAbsorberC_1_C = new G4Box("solidAbsorberC_1_C", AbsorberC_1_dx/2.0, AbsorberC_1_dy/2.0, 
                                     AbsorberC_1_dz/2.0); 
  G4Box *solidAbsorberC_1_1X0 = new G4Box("solidAbsorberC_1_1X0", AbsorberC_1_dx/2.0, AbsorberC_1_dy/2.0, 
                                     AbsorberC_1_dz_1X0/2.0);
  G4Tubs *solidAbsorberPlugW_3_1X0 = new G4Tubs("solidAbsorberPlugW_3_1X0", 0.0, AbsorberGr_1_dr, AbsorberC_1_dz_1X0/2.0, 0.0, 2.0*M_PI);                                   
  G4LogicalVolume *logicAbsorberPlugW_3_1X0 = new G4LogicalVolume(solidAbsorberPlugW_3_1X0, AbsorberPlugWMaterial, "logicAbsorberPlugW_3_1X0");                                    
  G4Tubs *solidAbsorberC_3_1X0 = new G4Tubs("solidAbsorberC3_1X0", 0.0, AbsorberGr_1_dr, AbsorberC_1_dz_1X0, 0.0, 2.0*M_PI);
  G4SubtractionSolid* solidAbsorberC_1X0_S = new G4SubtractionSolid("solidAbsorberC_1X0_S", solidAbsorberC_1_1X0,
                                                               solidAbsorberC_3_1X0, transformh);                                     
  G4LogicalVolume *logicAbsorberC_1_1X0 = new G4LogicalVolume(solidAbsorberC_1X0_S, AbsorberGraphiteMaterial, "logicAbsorberC_1_1X0"); 
  G4LogicalVolume *logicAbsorberC_C = new G4LogicalVolume(solidAbsorberC_1_C, environmentMaterial, "logicAbsorberC_C");

  for (int ii = 0; ii < 7; ++ii) {
//for (int ii = -3; ii < 4; ++ii) {
//      G4String aname = G4String("AbsorberGraphite_X0_") + std::to_string(ii+3);
//      new G4PVPlacement(0, G4ThreeVector(alpdx[ii], alpdy[ii], 0.0), logicAlpidTrk,
//                        vname, logicAlpidContainer, false, ii, lxs->OverlapTest);
//(2.0*ii+1)*AbsorberC_1_dz_1X0/2.0 - AbsorberC_1_dz_1X0/2.0      
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (ii - 5.0/2.0)*AbsorberC_1_dz_1X0 - AbsorberC_1_dz_1X0/2.0), logicAbsorberPlugW_3_1X0,
                     "WPlugGraphite_1X0", logicAbsorberC_C, false, ii, lxs->OverlapTest); 
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (ii - 5.0/2.0)*AbsorberC_1_dz_1X0 - AbsorberC_1_dz_1X0/2.0), logicAbsorberC_1_1X0,
                     "AbsorberGraphite_1X0", logicAbsorberC_C, false, ii, lxs->OverlapTest);
    }                                      
// Graphite Absorber 1 Container placement                      
  G4double AbsorberC_1_pos = lxs->AbsorberBeOH_2_pos + lxs->AbsorberBeO_2_dz/2.0 + AbsorberC_1_dz/2.0;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AbsorberC_1_pos), logicAbsorberC_C,
                     "AbsorberGraphite_Container", fLogicWorld, false, 0, lxs->OverlapTest);  
//  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AbsorberC_1_pos), logicAbsorberPlugW_3,
//                     "AbsoPlug3", fLogicWorld, false, 0, lxs->OverlapTest);

                                                                  
/* //Graphite Absorber Plug
  G4Tubs *solidAbsorberPlugW_3 = new G4Tubs("solidAbsorberPlugW_3", 0.0, AbsorberGr_1_dr, AbsorberC_1_dz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicAbsorberPlugW_3 = new G4LogicalVolume(solidAbsorberPlugW_3, AbsorberPlugWMaterial, "logicAbsorberPlugW_3");

  G4Tubs *solidAbsorberC_3 = new G4Tubs("solidAbsorberC3", 0.0, AbsorberGr_1_dr, AbsorberC_1_dz, 0.0, 2.0*M_PI);

  G4Box *solidAbsorberC_1 = new G4Box("solidAbsorberC_1", AbsorberC_1_dx/2.0, AbsorberC_1_dy/2.0, 
                                     AbsorberC_1_dz/2.0);
  G4SubtractionSolid* solidAbsorberC_S = new G4SubtractionSolid("solidAbsorberC_S", solidAbsorberC_1,
                                                               solidAbsorberC_3, transformh);
  G4LogicalVolume *logicAbsorberC_1 = new G4LogicalVolume(solidAbsorberC_S, AbsorberGraphiteMaterial, "logicAbsorberC_1");
  G4double AbsorberC_1_pos = AbsorberBeOH_2_pos + AbsorberBeO_2_dz/2.0 + AbsorberC_1_dz/2.0;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AbsorberC_1_pos), logicAbsorberC_1,
                     "AbsorberGraphite_1", fLogicWorld, false, 0, lxs->OverlapTest);  
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AbsorberC_1_pos), logicAbsorberPlugW_3,
                     "AbsoPlug3", fLogicWorld, false, 0, lxs->OverlapTest);*/

 //Graphite Absorber Wall
  G4double AbsorberC_Wall_dz = 1800.0 *mm;
  G4double AbsorberC_Wall_dr =3000.0 *mm;
  G4Tubs *solidAbsorberC_Wall = new G4Tubs("solidAbsorberC_Wall", 0.0, AbsorberC_Wall_dr, AbsorberC_Wall_dz/2.0, 0.0, 2.0*M_PI);
  
  G4LogicalVolume *logicAbsorberC_Wall = new G4LogicalVolume(solidAbsorberC_Wall, AbsorberGraphiteMaterial, "logicAbsorberC_Wall");

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->AbsorberC_Wall_pos), logicAbsorberC_Wall,
                     "AbsorberGraphite_Wall", fLogicWorld, false, 0, lxs->OverlapTest);
//   AddSegmentation();
}


void N6MuTracker::AddSegmentation()
{

}



G4AssemblyVolume* N6MuTracker::ConstructSupportAssembly()
{
//   N6SetUp *lxs = N6SetUp::Instance();
  G4AssemblyVolume *supportAssembly = new G4AssemblyVolume();
  return supportAssembly;
}



void N6MuTracker::CreateMaterial()
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
  G4Material *FR4_b = new G4Material("FR4_b",density=1.93*g/cm3,ncomponents=3);
  FR4_b->AddMaterial( Epoxy, fractionmass=0.34);
  FR4_b->AddMaterial( fiberglass, fractionmass=0.56);
  FR4_b->AddMaterial( Copper, fractionmass=0.1);
  
  //FR4 (Glass + Epoxy)
  G4Material* FR4 = new G4Material("FR4"  , density = 1.86*g/cm3, ncomponents=2);
  FR4->AddMaterial(fiberglass, fractionmass=0.528);
  FR4->AddMaterial(Epoxy, fractionmass=0.472);
  

}

