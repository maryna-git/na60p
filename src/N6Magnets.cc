//
/// \brief Implementation of the DetectorConstruction class
//

#include <algorithm>
#include <functional>
#include <fstream>

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4Paraboloid.hh"
#include "G4GenericTrap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "G4GeometryManager.hh"
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
#include "N6Magnets.hh"

/////////////////////////////////////////////////////////////////////////////////////
// MagnetAssembly


MEP48MagnetAssembly::MEP48MagnetAssembly(const G4String mtypename) :
  fMagnetType(mtypename)
{}



G4LogicalVolume* MEP48MagnetAssembly::GetFieldVolume(const G4String lvname)
{
  N6SetUp *lxs = N6SetUp::Instance();
  G4Material* EnvironmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Box *solidFlashMFieldVol = new G4Box(lvname + G4String("solid"), 1.0, 1.0, 1.0);
  G4LogicalVolume *logicFlashMFieldVol = new G4LogicalVolume(solidFlashMFieldVol, EnvironmentMaterial, lvname);
  return logicFlashMFieldVol;
}



void MEP48MagnetAssembly::ConstructMagnet()
{
  N6SetUp *lxs = N6SetUp::Instance();

  G4Material* magMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->MagnetMaterial);
  G4Material* wireMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TypMBWireMaterial);

  G4double IronPoleR = 0.5 *m;
  G4double IronPoleR2 = 0.6 *m;
  G4double IronPoleH = 0.3 *m;
  G4double CuPoleR = 0.85 *m;
  G4double CuPoleH = 0.2 *m;
  G4double IronPoleYpos = 0.2 *m;
  
  G4Paraboloid *solidIronPole = new  G4Paraboloid("solidIronPole", IronPoleH/2.0, IronPoleR, IronPoleR2);
  G4LogicalVolume *logicIronPole = new G4LogicalVolume(solidIronPole, magMaterial, "logicIronPole");
  
//  G4Tubs *solidIronPole = new G4Tubs("solidIronPole",  0.0,  IronPoleR, IronPoleH, 0.0, 2.0*M_PI);
//  G4LogicalVolume *logicIronPole = new G4LogicalVolume(solidIronPole, magMaterial, "logicIronPole");
  G4Tubs *solidCuPole = new G4Tubs("solidCuPole",  IronPoleR2,  CuPoleR, CuPoleH/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicCuPole = new G4LogicalVolume(solidCuPole, wireMaterial, "logicCuPole");
  
  G4double IronBoxX = 0.625 *m;
  G4double IronBoxY = 0.5 *m;
  G4double IronBoxZ = 1.2 *m;
  G4Box *solidFeBox = new G4Box("solidFeBox",  IronBoxX/2.0,  IronBoxY/2.0, IronBoxZ/2.0);
  G4LogicalVolume *logicFeBox = new G4LogicalVolume(solidFeBox, magMaterial, "logicFeBox");
  
  G4double IronTrapL = 1.0 *m;
  G4double IronTrapW = 0.8 *m;
  G4Trd *solidFeTrap = new G4Trd("solidFeTrap",  IronBoxZ/2.0,  IronTrapW/2.0,  IronBoxY/2.0,  IronBoxY/2.0,  IronTrapL/2.0);
  G4LogicalVolume *logicFeTrap = new G4LogicalVolume(solidFeTrap, magMaterial, "logicFeTrap");


 //Magnet Roof 
  G4double  IronRoofY = 0.9 *m;
  G4double  IronRoofL = 3.3 *m;
  G4double  IronRoofCutH = -0.4 *m;
//  G4double  IronRoofW = 2.0*IronRoofL*CuPoleR/(IronRoofL-CuPoleR); //0.85 *m;
  G4double  IronRoofW = 2.0*IronRoofL*(CuPoleR - IronBoxZ/2.0)/(IronRoofL-CuPoleR) + IronBoxZ;  
  G4Trd *solidFeTrapRoof0 = new G4Trd("solidFeTrapRoof0",  IronRoofW/2.0,  IronBoxZ/2.0,  IronRoofY/2.0,  IronRoofY/2.0,  IronRoofL/2.0);
  G4Tubs *solidRoofCut0 = new G4Tubs("solidRoofCut0",  CuPoleR,  2.0*CuPoleR, 2.0*IronRoofY, 0.0, M_PI);
  G4Transform3D trcut0(G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, 0.0, CuPoleR-IronRoofL/2.0));
  G4SubtractionSolid* solidFeTrapRoof1 = new G4SubtractionSolid("solidFeTrapRoof1", solidFeTrapRoof0, solidRoofCut0, trcut0);
  G4Transform3D trcut1(G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/6.0), 
                       G4ThreeVector(0.0, IronRoofCutH, CuPoleR-IronRoofL/2.0));
  G4SubtractionSolid* solidFeTrapRoof2 = new G4SubtractionSolid("solidFeTrapRoof2", solidFeTrapRoof1, solidRoofCut0, trcut1);
  G4Transform3D trcut2(G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), M_PI/4.0), 
                       G4ThreeVector(0.0, IronRoofCutH, IronRoofL/2.0-CuPoleR));
  G4SubtractionSolid* solidFeTrapRoof3 = new G4SubtractionSolid("solidFeTrapRoof3", solidFeTrapRoof2, solidRoofCut0, trcut2);

  G4Box *solidRoofCut1 = new G4Box("solidRoofCut1",  2*IronRoofW,  IronRoofW/2.0, IronRoofW/2.0);
  G4RotationMatrix rotcut3(G4ThreeVector(1.0, 0.0, 0.0), M_PI/4.0);
  G4double dalpha = M_PI/5.0;
  rotcut3.rotateZ(dalpha);
  G4double dy = 0.0*m;
  G4Transform3D trcut3(rotcut3, G4ThreeVector(0.0, -IronRoofW-IronRoofY/2.0 + dy, 0.0));// -(IronRoofW+IronRoofY)
  G4SubtractionSolid* solidFeTrapRoof4 = new G4SubtractionSolid("solidFeTrapRoof4", solidFeTrapRoof3, solidRoofCut1, trcut3);
  
  G4RotationMatrix rotcut4(G4ThreeVector(1.0, 0.0, 0.0), M_PI/4.0);
  rotcut4.rotateZ(-dalpha);
  G4Transform3D trcut4(rotcut4, G4ThreeVector(0.0, -IronRoofW-IronRoofY/2.0 + dy, 0.0));
  G4SubtractionSolid* solidFeTrapRoof = new G4SubtractionSolid("solidFeTrapRoof", solidFeTrapRoof4, solidRoofCut1, trcut4);

 
  G4LogicalVolume *logicFeTrapRoof = new G4LogicalVolume(solidFeTrapRoof, magMaterial, "logicFeTrapRoof");
  
//Collect everything to the assemblies
  G4AssemblyVolume  *magnetHalfAssembly = new G4AssemblyVolume();
  G4ThreeVector trmh(0.0, IronPoleYpos + IronPoleH/2.0, 0.0);
  magnetHalfAssembly->AddPlacedVolume(logicIronPole, trmh, new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), 3.0*0.5*M_PI));
  G4ThreeVector trmh2(0.0, IronPoleYpos + IronPoleH - CuPoleH/2.0, 0.0);
  magnetHalfAssembly->AddPlacedVolume(logicCuPole, trmh2, new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), 0.5*M_PI));

  G4ThreeVector trmh3(CuPoleR + IronTrapL + IronBoxX/2.0, IronBoxY/2.0, 0.0);
  magnetHalfAssembly->AddPlacedVolume(logicFeBox, trmh3,  new G4RotationMatrix());
  
  G4ThreeVector trmh4(CuPoleR + IronTrapL/2.0, IronBoxY/2.0, 0.0);
  magnetHalfAssembly->AddPlacedVolume(logicFeTrap, trmh4,  new G4RotationMatrix(G4ThreeVector(0.0, -1.0, 0.0), 0.5*M_PI));
  
  G4ThreeVector trmh5(IronRoofL/2.0 - CuPoleR, IronBoxY + IronRoofY/2.0, 0.0);
  magnetHalfAssembly->AddPlacedVolume(logicFeTrapRoof, trmh5,  new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), 0.5*M_PI));

  fMagnetAssembly = new G4AssemblyVolume();
  G4ThreeVector trmhassmbly(0.0, 0.0, 0.0);
   fMagnetAssembly->AddPlacedAssembly(magnetHalfAssembly, trmhassmbly, 0);
   fMagnetAssembly->AddPlacedAssembly(magnetHalfAssembly, trmhassmbly,
                                      new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), M_PI));
}



void MEP48MagnetAssembly::CostructSupport(G4LogicalVolume* lWorld, const G4ThreeVector& pos, const G4String& mname,
                                 G4bool rotate)
{
// Support
//  N6SetUp *lxs = N6SetUp::Instance();
}




