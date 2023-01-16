
#include <vector>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "N6SetUp.hh"


N6SetUp *N6SetUp::flxSetup = 0;


N6SetUp::N6SetUp()
{
  OverlapTest = true;

  EnvironmentMaterial = "Air20";
//  EnvironmentMaterial = "Galactic";
  FloorSurfaceYpos = -2478.0 *mm; //-1700.0 *mm;
  FloorY = 500.0 *mm;


  IPBoxX = 45.0 *cm;
  IPBoxY = 25.0 *cm;
//   IPContainerZ = 100.0 *cm;
  IPBoxZ = 60.0 *cm;
  IPBoxThickness = 5.0 *mm;

  BPipeR = 24.13 *mm;          // DN40;
  BPipeThickness = 1.651 *mm;  // seems to be the thinnest for DN40  //2.0 *mm;

  BeamPipeMaterial = "Aluminium";  // "StainlessSteel";  //"Iron";
  BeamPipeWindowMaterial = "Aluminium";

  BeamPipeVacuumMaterial = "XFELVacuum"; //"Galactic";
//   GTargetType = twire;
  GTargetType = tfoil;
  MagnetMaterial = "Iron";
  
  IPMagFieldY = -10000.0*gauss; //-16000.0*gauss;
  TypMBWireMaterial = "Copper";
  ScintCerenkovPhysics = false;
  
  BSMCaloX = 2005.0 *mm;
  BSMCaloY = 2005.0 *mm;
  BSMCaloNCellX = 401;
  BSMCaloNCellY = 401;

// Muon chambers Low-E setup
    ms0_pos = 293.9695 *cm;
    ms1_pos = 360.9695 *cm;
    ms2_pos = 736.965 *cm;
    ms3_pos = 799.965 *cm;    
    ms4_pos = 999.9335 *cm; 
    ms5_pos = 1039.9335 *cm;
/* Muon chambers High-E setup
    ms0_pos = 624. *cm;
    ms1_pos = 691. *cm;
    ms2_pos = 1067. *cm;
    ms3_pos = 1130. *cm;    
    ms4_pos = 1330. *cm; 
    ms5_pos = 1370 *cm;*/    
//Toroid
  TMagnetZpos = 547. *cm;//547.5 *cm; Low-E setup
//TMagnetZpos = 877. *cm;//877.5 *cm; High-E setup
  tMagRMax = 2971.0 *mm;
  tMagRMin = 1950.0 *mm;
  tMagnetCoilMaterial = "Aluminium";
//Dipole
  DipoleFieldB = 1.5 *tesla;
  DipoleFieldY = 400.0 *mm;
  DipoleFieldR = 500.0 *mm;
  
  DipoleMagnetX = 3300.0*mm;
  DipoleMagnetY = 3000.0 *mm;
  DipoleMagnetZ = 2000.0 *m;
  DipoleMagnetXpos = 0.0 *mm; // -50.0 *mm;
  DipoleMagnetYpos = 0.0 *mm;// - 1600.0 *mm;
  DipoleMagnetZpos = 0.0 *mm;// -650.0 *cm;

  DipoleMagnetThickness = 50.0 *cm;
  
  //Targets positions
  PbTarget_zpos_0 = -5.3 *cm;
  PbTarget_zpos_1 = -4.1 *cm;
  PbTarget_zpos_2 = -2.9 *cm;
  PbTarget_zpos_3 = -1.7 *cm;
  PbTarget_zpos_4 = -0.5 *cm;
  //Pythia Input
  PythiaInutZshift = PbTarget_zpos_2; // to start from the middle PB target
 
  //Absorbers Low_E setup
  AbsorberBeOH_1_pos = 71.5 *cm;
  AbsorberBeOH_2_pos = 124.0 *cm;//AbsorberBeOH_1_pos +AbsorberBeO_1_dz/2.0 + AbsorberBeO_2_dz/2.0; //1125 *mm;
  AbsorberBeO_1_dz = 33.0 *cm;
  AbsorberBeO_2_dz = 72.0 *cm;
  
  AbsorberC_Wall_pos =  900.0 *cm; //910.0 *cm;   Fluka has bug:wall overlaps with MS5

  //High-E setup
//  AbsorberBeOH_1_pos = 62.25*cm;// 
//  AbsorberBeOH_2_pos = 114.75*cm;/
//  AbsorberBeO_1_dz = 34.5 *cm;
//  AbsorberBeO_2_dz = 70.5 *cm;

//  AbsorberC_Wall_pos = 1230.0 *cm;
}



N6SetUp *N6SetUp::Instance() 
{
  if (!flxSetup) {
    flxSetup = new N6SetUp();
  }
  return flxSetup;
}




