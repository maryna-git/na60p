
#include "globals.hh"
// #include "G4Cache.hh"

#ifndef N6SETUP_h
#define N6SETUP_h 1

class N6SetUp
{
public: 
  static N6SetUp *Instance();
  ~N6SetUp();

protected:
  N6SetUp();
 

  static N6SetUp *flxSetup;

public: 
  enum tTargetType : G4int {tfoil, twire};

public: 
  G4bool OverlapTest; 
   
  G4String EnvironmentMaterial;
  G4double FloorSurfaceYpos;
  G4double FloorY;
  
  G4double CeilingSurfaceYpos;

  G4double IPBoxX;
  G4double IPBoxY;
  G4double IPContainerZ;
  G4double IPBoxZ;
  G4double IPBoxThickness;
   
  G4double BPipeR;
  G4double BPipeThickness;
  G4String BeamPipeMaterial;
  G4String BeamPipeWindowMaterial;
  G4String BeamPipeVacuumMaterial;


  G4String MagnetMaterial;
  G4String TypMBWireMaterial;
  
  G4double DumpMagFieldY;
  G4double IPMagFieldY;
  
  G4String  ScintCerenkovPhysics;
  G4String  GTargetType;
  G4double   BSMCaloX;
  G4double   BSMCaloY;
  G4double   BSMCaloNCellX;
  G4double   BSMCaloNCellY;
  
  G4double IPMagnetXpos;
  G4double IPMagnetYpos;
  G4double IPMagnetZpos;

  
  G4double ms0_pos;
  G4double ms1_pos;
  G4double ms2_pos;
  G4double ms3_pos;    
  G4double ms4_pos; 
  G4double ms5_pos;

  G4double TMagnetZpos;
  G4double tMagRMax;
  G4double tMagRMin;
  G4String tMagnetCoilMaterial;

  G4double DipoleFieldB;
  G4double DipoleFieldY;
  G4double DipoleFieldR;
  G4double DipoleMagnetThickness;
  G4double DipoleMagnetXpos;
  G4double DipoleMagnetYpos;
  G4double DipoleMagnetZpos;
  G4double DipoleMagnetX;
  G4double DipoleMagnetY;
  G4double DipoleMagnetZ;
  
  G4double PbTarget_zpos_0;
  G4double PbTarget_zpos_1;
  G4double PbTarget_zpos_2;
  G4double PbTarget_zpos_3;
  G4double PbTarget_zpos_4;
  G4double PythiaInutZshift;
  
  G4double AbsorberBeOH_1_pos;
  G4double AbsorberBeOH_2_pos;
  G4double AbsorberBeO_1_dz;
  G4double AbsorberBeO_2_dz;
  
  G4double AbsorberC_Wall_pos;
};

#endif
