
#ifndef LxBFields_h
#define LxBFields_h 1

#include <vector>

#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"



class TorBField: public G4MagneticField
{
public:
  TorBField() : G4MagneticField(), fPos(G4ThreeVector()),
                fBMax(0.0), fr0Min(0.0), fr0Max(0.0), fr1Min(0.0), fr1Max(0.0), fzMin(0.0), fzMax(0.0),
                ftfd(0.0), ftgrMin(1.0), ftgrMax(1.0) {};
  TorBField(const G4ThreeVector pos): G4MagneticField(), fPos(pos), fBMax(0.0),
                fr0Min(0.0), fr0Max(0.0), fr1Min(0.0), fr1Max(0.0), fzMin(0.0), fzMax(0.0),
                ftfd(0.0), ftgrMin(1.0), ftgrMax(1.0) {};
  ~TorBField() {};

  void SetParameters(const G4double bfmax, const G4double r0Min, const G4double r0Max, const G4double r1Min, const G4double r1Max,
                     const G4double zMin, const G4double zMax, const G4double tfd) {
                     fBMax = bfmax; fr0Min = r0Min; fr0Max = r0Max; fr1Min = r1Min; fr1Max = r1Max;
                     fzMin = zMin; fzMax = zMax; ftfd = tfd;
                     ftgrMax = (fr1Max - fr0Max)/(fzMax - fzMin);
                     ftgrMin = (fr1Min - fr0Min)/(fzMax - fzMin); }

  virtual void GetFieldValue(const G4double ppos[4], G4double *MagField) const;

  void DumpBFieldModel(const G4double phi) const;

protected:

  G4ThreeVector fPos;

  G4double fBMax;
  G4double fr0Min, fr0Max, fr1Min, fr1Max, fzMin, fzMax, ftfd;
  G4double ftgrMin, ftgrMax;

};



class DipoleBField: public G4MagneticField
{
public:
  DipoleBField() : G4MagneticField(), fPos(G4ThreeVector()),
                fBMax(0.0), fr0Min(0.0), fr0Max(0.0), fr1Min(0.0), fr1Max(0.0), fzMin(0.0), fzMax(0.0),
                ftfd(0.0), ftgrMin(1.0), ftgrMax(1.0) {};
  DipoleBField(const G4ThreeVector pos): G4MagneticField(), fPos(pos), fBMax(0.0),
                fr0Min(0.0), fr0Max(0.0), fr1Min(0.0), fr1Max(0.0), fzMin(0.0), fzMax(0.0),
                ftfd(0.0), ftgrMin(1.0), ftgrMax(1.0) {};
  ~DipoleBField() {};

  void SetParameters(const G4double bfmax, const G4double r0Min, const G4double r0Max, const G4double r1Min, const G4double r1Max,
                     const G4double zMin, const G4double zMax, const G4double tfd) {
                     fBMax = bfmax; fr0Min = r0Min; fr0Max = r0Max; fr1Min = r1Min; fr1Max = r1Max;
                     fzMin = zMin; fzMax = zMax; ftfd = tfd;
                     ftgrMax = (fr1Max - fr0Max)/(fzMax - fzMin);
                     ftgrMin = (fr1Min - fr0Min)/(fzMax - fzMin); }

  virtual void GetFieldValue(const G4double ppos[4], G4double *MagField) const;

  void DumpBFieldModel(const G4double phi) const;

protected:

  G4ThreeVector fPos;

  G4double fBMax;
  G4double fr0Min, fr0Max, fr1Min, fr1Max, fzMin, fzMax, ftfd;
  G4double ftgrMin, ftgrMax;

};


#endif
