//
/// \brief Implementation of the Magnetic field models classes
//

#include <algorithm>
#include <functional>
#include <fstream>

#include "G4UnitsTable.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4UserLimits.hh"
#include "G4TwoVector.hh"

#include "N6BFields.hh"




void TorBField::GetFieldValue(const G4double ppos[4], G4double *MagField) const
{
  G4ThreeVector pp(ppos[0], ppos[1], ppos[2]);
  pp -= fPos;
  G4TwoVector rv(pp[0], pp[1]);
//  G4TwoVector rv(pp[0], pp[1]);
  G4TwoVector bfxy = G4TwoVector(-pp[1], pp[0]).unit();
  G4double r_max = fr0Max + (pp[2] - fzMin) * ftgrMax;
  G4double r_min = fr0Min + (pp[2] - fzMin) * ftgrMin;

  G4double rr = rv.mag();
  G4double bf = 1.0 / ( (1.0 + exp((r_min-rr)/ftfd)) * (1.0 + exp((rr-r_max)/ftfd))) * 1.0/(rr + 0.01);

  bfxy *= bf * fBMax;
  MagField[0] = bfxy(0);
  MagField[1] = bfxy(1);
  MagField[2] = 0.0;

}



void TorBField::DumpBFieldModel(const G4double phi) const
{
  G4double pos[4] = {fPos.x(), fPos.y(), fPos.z()+fzMax, 0.0};
  G4double bfield[3];
  G4cout << "Radial point for the toroidal magnet:" << G4endl;
  G4int np = 200;
  for (G4int ii = 0; ii < np; ++ii) {
	G4double rr = 0.01*fr0Min + ii * 1.5*fr1Max/static_cast<G4double>(np-1);
	pos[0] = rr * cos(phi);
	pos[1] = rr * sin(phi);
	GetFieldValue(pos, bfield);
	G4ThreeVector bb(bfield[0], bfield[1], bfield[2]);
	bb *= 1.0/tesla;
	G4cout << "rr/B/Bx/By/Bz:  " << rr << "   " <<  bb.mag() << "  " << bb.x() << "  " << bb.y() << "  " << bb.z() << G4endl;
  }
  G4cout << "Radial point for the toroidal magnet done." << G4endl;
}


//////////////////////////////////////////////////////////////////////////////////
void DipoleBField::GetFieldValue(const G4double ppos[4], G4double *MagField) const
{
  G4ThreeVector pp(ppos[0], ppos[1], ppos[2]);
  pp -= fPos;
  G4TwoVector rv(pp[0], pp[1]);
  G4TwoVector bfxy = rv.orthogonal().unit();

  G4double rr = rv.mag();
  G4double bf = 1.0/(rr + 0.01);

  bfxy *= bf * fBMax;
  MagField[0] = bfxy(0);
  MagField[1] = 0.0;
  MagField[2] = 0.0;

}
