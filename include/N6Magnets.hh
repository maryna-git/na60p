

#ifndef N6Magnets_h
#define N6Magnets_h 1

#include <vector>

#include "G4ThreeVector.hh"


class N6MagnetAssembly
{
  public:
    N6MagnetAssembly() {};
    virtual ~N6MagnetAssembly() {};
    virtual G4LogicalVolume* GetFieldVolume(const G4String lvname) = 0;
    virtual G4AssemblyVolume* GetAssembly() = 0;
    virtual void CostructSupport(G4LogicalVolume* lWorld, const G4ThreeVector& pos, const G4String& mname,
                                 G4bool rotate = false) = 0;

};



class MEP48MagnetAssembly: public N6MagnetAssembly
{
  public:
    MEP48MagnetAssembly(const G4String mtypename);
    virtual ~MEP48MagnetAssembly() {if (fMagnetAssembly && !fMagnetAssembly->GetImprintsCount() ) {delete fMagnetAssembly;}};
    virtual G4AssemblyVolume* GetAssembly() { ConstructMagnet(); return fMagnetAssembly; };
    virtual G4LogicalVolume* GetFieldVolume(const G4String lvname);
    virtual void CostructSupport(G4LogicalVolume* lWorld, const G4ThreeVector& pos, const G4String& mname,
                                 G4bool rotate = false);

  protected:
    void ConstructMagnet();

  protected:
    G4String          fMagnetType;
    G4AssemblyVolume *fMagnetAssembly;
};

#endif

