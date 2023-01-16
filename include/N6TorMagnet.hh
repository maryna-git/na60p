
#ifndef N6TorMagnet_h
#define N6TorMagnet_h 1


class LxDetector;

class N6TorMagnet: public LxDetector
{
  public:
    N6TorMagnet(DetectorConstruction *detc = 0): LxDetector(detc), fFieldVolume(0), fCoilSideAngle(0), fCoilR0Max(0), fCoilR1Max(0) {};
    virtual ~N6TorMagnet() {};
    virtual void Construct();
    G4LogicalVolume* GetFieldVolume() const { return fFieldVolume; }
    void GetFieldVolumeParameters(G4double &r0max, G4double &r1max) const { r0max = fCoilR0Max; r1max = fCoilR1Max; }

  protected:
    void CreateMaterial();
    void AddSegmentation();
    G4AssemblyVolume* ConstructSupportAssembly();

  protected:
    G4LogicalVolume  *fFieldVolume;

    G4double fCoilSideAngle, fCoilR0Max, fCoilR1Max;

};



#endif
