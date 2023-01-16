
#ifndef N6SiTracker_h
#define N6SiTracker_h 1


class LxDetector;

class N6SiTracker: public LxDetector
{
  public:
    N6SiTracker(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~N6SiTracker() {};
    virtual void Construct();

  protected:
    void CreateMaterial();
    void AddSegmentation();
    G4AssemblyVolume* ConstructSupportAssembly();
    
    void ConstructMagnet(G4LogicalVolume *mfContainer, const G4ThreeVector mf);
};


#endif
