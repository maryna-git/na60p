
#ifndef N6MuTracker_h
#define N6MuTracker_h 1


class LxDetector;

class N6MuTracker: public LxDetector
{
  public:
    N6MuTracker(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~N6MuTracker() {};
    virtual void Construct();

  protected:
    void CreateMaterial();
    void AddSegmentation();
    G4AssemblyVolume* ConstructSupportAssembly();
};



#endif
