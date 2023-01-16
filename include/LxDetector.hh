
#ifndef LxDetector_h
#define LxDetector_h 1

class DetectorConstruction;
class G4AssemblyVolume;


class LxDetector
{
  public:
    LxDetector(DetectorConstruction *detc = 0): fDetector(detc) {};
    virtual ~LxDetector() {};
    virtual void Construct() = 0;
 
  protected:
    DetectorConstruction*  fDetector;
};


#endif
