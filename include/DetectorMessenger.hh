//
/// \brief Definition of the DetectorMessenger class
//

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;


class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction* );
   ~DetectorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    DetectorConstruction*      fDetector;

    G4UIdirectory*             fTestemDir;
    G4UIdirectory*             fDetDir;

    G4UIcmdWithAString*        fWorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* fWorldZCmd;
    G4UIcmdWithADoubleAndUnit* fWorldXYCmd;


    G4UIdirectory*  fDetBFieldDir;
    G4UIcommand*    fSetBFieldValue;
    G4UIcommand*    fSetBFieldDistrib;

};

#endif

