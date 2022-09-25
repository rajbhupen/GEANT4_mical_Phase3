
// $Id: micalDetectorMessenger.hh,v 1.6 2002/12/16 16:37:26 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef micalCal1SDMessenger_h
#define micalCal1SDMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "micalCal1SD.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"

class micalcal1SD;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalcal1SDMessenger: public G4UImessenger
{
  public:
    micalcal1SDMessenger(micalcal1SD* );
   ~micalcal1SDMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    micalcal1SD* micalcal1SDptr;
    
    G4UIdirectory*             micalDir;
    G4UIdirectory*             cal1SDDir;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

