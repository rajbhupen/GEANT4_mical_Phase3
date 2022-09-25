// $Id: micalDetectorMessenger.cc,v 1.9 2003/09/15 15:38:18 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "micalCal1SDMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalcal1SDMessenger::micalcal1SDMessenger(
                                           micalcal1SD* cal1SD)
:micalcal1SDptr(cal1SD) { 
  micalDir = new G4UIdirectory("/mical/");
  micalDir->SetGuidance("UI commands of this example");
  
  cal1SDDir = new G4UIdirectory("/mical/cal1SD/");
  cal1SDDir->SetGuidance("digi control");

  

  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalcal1SDMessenger::~micalcal1SDMessenger() {
  delete micalDir;
  delete cal1SDDir;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalcal1SDMessenger::SetNewValue(G4UIcommand* command,G4String newValue) { 
 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
