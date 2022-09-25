// $Id: micalDetectorMessenger.cc,v 1.9 2003/09/15 15:38:18 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "micalCal0SDMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalcal0SDMessenger::micalcal0SDMessenger(
                                           micalcal0SD* cal0SD)
  :micalcal0SDptr(cal0SD) { 
  micalDir = new G4UIdirectory("/mical/");
  micalDir->SetGuidance("UI commands of this example");
  
  cal0SDDir = new G4UIdirectory("/mical/cal0SD/");
  cal0SDDir->SetGuidance("digi control");

  CorrTimeSmrCmd = new G4UIcmdWithADoubleAndUnit("/mical/cal0SD/CorrTimeSmr",this);
  CorrTimeSmrCmd->SetGuidance("Set Correlated Time Smear");
  CorrTimeSmrCmd->SetParameterName("CorrTimeSmr",true, true);
  CorrTimeSmrCmd->SetDefaultUnit("ns");
  CorrTimeSmrCmd->SetDefaultValue(1.0);
  CorrTimeSmrCmd->SetUnitCategory("Time");

  UnCorrTimeSmrCmd = new G4UIcmdWithADoubleAndUnit("/mical/cal0SD/UnCorrTimeSmr",this);
  UnCorrTimeSmrCmd->SetGuidance("Set UnCorrelated Time Smear");
  UnCorrTimeSmrCmd->SetParameterName("UnCorrTimeSmr",true, true);
  UnCorrTimeSmrCmd->SetDefaultUnit("ns");
  UnCorrTimeSmrCmd->SetDefaultValue(1.0);
  UnCorrTimeSmrCmd->SetUnitCategory("Time");



  CorrIneffiCmd = new G4UIcmdWithADouble("/mical/cal0SD/CorrIneffiCmd",this);
  CorrIneffiCmd->SetGuidance("Set Correlated InEffi");
  CorrIneffiCmd->SetParameterName("CorrIneffi",true, true);
  //  CorrIneffiCmd->SetDefaultUnit("ns");
  CorrIneffiCmd->SetDefaultValue(0.1);
  // CorrIneffiCmd->SetUnitCategory("Time");


  UnCorrXIneffiCmd = new G4UIcmdWithADouble("/mical/cal0SD/UnCorrXIneffiCmd",this);
  UnCorrXIneffiCmd->SetGuidance("Set UnCorrelated X Ineffieciency in thye strip ");
  UnCorrXIneffiCmd->SetParameterName("UnCorrXIneffi",true, true);
  // UnCorrXIneffiCmd->SetDefaultUnit("ns");
  UnCorrXIneffiCmd->SetDefaultValue(1.0);
  // UnCorrXIneffiCmd->SetUnitCategory("Time");

  UnCorrYIneffiCmd = new G4UIcmdWithADouble("/mical/cal0SD/UnCorrYIneffiCmd",this);
  UnCorrYIneffiCmd->SetGuidance("Set UnCorrelated Y Ineffieciency in thye strip ");
  UnCorrYIneffiCmd->SetParameterName("UnCorrYIneffi",true, true);
  //  UnCorrYIneffiCmd->SetDefaultUnit("ns");
  UnCorrYIneffiCmd->SetDefaultValue(1.0);
  // UnCorrYIneffiCmd->SetUnitCategory("Time");

  


  
  TimeToDigiConvCmd = new G4UIcmdWithADouble("/mical/cal0SD/TimeToDigiConv",this);
  TimeToDigiConvCmd->SetGuidance("Set time to digi Conv Assuming Minimum scale of timing ~100 ps = 0.1 ns");
  TimeToDigiConvCmd->SetParameterName("TimeToDigiConv",true, true);
  TimeToDigiConvCmd->SetDefaultValue(0.1);

  SigSpeedCmd = new G4UIcmdWithADouble("/mical/cal0SD/SigSpeed",this);
  SigSpeedCmd->SetGuidance("Set signal speed in the strip in the units of ns/strip");
  SigSpeedCmd->SetParameterName("SigSpeed",true, true);
  SigSpeedCmd->SetDefaultValue(0.15);



  CorrNoiseCmd1 = new G4UIcmdWithADouble("/mical/cal0SD/CorrNoise1",this);
  CorrNoiseCmd1->SetGuidance("Set CorrNoise1 in the strip");
  CorrNoiseCmd1->SetParameterName("CorrNoise1",true, true);
  CorrNoiseCmd1->SetDefaultValue(0.);

  

  CorrNoiseCmd2 = new G4UIcmdWithADouble("/mical/cal0SD/CorrNoise2",this);
  CorrNoiseCmd2->SetGuidance("Set CorrNoise2 in the strip");
  CorrNoiseCmd2->SetParameterName("CorrNoise2",true, true);
  CorrNoiseCmd2->SetDefaultValue(0.);

  CorrNoiseCmd3 = new G4UIcmdWithADouble("/mical/cal0SD/CorrNoise3",this);
  CorrNoiseCmd3->SetGuidance("Set CorrNoise3 in the strip");
  CorrNoiseCmd3->SetParameterName("CorrNoise3",true, true);
  CorrNoiseCmd3->SetDefaultValue(0.);

  
  RootRandomCmd = new G4UIcmdWithAnInteger("/mical/cal0SD/RootRandom",this);
  RootRandomCmd->SetGuidance("To switch root random on or off");
  RootRandomCmd->SetParameterName("RootRandom",true,true);
  RootRandomCmd->SetDefaultValue(0);


  RandomNoiseCmd = new G4UIcmdWithAnInteger("/mical/cal0SD/RandomNoise",this);
  RandomNoiseCmd->SetGuidance("To set Random Noise in the strip");
  RandomNoiseCmd->SetParameterName("RandomNoise",true,true);
  RandomNoiseCmd->SetDefaultValue(0);

  

  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalcal0SDMessenger::~micalcal0SDMessenger() {
  delete micalDir;
  delete cal0SDDir;
  delete CorrTimeSmrCmd;    delete UnCorrTimeSmrCmd;
  delete CorrIneffiCmd;
  delete UnCorrYIneffiCmd;
  delete UnCorrXIneffiCmd;
  delete SigSpeedCmd;
  delete TimeToDigiConvCmd;
  delete RootRandomCmd;
  delete RandomNoiseCmd;
  delete CorrNoiseCmd1;
  delete CorrNoiseCmd2;
  delete CorrNoiseCmd3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalcal0SDMessenger::SetNewValue(G4UIcommand* command,G4String newValue) { 
  if( command == CorrTimeSmrCmd )
    { micalcal0SDptr->SetCorrTimeSmear(CorrTimeSmrCmd->GetNewDoubleValue(newValue));}
  if( command == UnCorrTimeSmrCmd )
    { micalcal0SDptr->SetUnCorrTimeSmear(UnCorrTimeSmrCmd->GetNewDoubleValue(newValue));}
  if( command == TimeToDigiConvCmd )
    { micalcal0SDptr->SetTimeToDigiConv(TimeToDigiConvCmd->GetNewDoubleValue(newValue));}
  if( command == SigSpeedCmd )
    { micalcal0SDptr->SetSignalSpeed(SigSpeedCmd->GetNewDoubleValue(newValue));}
  if( command == RootRandomCmd )
    { micalcal0SDptr->SetRootRandom(RootRandomCmd->GetNewIntValue(newValue));}

  if( command == RandomNoiseCmd )
    { micalcal0SDptr->SetRandomNoise(RandomNoiseCmd->GetNewIntValue(newValue));}



  if( command == CorrNoiseCmd1 )
    { micalcal0SDptr->SetCorrNoise1(CorrNoiseCmd1->GetNewDoubleValue(newValue));}

  if( command == CorrNoiseCmd2 )
    { micalcal0SDptr->SetCorrNoise2(CorrNoiseCmd2->GetNewDoubleValue(newValue));}

  if( command == CorrNoiseCmd3 )
    { micalcal0SDptr->SetCorrNoise3(CorrNoiseCmd3->GetNewDoubleValue(newValue));}


  if( command == CorrIneffiCmd )
    { micalcal0SDptr->SetCorrInefficiency(CorrIneffiCmd->GetNewDoubleValue(newValue));}


  if( command == UnCorrXIneffiCmd )
    { micalcal0SDptr->SetUnCorrXInefficiency(UnCorrXIneffiCmd->GetNewDoubleValue(newValue));}

  if( command == UnCorrYIneffiCmd )
    { micalcal0SDptr->SetUnCorrYInefficiency(UnCorrYIneffiCmd->GetNewDoubleValue(newValue));}
  
  

 





 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
