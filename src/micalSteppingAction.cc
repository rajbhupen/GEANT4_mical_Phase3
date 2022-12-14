// $Id: micalSteppingAction.cc,v 1.9 2005/02/02 17:11:11 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "micalSteppingAction.hh"

#include "micalDetectorConstruction.hh"
#include "micalEventAction.hh"
#include "MultiSimAnalysis.hh"
#include "G4Track.hh"

////#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalSteppingAction::micalSteppingAction(micalDetectorConstruction* det,
                                         micalEventAction* evt, MultiSimAnalysis* panalysis)
  :detector(det), eventaction(evt), pAnalysis(panalysis)					 
{
  // ical0Field = new micalElectroMagneticField();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalSteppingAction::~micalSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalSteppingAction::UserSteppingAction(const G4Step* aStep) {
  // cout<<"   void micalSteppingAction::UserSteppingAction(const G4Step* aStep) {        "<<endl;
  // micalElectroMagneticField* ical0Field = micalElectroMagneticField::EMFPointer;

  const G4Track* track = aStep->GetTrack();
  G4VPhysicalVolume* volume = track->GetVolume();
  
  G4Material* material = track->GetMaterial();
  //  G4Material* material2 = volume->GetLogicalVolume()->GetMaterial();

  //  if (material->GetName() !="Iron") G4cout <<"material "<<material->GetName()<<" "<<material2->GetName()<<" "<<" "<<detector->GetActiveMaterial()->GetName()<<" "<<aStep->GetTotalEnergyDeposit()<<G4endl;

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();
  //G4double edep = aStep->GetDeltaEnergy();

  G4String namex = volume->GetLogicalVolume()->GetName(); // material->GetName();
  // G4TouchableHistory* theTouchable = (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );
  // cout<<"XXXXXXXXXXXXXXXXXXXXXXX"<<endl;






G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
G4double TotalEdep = aStep->GetTotalEnergyDeposit()/keV;
G4int PdgId = track->GetDefinition()->GetPDGEncoding();

// if (abs(PdgId)!=13) { aStep->GetTrack()->SetTrackStatus(fStopAndKill);} //Manually kill a track
 
G4VPhysicalVolume* pv = preStepPoint->GetPhysicalVolume();
G4String name=pv->GetName();

G4VPhysicalVolume* pv1 = preStepPoint->GetPhysicalVolume();
G4String name1=pv1->GetName();
 
const G4VProcess* tpr2 = postStepPoint->GetProcessDefinedStep();
const G4VProcess* tpr3 = track->GetCreatorProcess();
 
// cout<<aStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() <<endl;

 G4int trkid = track->GetTrackID();
G4int parentId = track->GetParentID();
G4double edepNonIon = aStep->GetNonIonizingEnergyDeposit()/keV;
// cout<<aStep->GetPreStepPoint()->GetPosition();
// cout<<aStep->GetPostStepPoint()->GetPosition();
track->GetVertexKineticEnergy()/keV;
// cout<<track->GetLogicalVolumeAtVertex()->GetName();
// cout<<track->GetTotalEnergy()/keV;

// cout<<"name "<<name<<endl;





























  
  // for(int ij=0; ij<14; ij++) {
  //   cout<<"ij "<<ij<<" physiName "<<theTouchable->GetVolume(ij)->GetName()<<" material = "<<theTouchable->GetVolume(ij)->GetLogicalVolume()->GetMaterial()->GetName()<<endl;
  // }
  // cout<<"XXXXXXXXXXXXXXXXXXXXXXX"<<endl;
	//	if(volume->GetLogicalVolume()->GetName() == "IRLAYElog" && abs(track->GetDefinition()->GetPDGEncoding())==13) {
	//	if (volume->GetLogicalVolume()->GetName() == "GASRlog" && abs(track->GetDefinition()->GetPDGEncoding())==13) {
	if (abs(track->GetDefinition()->GetPDGEncoding())==13) {
		
		G4ThreeVector glbpos = 0.5*(aStep->GetPreStepPoint()->GetPosition() + aStep->GetPostStepPoint()->GetPosition());
		//		double magpos[3];
		//		magpos[0] = glbpos.x();
		//		magpos[1] = glbpos.y();
		//		magpos[2] = glbpos.z();
    
    // // cout<<"pos "<<glbpos<<endl;
		//		double magfld[6];
		//		ical0Field->GetFieldValue(magpos,magfld);
		//		cout<<"step:glbpos "<<glbpos<<" "<<volume->GetLogicalVolume()->GetName()<<" "<<volume->GetLogicalVolume()->GetMaterial()->GetName()<<" "<<edep<<" "<<track->GetDefinition()->GetPDGEncoding()<<endl;//" "<<magfld[0]/tesla<<" "<<magfld[1]/tesla<<endl;
    
	} // if(volume->GetLogicalVolume()->GetName() == "IRLAYElog")


  if (namex=="GASRlog") {
    pAnalysis->pdedz[0]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="QURZlog") {
    pAnalysis->pdedz[1]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="COATlog") {
    pAnalysis->pdedz[2]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="MYLARlog") {
    pAnalysis->pdedz[3]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="CUPLlog") {
    pAnalysis->pdedz[4]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="HoneyComblog") {
    pAnalysis->pdedz[5]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="ALlog") {
    pAnalysis->pdedz[6]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="G10Clog") {
    pAnalysis->pdedz[7]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="LAYElog") {
    pAnalysis->pdedz[8]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="IRONlog") {
    pAnalysis->pdedz[9]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="IRMODUlog") {
    pAnalysis->pdedz[10]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="INOMlog") {
    pAnalysis->pdedz[11]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (namex=="World") {
    pAnalysis->pdedz[12]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (strstr(namex,"COILlog")) {
    pAnalysis->pdedz[13]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  } else if (strstr(namex,"SPACER") || strstr(namex,"G4_Fe IRLAYElog") || strstr(namex,"G4_Fe")) {
    pAnalysis->pdedz[14]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
    // cout<<"pAnalysis->pdedz[0-14]->Fill();"<<endl;
  } else {
    //    G4cout <<"pox "<<aStep->GetPreStepPoint()->GetPosition()<<" "<<aStep->GetPreStepPoint()->GetPosition().z()/cm<<" "<<material->GetName()<<" "<<volume->GetLogicalVolume()->GetName()<<G4endl;
  }
  pAnalysis->pdedz[15]->Fill(aStep->GetPreStepPoint()->GetPosition().z()/cm, edep);
  // cout<<"pAnalysis->pdedz[15]->Fill()="<< aStep->GetPreStepPoint()->GetPosition().z()/cm <<endl;
  if (edep >0 && (material->GetName() =="rpcgas")) { // GASRphy  250709 ||  volume->GetName() =="quartz")) { 
    if (volume->GetMotherLogical() !=0) { 
      //      G4cout <<"xxxxxxx "<<volume->GetName()<<" "<<volume->GetMultiplicity()<<" "<<volume->GetCopyNo()<<" x "<<volume->GetMotherLogical()->GetName()<<" "<<volume->GetLogicalVolume()->GetName()<<" y "<<volume->GetLogicalVolume()->GetMaterial()->GetName()<<" "<<edep<<" "<<aStep->GetPreStepPoint()->GetPosition()<<" "<<edep<<G4endl;
      
    } else {
      G4cout <<"volname "<<volume->GetName()<<" "<<volume->GetMultiplicity()<<" "<<volume->GetCopyNo()<<" x "<<volume->GetLogicalVolume()->GetName()<<" y "<<volume->GetLogicalVolume()->GetMaterial()->GetName()<<" "<<aStep->GetPreStepPoint()->GetPosition()<<" "<<edep<<G4endl;
      
    }
  }
  
  G4double stepl = 0.;
  if (track->GetDefinition()->GetPDGCharge() != 0.) {
    stepl = aStep->GetStepLength();
  }

  float ds = 0.0;
  if ( (track->GetDefinition()->GetPDGEncoding() == 13 && material->GetName() == "G4_Fe") || 
       (track->GetDefinition()->GetPDGEncoding() == -13 && material->GetName() == "G4_Fe")) {
    ds = stepl;
    eventaction->rang += ds; // meghna
  }

  //  cout<<volume->GetName()<<" "<< detector->Getcal0()->GetName()<<endl;

  
  if (volume->GetName() == detector->Getcal0()->GetName()) { // GMA14 gdml

    //   cout<<volume->GetName()<<" "<< detector->Getcal0()->GetName()<<endl;
    eventaction->Addcal0(edep,stepl);
  } else {
    eventaction->AddGap(edep,stepl);
  }


  if (volume->GetName() == detector->Getcal1()->GetName()) {
      cout<<volume->GetName()<<" "<< detector->Getcal1()->GetName()<<endl;
    
    eventaction->Addcal1(edep,stepl);















    
  } 

  if (volume->GetName() == detector->Getcal2()->GetName()) {

    //  cout<<volume->GetName()<<" "<< detector->Getcal2()->GetName()<<endl;
    eventaction->Addcal2(edep,stepl);
  } 




  
  
  //  if (volume == detector->Getcal1()) eventaction->AddGap(edep,stepl);



  
       
 // save the random number seed of this event, under condition
 //// if(condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



