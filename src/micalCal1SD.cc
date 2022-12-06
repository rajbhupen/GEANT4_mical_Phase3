//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B2TrackerSD.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B2TrackerSD.cc
/// \brief Implementation of the B2TrackerSD class
//refer:https://www.slac.stanford.edu/xorg/geant4/KISTI2019/HandsOn3/#ex1s6
#include "G4VSensitiveDetector.hh"
#include "micalCal1SD.hh" //
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RunManager.hh"//..
#include "micalRunAction.hh"//..
#include "micalDetectorParameterDef.hh"
#include <iostream>
#include <fstream>
#include "TRandom3.h"
using namespace std; 


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//B1TrackerSD::B1TrackerSD(B1RunAction* runAction)
//: G4UserEventAction(),
//fRunAction(runAction)

//{} 
/////////////////////////////////////////////////////

micalcal1SD::micalcal1SD(const G4String name)
  : G4VSensitiveDetector(name),
    cal1Collection(NULL),
    // SWidth(0),
    // fEnvelopeBox(0),
    numberInCell(20000),
    Counter(0),
    InCell(0)
{

  cout<<"cal1sd constructor"<<endl;
  G4String HCname;
  paradef = micalDetectorParameterDef::AnPointer;
  collectionName.insert(HCname="cal1Collect");
  cal1SDMessenger = new micalcal1SDMessenger(this);  
  pAnalysis = MultiSimAnalysis::AnPointer;
  //	CmvStrip_pointer = new CmvStrip_Manager();
  //	SipmHit_pointer = new SipmHit_Manager();
  //	CmvHit_pointer = new CmvHit_Manager();
  SetPhotonSpeed(162.0);
  SetCMVadctons(0.1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalcal1SD::~micalcal1SD() 
{
  cout<<"cal1sd distructor"<<endl;
  
  //delete fMessenger;
  // for (unsigned ij=0; ij<CmvStrip_pointer->CmvStrip_list.size(); ij++) {
  //   if (CmvStrip_pointer->CmvStrip_list[ij]) {
  //     //  cout <<"ij "<< ij<<" "<<CmvStrip_pointer->CmvStrip_list.size()<<endl;
  //     delete CmvStrip_pointer->CmvStrip_list[ij]; CmvStrip_pointer->CmvStrip_list[ij]=0;
  //   }
  // }

  // CmvStrip_pointer->CmvStrip_list.clear();
  // if (CmvStrip_pointer) {delete CmvStrip_pointer; CmvStrip_pointer=0;}

  // for (unsigned ij=0; ij<SipmHit_pointer->SipmHit_list.size(); ij++) {
  //   if (SipmHit_pointer->SipmHit_list[ij]) {
  //     //  cout <<"ij "<< ij<<" "<<SipmHit_pointer->SipmHit_list.size()<<endl;
  //     delete SipmHit_pointer->SipmHit_list[ij]; SipmHit_pointer->SipmHit_list[ij]=0;
  //   }
  // }

  // SipmHit_pointer->SipmHit_list.clear();
  // if (SipmHit_pointer) {delete SipmHit_pointer; SipmHit_pointer=0;}


	
}
//.......................................................................

void micalcal1SD::Initialize(G4HCofThisEvent* hce)
{
  cout<<"micalcal1SD::Initialize start"<<endl;
  //InCell =0;

 
  CmvStrip_pointer = new CmvStrip_Manager();
  SipmHit_pointer = new SipmHit_Manager();
  
  CmvStrip_pointer->CmvStrip_list.clear();
  SipmHit_pointer->SipmHit_list.clear();
  //	CmvHit_pointer->CmvHit_list.clear();
  
  for(int op=0; op<3;op++) {
    partopscint[op] = paradef->partopscint[op];
  }
  
  AirGapScintTop= paradef->AirGapScintTop;


  int jmax;
  for(int i =0;i<4;i++){

    if(i==0){
      jmax=4;
    }
    else{
      jmax=3;
    }
    for(int j=0;j<jmax;j++){
      for(int k=0;k<3;k++){
	PhyVolGlPos[i][j][k] = paradef->ScintLayGPos[i][j][k];

      }//k
  
    }//j

 
  }//i
	

  NoScntStrpTop = paradef->GetNoScntStrpTop();//88
  NoScntStrpSide = paradef->GetNoScntStrpSide();//40
  Sipm_Pedestal = paradef->GetSipm_Pedestal();//100
  Cmv_Threshold = paradef->GetCmv_Threshold();//120
	      
		
  //..................................................................................................	
	
  // Create hits collection
  
  cal1Collection 
    = new micalcal1HitsCollection(SensitiveDetectorName, collectionName[0]); 
  
  //cout<<"SensitiveDetectorName" <<SensitiveDetectorName<<endl;
  //cout<<"collectionName[0]" <<collectionName[0]<<endl;
  // Add this collection in hce
  
  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, cal1Collection ); 
  
  paradef =  micalDetectorParameterDef::AnPointer;
  
  
  cout<<"micalcal1SD::Initialize ends 	"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool micalcal1SD::ProcessHits(G4Step* aStep,  
                                G4TouchableHistory*)  //called everytime when the particle hts the sensitive dete
{  
  
  cout<<endl<<"micalcal1SD::ProcessHits start"<<endl;
  
  G4double edep = aStep->GetTotalEnergyDeposit()/keV;

  G4TouchableHistory* theTouchable = (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );
  
  int pdgid = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  //to  print copy number of volumes
   int level = theTouchable->GetHistoryDepth();
  
  cout<<"level="<<level<<endl; //level=4
  for(int ij=0; ij<level+1; ij++) { //level=4
  
  cout<<ij<<" volname "<<theTouchable->GetVolume(ij)->GetName()
  <<" copyNo. "<<theTouchable->GetCopyNumber(ij) <<"edep "<< edep << "pid "<<pdgid<< "nInLA:"<<theTouchable->GetCopyNumber(1)<<endl;
  }
  
  G4ThreeVector parmom = aStep->GetTrack()->GetMomentum();
  
  //  double momentum= parmom.mag();
  //  double polang	= parmom.theta();
  //  double aziang	= parmom.phi();
  //  cout<<"pid "<<aStep->GetTrack()->GetDefinition()->GetPDGEncoding()<<endl;
  //  if (abs(aStep->GetTrack()->GetDefinition()->GetPDGEncoding()) !=13) return false;
  //cout<<"edep "<< edep<<endl;
    if (edep==0.)  return false;

  double plen = aStep->GetStepLength();//pathlength
  cout<<"before birks "<< edep<<endl;

    edep = 0.973*edep/(1+(0.001*0.126*edep/plen));//birks law kB 0.126mm/MeV polystyrene 
  //edep =  gRandom->Poisson(edep);
  cout<<"after birks "<< edep<<endl;

    
  cout<<"pid "<<aStep->GetTrack()->GetDefinition()->GetPDGEncoding()<<endl;  
  //We are simulating a  detector that will trigger only if some energy has been deposited (i.e. via ionization), 
  //for example if a neutron passes through the detector (without making interactions) its passage should not be recorded.
  // Check the energy deposited in the step, if zero do not do anything.
  
  
  G4ThreeVector glbpos =0.5*(aStep->GetPreStepPoint()->GetPosition() + aStep->GetPostStepPoint()->GetPosition()); //avg distance betwen the two

 
 



  
  G4ThreeVector tmp;
  
  //  tmp = (1/m)*glbpos;
  // tmp.y() = (1/m)*glbpos.y();
  // tmp.z() = (1/m)*glbpos.z();
  
  //cout<<"tmpx" << tmp.x() <<endl;
  //cout<<"tmpy" << tmp.y() <<endl;
  //cout<<"tmpz" << tmp.z() <<endl;
  /*
  G4int nScntStrp = theTouchable->GetCopyNumber(0);
  G4int nInLA = theTouchable->GetCopyNumber(1) ;
  //  G4int tmplay = theTouchable->GetCopyNumber(2) ;
  G4String loc = theTouchable->GetVolume(1)->GetName();
  
  // cout<<"Layer No: "<<nInLA<<" presteppoint: "<<aStep->GetPreStepPoint()->GetPosition()<<" poststeppoint: "<<aStep->GetPostStepPoint()->GetPosition()<<endl;
  G4int loc_no=-1;  //Number starting from 1, not from ZERO



  if(strstr(loc,"physiTopScint_1cm") || strstr(loc,"physiTopScint_2cm")){
    loc_no = 1;
  }

  else if (strstr(loc,"physiSideScint_L")){
    loc_no = 2;
  }

    else if (strstr(loc,"physiSideScint_R")){
    loc_no = 3;
  }

    else if (strstr(loc,"physiSideScint_D")){
    loc_no = 4;
  }

  //5 is reserved for front wall
      else if (strstr(loc,"physiSideScint_miniL")){
    loc_no = 6;
  }

    else if (strstr(loc,"physiSideScint_miniR")){
    loc_no = 7;
  }

*/
  //Updated for Tile://18062022 RSA

 G4int nInTile = theTouchable->GetCopyNumber(0);
G4int TileNo = theTouchable->GetCopyNumber(2);
 G4int nScntStrp = nInTile +8*TileNo; //8 EPS  makeup one Tile
 G4int nInLA = theTouchable->GetCopyNumber(3) ;
  //  G4int tmplay = theTouchable->GetCopyNumber(2) ;
  G4String loc = theTouchable->GetVolume(3)->GetName();
  
  // cout<<"Layer No: "<<nInLA<<" presteppoint: "<<aStep->GetPreStepPoint()->GetPosition()<<" poststeppoint: "<<aStep->GetPostStepPoint()->GetPosition()<<endl;
  G4int loc_no=-1;  //Number starting from 1, not from ZERO



  if(strstr(loc,"physiTopLay_1cm") || strstr(loc,"physiTopLay_2cm")){
    loc_no = 1;
  }

  else if (strstr(loc,"physiSideScint_L")){
    loc_no = 2;
  }

    else if (strstr(loc,"physiSideScint_R")){
    loc_no = 3;
  }

    else if (strstr(loc,"physiSideScint_B")){
    loc_no = 4;
  }

  //5 is reserved for front wall
      else if (strstr(loc,"physiSideScint_miniL")){
    loc_no = 6;
  }

    else if (strstr(loc,"physiSideScint_miniR")){
    loc_no = 7;
  }

  

  
  //  cout<<"pdg "<<pdgid<<" loc " <<loc_no<<" "<< loc<<" nInLA " << nInLA<<" "<<tmplay<<" nScntStrp " << nScntStrp<<endl;
  //cout<<"loc_no " << loc_no<<endl;
  
  G4double atime = aStep->GetPreStepPoint()->GetGlobalTime()/(ns);
  G4ThreeVector localpos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(glbpos);
  cout<<" Cal1SD  -------"<<atime<<" "<<glbpos<<" "<<localpos<<endl;
  // double tra_LPosy = localpos.y() + partopscint[2]; // shift of origin 
  //To reduce the usage of memory..................................
  // unsigned int ScntStrpid;   // 16 bits  declared in hh   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
  
  unsigned int IdSiPM=0;
  
  IdSiPM = loc_no;     //occupies 3 bits for CMVD number number (if we have front CMVD need 3 bits)
  // cout<<"Process hits"<<IdSiPM	<<endl;						 
  IdSiPM<<=2;//
  // cout<<IdSiPM	<<endl;	
  IdSiPM +=nInLA;      //occupies 2 bits for layer
  // cout<<IdSiPM	<<endl;	
  IdSiPM<<=7;
  // cout<<IdSiPM	<<endl;	
  IdSiPM +=nScntStrp;// 7 bits space strip
  // cout<<IdSiPM	<<endl;	
  IdSiPM <<=2; // Shift just to be compatible with Digi sample, which will include SiPM
  cout<<"New Hit "<<Counter<<" "<<IdSiPM<<" "<<edep<<" "<<loc_no<<" "<<nInLA<<" "<<nScntStrp<<endl;


  int oldCellId = -1;
  for (int ij=0; ij<InCell; ij++) {
    // cout<<"ij "<<ij<<endl;
    if (IdSiPM ==CellDetID[ij]) {oldCellId = ij;//cout<<"oldCellId"<<oldCellId<<endl;
    }
  }
      
  if (oldCellId ==-1 && InCell <numberInCell -1 ) {

 	

	
    Counter++;
    micalcal1Hit* newHit = new micalcal1Hit();
	
    newHit->SetHitId(IdSiPM);
      
    //	pdgid = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    //newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
    newHit->SetpdgId(pdgid);
    newHit->SetEdep(edep);
    newHit->SetTime(atime);
    newHit->SetPos(glbpos);
   newHit->SetLocalPos(localpos);

	
   //  newHit->SetLocalXPos(localpos.x());
   // newHit->SetLocalYPos(localpos.y());		     
   //  newHit->SetLocalZPos(localpos.z());	
   newHit->SetMom( aStep->GetTrack()->GetMomentum());//Get momentum from prestep
      
    InCell = cal1Collection->insert( newHit );
    //  cout<<"InCell "<<InCell<<endl;
    CellDetID[InCell-1] = IdSiPM;      //?
    //  cout<<"CellDetID[InCell-1] "<<CellDetID[InCell-1]<<endl;


  }
  //.......................................................................................................
  cout<<"..ProcessHits ends..."<<endl;


  if (oldCellId >=0) {
    (*cal1Collection)[oldCellId]->AddEdep(edep,glbpos,localpos);
    cout<<"addedep "<<endl;
    cout<<(*cal1Collection)[oldCellId]->GetPos()<< " "<< (*cal1Collection)[oldCellId]->GetLocalPos()<<" "<<(*cal1Collection)[oldCellId]->GetEdep()<<endl;
    if (atime <(*cal1Collection)[oldCellId]->GetTime()) {

      (*cal1Collection)[oldCellId]->SetTime(atime);
    }
  }
    







	
  return true;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalcal1SD::EndOfEvent(G4HCofThisEvent*)
{
  //  double sigspeed = paradef->sigspeed;//16.3cm/ns 
  cout<<"micalcal1SD::EndOfEvent START "<<endl;
  InCell = 0;



  cout<<"hello "<<SipmHit_pointer->SipmHit_list.size()<<" "<<CmvStrip_pointer->CmvStrip_list.size()<<endl;







  
  //	G4float ScintHitGPos[3];
  //	G4float PhyVolGlPos[3];
  if (pAnalysis->InputOutput ==3 || pAnalysis->InputOutput ==4) {
    pAnalysis->inputRootFile->cd();
    gDirectory->pwd();
    if(pAnalysis->FirstEvt+pAnalysis->ievent< pAnalysis->inputEventTree->GetEntries()){ 
      //      pAnalysis->inputEventTree->GetEntry(pAnalysis->FirstEvt+pAnalysis->ievent++); will be done at cal1SD
    } else{
      cout<<"\n Error: Event no. greater than total no. of entries in the input file. \n";
      exit(1);
    }
    for(unsigned int rr1=0; rr1<cal1Collection->entries(); rr1++) {
      cout<<rr1<<" ";
      (*cal1Collection)[rr1]->Print();
    }
		
    cout <<"siminput "<< pAnalysis->cmv_nsimhit<<endl;
    cout<<"Before loop: "<<cal1Collection->entries()<<endl;
    for(unsigned ij=0;ij<pAnalysis->cmv_nsimhit;ij++) {
      micalcal1Hit* newHit = new micalcal1Hit();
          G4ThreeVector mom(pAnalysis->cmv_simpx[ij],pAnalysis->cmv_simpy[ij],pAnalysis->cmv_simpz[ij]);
      //      mom.setMag(pAnalysis->cmv_simmom[ij]);
      //      mom.setTheta(pAnalysis->cmv_simthe[ij]);
      //      mom.setPhi(pAnalysis->cmv_simphi[ij]);
      G4ThreeVector pos(pAnalysis->cmv_simposx[ij],pAnalysis->cmv_simposy[ij],pAnalysis->cmv_simposz[ij]);
      newHit->SetHitId(pAnalysis->cmv_detid[ij]);
      newHit->SetpdgId(pAnalysis->cmv_simpdgid[ij]);
      newHit->SetEdep( pAnalysis->cmv_simenr[ij] );
      newHit->SetTime( pAnalysis->cmv_simtime[ij] );
			
      newHit->SetPos( pos );
      newHit->SetMom(mom);
	  G4ThreeVector localpos(pAnalysis->cmv_simlocx[ij],pAnalysis->cmv_simlocy[ij],pAnalysis->cmv_simlocz[ij]);		
	  newHit->SetLocalPos( localpos );
	  //      newHit->SetLocalXPos(pAnalysis->cmv_simlocx[ij]);
	  //      newHit->SetLocalYPos(pAnalysis->cmv_simlocy[ij]);
	  //        newHit->SetLocalZPos(pAnalysis->cmv_simlocz[ij]);
      cout<<"ij "<<ij<<" "<<pos<<endl;
      // newHit->Print();
      // cout <<"newhits "<< newHit->GetTime()<<endl;
      cal1Collection->insert( newHit );
    }
    cout<<"cal1Collection->size "<<cal1Collection->entries()<<endl;
    for(unsigned int rr1=0; rr1<cal1Collection->entries(); rr1++) {
      cout<<rr1<<" ";
      (*cal1Collection)[rr1]->Print();
    }
    pAnalysis->inputRootFile->cd();
  }

  if (pAnalysis->InputOutput <=4) {
    if (pAnalysis->InputOutput==2 || pAnalysis->InputOutput==0 ||  pAnalysis->InputOutput==1   ) { //gen to sim. 0 added on 17022022
  gDirectory->pwd();
      pAnalysis->pRootFile->cd();
  gDirectory->pwd();
      pAnalysis->cmv_nsimhit = cal1Collection->entries();
      cout<<"pAnalysis->InputOutput==2 or ==0 "<<cal1Collection->entries()<<" "<<pAnalysis->cmv_nsimhit<<" >  " <<pAnalysis->cmv_nsimhtmx<<endl;
      if (pAnalysis->cmv_nsimhit >pAnalysis->cmv_nsimhtmx) pAnalysis->cmv_nsimhit =pAnalysis->cmv_nsimhtmx;
      for (int ij=0; ij< (int)cal1Collection->entries() && ij<(int)pAnalysis->cmv_nsimhit; ij++) {
        pAnalysis->cmv_detid[ij] =  (*cal1Collection)[ij]->GetHitId();// id shifted by 2 bits no sipm no. added..

	cout<< "hit id  " <<  pAnalysis->cmv_detid[ij]<<endl;
        pAnalysis->cmv_simpdgid[ij] =  (*cal1Collection)[ij]->GetpdgId();
        pAnalysis->cmv_simtime[ij] = (*cal1Collection)[ij]->GetTime();
        pAnalysis->cmv_simenr[ij] = (*cal1Collection)[ij]->GetEdep();
        
        G4ThreeVector posvec1 = (*cal1Collection)[ij]->GetPos();
        pAnalysis->cmv_simposx[ij] = posvec1.x();
        pAnalysis->cmv_simposy[ij] = posvec1.y();
        pAnalysis->cmv_simposz[ij] = posvec1.z();
        
        G4ThreeVector momvec = (*cal1Collection)[ij]->GetMom();
        pAnalysis->cmv_simpx[ij] = momvec.mag();
        pAnalysis->cmv_simpy[ij] = momvec.theta();
        pAnalysis->cmv_simpz[ij] = momvec.phi();
        
        pAnalysis->cmv_simlocx[ij] = (*cal1Collection)[ij]->GetLocalXPos();
        pAnalysis->cmv_simlocy[ij] = (*cal1Collection)[ij]->GetLocalYPos();
	   pAnalysis->cmv_simlocz[ij] = (*cal1Collection)[ij]->GetLocalZPos();
        
        if (ij >= (int)pAnalysis->cmv_nsimhtmx) break; ; //redundant

 if(debug)	cout<<"sim data stored....."<<endl;
      }
      //      pAnalysis->pEventTree->Fill();
     if(debug)  cout<<"io: "<<pAnalysis->InputOutput<<endl;
    }
         if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput==1 || pAnalysis->InputOutput==3 || pAnalysis->InputOutput==4  ) { //else of if (pAnalysis->InputOutput==2), 0, 1, 3 & 4 (Sim to Digi step) //earlier only else{ 17032022
      //Digitisation
    // else{
       //  pAnalysis->pRootFile->cd();			
        cout<<"1cal1Collection->entries() "<<cal1Collection->entries()<<endl;
      for (unsigned int ij=0; ij<cal1Collection->entries(); ij++) {
	//				double edep = (*cal1Collection)[ij]->GetEdep();
	G4ThreeVector posvec2 = (*cal1Collection)[ij]->GetPos();
	unsigned long detid = (*cal1Collection)[ij]->GetHitId();
	//	cout<<"Check 1"<<  detid;// This id is just shifted by 2 bits, sipm no. not added
	CmvStrip* strip = new CmvStrip(); //GMA Memory leakages ??

		        
	strip->SetId(detid);
	//
	 unsigned long  Id = detid;
	 Id>>=9;
				
	 int laynoo = Id%4;

	 Id = detid;
	 Id>>=11;
	 int locno = Id%8;
 if(debug)	 cout<<"Fill:: "<<laynoo<< " "<<locno<< endl;
	 // if(laynoo==0 && locno==2 ){
	 //    pAnalysis->hist33->Fill(posvec2.x());
	 //    //	  pAnalysis->hist11->Fill(posvec2.z());
	 // }
	
	//	cout <<"detid "<<detid<<endl;
	strip->SetpdgId((*cal1Collection)[ij]->GetpdgId());
	strip->SetXPos(posvec2.x());
	strip->SetYPos(posvec2.y());
	strip->SetZPos(posvec2.z());
	strip->SetXLocPos((*cal1Collection)[ij]->GetLocalXPos());
	strip->SetYLocPos((*cal1Collection)[ij]->GetLocalYPos());
	strip->SetZLocPos((*cal1Collection)[ij]->GetLocalZPos());
	strip->SetTime((*cal1Collection)[ij]->GetTime());
	strip->SetPulse((*cal1Collection)[ij]->GetEdep());
	G4ThreeVector momvec2 = (*cal1Collection)[ij]->GetMom();
	strip->SetSimMom(momvec2.mag());
	strip->SetSimThe(momvec2.theta());
	strip->SetSimPhi(momvec2.phi());
	if(debug)	cout<<"cmv: strip:"<<endl;
	 strip->Print();
	if(debug)	cout<<endl;
	CmvStrip_pointer->CmvStrip_list.push_back(strip);
	for (int jk=0; jk<4; jk++) {
	  if(debug)	  cout<<"Strip to SiPM "<<endl;
	  SipmHit* sipmht = new SipmHit(strip, jk); //GMA Memory leakage ??
	   if(debug) sipmht->Print();
					
	  int isipmid = sipmht->GetId();
	  if(debug) cout<<"isipmid "<<isipmid<<endl;
	  int oldid=-1;
	  for (unsigned int kl=0; kl<SipmHit_pointer->SipmHit_list.size(); kl++) {
	    if (isipmid==SipmHit_pointer->SipmHit_list[kl]->GetId()) {
	      if(debug) cout<<"check Update"<<endl;
	      SipmHit_pointer->SipmHit_list[kl]->Update(sipmht->GetPulse(), sipmht->GetTime());
	      oldid = kl;
	      break;
	    }
	  }
	  if (oldid<0 && sipmht->GetPulse() > 0) {//0.2pC
	    SipmHit_pointer->SipmHit_list.push_back(sipmht);
	  }
	  //	delete sipmht;	
	}
	//	cout <<"ijcmstrip "<< ij<<" "<<(*cal1Collection)[ij]->GetEdep()<<" "<<SipmHit_pointer->SipmHit_list.size()<<endl;
	//delete strip;
      } // for (int ij=0; ij<cal1Collection->entries(); ij++)
			
    


      /*
      // Addition of noise:
      double noise,tim;
      unsigned int id=0;
      int cmax,bmax;
      
      for(int ia=1; ia<5; ia++){//locno
	id++;
	id<<=2;
				
	cmax = (ia==1) ? NoScntStrpTop : NoScntStrpSide; //88:40
	bmax = (ia==1) ? 4 : 3;
	
	for(int jb =0; jb<bmax; jb++) {//layerno
	  if(jb==0){ id+=0;} else{id++;}
	  id<<=7;
	  
	  for(int kc =0; kc<cmax; kc++){//stripno
	    if(kc==0){id+=0;} else{	id++; }
	    id<<=2;
	    
	    for(int med=0; med<4; med++){//sipmno
	      if(med==0){id+=0;} else{id++;	} 
	      //  cout<<"id " <<id<<"  "<<ia<<"  "<<jb<<"  "<<kc<<"  " <<med<<endl;
	      SipmHit* hit = new SipmHit();
	    
	      noise = pAnalysis->noise_hist[1][2]->GetRandom();
	      tim = 1000*(gRandom->Uniform());
	      //tim = 500*(2*gRandom->Uniform()-1);
	      //	cout<<"noise " <<noise<<" time "<<tim<<" "<<Sipm_Pedestal<<" "<<Sipm_Pedestal+100*noise<<" Getpulse "<<hit->GetPulse()<<endl;
	      hit->SetId(id);
	      
	      hit->SetPulse(max(0.0, Sipm_Pedestal +100*noise)) ;//pedestal = 100 in 0.01pC
	      hit->SetTime(tim);
	      //	cout<<"noisepulse noisetime "<<hit->GetPulse()<<" "<<hit->GetTime()<<endl;
	      int isipmid = hit->GetId();
	      int oldid=-1;
	      for (unsigned int jk=0; jk<SipmHit_pointer->SipmHit_list.size(); jk++) {
		if (isipmid==SipmHit_pointer->SipmHit_list[jk]->GetId()) {
		  
		  cout<<"Adding Noise in the hits: "<<isipmid<<" "<<SipmHit_pointer->SipmHit_list[jk]->GetPulse()<<" "<<SipmHit_pointer->SipmHit_list[jk]->GetTime()<<endl;
		  SipmHit_pointer->SipmHit_list[jk]->Update(hit->GetPulse(), hit->GetTime());
		  cout<<"After Update: "<<SipmHit_pointer->SipmHit_list[jk]->GetPulse()<<" "<<SipmHit_pointer->SipmHit_list[jk]->GetTime()<<endl;
		  oldid = jk;
		  break;
		}
	      }
	      if (oldid<0 && hit->GetPulse() > Cmv_Threshold ) { //120
		cout<<"New hit created only due to Noise: "<<endl;

		
		if(ia==1){

		  hit->SetXPos( PhyVolGlPos[ia-1][jb][0] - 0.5*((2*NoScntStrpTop*partopscint[0])+((NoScntStrpTop+1)*AirGapScintTop)) +(kc+1)*(AirGapScintTop)+(2*kc+1)*partopscint[0]);
		  hit->SetXLocPos( - 0.5*((2*NoScntStrpTop*partopscint[0])+((NoScntStrpTop+1)*AirGapScintTop)) +(kc+1)*(AirGapScintTop)+(2*kc+1)*partopscint[0]);

		  hit->SetYPos(PhyVolGlPos[ia-1][jb][1]);
		  hit->SetYLocPos(0);
		    
		  hit->SetZPos(PhyVolGlPos[ia-1][jb][2]);
	          hit->SetZLocPos(0);
		}//	if(ia==1){



		else	if(ia==2 || ia==3){

		  hit->SetXPos(PhyVolGlPos[ia-1][jb][0]);
        
		  hit->SetXLocPos( - 0.5*((2* NoScntStrpSide*partopscint[0])+((NoScntStrpSide+1)*AirGapScintTop)) +(kc+1)*(AirGapScintTop)+(2*kc+1)*partopscint[0]);
	
		  hit->SetYPos(PhyVolGlPos[ia-1][jb][1]);
		  hit->SetYLocPos(0);

		  hit->SetZPos( PhyVolGlPos[ia-1][jb][2] - 0.5*((2*NoScntStrpSide*partopscint[0])+((NoScntStrpSide+1)*AirGapScintTop)) +(kc+1)*(AirGapScintTop)+(2*kc+1)*partopscint[0]);
		  hit->SetZLocPos(0);


		}//	if(ia==2 || ia==3){


		else	if(ia==4){

		  hit->SetXPos(PhyVolGlPos[ia-1][jb][0]);
		  
		  hit->SetXLocPos(-0.5*((2*NoScntStrpSide*partopscint[0])+((NoScntStrpSide+1)*AirGapScintTop)) +( kc+1)*(AirGapScintTop)+(2*kc+1)*partopscint[0]);


		  hit->SetYLocPos(0);
		  hit->SetYPos(PhyVolGlPos[ia-1][jb][1]);

		  hit->SetZPos( PhyVolGlPos[ia-1][jb][2] -  0.5*((2*NoScntStrpSide*partopscint[0])+((NoScntStrpSide+1)*AirGapScintTop)) +( kc+1)*(AirGapScintTop)+(2*kc+1)*partopscint[0]);
		  hit->SetZLocPos(0);
		  
		}


		  

		
		//	cout<<"oldid=-1  "<<hit->GetPulse()<<" "<<hit->GetTime()<<endl;
        
		SipmHit_pointer->SipmHit_list.push_back(hit);
		hit->Print();

		cout<<hit->GetId()<<" "<<hit->GetStripId()<<" "<<hit->GetPlane()<<" "<<hit->GetLayer()<<" "<<hit->GetSiPM()<<" "<<
		  hit->GetStrip()<<" "<<endl;		
	      }
	      else{ delete hit; } // added to save memory since these hits are not useful
	     
	    }//sipm no
	    id-=3;
	    id>>=2;
	  }//stripno
	  id-=(cmax-1);
	  id>>=7;
	  
	  //cout<<endl;
	}//layerno
	id-=(bmax-1);
	id>>=2;
	//cout<<endl<<endl;
      }//locno
      
              */
    }
     /*
    //delete the entry with pulse(signal + noise) < 12
    cout<<"SipmHit_pointer->SipmHit_list.size() "<<SipmHit_pointer->SipmHit_list.size()<<endl;
    for(int ij=0;ij<SipmHit_pointer->SipmHit_list.size();ij++){
	   
      cout<<SipmHit_pointer->SipmHit_list[ij]->GetId()<<"  "<<SipmHit_pointer->SipmHit_list[ij]->GetPulse()<<endl;
      if(SipmHit_pointer->SipmHit_list[ij]->GetPulse()<12){
	cout<<"/delete the entry with pulse(signal + noise) < 12 "<<endl;
	cout<<SipmHit_pointer->SipmHit_list[ij]->GetId()<<endl;
	SipmHit_pointer->SipmHit_list.erase(SipmHit_pointer->SipmHit_list.begin()+ij);
	       
	ij--;
 
      }
    }

    cout<<"After  SipmHit_pointer->SipmHit_list.size() "<<SipmHit_pointer->SipmHit_list.size()<<endl;
     */

	   if(debug)  cout<<"After  SipmHit_pointer->SipmHit_list.size() "<<SipmHit_pointer->SipmHit_list.size()<<endl;
	    if (pAnalysis->InputOutput==1 || pAnalysis->InputOutput==4  ||   pAnalysis->InputOutput==0 ) { // 0 added on 17022022
      pAnalysis->pRootFile->cd();
      pAnalysis->cmv_ndigihit = SipmHit_pointer->SipmHit_list.size();
     if(debug)  cout<<"  SipmHit_pointer->SipmHit_list.size():  "<<  SipmHit_pointer->SipmHit_list.size()<<endl;
     if(debug)  cout <<"  pAnalysis->cmv_ndigihit "<<  pAnalysis->cmv_ndigihit<<endl;
	    
      if (pAnalysis->cmv_ndigihit >pAnalysis->cmv_ndigihtmx) pAnalysis->cmv_ndigihit =pAnalysis->cmv_ndigihtmx;
      for (unsigned ij=0; ij<SipmHit_pointer->SipmHit_list.size() && ij<pAnalysis->cmv_ndigihit; ij++) {
	pAnalysis->cmv_digipdgid[ij] =SipmHit_pointer->SipmHit_list[ij]->GetpdgId();

 if(debug)	cout<<"sipmid "<<SipmHit_pointer->SipmHit_list[ij]->GetId()<<endl;
        pAnalysis->cmv_sipmid[ij] =SipmHit_pointer->SipmHit_list[ij]->GetId();
 if(debug)	cout<< pAnalysis->cmv_sipmid[ij]<<endl;
	pAnalysis->cmv_digitimpul[ij] =SipmHit_pointer->SipmHit_list[ij]->GetTimePulse();
	pAnalysis->cmv_digitime[ij] =SipmHit_pointer->SipmHit_list[ij]->GetTime();
	pAnalysis->cmv_digipul[ij] =SipmHit_pointer->SipmHit_list[ij]->GetPulse();
	pAnalysis->cmv_digiposx[ij] =SipmHit_pointer->SipmHit_list[ij]->GetXPos();
	pAnalysis->cmv_digiposy[ij] =SipmHit_pointer->SipmHit_list[ij]->GetYPos();
	pAnalysis->cmv_digiposz[ij] =SipmHit_pointer->SipmHit_list[ij]->GetZPos();
	pAnalysis->cmv_digimom[ij] =SipmHit_pointer->SipmHit_list[ij]->GetSimMom();
	pAnalysis->cmv_digithe[ij] =SipmHit_pointer->SipmHit_list[ij]->GetSimThe();
	pAnalysis->cmv_digiphi[ij] =SipmHit_pointer->SipmHit_list[ij]->GetSimPhi();
	pAnalysis->cmv_digilocx[ij] =SipmHit_pointer->SipmHit_list[ij]->GetXLocPos();
	pAnalysis->cmv_digilocy[ij] =SipmHit_pointer->SipmHit_list[ij]->GetYLocPos();
	pAnalysis->cmv_digilocz[ij] =SipmHit_pointer->SipmHit_list[ij]->GetZLocPos();
	
	if (ij >=pAnalysis->cmv_ndigihtmx) break; //redundant
      }
      //			pAnalysis->pEventTree->Fill();
    }
  } else { //if (pAnalysis->InputOutput <=4) : Red digin
    pAnalysis->inputRootFile->cd();
		
    // if(pAnalysis->FirstEvt+pAnalysis->ievent<= pAnalysis->inputEventTree->GetEntries()){    
    //   //      pAnalysis->inputEventTree->GetEntry(pAnalysis->FirstEvt+pAnalysis->ievent++); //Will be done at cal1SD
    // } else {
    //   cout<<"\n Error: Event no. greater than total no. of entries in the input file.\n";
    //   exit(1);
    // }
    cout<<"pAnalysis->cmv_ndigihit "<<pAnalysis->cmv_ndigihit<<endl;
    for(unsigned ij=0;ij<pAnalysis->cmv_ndigihit;ij++) {
      //      unsigned istrp = pAnalysis->stripid[ij];
      SipmHit*  sipmht = new SipmHit(); //VALGRIND
      sipmht->SetId(pAnalysis->cmv_sipmid[ij]);
      sipmht->SetpdgId(pAnalysis->cmv_digipdgid[ij]);
      sipmht->SetTimePulse(pAnalysis->cmv_digitimpul[ij]);
      sipmht->SetXPos(pAnalysis->cmv_digiposx[ij]);
      sipmht->SetYPos(pAnalysis->cmv_digiposy[ij]);
      sipmht->SetZPos(pAnalysis->cmv_digiposz[ij]);
      sipmht->SetSimMom(pAnalysis->cmv_digimom[ij]);
      sipmht->SetSimThe(pAnalysis->cmv_digithe[ij]);
      sipmht->SetSimPhi(pAnalysis->cmv_digiphi[ij]);
      sipmht->SetXLocPos(pAnalysis->cmv_digilocx[ij]);
      sipmht->SetYLocPos(pAnalysis->cmv_digilocy[ij]);
    sipmht->SetZLocPos(pAnalysis->cmv_digilocz[ij]);
      int isipmid = sipmht->GetId();
      int oldid=-1;
      cout<<"check11 "<<SipmHit_pointer->SipmHit_list.size()<<endl;
      for (unsigned int jk=0; jk<SipmHit_pointer->SipmHit_list.size(); jk++) {
	cout<<"jk "<<jk<<endl;
	if (isipmid==SipmHit_pointer->SipmHit_list[jk]->GetId()) {
	  SipmHit_pointer->SipmHit_list[jk]->Update(sipmht->GetPulse(), sipmht->GetTime());
	  oldid = jk;
	  break;
	}
      }
      if (oldid<0 && sipmht->GetPulse() > 0.0001) {
	//	sipmht->Print();
	SipmHit_pointer->SipmHit_list.push_back(sipmht);
      }
      //	delete sipmht;	
    }
    pAnalysis->pRootFile->cd();
  } //if (pAnalysis->InputOutput <=4)

  /*
  // Convert SiPM hit to Cmv Hit
  SipmHit* foursipm[4]={0}; //GMA memory leakage ?
  if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput ==3 || pAnalysis->InputOutput==5) {
  for (int ij=0; ij<SipmHit_pointer->SipmHit_list.size(); ij++) {
  int tmpstripid = -1;
  int tmpside = -1; //Used this to find global position of the layer
  if (!(SipmHit_pointer->SipmHit_list[ij]->isUsed())) {
  tmpstripid = SipmHit_pointer->SipmHit_list[ij]->GetStripId();
  foursipm[SipmHit_pointer->SipmHit_list[ij]->GetSiPM()] = SipmHit_pointer->SipmHit_list[ij];
  SipmHit_pointer->SipmHit_list[ij]->SetUsed(true);
				
  tmpside = SipmHit_pointer->SipmHit_list[ij]->GetPlane()-1; //We had added 1 while storing it.

  //Look for all other SiPM of same strip
  for (int jk=ij+1; jk<SipmHit_pointer->SipmHit_list.size(); jk++) {
  if (!(SipmHit_pointer->SipmHit_list[jk]->isUsed())) {	
  int tmpstripid2 = SipmHit_pointer->SipmHit_list[jk]->GetStripId();
  if (tmpstripid !=tmpstripid2) continue;
  foursipm[SipmHit_pointer->SipmHit_list[jk]->GetSiPM()] = SipmHit_pointer->SipmHit_list[jk];
  SipmHit_pointer->SipmHit_list[jk]->SetUsed(true);
  }
  }
  }
  if (tmpside>=0) {
  CmvHit* tmpcmvHit = new CmvHit(foursipm[0], foursipm[1], foursipm[2], foursipm[3], 	PhyVolGlPos[tmpside][0]);
				
  }
  }
  */
	
  //   G4int nofHits = cal1Collection->entries();
  //   //G4cout << G4endl
  //   //<< "-------->Hits Collection: in this event there are " << nofHits 
  //   //<< " hits in the tracker chambers: " << G4endl;
  //   //for ( G4int i=0; i<nofHits; i++ ) (*cal1Collection)[i]->Print();
  //   cout<<"nofHits::"<<nofHits<<endl;
  //   pAnalysis->CMV_nsimhit = cal1Collection->entries();
  //   //cout<<"check-1"<<endl;
  //   for (int ij=0; ij<cal1Collection->entries(); ij++) {
  //     // cout<<"ij"<<ij<<endl;
  //     unsigned int SiPMId = (*cal1Collection)[ij]->GetHitId();
  //     G4ThreeVector posvec = (*cal1Collection)[ij]->GetPos();  // get glb poistion
  //     int pdgid = (*cal1Collection)[ij]->GetpdgId();
  //     double atime = (*cal1Collection)[ij]->GetTime();          //detid is stored in Sethitid
  //     double Edep = (*cal1Collection)[ij]->GetEdep(); 
  //     G4ThreeVector localpos = (*cal1Collection)[ij]->GetLocalPos();
    
    
  //     cout<<"ij "<<ij<<" Edep "<<Edep<<" pdgid "<<pdgid<<" posvec "<<posvec<<" SiPMId  "<<SiPMId <<endl;
  //     //cout<<"check1"<<endl;
    
  //     pAnalysis->CMV_detid[ij] = SiPMId;
  //     pAnalysis->CMV_simpdgid[ij] =  pdgid;
  //     pAnalysis->CMV_simtime[ij] = atime;
  //     pAnalysis->CMV_simenr[ij] = Edep;
  //     pAnalysis->CMV_simposx[ij] = posvec.x();
  //     pAnalysis->CMV_simposy[ij] = posvec.y();
  //     pAnalysis->CMV_simposz[ij] = posvec.z();
    
    
  //     int LayerNo, loc_no, ScntStrpNo, SiPMNo;
  //     //cout<<"ScntStrpid"<<ScntStrpid<<endl;
    
  //     SiPMNo = SiPMId%4; //2^2
  //     SiPMId>>=2;
  //     ScntStrpNo= SiPMId%128; //2^7
  //     // cout<<"ScntStrpNo"<<ScntStrpNo<<endl;
  //     SiPMId>>=7;
  //     LayerNo = SiPMId%4;
  //     // cout<<"LayerNo"<<LayerNo<<endl;
  //     SiPMId>>=2;      
  //     loc_no = SiPMId;
  //     //cout<<"Loc_No"<<loc_no<<endl;
    
  //     pAnalysis->CMV_digiSiPMNo[ij] =SiPMNo;
  //     pAnalysis->CMV_digiScntStrpNo[ij] = ScntStrpNo;
  //     pAnalysis->CMV_digiLayerNo[ij] =  LayerNo; 
  //     pAnalysis->CMV_digiLocNo[ij] =	loc_no;
    
  //     cout<<"local position "<<localpos<<endl;
  //     cout<<" SiPMNo "<<SiPMNo <<" ScntStrpNo  "<<ScntStrpNo <<" LayerNo  "<<LayerNo <<" Loc_No "<<loc_no<<endl;
   
  //     double sigpos;
  //     if(loc_no==3) {
  //       sigpos=localpos.x();
  //       // cout<<"signal position "<<sigpos<<endl;
  //     } else {
  //       sigpos = localpos.y();
  //       // cout<<"signal position "<<sigpos<<endl;
  //     }
    
  //     double sigtim = atime +  (partopscint[1]+ sigpos)/sigspeed + gRandom->Gaus(0,1)*ns;
  //     double  sigtimdash = atime + (partopscint[1]- sigpos)/sigspeed + gRandom->Gaus(0,1)*ns ;
    
  //     if(SiPMNo ==0 || SiPMNo==1) {
  //       cout<<atime<<" signal time "<<sigtim/ns<<endl;
  //       pAnalysis->CMV_digisigtim[ij]= sigtim;
  //     } else {
  //       cout<<atime<<" signal time "<<sigtimdash/ns<<endl;
  //       pAnalysis->CMV_digisigtim[ij]=sigtimdash; //GMA it was an error
  //     }
    
  //     if ( loc_no == 1){
  //       PhyVolGlPos[0]= paradef->Phys_TopScint_GPos[LayerNo][0];
  //       PhyVolGlPos[1]=paradef->Phys_TopScint_GPos[LayerNo][1];
  //       PhyVolGlPos[2]= paradef->Phys_TopScint_GPos[LayerNo][2];
  //     } else if (loc_no == 2){
  //       PhyVolGlPos[0]= paradef->Phys_SideScint_L_GPos[LayerNo][0];
  //       PhyVolGlPos[1]=paradef->Phys_SideScint_L_GPos[LayerNo][1];
  //       PhyVolGlPos[2]= paradef->Phys_SideScint_L_GPos[LayerNo][2];
  //     } else if (loc_no == 3){
  //       PhyVolGlPos[0]= paradef->Phys_SideScint_R_GPos[LayerNo][0];
  //       PhyVolGlPos[1]=paradef->Phys_SideScint_R_GPos[LayerNo][1];
  //       PhyVolGlPos[2]= paradef->Phys_SideScint_R_GPos[LayerNo][2];
  //     } else if (loc_no == 4){
  //       PhyVolGlPos[1]= paradef->Phys_SideScint_D_GPos[LayerNo][1];
  //       PhyVolGlPos[0]=paradef->Phys_SideScint_D_GPos[LayerNo][0];
  //       PhyVolGlPos[2]= paradef->Phys_SideScint_D_GPos[LayerNo][2];
  //     }
    
  //     //cout<<"loc_no " << loc_no<<endl;
    
    
  //     // cout<<" ---------"<<PhyVolGlPos[0]<<endl;
    
  //     if(loc_no == 1){
      
  //       ScintHitGPos[0] = PhyVolGlPos[0] - 0.5*((88*2*partopscint[0])+(89*AirGapScintTop)) +( ScntStrpNo+1)*(AirGapScintTop)+(2*ScntStrpNo+1)*partopscint[0];
  //       ScintHitGPos[1] =PhyVolGlPos[1];
  //       ScintHitGPos[2] = PhyVolGlPos[2];
  //       // pAnalysis->ScntStrpXPos[LayerNo].push_back(ScintHitGPos[0]);
  //       // pAnalysis->ScntStrpNo[LayerNo].push_back(ScntStrpNo);
      
  //       cout<<" ---------"<<ScintHitGPos[0]<<endl;
  //     } else if(loc_no == 2 ||loc_no == 3 ){
  //       ScintHitGPos[0] = PhyVolGlPos[0] ;
  //       ScintHitGPos[1]= PhyVolGlPos[1];
  //       ScintHitGPos[2] = PhyVolGlPos[2]-0.5*((40*2*partopscint[0])+(41*AirGapScintTop))+( ScntStrpNo+1)*(AirGapScintTop)+(2*ScntStrpNo+1)*partopscint[0];
      
  //     }else if(loc_no == 4){
  //       ScintHitGPos[0]=PhyVolGlPos[0];
  //       ScintHitGPos[1] = PhyVolGlPos[1] ;
  //       ScintHitGPos[2] = PhyVolGlPos[2]-0.5*((40*2*partopscint[0])+(41*AirGapScintTop))+( ScntStrpNo+1)*(AirGapScintTop)+(2*ScntStrpNo+1)*partopscint[0];
      
  //     }
 
  //     pAnalysis->CMV_digiposx[ij] =	ScintHitGPos[0]; 
  //     pAnalysis->CMV_digiposy[ij] =	ScintHitGPos[1]; 
  //     pAnalysis->CMV_digiposz[ij] =	ScintHitGPos[2]; 	
  //   }//for ij
  
  cout<<"micalcal1SD::EndOfEvent ends "<<endl;


 

 
}

void micalcal1SD::SetPhotonSpeed(G4double val) {
  cout<<"void micalcal1SD::SetPhotonSpeed(G4double "<<val<<")"<<endl;
  // cout<<"SignalSpeed = "<<val<<endl;
  // cout<<"...}"<<endl;
  //  SignalSpeed = val;
  pAnalysis->SetPhotonSpeedVal(val);
}

void micalcal1SD::SetCMVadctons(G4double val) {
  cout<<"void micalcal1SD::SetCMVadctons(G4double "<<val<<")"<<endl;
  // cout<<"SignalSpeed = "<<val<<endl;
  // cout<<"...}"<<endl;
  //  SignalSpeed = val;
  pAnalysis->SetCMVadctons(val);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
