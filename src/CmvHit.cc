//12/07/09 Put position of hit in this class from the strip information
// Where should we use mapping ? Electronic channel to geometrical channels 

#include "CmvHit.h"
//#include "Validity/VldContext.h"
#include <cmath>
#include "TMath.h"
#include <iostream>
using namespace std;

CmvHit::CmvHit() :
  fLeMiSipm(0), fLePlSipm(0), fRiMiSipm(0), fRiPlSipm(0),
  fStripId(-100),
  fXTruePosX(-1000), fYTruePosY(-1000), fZTruePosZ(-1000),
  fXRecoPosX(-1000), fYRecoPosY(-1000), fZRecoPosZ(-1000),
  fLeTime (-1000),  fRiTime(-1000),
  fXPosErr(-1000),  fYPosErr(-1000), fZPosErr(-1000),
  fMomentum(-1000), fTheta(-1000),  fPhi(-1000), PhotonSpeed(162.0), layergap(0.02),
  isUsed(0)
{; }

CmvHit::CmvHit(CmvHit* fd) {

  fLeMiSipm = fd->fLeMiSipm;   //Strip for L-Y-SiPM 
  fLePlSipm = fd->fLePlSipm;   //Strip for L+Y-SiPM
  fRiMiSipm = fd->fRiMiSipm;   //Strip for R-Y-SiPM 
  fRiPlSipm = fd->fRiPlSipm;   //Strip for R+Y-SiPM

  fStripId  = fd->fStripId; // All four SiPM should have come fro mthe same StripID

  fXTruePosX = fd->fXTruePosX;      //Global X-position from Simulation
  fYTruePosY = fd->fYTruePosY;      //Global Y-position from Simulation
  fZTruePosZ = fd->fZTruePosZ;      //Global Z-position from Simulation

  fXRecoPosX = fd->fXRecoPosX;      //Global X-position from Reconstruction
  fYRecoPosY = fd->fYRecoPosY;      //Global Y-position from Reconstruction
  fZRecoPosZ = fd->fZRecoPosZ;      //Global Z-position from Reconstruction

  fLeTime = fd->fLeTime;      // Timing of Left-SiPMs (ns)
  fRiTime = fd->fRiTime;      // Timing of Right-SiPMs (ns)

  fLePulse = fd->fLePulse;      // Timing of Left-SiPMs (ns)
  fRiPulse = fd->fRiPulse;      // Timing of Right-SiPMs (ns)
	
  fXPosErr = fd->fXPosErr; //Error in X-position (m)
  fYPosErr = fd->fYPosErr; //Error in Y-position (m)
  fZPosErr = fd->fZPosErr; //Error in Y-position (m)
	
  fMomentum = fd->fMomentum; //Momentum of track which is behind this hit (use only one SiPM)
  fTheta = fd->fTheta;  //Theta
  fPhi = fd->fPhi;    // Phi
  PhotonSpeed = fd->PhotonSpeed;
  layergap = fd->layergap;
  isUsed    = fd->isUsed;
	
}

CmvHit::CmvHit(SipmHit* L0,  SipmHit* L1, SipmHit* R0, SipmHit* R1, double* posvec) :
  fLeMiSipm(0), fLePlSipm(0), fRiMiSipm(0), fRiPlSipm(0),
  fStripId(-100),
  fXTruePosX(-1000), fYTruePosY(-1000), fZTruePosZ(-1000),
  fXRecoPosX(-1000), fYRecoPosY(-1000), fZRecoPosZ(-1000),
  fLeTime (-1000),  fRiTime(-1000),
  fLePulse (-1000),  fRiPulse(-1000),
  fXPosErr(-1000),  fYPosErr(-1000), fZPosErr(-1000),
  fMomentum(-1000), fTheta(-1000),  fPhi(-1000), PhotonSpeed(162.0), layergap(0.02), isUsed(0)
{

  paradef = micalDetectorParameterDef::AnPointer; //AAR:
  pAnalysis = MultiSimAnalysis::AnPointer;
  PhotonSpeed = pAnalysis->GetPhotonSpeedVal();
  double CMVadctons = pAnalysis->GetCMVadctons();
  //	CMVadctons = 1;
  double pulsx[4]={-1,-1,-1,-1};
  double timex[4]={-1,-1,-1,-1};
  int index[4]={-1,-1,-1,-1}; 
  int nvalidHit=0;
  double pos[3]={0,0,0};

  if(debug)	cout<<"debug "<<posvec[2]<<endl;
	  
	
  //	cout<< L0->GetStripId()<<" "<< L1->GetStripId()<<" "<< R0->GetStripId()<<" "<< R1->GetStripId()<<" "<<endl;
	
  fStripId=-1;
  if (L0 && L0->GetStripId()) {//##
    if(debug) cout<<"L0 "<< L0->GetStrip()<<endl;
    fLeMiSipm = L0;
    int id=L0->GetSiPM();
    pulsx[id] = L0->GetPulse();
    timex[id] =  L0->GetTime();
	
    index[id] = id;
    if (fStripId<0) {
      fStripId = L0->GetStripId(); //((L0->GetStripId()>>2)<<2);
      fMomentum = L0->GetSimMom();
      fTheta = L0->GetSimThe();
      fPhi = L0->GetSimPhi();
    }

    if(debug)	cout<<"zpos before add"<<L0->GetZPos()<<" "<<L0->GetStripId()<<endl;
    pos[0] += L0->GetXPos();
    pos[1] += L0->GetYPos();
    pos[2] += L0->GetZPos();
    if(debug)	cout<<"L0"<<endl;
    for(int ii=0;ii<3;ii++){	if(debug)	cout<<""<<pos[ii]<<" ";}
    if(debug)	cout<<endl;
    nvalidHit++;
  }
	
  if (L1 && L1->GetStripId()) {

    if(debug)  cout<<"L1 "<< L1->GetStrip()<<endl;
    fLePlSipm = L1;
    int id=L1->GetSiPM();
    pulsx[id] = L1->GetPulse();
    timex[id] =  L1->GetTime();
    index[id]=id;
    if (fStripId<0) {
      fStripId = L1->GetStripId(); //((L1->GetStripId()>>2)<<2);
      fMomentum = L1->GetSimMom();
      fTheta = L1->GetSimThe();
      fPhi = L1->GetSimPhi();
    }

    if(debug)	cout<<"zpos before add"<<L1->GetZPos()<<" "<<L1->GetStripId()<<endl;
    pos[0] += L1->GetXPos();
    pos[1] += L1->GetYPos();
    pos[2] += L1->GetZPos();


    if(debug)	cout<<"L1"<<endl;
    for(int ii=0;ii<3;ii++){	if(debug)	cout<<""<<pos[ii]<<" ";}
    if(debug)	cout<<endl;

    nvalidHit++;
  }
	
  if (R0 && R0->GetStrip()) {
    if(debug)  cout<<"R0 "<< R0->GetStripId()<<endl;
    fRiMiSipm=R0;
    int id=R0->GetSiPM();
    pulsx[id] = R0->GetPulse();
    timex[id] =  R0->GetTime();
    index[id]=id;
    if (fStripId<0) {
      fStripId = R0->GetStripId(); //((R0->GetStripId()>>2)<<2);
      fMomentum = R0->GetSimMom();
      fTheta = R0->GetSimThe();
      fPhi = R0->GetSimPhi();
    }

    if(debug)	cout<<"zpos before add"<<R0->GetZPos()<<" "<<R0->GetStripId()<<endl;
    pos[0] += R0->GetXPos();
    pos[1] += R0->GetYPos();
    pos[2] += R0->GetZPos();

    if(debug)	cout<<"R0"<<endl;
    for(int ii=0;ii<3;ii++){	cout<<""<<pos[ii]<<" ";}
    if(debug)	cout<<endl;

    nvalidHit++;
  }
	
  if (R1 && R1->GetStripId()) {

    if(debug)  cout<<"R1 "<< R1->GetStrip()<<endl;
    fRiPlSipm = R1;
    int id=R1->GetSiPM();
    pulsx[id] = R1->GetPulse();
    timex[id] =  R1->GetTime();
    index[id]=id;
    if (fStripId<0) {
      fStripId = R1->GetStripId(); //((R1->GetStripId()>>2)<<2);
      fMomentum = R1->GetSimMom();
      fTheta = R1->GetSimThe();
      fPhi = R1->GetSimPhi();
    }

    if(debug)	cout<<"zpos before add"<<R1->GetZPos()<<" "<<R1->GetStripId()<<endl;
    pos[0] += R1->GetXPos();
    pos[1] += R1->GetYPos();
    pos[2] += R1->GetZPos();


    if(debug)	cout<<"R1"<<endl;
    for(int ii=0;ii<3;ii++){	if(debug)	cout<<""<<pos[ii]<<" ";}
    if(debug)	cout<<endl;

    nvalidHit++;
  }
  double flefact = 1;
  fLePulse=-10;
  fLeTime=-1000;
  if (index[0]>=0 && index[1]>=0) {
    fLePulse = 0.5*(pulsx[0]+pulsx[1]); fLeTime = 0.5*(timex[0]+timex[1]);
    flefact = 0.707;
  } else if (index[0]>=0) {
    fLePulse = pulsx[0]; fLeTime = timex[0];
		
  } else if (index[1]>=0) {
    fLePulse = pulsx[1]; fLeTime = timex[1];
  }

  double frifact= 1;
  fRiPulse=-10;
  fRiTime=-1000;
  if (index[2]>=0 && index[3]>=0) {
    fRiPulse = 0.5*(pulsx[2]+pulsx[3]); fRiTime = 0.5*(timex[2]+timex[3]);
    frifact = 0.707;
  } else if (index[2]>=0) {
    fRiPulse = pulsx[2]; fRiTime = timex[2];
  } else if (index[3]>=0) {
    fRiPulse = pulsx[3]; fRiTime = timex[3];
  }
  if (nvalidHit<1) {nvalidHit=1;}
  fXTruePosX =pos[0]/nvalidHit;
  fYTruePosY =pos[1]/nvalidHit;
  fZTruePosZ =pos[2]/nvalidHit;


  if(debug)	cout<<"true pos debug "<<endl;
		
  if(debug)	cout<<fZTruePosZ<<" "<<nvalidHit<<endl;
	
  if(debug)	cout<<"iplane"<<fStripId<<endl;
  int iplane = (fStripId>>9)&0x07;
  int ilayer =  (fStripId>>7)&0x03;
  int istrip = fStripId&0x7f;
  int NoofEPSinTile = paradef->GetNoofEPSinTile();
  int tileno= istrip/NoofEPSinTile;

  double scintwidth =  paradef->partopscint[0];
  double scintheight  =  paradef->partopscint[2];//0.5cm halflength

  double TopPlaneHalfLength = paradef->TopPlaneHalfLength;//#R
  double SidePlaneHalfLength = paradef->SidePlaneHalfLength;//#R
  double  SideSmallPlaneHalfLength = paradef->GetSideSmallPlaneHalfLength();
  double GapBtwTiles = paradef->GetGapBtwTiles();
  if(debug)	cout <<"pulse "<< fLePulse <<" "<<	fLeTime<<" "<<fRiPulse <<" "<<	fRiTime<<" "<<iplane<<" "<<istrip<<" "<<scintwidth<<endl;


  bool ispos = (timex[0]>=0 || timex[1]>=0 ) && (timex[2]>=0 || timex[3]>=0);


	
  if(iplane == 1){
    if(debug)	  cout<<"iplane is 1"<<endl;
    if(ilayer==0 || ilayer==2){
      scintheight = 9;//2 cm 1 cm 2 cm 1cm
    }
	 
    //	fXRecoPosX = posvec[0] -TopPlaneHalfLength +(2* istrip)*(paradef->AirGapScintTop)+(2*istrip+1)*scintwidth;


	  
    fXRecoPosX = posvec[0] -TopPlaneHalfLength  +tileno*2*GapBtwTiles + (NoofEPSinTile-1)*tileno* (2*paradef->AirGapScintTop) + (istrip-(NoofEPSinTile*tileno))*(2*paradef->AirGapScintTop)  +(2*istrip+1)*scintwidth;

	  
    fYRecoPosY = posvec[1] +0.5* PhotonSpeed*(fLeTime - fRiTime)*CMVadctons; //Check sign of second term #R

    if(debug) 	cout <<"fYRecoPosY "<<fYTruePosY<<" "<<fYRecoPosY<<" "<< posvec[1]<<" "<< PhotonSpeed<<" "<<fLeTime <<" "<< fRiTime<<endl;

    if(debug)	cout <<"fXRecoPosX "<<fXTruePosX<<" "<<fXRecoPosX<<" "<< posvec[0]<<" "<< TopPlaneHalfLength<<" "<<tileno <<" "<< istrip<<endl<<" "<<NoofEPSinTile<<" " <<GapBtwTiles<<" "<<paradef->AirGapScintTop<<endl;
		
		
    fZRecoPosZ = posvec[2]; // + (ilayer - 1.5)*layergap;
    if(debug)	cout<<"reco pos z:"<<fZRecoPosZ<<" "<<posvec[2]<<endl;
    fXPosErr = 2*scintwidth/sqrt(12); // 2/sqrt(12)
    if(ispos){
      double err1 = flefact*0.707* exp(1.151-0.000155*(fLePulse));
      double err2 = frifact*0.707* exp(1.151-0.000155*(fRiPulse));
      fYPosErr = pow(err1*err1+err2*err2,0.5);
    }else{
		  
		  
		  
      fYPosErr =  -100; //1p.e. = 0.2pC pulse in 10 fC
    }
    fZPosErr = 2*scintheight/sqrt(12);	//2/sqrt(12)																								
		
  } else if(iplane == 2 || iplane == 3 ){ //Sides
    fXRecoPosX = posvec[0]; // + (ilayer - 1)*layergap;
    fYRecoPosY = posvec[1] +0.5* PhotonSpeed*(fLeTime - fRiTime)*CMVadctons; //Check sign of second term #R

	  
    fZRecoPosZ = posvec[2] -SidePlaneHalfLength  +tileno*2*GapBtwTiles + (NoofEPSinTile-1)*tileno* (2*paradef->AirGapScintTop) + (istrip-(NoofEPSinTile*tileno))*(2*paradef->AirGapScintTop)  +(2*istrip+1)*scintwidth;

	  
    // fZRecoPosZ = posvec[2]-SidePlaneHalfLength+( 2*istrip+1)*(paradef->AirGapScintTop)+(2*istrip+1)*scintwidth;
    if(debug)  cout<<"LeftposZ:"<<"  "<<fZRecoPosZ<<"  "<<SidePlaneHalfLength<<"   "<<istrip<<"    "<<paradef->AirGapScintTop<<"   "<<scintwidth<<endl;
    fXPosErr =2*scintheight/sqrt(12) ; 
    if(ispos){
      double err1 = flefact*0.707* exp(1.151-0.000155*(fLePulse));
      double err2 = frifact*0.707* exp(1.151-0.000155*(fRiPulse));
      fYPosErr = pow(err1*err1+err2*err2,0.5);
    }else{
	    
	    
	    
      fYPosErr =  -100; //1p.e. = 0.2pC pulse in 10 fC
    }
	  
	  
	
    fZPosErr =2* scintwidth*sqrt(12); // 
	  
  } else if(iplane == 4){ //Back

    if(debug)  cout<<"Back_uniform_smr_err"<<fLeTime<<" "<<fRiTime<<" "<<fLeTime - fRiTime<<endl;
    if(debug)  cout<<posvec[0]<<" "<<posvec[1]<<" "<<posvec[2]<<endl;
    fXRecoPosX = posvec[0] +0.5* PhotonSpeed*(fLeTime - fRiTime)*CMVadctons; //Check sign of second term; #R
	
    fYRecoPosY = posvec[1]; // + (ilayer - 1)*layergap;
    //	  fZRecoPosZ = posvec[2]-SidePlaneHalfLength+( 2*istrip+1)*(paradef->AirGapScintTop)+(2*istrip+1)*scintwidth;
          
    fZRecoPosZ = posvec[2] -SidePlaneHalfLength  +tileno*2*GapBtwTiles + (NoofEPSinTile-1)*tileno* (2*paradef->AirGapScintTop) + (istrip-(NoofEPSinTile*tileno))*(2*paradef->AirGapScintTop)  +(2*istrip+1)*scintwidth;
    if(ispos){
      double err1 = flefact*0.707* exp(1.151-0.000155*(fLePulse));
      double err2 = frifact*0.707* exp(1.151-0.000155*(fRiPulse));
      fXPosErr = pow(err1*err1+err2*err2,0.5);
    }else{
	  
	  
	  
      fXPosErr =  -100; //1p.e. = 0.2pC pulse in 10 fC
    }
	
    fYPosErr = 2*scintheight/sqrt(12);//
    fZPosErr =2* scintwidth/sqrt(12); // 
  }

  else if(iplane==6 || iplane==7){ //miniLeft miniRight
    scintheight = 9;//2 cm
    if(debug)	  cout<<"mini_uniform_smr_err"<<fLeTime<<" "<<fRiTime<<" "<<fLeTime - fRiTime<<endl;
    if(debug)  cout<<posvec[0]<<" "<<posvec[1]<<" "<<posvec[2]<<endl;
    fZRecoPosZ = posvec[2] +0.5* PhotonSpeed*(fLeTime - fRiTime)*CMVadctons; //Check sign of second term; #R
	
    fXRecoPosX = posvec[0]; // + (ilayer - 1)*layergap;
    //  fYRecoPosY = posvec[1]-SideSmallPlaneHalfLength+( 2*istrip+1)*(paradef->AirGapScintTop)+(2*istrip+1)*scintwidth;
	    
    fYRecoPosY = posvec[1] -SideSmallPlaneHalfLength  +tileno*2*GapBtwTiles + (NoofEPSinTile-1)*tileno* (2*paradef->AirGapScintTop) + (istrip-(NoofEPSinTile*tileno))*(2*paradef->AirGapScintTop)  +(2*istrip+1)*scintwidth;


    if(debug)	cout <<"fYRecoPosY "<<fYTruePosY<<" "<<fYRecoPosY<<" "<< posvec[1]<<" "<< SideSmallPlaneHalfLength<<" "<<tileno <<" "<< istrip<<" "<<NoofEPSinTile<<" " <<GapBtwTiles<<" "<<paradef->AirGapScintTop<<endl;





    if(ispos){
      double err1 = flefact*0.707* exp(1.151-0.000155*(fLePulse));
      double err2 = frifact*0.707* exp(1.151-0.000155*(fRiPulse));
      fZPosErr = pow(err1*err1+err2*err2,0.5);
    }else{
	  
	  
	  
      fXPosErr =  -100; //1p.e. = 0.2pC pulse in 10 fC
    }
	
    fXPosErr = 2*scintheight/sqrt(12);//
    fYPosErr =2* scintwidth/sqrt(12); // 
  }





	

  isUsed = false;	
}

CmvHit::~CmvHit()
{
}

// bool CmvHit::isIdentical(CmvHit* hit) {

//   if ((hit->GetXStripNum() == GetXStripNum()) &&
//       (hit->GetYStripNum() == GetYStripNum())) {
//     return true;
//   }
//   return false;
// }

void CmvHit::Print() {
  //	cout<<"----------------------------------------------------------------------"<<endl;
  //	cout<<"Hit combination "<<endl;
  cout<< "CmvHits():" 
    // <<std::setw(4) <<jk <<" "
      << "Idd="<<std::setw(3)<<GetStripId()
      << " Plane "<< std::setw(3)<<  GetPlane()
      << " Layer "<< std::setw(3)<<  GetLayer()
      << " Strip No "	<< std::setw(3)<<  GetStrip()
      << " X_Reco="   <<std::setw(8)<<  GetRecoPosX() 
      << " Y_Reco=" <<std::setw(8)<<  GetRecoPosY()
      << " Z_reco=" <<std::setw(8)<<  GetRecoPosZ()
      << " X_Pos=" <<std::setw(8)<<  GetTruePosX()
      << " Y_Pos=" <<std::setw(8)<<  GetTruePosY()
      << " Z_Pos=" <<std::setw(8)<<  GetTruePosZ()
      << " chg="   <<std::setw(8)<<  GetPulse()
     

      << endl;
  //  cout<<"......................................................................"<<endl;
}

// void CmvHit::Print() {
//   cout<<"----------------------------------------------------------------------"<<endl;



//   cout<<"......................................................................"<<endl;


// }
