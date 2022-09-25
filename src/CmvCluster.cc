
#include "CmvCluster.h"
#include "CmvHit.h"
//#include "Validity/VldContext.h"
#include <cmath>
#include "TMath.h"
#include <iostream>
using namespace std;
/*
CmvCluster::CmvCluster() :
  clustersize(-1),  isUsed(0),
  plane(-1), layer(-1), strip(-1),
  fXTruePosX(-1000), fYTruePosY(-1000), fZTruePosZ(-1000),
  fXRecoPosX(-1000), fYRecoPosY(-1000), fZRecoPosZ(-1000),
  fLeTime (-1000),  fRiTime(-1000),
  fXPosErr(-1000),  fYPosErr(-1000), fZPosErr(-1000),
  fMomentum(-1000), fTheta(-1000),  fPhi(-1000), PhotonSpeed(162.0), layergap(0.02)
{; }

CmvCluster::CmvCluster(CmvCluster* fd) {
  plane = fd->plane;
  layer = fd->layer;
  strip = fd->strip;
  clustersize = fd->clustersize;
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
	
  fMomentum = fd->fMomentum; //Momentum of track which is behind this cluster (use only one SiPM)
  fTheta = fd->fTheta;  //Theta
  fPhi = fd->fPhi;    // Phi
  PhotonSpeed = fd->PhotonSpeed;
  layergap = fd->layergap;

  isUsed    = fd->isUsed;
}



CmvCluster::CmvCluster(CmvHit* L0, CmvHit* L1):
  clustersize(-1),
  plane(-1), layer(-1), strip(-1),
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
  // double CMVadctons = pAnalysis->GetCMVadctons();
  // double pos[3]={0,0,0};
   double pulse[2]={0};



  
  if (L0 && L1 &&  abs(L0->GetStrip() -  L1->GetStrip())==1 ) {//##
    cout<<"Found Two neighouring strip hits "<<endl;
    cout<<"For L0"<<endl;
    plane = L0->GetPlane();
    layer = L0->GetLayer();
    strip = L0->GetStrip();
    clustersize = 2;

    
    int plane2 = L1->GetPlane();
    int layer2 = L1->GetLayer();
    int strip2 = L1->GetStrip();
    
    cout<<" plane1 "<<plane<<" layer1 "<<layer<<" strip1 "<<strip<<endl;

 
    cout<<"L0 Pos: "<<L0->GetRecoPosX()<<" "<<L0->GetRecoPosY()<<" "<<L0->GetRecoPosZ()<<" "<<endl;
  
    cout<<"L0 Pos Err "<<L0->GetPosXErr()<<" "<<L0->GetPosYErr()<<" "<<L0->GetPosYErr()<<endl;
    cout<<"L0 Pulse: "<<L0->GetLePulse()<<" "<<L0->GetRiPulse()<<endl;
  
	  
     pulse[0]=(L0->GetLePulse() + L0->GetRiPulse())/2;
     pulse[1]=(L1->GetLePulse() + L1->GetRiPulse())/2;

    //  cout<<"pulse "<<pulse[0]<<" "<<pulse[1]<<endl;
    
    cout<<"clustersize "<<clustersize<<endl;

    cout<<endl<<"For L1 "<<endl;

    cout<<" plane2 "<<plane2<<" layer2 "<<layer2<<" strip2 "<<strip2<<endl;
    cout<<"L1 Pos: "<<L1->GetRecoPosX()<<" "<<L1->GetRecoPosY()<<" "<<L1->GetRecoPosZ()<<" "<<endl;
  cout<<"L1 Pos Err "<<L1->GetPosXErr()<<" "<<L1->GetPosYErr()<<" "<<L1->GetPosYErr()<<endl;
    
    cout<<"L1 Pulse: "<<L1->GetLePulse()<<" "<<L1->GetRiPulse()<<endl;











     
     
 // refer louis lyons stat book sec. 10.10:

    
    double sigmax = pow(L0->GetPosXErr(),-2) + pow(L1->GetPosXErr(),-2); //=1/sigma^2
    double sigmay = pow(L0->GetPosYErr(),-2) + pow(L1->GetPosYErr(),-2);
    double sigmaz = pow(L0->GetPosZErr(),-2) + pow(L1->GetPosZErr(),-2);
	
    fXRecoPosX = ( (L0->GetRecoPosX()/( L0->GetPosXErr()* L0->GetPosXErr())) + ( L1->GetRecoPosX()/( L1->GetPosXErr()* L1->GetPosXErr()))) /sigmax;	
         fXTruePosX = (pulse[0]* L0->GetTruePosX() + pulse[1]* L1->GetTruePosX() )/(pulse[0]+pulse[1]);
    //       fXTruePosX = ( L0->GetTruePosX() + L1->GetTruePosX())/2;

    fYRecoPosY = (( L0->GetRecoPosY()/( L0->GetPosYErr()* L0->GetPosYErr())) + ( L1->GetRecoPosY()/( L1->GetPosYErr()* L1->GetPosYErr()))) /sigmay;	
        fYTruePosY = (pulse[0]* L0->GetTruePosY() +pulse[1]* L1->GetTruePosY())/(pulse[0]+pulse[1]);

 

    fZRecoPosZ = (( L0->GetRecoPosZ()/( L0->GetPosZErr()* L0->GetPosZErr())) + ( L1->GetRecoPosZ()/( L1->GetPosZErr()* L1->GetPosZErr()))) /sigmaz;

     fZTruePosZ = (pulse[0]* L0->GetTruePosZ() +pulse[1]* L1->GetTruePosZ())/(pulse[0]+pulse[1]);
     //    fZTruePosZ = ( L0->GetTruePosZ() + L1->GetTruePosZ())/2;
	
    fXPosErr = (L0->GetPosXErr()* L1->GetPosXErr())/pow((L0->GetPosXErr()*L0->GetPosXErr())+(L1->GetPosXErr()*L1->GetPosXErr()),0.5); 
    fYPosErr = (L0->GetPosYErr()* L1->GetPosYErr())/pow((L0->GetPosYErr()*L0->GetPosYErr())+(L1->GetPosYErr()*L1->GetPosYErr()),0.5);
    fZPosErr = (L0->GetPosZErr()* L1->GetPosZErr())/pow((L0->GetPosZErr()*L0->GetPosZErr())+(L1->GetPosZErr()*L1->GetPosZErr()),0.5); 



    fLeTime = (L0->GetLeTime() +   L1->GetLeTime())/2 ;      // Timing of Left-SiPMs (ns)
    fRiTime = (L0->GetRiTime() +  L1->GetRiTime())/2  ;      // Timing of Right-SiPMs (ns)

    fLePulse = (L0->GetLePulse() +    L1->GetLePulse())/2;      // Timing of Left-SiPMs (ns)
    fRiPulse = (L0->GetRiPulse() +   L1->GetLePulse())/2;


    //true pos energy weight only along 





    


    
  }//    	if (L0 && L1 &&  L0->GetStrip() &&  L1->GetStrip() ) {//##



  else  if (L0 ) {//##//&& L0->GetStrip()
    cout<<"Only L0 "<<endl;
    plane = L0->GetPlane();
    layer = L0->GetLayer();
    strip = L0->GetStrip();
    clustersize = 1;
	   
    cout<<" plane "<<plane<<" layer "<<layer<<" strip "<<strip<<endl;

  
      
    fXRecoPosX = L0->GetRecoPosX();
    fXTruePosX = L0->GetTruePosX();
	    
    fYRecoPosY = L0->GetRecoPosY();
    fYTruePosY = L0->GetTruePosY();
	
    fZTruePosZ = L0->GetTruePosZ();
    fZRecoPosZ = L0->GetRecoPosZ();


    fLeTime = L0->GetLeTime();      // Timing of Left-SiPMs (ns)
    fRiTime = L0->GetRiTime();      // Timing of Right-SiPMs (ns)

    fLePulse = L0->GetLePulse();      // Timing of Left-SiPMs (ns)
    fRiPulse = L0->GetRiPulse();      // Timing of Right-SiPMs (ns)



	
    fXPosErr = L0->GetPosXErr(); //Error in X-position (m)
    fYPosErr = L0->GetPosYErr(); //Error in Y-position (m)
    fZPosErr = L0->GetPosZErr(); //Error in Y-position (m)



      
	  

  }//    	if (L0 && L1 &&  L0->GetStrip() &&  L1->GetStrip() ) {//##


  cout<<" plane "<<plane<<" layer "<<layer<<" strip "<<strip<<" clustersize" <<clustersize<<endl;

  cout<<"TruePos "<<fXTruePosX<<" "<<fYTruePosY<<" "<<fZTruePosZ<<endl;
  cout<<"RecoPos "<<fXRecoPosX<<" "<<fYRecoPosY<<" "<<fZRecoPosZ<<endl;

  isUsed = false;


		
}

//






void CmvCluster::CombineClusts(CmvCluster* C2)

{

  paradef = micalDetectorParameterDef::AnPointer; //AAR:
  pAnalysis = MultiSimAnalysis::AnPointer;
  PhotonSpeed = pAnalysis->GetPhotonSpeedVal();
  double CMVadctons = pAnalysis->GetCMVadctons();
  double pos[3]={0,0,0};
  double pulse[2]={0};



  cout<<" plane "<<plane<<" layer "<<layer<<" strip "<<strip<<" clustersize" <<clustersize<<endl;

  cout<<"TruePos "<<fXTruePosX<<" "<<fYTruePosY<<" "<<fZTruePosZ<<endl;
  cout<<"RecoPos "<<fXRecoPosX<<" "<<fYRecoPosY<<" "<<fZRecoPosZ<<endl;



  int plane2 = C2->GetPlane();
  int layer2 = C2->GetLayer();
  int strip2 = C2->GetStrip();
  int clustersize2 =0 ;
  cout<<" plane2 "<<plane2<<" layer2 "<<layer2<<" strip2 "<<strip2<<" clustersize2" <<clustersize2<<endl;



  double  ClustRecoPosX = C2->GetRecoPosX();
  double ClustTruePosX = C2->GetTruePosX();
	    
  double ClustRecoPosY = C2->GetRecoPosY();
  double ClustTruePosY = C2->GetTruePosY();
	
  double ClustTruePosZ = C2->GetTruePosZ();
  double ClustRecoPosZ = C2->GetRecoPosZ();

  double ClustRecoPosXErr = C2->GetPosXErr();
  double ClustRecoPosYErr = C2->GetPosYErr();
  double ClustRecoPosZErr = C2->GetPosZErr();


  double   ClustfLeTime = C2->GetLeTime();     
  double  ClustfRiTime = C2->GetRiTime();     

  double  ClustfLePulse = C2->GetLePulse();   
  double  ClustfRiPulse = C2->GetRiPulse();   

  double ClustPulse_avg = (ClustfLePulse+ClustfRiPulse)/2;
  
  // if averaging use:
  //  fXTruePosX = fXTruePosX*2;
  // fXTruePosX += ClustTruePosX;
  // fXTruePosX = fXTruePosX/3;


  //   fYTruePosY = fYTruePosY*2;
  // fYTruePosY += ClustTruePosY;
  // fYTruePosY = fYTruePosY/3;
    
  // fZTruePosZ = fZTruePosZ*2;
  // fZTruePosZ += ClustTruePosZ;
  // fZTruePosZ = fZTruePosZ/3;

  

  // if using energy weighted position use:
  double pulse_avg = (fLePulse+fRiPulse)/2;
  fXTruePosX = fXTruePosX*pulse_avg;
  fXTruePosX+= ClustPulse_avg* ClustTruePosX;
  fXTruePosX = fXTruePosX/(pulse_avg + ClustPulse_avg);

  fYTruePosY = fYTruePosY*pulse_avg;
  fYTruePosY+= ClustPulse_avg* ClustTruePosY;
  fYTruePosY = fYTruePosY/(pulse_avg + ClustPulse_avg);

  fZTruePosZ = fZTruePosZ*pulse_avg;
  fZTruePosZ+= ClustPulse_avg* ClustTruePosZ;
  fZTruePosZ = fZTruePosZ/(pulse_avg + ClustPulse_avg);


  //

  
  
  fXRecoPosX = fXRecoPosX*pow(fXPosErr,-2);
  fXRecoPosX += ClustRecoPosX/(ClustRecoPosXErr*ClustRecoPosXErr);
  fXRecoPosX = fXRecoPosX/(pow(fXPosErr,-2) +pow(ClustRecoPosXErr,-2) );




  fYRecoPosY = fYRecoPosY*pow(fYPosErr,-2);
  fYRecoPosY += ClustRecoPosY/(ClustRecoPosYErr*ClustRecoPosYErr);
  fYRecoPosY = fYRecoPosY/(pow(fYPosErr,-2) +pow(ClustRecoPosYErr,-2) );




  fZRecoPosZ = fZRecoPosZ*pow(fZPosErr,-2);
  fZRecoPosZ += ClustRecoPosZ/(ClustRecoPosZErr*ClustRecoPosZErr);
  fZRecoPosZ = fZRecoPosZ/(pow(fZPosErr,-2) +pow(ClustRecoPosZErr,-2) );

  fXPosErr = pow(fXPosErr,-2);
  fXPosErr += pow(ClustRecoPosXErr,-2) ;
fXPosErr = pow(fXPosErr,-0.5);


  fYPosErr = pow(fYPosErr,-2);
  fYPosErr += pow(ClustRecoPosYErr,-2);
fYPosErr = pow(fYPosErr,-0.5);

  fZPosErr = pow(fZPosErr,-2);
  fZPosErr += pow(ClustRecoPosZErr,-2 );
fZPosErr = pow(fZPosErr,-0.5);


 fLeTime = fLeTime*2;
  fLeTime +=  ClustfLeTime;
  fLeTime = fLeTime/3;
  
  fRiTime = fRiTime*2;
  fRiTime +=  ClustfRiTime;     
  fRiTime = fRiTime/3;

    fLePulse = fLePulse*2;
  fLePulse +=  ClustfLePulse;
  fLePulse = fLePulse/3;

    fRiPulse = fRiPulse*2;
  fRiPulse +=  ClustfRiPulse;   
  fRiPulse = fRiPulse/3;
  clustersize=3;
  isUsed = false;
		
}


















//previous
CmvCluster::CmvCluster(CmvCluster* C0, CmvCluster* C1):
  clustersize(-1),
  plane(-1), layer(-1), strip(-1),
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
  double pos[3]={0,0,0};
  double pulse[2]={0};
  if (C0 && C1 &&  C0->GetStrip() &&  C1->GetStrip() ) {//##
    plane = C0->GetPlane();
    layer = C0->GetLayer();
    strip = C0->GetStrip();

    clustersize = 2;
    cout<<" plane "<<plane<<" layer "<<layer<<" strip "<<strip<<endl;


    cout<<C0->GetRecoPosX()<<" "<<C0->GetRecoPosY()<<" "<<C0->GetRecoPosZ()<<" "<<endl;
    cout<<C1->GetRecoPosX()<<" "<<C1->GetRecoPosY()<<" "<<C1->GetRecoPosZ()<<" "<<endl;


	  
    pulse[0]=(C0->GetLePulse() + C0->GetRiPulse())/2;
    pulse[1]=(C1->GetLePulse() + C1->GetRiPulse())/2;

    cout<<"pulse "<<pulse[0]<<" "<<pulse[1]<<endl;

    // if(C0->GetPlane()==C1->GetPlane()==1){
    //     cout<<"case2 "<< C0->GetPlane()<<" "<< C1->GetPlane()<<endl;
    fXRecoPosX =(( pulse[0]* C0->GetRecoPosX()) + ( pulse[1]* C1->GetRecoPosX()))/(pulse[0]+pulse[1]);
    fXTruePosX =(( pulse[0]* C0->GetTruePosX()) + ( pulse[1]* C1->GetTruePosX()))/(pulse[0]+pulse[1]);
	    
    fYRecoPosY =(( pulse[0]* C0->GetRecoPosY()) + ( pulse[1]* C1->GetRecoPosY()))/(pulse[0]+pulse[1]);
    fYTruePosY =(( pulse[0]* C0->GetTruePosY()) + ( pulse[1]* C1->GetTruePosY()))/(pulse[0]+pulse[1]);
	
    fZTruePosZ =(( pulse[0]* C0->GetTruePosZ()) + ( pulse[1]* C1->GetTruePosZ()))/(pulse[0]+pulse[1]);
    fZRecoPosZ =(( pulse[0]* C0->GetRecoPosZ()) + ( pulse[1]* C1->GetRecoPosZ()))/(pulse[0]+pulse[1]);
	    
    //   }//    if(C0->GetPlane()==C1->GetPlane()==1){



    // else  if(C0->GetPlane()==C1->GetPlane()==2 ||C0->GetPlane()==C1->GetPlane()==3 ){
    //     cout<<"case1 "<< C0->GetPlane()<<" "<< C1->GetPlane()<<endl;
    //     fXRecoPosX =( C0->GetRecoPosX() +  C1->GetRecoPosX())/clustersize;
    //     fXTruePosX =( ( pulse[0]* C0->GetTruePosX()) + ( pulse[1]* C1->GetTruePosX()))/(pulse[0]+pulse[1]);
	  
    //     fYRecoPosY =(( pulse[0]* C0->GetRecoPosY()) + ( pulse[1]* C1->GetRecoPosY()))/(pulse[0]+pulse[1]);
    //     fYTruePosY =(( pulse[0]* C0->GetTruePosY()) + ( pulse[1]* C1->GetTruePosY()))/(pulse[0]+pulse[1]);
	  
    //     fZRecoPosZ =(( pulse[0]* C0->GetRecoPosZ()) + ( pulse[1]* C1->GetRecoPosZ()))/(pulse[0]+pulse[1]);
    //     fZTruePosZ =(( pulse[0]* C0->GetTruePosZ()) + ( pulse[1]* C1->GetTruePosZ()))/(pulse[0]+pulse[1]);

    //   }//    else  if(C0->GetPlane()==C1->GetPlane()==2 ||C0->GetPlane()==C1->GetPlane()==3 ){




    //   else   if(C0->GetPlane()==4 || C0->GetPlane()==4){
    //     cout<<"case4 "<< C0->GetPlane()<<" "<< C1->GetPlane()<<endl;
    //    fXRecoPosX =(( pulse[0]* C0->GetRecoPosX()) + ( pulse[1]* C1->GetRecoPosX()))/(pulse[0]+pulse[1]);
      
    
    //     fXTruePosX =(( pulse[0]* C0->GetTruePosX()) + ( pulse[1]* C1->GetTruePosX()))/(pulse[0]+pulse[1]);
	  
    //         fYRecoPosY =(C0->GetRecoPosY()+ C1->GetRecoPosY())/clustersize;
    //     fYTruePosY =(( pulse[0]* C0->GetTruePosY()) + ( pulse[1]* C1->GetTruePosY()))/(pulse[0]+pulse[1]);
	  
    //     fZRecoPosZ =(( pulse[0]* C0->GetRecoPosZ()) + ( pulse[1]* C1->GetRecoPosZ()))/(pulse[0]+pulse[1]);
    //     fZTruePosZ =(( pulse[0]* C0->GetTruePosZ()) + ( pulse[1]* C1->GetTruePosZ()))/(pulse[0]+pulse[1]);

    //   }//   else   if(C0->GetPlane()==2 && C1->GetPlane()==2 ||C0->GetPlane()==3 && C1->GetPlane()==3 ){











  


  }//    if (C0 && C1 &&  C0->GetStrip() &&  C1->GetStrip() ) {//##





  else  if (C0 && C0->GetStrip() ) {//##
    cout<<"Only C0 "<<endl;
    plane = C0->GetPlane();
    layer = C0->GetLayer();
    strip = C0->GetStrip();
    clustersize = -1;
	   
    cout<<" plane "<<plane<<" layer "<<layer<<" strip "<<strip<<endl;

  
      
    fXRecoPosX = C0->GetRecoPosX();
    fXTruePosX = C0->GetTruePosX();
	    
    fYRecoPosY = C0->GetRecoPosY();
    fYTruePosY = C0->GetTruePosY();
	
    fZTruePosZ = C0->GetTruePosZ();
    fZRecoPosZ = C0->GetRecoPosZ();

    fLeTime = C0->GetLeTime();      // Timing of Left-SiPMs (ns)
    fRiTime = C0->GetRiTime();      // Timing of Right-SiPMs (ns)

    fLePulse = C0->GetLePulse();      // Timing of Left-SiPMs (ns)
    fRiPulse = C0->GetRiPulse();      // Timing of Right-SiPMs (ns)

  	
    fXPosErr = C0->GetPosXErr(); //Error in X-position (m)
    fYPosErr = C0->GetPosYErr(); //Error in Y-position (m)
    fZPosErr = C0->GetPosZErr(); //Error in Y-position (m)


   

  
	  

  }//    	if (C0 && L1 &&  C0->GetStrip() &&  L1->GetStrip() ) {//##













  

  cout<<" plane "<<plane<<" layer "<<layer<<" strip "<<strip<<" clustersize" <<clustersize<<endl;

  cout<<"TruePos "<<fXTruePosX<<" "<<fYTruePosY<<" "<<fZTruePosZ<<endl;
  cout<<"RecoPos "<<fXRecoPosX<<" "<<fYRecoPosY<<" "<<fZRecoPosZ<<endl;

  isUsed = false;
		
}











CmvCluster::~CmvCluster()
{
}

// bool CmvCluster::isIdentical(CmvCluster* cluster) {

//   if ((cluster->GetXStripNum() == GetXStripNum()) &&
//       (cluster->GetYStripNum() == GetYStripNum())) {
//     return true;
//   }
//   return false;
// }

void CmvCluster::Print() {
  //	cout<<"----------------------------------------------------------------------"<<endl;
  //	cout<<"Cluster combination "<<endl;
  cout<< "CmvClusters():" 
    // <<std::setw(4) <<jk <<" "
    // << "Idd="<<std::setw(3)<<GetStripId()
      << "Plane "<< std::setw(3)<<  GetPlane()
      << std::setw(3)<<  GetLayer()
      << std::setw(4)<<  GetStrip()
       << std::setw(4)<<  GetClustersize()
      << " X_Reco="   <<std::setw(8)<<  GetRecoPosX() 
      << " Y_Reco=" <<std::setw(8)<<  GetRecoPosY()
      << " Z_reco=" <<std::setw(8)<<  GetRecoPosZ()
      << " X_Pos=" <<std::setw(8)<<  GetTruePosX()
      << " Y_Pos=" <<std::setw(8)<<  GetTruePosY()
      << " Z_Pos=" <<std::setw(8)<<  GetTruePosZ()
    // << " chg="   <<std::setw(8)<<  GetPulse()
      << " time="  <<std::setw(8)<<  GetLeTime()
      << endl;
  //  cout<<"......................................................................"<<endl;
}

// void CmvCluster::Print() {
//   cout<<"----------------------------------------------------------------------"<<endl;



//   cout<<"......................................................................"<<endl;


// }

*/

CmvCluster::CmvCluster(CmvHit* hit) :
  fPlane(-1),fLayer(-1),
 fStrip(-1), fBegStrip(-1),fEndStrip(-1),

  fClusterSize(0), isUsed(0),
  
  fTruePosX(-1000), fBegTruePosX(-1000), fEndTruePosX(-1000),
  fTruePosY(-1000), fBegTruePosY(-1000), fEndTruePosY(-1000),
  fTruePosZ(-1000), fBegTruePosZ(-1000), fEndTruePosZ(-1000),

  fRecoPosX(-1000), fBegRecoPosX(-1000), fEndRecoPosX(-1000),
  fRecoPosY(-1000), fBegRecoPosY(-1000), fEndRecoPosY(-1000),
  fRecoPosZ(-1000), fBegRecoPosZ(-1000), fEndRecoPosZ(-1000),

  fWRecoPosX(-1000), fWRecoPosY(-1000),fWRecoPosZ(-1000),
  
  // fLeTime (-1000), fBegLeTime (-1000), fEndLeTime (-1000),
  // fRiTime(-1000), fBegRiTime(-1000), fEndRiTime(-1000),

  fPulse (-1000), fBegPulse (-1000), fEndPulse (-1000),


  
  fXPosErr(-1000),  fYPosErr(-1000), fZPosErr(-1000),
  fMomentum(-1000), fTheta(-1000),  fPhi(-1000), PhotonSpeed(162.0), layergap(0.02)


{
  paradef = micalDetectorParameterDef::AnPointer; //AAR: 

  this->AddHits(hit);
}


CmvCluster::~CmvCluster()
{
  HitsInCmvCluster.clear();
}


void CmvCluster::AddHits(CmvHit* hit) {
  if(HitsInCmvCluster.size()==0) {
    HitsInCmvCluster.push_back(hit);
    fPlane=hit->GetPlane();
    fLayer = hit->GetLayer();

    fStrip = hit->GetStrip();
    fBegStrip = hit->GetStrip();
    fEndStrip = hit->GetStrip();    
    
    isUsed    = hit->isUsed;
    
    fTruePosX = hit->GetTruePosX();   fBegTruePosX = hit->GetTruePosX();  fEndTruePosX = hit->GetTruePosX();    
    fTruePosY = hit->GetTruePosY();   fBegTruePosY = hit->GetTruePosY();  fEndTruePosY = hit->GetTruePosY(); 
    fTruePosZ = hit->GetTruePosZ();   fBegTruePosZ = hit->GetTruePosZ();  fEndTruePosZ = hit->GetTruePosZ();
 


    fRecoPosX = hit->GetRecoPosX();   fBegRecoPosX = hit->GetRecoPosX();  fEndRecoPosX = hit->GetRecoPosX();    
    fRecoPosY = hit->GetRecoPosY();   fBegRecoPosY = hit->GetRecoPosY();  fEndRecoPosY = hit->GetRecoPosY(); 
    fRecoPosZ = hit->GetRecoPosZ();   fBegRecoPosZ = hit->GetRecoPosZ();  fEndRecoPosZ = hit->GetRecoPosZ();


    fWRecoPosX = hit->GetRecoPosX();
    fWRecoPosY = hit->GetRecoPosY();
    fWRecoPosZ = hit->GetRecoPosZ();
    
  // fLeTime = hit->fLeTime; fBegLeTime = hit->fLeTime; fEndLeTime = hit->fLeTime; 
  // fRiTime = hit->fRiTime; fBegRiTime = hit->fRiTime; fEndRiTime = hit->fRiTime;

  
    fBegPulse = hit->GetPulse();
    fEndPulse = hit->GetPulse();

  
  fXPosErr = hit->GetPosXErr();
  fYPosErr = hit->GetPosYErr(); 
  fZPosErr = hit->GetPosZErr(); 
	
  fMomentum = hit->GetGenMom(); 
  fTheta = hit->GetGenThe();  
  fPhi = hit->GetGenPhi();   
  PhotonSpeed = hit->PhotonSpeed;
  layergap = hit->layergap;

  fPulse = 0;

  } else {
    if(this->ContainsHit(hit)==true) {return;}
    HitsInCmvCluster.push_back(hit);   
    
    ////  if (hit->GetXYPosErr()<100) {
      if(hit->GetStrip()<fBegStrip) fBegStrip=hit->GetStrip();
      if(hit->GetStrip()>fEndStrip) fEndStrip=hit->GetStrip();

      if(hit->GetTruePosX()<fBegTruePosX) fBegTruePosX=hit->GetTruePosX();
      if(hit->GetTruePosY()<fBegTruePosY) fBegTruePosY=hit->GetTruePosY();
      if(hit->GetTruePosZ()<fBegTruePosZ) fBegTruePosZ=hit->GetTruePosZ();
 
      if(hit->GetRecoPosX()<fBegRecoPosX) fBegRecoPosX=hit->GetRecoPosX();
      if(hit->GetRecoPosY()<fBegRecoPosY) fBegRecoPosY=hit->GetRecoPosY();
      if(hit->GetRecoPosZ()<fBegRecoPosZ) fBegRecoPosZ=hit->GetRecoPosZ();

     if(hit->GetPulse()<fBegPulse) fBegPulse=hit->GetPulse();
      if(hit->GetPulse()>fEndPulse) fEndPulse=hit->GetPulse();

 

      //just Avg:
    fTruePosX = ( fClusterSize*fTruePosX + hit->GetTruePosX() )/(fClusterSize + 1); 
    fTruePosY = ( fClusterSize*fTruePosY + hit->GetTruePosY() )/(fClusterSize + 1);
    fTruePosZ = ( fClusterSize*fTruePosZ + hit->GetTruePosZ() )/(fClusterSize + 1); 

    fRecoPosX = ( fClusterSize*fRecoPosX + hit->GetRecoPosX() )/(fClusterSize + 1); 
    fRecoPosY = ( fClusterSize*fRecoPosY + hit->GetRecoPosY() )/(fClusterSize + 1);
    fRecoPosZ = ( fClusterSize*fRecoPosZ + hit->GetRecoPosZ() )/(fClusterSize + 1); 

    cout<<".. fPulse fWRecoPosX hit->GetRecoPosX() hit->GetPulse() "<<  fPulse<<" "<< fWRecoPosX<<" "<< hit->GetRecoPosX()<<" "<< hit->GetPulse()<<endl;
    //energy weighted average
    fWRecoPosX = ( fPulse*fWRecoPosX + hit->GetRecoPosX()*hit->GetPulse() )/(fPulse+hit->GetPulse() );
    fWRecoPosY = ( fPulse*fWRecoPosY + hit->GetRecoPosY()*hit->GetPulse() )/(fPulse+hit->GetPulse() );
    fWRecoPosZ = ( fPulse*fWRecoPosZ + hit->GetRecoPosZ()*hit->GetPulse() )/(fPulse+hit->GetPulse() );
      
  }
  

    
  fPulse += hit->GetPulse();


    cout<<" fPulse fWRecoPosX hit->GetRecoPosX() hit->GetPulse() "<<  fPulse<<" "<< fWRecoPosX<<" "<< hit->GetRecoPosX()<<" "<< hit->GetPulse()<<endl;


  fClusterSize++;
  return;
}

bool CmvCluster::ContainsHit(CmvHit* hit) {
  for(unsigned int ij=0; ij<HitsInCmvCluster.size(); ++ij) {
    if(hit==HitsInCmvCluster[ij]) {return true;}
  }
  return false;
}



CmvHit* CmvCluster::GetHit(unsigned int ij) const {
  if(ij<HitsInCmvCluster.size()) {return HitsInCmvCluster[ij];}
  else {return 0;}

}

void CmvCluster::Print() {

  cout<<" CmvCluster():"
      <<" Plane = "<<fPlane
      <<" Layer = "<<fLayer
      <<" Beg Hit = "<< std::setw(8)<<fBegStrip 
      <<" End Hit = "<< std::setw(8)<<fEndStrip 
      <<" ClusterNo = "<< std::setw(8)<<fStrip 
      <<" Beg Pulse = "<< std::setw(8)<<fBegPulse 
      <<" End Pulse = "<< std::setw(8)<<fEndPulse 
      <<" Pulse = "<< std::setw(8)<<fPulse 
        
      <<" BegTruePos = ("<<fBegTruePosX<<","<<fBegTruePosY<<","<<fBegTruePosZ<<")"
      <<" EndTruePos = ("<<fEndTruePosX<<","<<fEndTruePosY<<","<<fEndTruePosZ<<")"
      <<" TruePos = ("<<fTruePosX<<","<<fTruePosY<<","<<fTruePosZ<<")"

      <<" BegRecoPos = ("<<fBegRecoPosX<<","<<fBegRecoPosY<<","<<fBegRecoPosZ<<")"
      <<" EndRecoPos = ("<<fEndRecoPosX<<","<<fEndRecoPosY<<","<<fEndRecoPosZ<<")"
      <<" RecoPos = ("<<fRecoPosX<<","<<fRecoPosY<<","<<fRecoPosZ<<")"
    
      <<" WRecoPos = ("<<fWRecoPosX<<","<<fWRecoPosY<<","<<fWRecoPosZ<<")"

      <<" ClusterSize = "<< std::setw(8)<<fClusterSize
      <<endl;
  // for(unsigned int ij=0; ij<HitsInCmvCluster.size(); ij++) {
  //   HitsInCmvCluster[ij]->Print();
  // }
 
}
