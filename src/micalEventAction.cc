//
// $Id: micalEventAction.cc,v 1.24 2005/05/30 14:24:31 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
//GMAA Store begin and end hit position and muo energy at those points for the resolution of reconstruction algorithms, which is much better than true muon energy resolution

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "micalEventAction.hh"
#include "micalEventActionMessenger.hh"

#include "micalCal0Hit.hh"
#include "micalCal1Hit.hh" //cmv
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"

#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "micalElectroMagneticField.hh"
#include "Randomize.hh"
#include <iomanip>
#include <utility>
#include "G4SDManager.hh"
#include "G4Trajectory.hh"
//#include <iomanip.h>

#include "vect_manager.h"
//#include "InoPatternRecognition.h" //GMA14
#include "InoTrackFinder.h"
#include "InoTrackFitAlg.h"
#include "InoVertex.h"
#include "TString.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TMath.h"
//#include "TVectorD.h"
//#include "TMatrixTBase.h"
#include "TMatrixDEigen.h"
#include "StraightLineFit.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool micalEventAction::ClosDistbtwLineEdge(double* Line,double* Plane, std::vector<double>& edge, double* GOnLine, double* GOnedge){
  cout<<"void micalEventAction::ClosDistbtwLineEdge(double* Line, std::vector<double>& edge, double* GOnLine, double* GOnedge){"<<endl;
  
  bool ok;	
  double b = Line[3]*Plane[3] + Line[4]*Plane[4] + Line[5]*Plane[5];
  cout<<Line[3]<<" "<<Plane[3]<<" "<<Line[4]<<" "<<Plane[4]<<" "<<Line[5]<<" "<<Plane[5]<<" "<<b<<endl;
  // ok= (fabs(b) > 1e-10) ? 1 : 0;//
  
  ok= ((1.0*b)> 1e-10) ? 1 : 0;//#we only do intersection with those planes which muon passed first..
  if(ok==1) {
  
 double rveccloseedge[3] ={ -Line[0]+edge[0],-Line[1]+edge[1],-Line[2]+edge[2] };
	    cout<<"rveccloseedge: "<<rveccloseedge[0]<<" "<<rveccloseedge[1]<<" "<<rveccloseedge[2]<<endl;
 
	    double A1dotA2 = Line[3]*edge[3]+Line[4]*edge[4]+Line[5]*edge[5];
	    double A1dotrvec =Line[3]*rveccloseedge[0]+Line[4]*rveccloseedge[1]+Line[5]*rveccloseedge[2];
	    double A2dotrvec = edge[3]*rveccloseedge[0]+edge[4]*rveccloseedge[1]+edge[5]*rveccloseedge[2];  
 
	    double lambda = (A1dotrvec - A2dotrvec*(A1dotA2))/(1-pow(A1dotA2,2));
	    double mu = lambda*A1dotA2 - A2dotrvec;
	       
	    //	    cout<<mu<<" "<<lambda<<endl;
	    GOnLine[0] = {Line[0]+lambda*Line[3]};
	    GOnLine[1] = {Line[1]+lambda*Line[4]};
	    GOnLine[2] = {Line[2]+lambda*Line[5]};
	    
	    // cout<<edge[0]<<" "<<mu<<" "<<edge[3]<<" "<<miniedge<<" "<<edge[0]+mu*edge[3]<<endl;
	    GOnedge[0] = {edge[0]+mu*edge[3]};
	    GOnedge[1] = {edge[1]+mu*edge[4]};
	    GOnedge[2] = {edge[2]+mu*edge[5]};  
  }
  else{
    GOnLine[0]=1000000; GOnLine[1]=1000000; GOnLine[2]=1000000;
    GOnedge[0]=1000000; GOnedge[1]=1000000; GOnedge[2]=1000000;

  }
  return ok;  

}
 //cmv straightline extrapolation

bool micalEventAction::LinePlaneInt(double* Line, double* Plane, double* Point){
  //  cout<<"bool micalEventAction::LinePlaneInt(double* Line, double* Plane, double* Point){"<<endl;
  //	G4double Dist;
  //	G4double a, b;
  bool ok;
	
  double b = Line[3]*Plane[3] + Line[4]*Plane[4] + Line[5]*Plane[5];
  //  cout<<Line[3]<<" "<<Plane[3]<<" "<<Line[4]<<" "<<Plane[4]<<" "<<Line[5]<<" "<<Plane[5]<<" "<<b<<endl;
  // ok= (fabs(b) > 1e-10) ? 1 : 0;//
  
  ok= ((1.0*b)> 1e-10) ? 1 : 0;//#we only do intersection with those planes which muon passed first..
  if(ok==1) {
    double a=(Plane[0]-Line[0])*Plane[3] +
      (Plane[1]-Line[1])*Plane[4] +
      (Plane[2]-Line[2])*Plane[5];
    //    cout<<"a "<<a<<endl;
    G4double Dist = a/b;
    //    cout<<"dist: "<<Dist<<endl;
    
    Point[0] = Line[0] + Line[3]*Dist;
    Point[1] = Line[1] + Line[4]*Dist;
    Point[2] = Line[2] + Line[5]*Dist;
    //    cout<<"Point "<<Point[0]<<" "<<Point[1]<<" "<<Point[2]<<endl;
  } else {
    //    cout<<"Setting Point =-100000 "<< endl;
    Point[0]=-1000000; Point[1]=-1000000; Point[2]=-1000000;
  }
  return ok;
}
//




void micalEventAction::CMVD_Extrapolation(){   

  
  cout<<"......................CMVD straight line extrapolation........................."<<endl;

  cout<<"Store in CmvClusterBank"<<endl;
  cout<<"CmvCluster_pointer->CmvCluster_list.size() "<<CmvCluster_pointer->CmvCluster_list.size()<<endl;
  for(int pl=0;pl<7;pl++){
    for(int la=0;la<4;la++){

      CmvClusterBank[pl][la].clear();


    }



  }
  for (unsigned int ix=0; ix<CmvCluster_pointer->CmvCluster_list.size(); ix++) {
    CmvCluster* tmpcluster = CmvCluster_pointer->CmvCluster_list[ix];
    tmpcluster->Print();
    CmvClusterBank[tmpcluster->GetPlane()-1][tmpcluster->GetLayer()].push_back(tmpcluster);
  }

  
  //  CmvLayExtra_pointer = new CmvLayExtra_Manager();
  //  CmvLayExtra_pointer->CmvLayExtra_list.clear();


  int counter=0;

  
  int ijmax=0;
  cout<<"check ab "<<pAnalysis->ntrkt<<endl;
  for (unsigned jk=0; jk<pAnalysis->ntrkt ; jk++) {

    pAnalysis->CMVDnLayer =pAnalysis->cmvdnlayermx;
    cout<<"ntrkt CMVDnLayer "<<pAnalysis->ntrkt<<" "<<pAnalysis->CMVDnLayer<< endl;
  
    //		double    momvx = pAnalysis->momvx[jk];
    //This theta phi represents a track going downward. The dirvector using this theta phi matched with generated direction vector (in PGA after...)
    double theta =pAnalysis->thevx[jk];
    double phi = pAnalysis->phivx[jk];
    double posx = pAnalysis->posxvx[jk]*m; //these were stored in metre and convert in mm unit
    double posy =pAnalysis->posyvx[jk]*m;
    double posz= pAnalysis->poszvx[jk]*mm;//already stored in mm only


    double atimslope = pAnalysis->atimslope[jk];
    double atiminter = pAnalysis->atiminter[jk];


    cout<<"atimslope "<<atimslope<<"atiminter "<<atiminter<<endl;
    cout<<"...jk..... "<<jk<<" "<<posx<<" "<<posy<<" "<<posz<<" "<<theta<<" "<<phi<<endl;
    //err calculation:
    
    double therr = pAnalysis->therr[jk];
    double pherr = pAnalysis->pherr[jk];
    // these are errors i.e. just sigma and are stored in meters..
    double xxerr = pAnalysis->xxerr[jk]*m;
    double yyerr = pAnalysis->yyerr[jk]*m;
    double txerr = pAnalysis->txerr[jk];// x = a + b z , err in b (no units)
    double tyerr = pAnalysis->tyerr[jk];// y = a + b z , err in b (no units)
    cout<<"error in positions and angles  "<<xxerr<<" "<<yyerr<<" "<<txerr<<" "<<tyerr<<" "<<therr<<" "<<pherr<<endl;
    
    double xxtxerr = pAnalysis->xxtxerr[jk]*m;
    double xxtyerr = pAnalysis->xxtyerr[jk]*m;
    double yytyerr = pAnalysis->yytyerr[jk]*m;
    double yytxerr = pAnalysis->yytxerr[jk]*m;
    double txtyerr = pAnalysis->txtyerr[jk];
    cout<<"xxtxerr "<<xxtxerr<<"yytyerr "<<yytyerr<<endl;
    //err cal.

    //    pAnalysis->chisq[jk]=-1.0; pAnalysis->chisq2[jk]=-1.0;
    //    pAnalysis->posxvx[jk]=0.0;  pAnalysis->posyvx[jk]=0.0;  pAnalysis->poszvx[jk]=0.0;
   


    //convert theta and phi to dxdz and dydz
    
    double dxdz = tan(theta)*cos(phi);   
    double dydz = tan(theta)*sin(phi);
    
    cout<< "Downward theta "<< theta<<" phi  "<< phi<< " dxdz "<< dxdz<<" dydz "<<dydz<<endl;


    double PI = acos(-1.0);
	
    //This gives a downward going track
    // G4ThreeVector dirVector1(0,0,1);       
    // dirVector1.setTheta(theta);
    // dirVector1.setPhi(phi);
    //We need to do extrapolation for those cmvd plane which a muon have passed before entering RPC stack. Thus we need to reverse sign of dirvector(parity) thus we add Pi to phi and do theta-pi..


    double dxdz_upward = tan(PI-theta)*cos(phi+PI);   //modified on 10.01.2022
    double dydz_upward = tan(PI-theta)*sin(phi+PI);
    //slope doesnt change on inverting all xyz to negative
    cout<< "Upward  theta "<< PI-theta<<" phi  "<< phi+PI<< " dxdz "<< dxdz_upward<<" dydz "<<dydz_upward<<endl;



    
    G4ThreeVector dirVector(0,0,1);  //     
    dirVector.setTheta(PI-theta);// this works atleast for SL
    dirVector.setPhi(phi+PI);

    //   dirVector.setTheta(theta);
    //  dirVector.setPhi(phi);
    

    cout<<"dirvector "<<dirVector<<" "<<endl;
    //thus this dirvector represents a track going backward towards the vertex/gen point.
    G4double G_Point[3]={0};
    G4double Ptxerr, Ptyerr, Ptzerr =0;
    G4double Line[6]={posx,posy,posz,dirVector.x(),dirVector.y(),dirVector.z()};
    
    cout<<"Line: "<<Line[0]<<" "<<Line[1]<<" "<<Line[2]<<" "<<Line[3]<<" "<<Line[4]<<" "<<Line[5]<<endl;



    // for (unsigned int ij=0; ij<CmvHit_pointer->CmvHit_list.size(); ij++) {

      
    //   CmvHit_pointer->CmvHit_list[ij]->Print();
    // }
    // //...
 

    unsigned int layid = 0;
    double ellip_diff[4][4];
    double erralgstrplen[4][4];
    double  xhat=0,yhat=0,zhat=0; //area unit vector of planes

    double layhalflength;
    for(int ijk=0;ijk<7;ijk++){//loc_no loop 0:Top, 1:Left 2:Right 3:Back 4:Front 5: miniLeft 6:miniRight //4
      if(ijk==4){continue;} //Reserved for Front Wall
      cout<<"locno: "<<ijk<<endl;
      if(ijk==0){ijmax=4;xhat=0;yhat=0;zhat=1;layhalflength = paradef->GetTopPlaneHalfLength();}//top
      else if(ijk==1 ){ijmax=3;xhat=-1;yhat=0;zhat=0; layhalflength = paradef->GetSidePlaneHalfLength(); }//left
      else if(ijk==2){ijmax=3;xhat=1;yhat=0;zhat=0;layhalflength = paradef->GetSidePlaneHalfLength(); }//right
      else if(ijk==3){ijmax=3;xhat=0;yhat=1;zhat=0;layhalflength = paradef->GetSidePlaneHalfLength();}//back
      else if(ijk==4){ijmax=3;xhat=0;yhat=-1;zhat=0;layhalflength = paradef->GetSidePlaneHalfLength();}//Front
      else if(ijk==5){ijmax=3;xhat=-1;yhat=0;zhat=0;layhalflength = paradef->GetSideSmallPlaneHalfLength();}//miniLeft
      else if(ijk==6){ijmax=3;xhat=1;yhat=0;zhat=0;layhalflength = paradef->GetSideSmallPlaneHalfLength();}//miniRight

      for(int ij=0;ij<ijmax;ij++){//layer loop
	double diffmx=1000000;
	cout<<" Layer No. "<<ij<<endl;
	//	cout<<ijk<<" "<<ij<<endl;
	G4double Plane[6]={PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2],xhat,yhat,zhat};

        cout<<"Plane: "<<Plane[0]<<" "<<Plane[1]<<" "<<Plane[2]<<" "<<Plane[3]<<" "<<Plane[4]<<" "<<Plane[5]<<endl;
	

        bool pl2 = LinePlaneInt (Line, Plane, G_Point);

	double layhalfbreadth = paradef->partopscint[1];     //2300 for left right
	if(ijk==3){
	  layhalfbreadth = paradef->partopscint[1]+50;//4.7 m    

	}
	if(ijk==0){
	  layhalfbreadth = paradef->partopscint[1]-50;//4.5 m    

	}

	if(ijk==5 || ijk==6){
	  layhalfbreadth = 1000;

	}



		
	//	cout<<"Loc_no: "<<ijk<<" layhalflength "<<layhalflength<<"layhalfbreadth  "<<layhalfbreadth<<endl;

	// cout<<"Point: "<<G_Point[0]<<" "<<G_Point[1]<<" "<<G_Point[2]<<endl;
	vector <double> edge[4];
    
	if(ijk==0){

	  edge[0] =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]+layhalfbreadth, PhyVolGlPos[ijk][ij][2],-1,0,0}; //backside edge
	  edge[1] =	      {PhyVolGlPos[ijk][ij][0]-layhalflength, PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2],0,-1,0}; //leftside edge
	  edge[2] =	      {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]-layhalfbreadth, PhyVolGlPos[ijk][ij][2],1,0,0}; //frontside edge
	  edge[3] =	    {PhyVolGlPos[ijk][ij][0]+layhalflength, PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2],0,1,0}; //rightside edge


	}
	      
	if(ijk==1){

	  edge[0]    =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]+layhalflength,0,1,0}; //topside edge
	  edge[1] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]+layhalfbreadth, PhyVolGlPos[ijk][ij][2],0,0,-1}; //backside edge
	  edge[2] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]-layhalflength,0,-1,0}; //bottomside edge
	  edge[3] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]-layhalfbreadth, PhyVolGlPos[ijk][ij][2],0,0,1} ; //frontside edge
		
		
	}
	      
	if(ijk==2){
		
	  edge[0]  =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]+layhalflength,0,-1,0}; //topside edge
	  edge[1] =     {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]+layhalfbreadth, PhyVolGlPos[ijk][ij][2],0,0,-1}; //backside edge
	  edge[2] =    {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]-layhalflength,0,1,0}; //bottomside edge
	  edge[3] =   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]-layhalfbreadth, PhyVolGlPos[ijk][ij][2],0,0,1}; //frontside edge
		
		  
	}
	      
	if(ijk==3){
		
	  edge[0]   =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]+layhalflength,1,0,0}; //topside edge
	  edge[1] =	{PhyVolGlPos[ijk][ij][0]+layhalfbreadth, PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2],0,0,-1}; //rightside edge
	  edge[2] =	  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]-layhalflength,-1,0,0}; //bottomside edge
	  edge[3] =    {PhyVolGlPos[ijk][ij][0]-layhalfbreadth, PhyVolGlPos[ijk][ij][1],PhyVolGlPos[ijk][ij][2],0,0,1}; //leftside edge
		
	       
	}
	//miniLeft
	if(ijk==5){

	  edge[0]    =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]+layhalfbreadth,0,1,0}; //topside edge
	  edge[1] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]+layhalflength, PhyVolGlPos[ijk][ij][2],0,0,-1}; //backside edge
	  edge[2] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]-layhalfbreadth,0,-1,0}; //bottomside edge
	  edge[3] =	   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]-layhalflength, PhyVolGlPos[ijk][ij][2],0,0,1} ; //frontside edge
		
		
	}
	      
	if(ijk==6){

	  edge[0]  =  {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]+layhalfbreadth,0,-1,0}; //topside edge
	  edge[1] =     {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]+layhalflength, PhyVolGlPos[ijk][ij][2],0,0,-1}; //backside edge
	  edge[2] =    {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1], PhyVolGlPos[ijk][ij][2]-layhalfbreadth,0,1,0}; //bottomside edge
	  edge[3] =   {PhyVolGlPos[ijk][ij][0], PhyVolGlPos[ijk][ij][1]-layhalflength, PhyVolGlPos[ijk][ij][2],0,0,1}; //frontside edge
		
		  
	}

	// CmvLayExtra* layexp = new CmvLayExtra(); //GMA Memory leakages ??

	     layid = ijk+1;
	     layid<<=2;
	     layid+=ij;

	     cout<<"layid: "<<layid<<endl;
	//     layexp->SetId(layid);
	 
        if(pl2){
        
	  cout<<"---Intersection with Plane found---"<<endl;
	 
	  bool isInside= false;
	  int delta = 0;// 100cm
	  //  cout<<"checkkk: "<<G_Point[2]<<" "<<PhyVolGlPos[ijk][ij][2]<<endl;
	    
	  switch (ijk){
	    //  cout<<"checkkk: "<<G_Point[2]<<PhyVolGlPos[ijk][ij][2]<<endl;
          case 0: isInside = ( (G_Point[0])<  (PhyVolGlPos[ijk][ij][0]+layhalflength+delta) &&  (G_Point[0])>  (PhyVolGlPos[ijk][ij][0]-layhalflength-delta) &&  (G_Point[1])<  (PhyVolGlPos[ijk][ij][1]+layhalfbreadth+delta) &&  (G_Point[1])>  (PhyVolGlPos[ijk][ij][1]-layhalfbreadth-delta)    && abs(G_Point[2]-PhyVolGlPos[ijk][ij][2])<1.e-10 );
	    // cout<<"case 0"<<endl;
	    //	    cout<<fixed<<G_Point[2]<<" "<<fixed<<PhyVolGlPos[ijk][ij][2]<<endl;
	    // cout<<"0"<<PhyVolGlPos[ijk][ij][0]<<" "<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<endl;
	    
	    // cout<<"   (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    "<<    (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   "<<     (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][0])+layhalflength   "<<    (PhyVolGlPos[ijk][ij][0])+layhalflength   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][0])-layhalflength     "<<     (PhyVolGlPos[ijk][ij][0])-layhalflength   <<endl;

   
	    
            break;
	    
          case 1: isInside = (  (G_Point[1])<  (PhyVolGlPos[ijk][ij][1]+layhalfbreadth+delta) &&  (G_Point[1])>  (PhyVolGlPos[ijk][ij][1]-layhalfbreadth-delta) &&  (G_Point[2])<  (PhyVolGlPos[ijk][ij][2]+layhalflength+delta) &&  (G_Point[2])>  (PhyVolGlPos[ijk][ij][2]-layhalflength-delta)       &&     abs(G_Point[0]-PhyVolGlPos[ijk][ij][0])<1.e-10  );
	    //	    cout<<"case 1"<<endl;
	    
	    //	    cout<<"1"<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<endl;
	    
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    "<<    (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   "<<     (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][2])+layhalflength   "<<    (PhyVolGlPos[ijk][ij][2])+layhalflength   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][2])-layhalflength     "<<     (PhyVolGlPos[ijk][ij][2])-layhalflength   <<endl;

   
	    
	    break;
	    
	  case 2: isInside = ( (G_Point[1])<  (PhyVolGlPos[ijk][ij][1]+layhalfbreadth+delta) &&  (G_Point[1])>  (PhyVolGlPos[ijk][ij][1]-layhalfbreadth-delta) &&  (G_Point[2])<  (PhyVolGlPos[ijk][ij][2]+layhalflength+delta) &&  (G_Point[2])>  (PhyVolGlPos[ijk][ij][2]-layhalflength-delta)   &&   abs(G_Point[0] -  PhyVolGlPos[ijk][ij][0])<1.e-10  );
	    //	    cout<<"case 2"<<endl;

	    //	    cout<<"2 "<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<endl;
	    
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    "<<    (PhyVolGlPos[ijk][ij][1])+layhalfbreadth    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   "<<     (PhyVolGlPos[ijk][ij][1])-layhalfbreadth   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][2])+layhalflength   "<<    (PhyVolGlPos[ijk][ij][2])+layhalflength   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][2])-layhalflength     "<<     (PhyVolGlPos[ijk][ij][2])-layhalflength   <<endl;

	    break;
	    
	  case 3: isInside = (  abs(G_Point[0]-PhyVolGlPos[ijk][ij][0])<layhalfbreadth+delta  /* (G_Point[0])<  (PhyVolGlPos[ijk][ij][0]+layhalfbreadth+delta) &&  (G_Point[0])>  (PhyVolGlPos[ijk][ij][0]-layhalfbreadth-delta) &&  (G_Point[2])<  (PhyVolGlPos[ijk][ij][2]+layhalflength+delta) &&  (G_Point[2])> (PhyVolGlPos[ijk][ij][2]-layhalflength-delta)     &&  abs(G_Point[1] -PhyVolGlPos[ijk][ij][1])<1.e-10*/ );
	    //	    cout<<"case 3:: "<<endl;


	    //	    cout<<"3"<<PhyVolGlPos[ijk][ij][0]<<" "<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<" "<<G_Point[0]<<endl;
	    
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][0])+layhalfbreadth    "<<    (PhyVolGlPos[ijk][ij][0])+layhalfbreadth    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][0])-layhalfbreadth   "<<     (PhyVolGlPos[ijk][ij][0])-layhalfbreadth   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][2])+layhalflength   "<<    (PhyVolGlPos[ijk][ij][2])+layhalflength   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][2])-layhalflength     "<<     (PhyVolGlPos[ijk][ij][2])-layhalflength   <<endl;





	    
	    break;



	  case 5: isInside = (  (G_Point[1])<  (PhyVolGlPos[ijk][ij][1]+layhalflength+delta) &&  (G_Point[1])>  (PhyVolGlPos[ijk][ij][1]-layhalflength-delta) &&  (G_Point[2])<  (PhyVolGlPos[ijk][ij][2]+layhalfbreadth+delta) &&  (G_Point[2])>  (PhyVolGlPos[ijk][ij][2]-layhalfbreadth-delta)  &&     abs(G_Point[0]-  PhyVolGlPos[ijk][ij][0])<1.e-10  );
	    //	    cout<<"case 1"<<endl;
	    
	    //	    cout<<"1"<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<endl;
	    
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][1])+layhalflength    "<<    (PhyVolGlPos[ijk][ij][1])+layhalflength    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][1])-layhalflength   "<<     (PhyVolGlPos[ijk][ij][1])-layhalflength   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][2])+layhalfbreadth   "<<    (PhyVolGlPos[ijk][ij][2])+layhalfbreadth   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][2])-layhalfbreadth     "<<     (PhyVolGlPos[ijk][ij][2])-layhalfbreadth   <<endl;

   
	    
	    break;



	  case 6: isInside = (  (G_Point[1])<  (PhyVolGlPos[ijk][ij][1]+layhalflength+delta) &&  (G_Point[1])>  (PhyVolGlPos[ijk][ij][1]-layhalflength-delta) &&  (G_Point[2])<  (PhyVolGlPos[ijk][ij][2]+layhalfbreadth+delta) &&  (G_Point[2])>  (PhyVolGlPos[ijk][ij][2]-layhalfbreadth-delta)       &&     abs(G_Point[0]-  PhyVolGlPos[ijk][ij][0])<1.e-10  );
	    //	    cout<<"case 1"<<endl;
	    
	    //	    cout<<"1"<<PhyVolGlPos[ijk][ij][1]<<" "<<PhyVolGlPos[ijk][ij][2]<<" "<<layhalflength<<" "<<layhalfbreadth<<endl;
	    
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][1])+layhalflength    "<<    (PhyVolGlPos[ijk][ij][1])+layhalflength    <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][1])-layhalflength   "<<     (PhyVolGlPos[ijk][ij][1])-layhalflength   <<endl;
	    //	    cout<<"   (PhyVolGlPos[ijk][ij][2])+layhalfbreadth   "<<    (PhyVolGlPos[ijk][ij][2])+layhalfbreadth   <<endl;
	    //	    cout<<"    (PhyVolGlPos[ijk][ij][2])-layhalfbreadth     "<<     (PhyVolGlPos[ijk][ij][2])-layhalfbreadth   <<endl;

   
	    
	    break;

	    
          default :isInside=false;
	    //	    cout<<"default case"<<endl;
	    break;
          }

	  double trg =0;
	  	  cout<<"isinside: "<<isInside<<endl;

	
	  // bool hitpresent = false;

	  // layexp->SetExtXPos(G_Point[0]);
	  // layexp->SetExtYPos(G_Point[1]);
	  // layexp->SetExtZPos(G_Point[2]);
	  // layexp->SetUsed(true);
	  // layexp->Print();

	  G4ThreeVector extposvec;
	  G4ThreeVector tmphtvec;
	  G4ThreeVector truehtvec;    

          if(isInside) {
	    
	    cout<<"Extrapolation inside the detector boundary"<<endl;


	    
	    //error calculation
	    if(ijk==0){ //top wall
	      Ptxerr = sqrt((xxerr*xxerr)+pow((G_Point[2]-posz),2)*txerr*txerr+2*(G_Point[2]-posz)*xxtxerr);
	      Ptyerr = sqrt((yyerr*yyerr)+pow((G_Point[2]-posz),2)*tyerr*tyerr+ 2*(G_Point[2]-posz)*yytyerr);
	      Ptzerr = 0;
		    
	    }
	    else if(ijk==3) {//for back wall
	      Ptxerr =sqrt((xxerr*xxerr)+pow((G_Point[1]-posy)*txerr/dydz,2)+ pow(((G_Point[1]-posy)*dxdz*tyerr)/(dydz*dydz),2)+       pow((dxdz*yyerr)/dydz,2)            + 2*((G_Point[1]-posy)*xxtxerr)/dydz-  2*pow((G_Point[1]-posy),2)*dxdz*txtyerr/pow(dydz,3) - 2*(G_Point[1]-posy)*dxdz*xxtyerr/pow(dydz,2)    - 2*(dxdz*xxtyerr)/dydz - 2*   (G_Point[1]-posy)*dxdz*yytxerr/pow(dydz,2) +2*dxdz*dxdz* (G_Point[1]-posy)  *yytyerr/pow(dydz,3)  );
	      Ptyerr =0 ;
	      Ptzerr =sqrt(  pow((G_Point[1]-posy)*tyerr/(dydz*dydz),2) + pow((yyerr/dydz),2) + 2*(G_Point[1]-posy)*yytyerr/pow(dydz,3)   ) ;
	    }
		  
	    else { //for left right walls
	      Ptxerr =0 ;
	      Ptyerr =sqrt((yyerr*yyerr)+pow((G_Point[0]-posx)*tyerr/dxdz,2)+ pow(((G_Point[0]-posx)*dydz*txerr)/(dxdz*dxdz),2)+pow((dydz*xxerr)/dxdz,2)          + 2*((G_Point[0]-posx)*yytyerr)/dxdz- 2*(dydz*xxtyerr)/dxdz -  2*pow((G_Point[0]-posx),2)*dydz*txtyerr/pow(dxdz,3)-2*((G_Point[0]-posx)*dydz*xxtyerr)/(dxdz*dxdz)+2*(G_Point[0]-posx)*dydz*dydz*xxtxerr/pow(dxdz,3) - 2*(G_Point[0]-posx)*dydz*yytxerr/pow(dxdz,2));
	      Ptzerr =sqrt(pow((G_Point[0]-posx)*txerr/(dxdz*dxdz),2) +pow((xxerr/dxdz),2) +2* ((G_Point[0]-posx)*xxtxerr)/pow(dxdz,3)) ;
	    }
	    // cout<<"Ptxerr "<<Ptxerr<<" Ptyerr "<<Ptyerr<<"Ptzerr "<<Ptzerr<<endl;

		  

	   
	    // cout<<"Point inside the detector boundary"<<endl;

	   


	    extposvec.setX(G_Point[0]);
	    extposvec.setY(G_Point[1]);
	    extposvec.setZ(G_Point[2]);
 
	    
	     cout<<"Point:"<<extposvec<<endl;
	     cout<<" CmvClusterBank[ijk][ij].size() "<<CmvClusterBank[ijk][ij].size()<<endl;
	    if(CmvClusterBank[ijk][ij].size()>0){
	      cout<<"..check.."<<endl;
            for (unsigned int ix=0; ix<CmvClusterBank[ijk][ij].size(); ix++) {

	      cout<<"## "<<CmvClusterBank[ijk][ij][ix]->GetPlane()-1<<" "<<CmvClusterBank[ijk][ij][ix]->GetLayer()<<" "<<ijk<<" "<<ij<<" "<<ix<<endl;
	      
                
        
		CmvClusterBank[ijk][ij][ix]->Print();


	     	
		//		cout<<"CmvClusterBank[ijk][ij][ix] inside loop "<<CmvClusterBank[ijk][ij][ix]->GetPlane()<<" "<<CmvClusterBank[ijk][ij][ix]->GetLayer()<<endl;
		tmphtvec.setX(CmvClusterBank[ijk][ij][ix]->GetRecoPosX());
		tmphtvec.setY(CmvClusterBank[ijk][ij][ix]->GetRecoPosY());
		tmphtvec.setZ(CmvClusterBank[ijk][ij][ix]->GetRecoPosZ());

		//		cout<<"recohtvec: "<<tmphtvec<<endl;
        	double difx = (extposvec-tmphtvec).mag();

		truehtvec.setX(CmvClusterBank[ijk][ij][ix]->GetTruePosX());
		truehtvec.setY(CmvClusterBank[ijk][ij][ix]->GetTruePosY());
		truehtvec.setZ(CmvClusterBank[ijk][ij][ix]->GetTruePosZ());
		//		cout<<" truehtvvec: "<<truehtvec<<endl;
		//Store the differences in three


		//For multiple clusters in a layer consider one with minimum difx
		if (difx < diffmx) {

		  diffmx = difx;
		  //		  cout<<" Diff "<<difx<<endl;
                
	
     		

		  cout<<" Extrapolated position "<<G_Point[0]<<" "<<G_Point[1]<<" "<<G_Point[2]<<endl;
		  cout<<" Reconstructed position "<<tmphtvec.x()<<" "<<tmphtvec.y()<<" "<<tmphtvec.z()<<endl;
		  cout<<" True Position "<<truehtvec.x()<<" "<<truehtvec.y()<<" "<<truehtvec.z()<<endl;

		  cout<<" Reco-Extrap: "<<  tmphtvec.x()-G_Point[0]<<" "<<  tmphtvec.y()-G_Point[1]<<" "    <<tmphtvec.z()-G_Point[2]<<endl;
		  cout<<" True-Extra: "<<truehtvec.x()-G_Point[0]<<" " << truehtvec.y()-G_Point[1]<<" " << truehtvec.z()-G_Point[2] <<endl;
    


		  
		  pAnalysis->CMVDTruePosX[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosX();
		  pAnalysis->CMVDTruePosY[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosY();
		  pAnalysis->CMVDTruePosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosZ();

		  
		  pAnalysis->CMVDRecoPosX[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosX();
		  pAnalysis->CMVDRecoPosY[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosY();
		  pAnalysis->CMVDRecoPosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosZ();


		  pAnalysis->CMVDWRecoPosX[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosX();
		  pAnalysis->CMVDWRecoPosY[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosY();
		  pAnalysis->CMVDWRecoPosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosZ();


		   
		  pAnalysis->CMVDExpnHit[layid] = 3;		  
		 		  
	  
		}//difx
		
	      

		
            }//for
	    }//if
	   
	    else{
	      pAnalysis->CMVDExpnHit[layid] = 2;
	    }
          }//isinside


	  else{
	    cout<<"..Extrapolated point outside.."<<endl;

	    if(CmvClusterBank[ijk][ij].size()>0){

	      cout<<"But Hit Present"<<endl;
	      for (unsigned int ix=0; ix<CmvClusterBank[ijk][ij].size(); ix++) {
		
		//		cout<<"## "<<CmvClusterBank[ijk][ij][ix]->GetPlane()-1<<" "<<CmvClusterBank[ijk][ij][ix]->GetLayer()<<endl;
	      
                
        
		CmvClusterBank[ijk][ij][ix]->Print();


	     	
		//		cout<<"CmvClusterBank[ijk][ij] inside loop "<<CmvClusterBank[ijk][ij][ix]->GetPlane()<<" "<<CmvClusterBank[ijk][ij][ix]->GetLayer()<<endl;
		tmphtvec.setX(CmvClusterBank[ijk][ij][ix]->GetRecoPosX());
		tmphtvec.setY(CmvClusterBank[ijk][ij][ix]->GetRecoPosY());
		tmphtvec.setZ(CmvClusterBank[ijk][ij][ix]->GetRecoPosZ());

		//		cout<<"recohtvec: "<<tmphtvec<<endl;
        	 double difx = (extposvec-tmphtvec).mag();

		truehtvec.setX(CmvClusterBank[ijk][ij][ix]->GetTruePosX());
		truehtvec.setY(CmvClusterBank[ijk][ij][ix]->GetTruePosY());
		truehtvec.setZ(CmvClusterBank[ijk][ij][ix]->GetTruePosZ());
		cout<<" truehtvvec: "<<truehtvec<<endl;
		//Store the differences in three


		//For multiple clusters in a layer consider one with minimum difx
		if (difx < diffmx) {

		  diffmx = difx;
		  //cout<<" Diff "<<difx<<endl;
                
	
     		

		  cout<<" Extrapolated position "<<G_Point[0]<<" "<<G_Point[1]<<" "<<G_Point[2]<<endl;
		  cout<<" Reconstructed position "<<tmphtvec.x()<<" "<<tmphtvec.y()<<" "<<tmphtvec.z()<<endl;
		  cout<<" True Position "<<truehtvec.x()<<" "<<truehtvec.x()<<" "<<truehtvec.y()<<" "<<truehtvec.z()<<endl;

		  cout<<"Reco-Extrap: "<<  tmphtvec.x()-G_Point[0]<<" "<<  tmphtvec.y()-G_Point[1]<<" "    <<tmphtvec.z()-G_Point[2]<<endl;
		  cout<<"True-Extra: "<<truehtvec.x()-G_Point[0]<<" " << truehtvec.y()-G_Point[1]<<" " << truehtvec.z()-G_Point[2] <<endl;
    


		  
		  pAnalysis->CMVDTruePosX[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosX();
		  pAnalysis->CMVDTruePosY[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosY();
		  pAnalysis->CMVDTruePosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetTruePosZ();

		  
		  pAnalysis->CMVDRecoPosX[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosX();
		  pAnalysis->CMVDRecoPosY[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosY();
		  pAnalysis->CMVDRecoPosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetRecoPosZ();


		  pAnalysis->CMVDWRecoPosX[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosX();
		  pAnalysis->CMVDWRecoPosY[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosY();
		  pAnalysis->CMVDWRecoPosZ[layid] = CmvClusterBank[ijk][ij][ix]->GetWRecoPosZ();


		   
		  pAnalysis->CMVDExpnHit[layid] = 1;		  
		 		  
	  
		}//difx
		
	      

		
	      }
	    }
	    else{
	    pAnalysis->CMVDExpnHit[layid] = 0;// exp point is outside and no hit is there
	    }

	  }


		  pAnalysis->CMVDExpPosX[layid] = G_Point[0];
		  pAnalysis->CMVDExpPosY[layid] = G_Point[1];
		  pAnalysis->CMVDExpPosZ[layid] = G_Point[2];

        

	}
	else{
	  cout<<"..No Line PLane Intersection.."<<endl;




		  
		  pAnalysis->CMVDTruePosX[layid] = -10000;
		  pAnalysis->CMVDTruePosY[layid] = -10000;
		  pAnalysis->CMVDTruePosZ[layid] = -10000;

		  
		  pAnalysis->CMVDRecoPosX[layid] = -10000;
		  pAnalysis->CMVDRecoPosY[layid] = -10000;
		  pAnalysis->CMVDRecoPosZ[layid] = -10000;


		  pAnalysis->CMVDWRecoPosX[layid] = -10000;
		  pAnalysis->CMVDWRecoPosY[layid] = -10000;
		  pAnalysis->CMVDWRecoPosZ[layid] = -10000;



		  pAnalysis->CMVDExpPosX[layid] = -10000;
		  pAnalysis->CMVDExpPosY[layid] = -10000;
		  pAnalysis->CMVDExpPosZ[layid] = -10000;


		   
		  pAnalysis->CMVDExpnHit[layid] = -10000;		  
		 		  

	}

	  
      }//layer loop
  
  
    }//loc_no loop

  }


	
}
micalEventAction::micalEventAction()
  :drawFlag("all"),printModulo(1000),eventMessenger(0)
   //  :drawFlag("charged"),printModulo(1),eventMessenger(0)
{
  eventMessenger = new micalEventActionMessenger(this);
  cal0CollID = -1;
  cal1CollID = -1;//cmv
  nevent = 0;
  rang = 0;
	
  //VtxPlane[]={0};
  //EndPlane[]={0};
  //  cal1CollID = -1;  //  cal2CollID = -1;  //  c0 =0;
	
  //	if (iflag <0) {
  //		iflag = 1;
  paradef = micalDetectorParameterDef::AnPointer;
		
  for(int op=0; op<3;op++) {
    partopscint[op] = paradef->partopscint[op];
  }
	
  AirGapScintTop= paradef->AirGapScintTop;
	
       
  int jmax;
  for(int i =0;i<7;i++){//<4
   
    jmax =  (i==0)? 4:3;
    // if(i==0){
    //   jmax=4;
    // }
    // else{
    //   jmax=3;
    // }
    for(int j=0;j<jmax;j++){
      for(int k=0;k<3;k++){
	PhyVolGlPos[i][j][k] = paradef->ScintLayGPos[i][j][k];
	// cout<< PhyVolGlPos[i][j][k]<<" ";
      }//k
      // cout<<endl;
    }//j
    // cout<<endl<<endl;
  }//i
	
  NoScntStrpTop = paradef->GetNoScntStrpTop();//88
  NoScntStrpSide = paradef->GetNoScntStrpSide();//40
  NoScntStrpSideSmallay= paradef->GetNoScntStrpSideSmallay();	 //8
  SidePlaneHalfLength = paradef->GetSidePlaneHalfLength();
  SideSmallPlaneHalfLength = paradef->GetSideSmallPlaneHalfLength();
  ScntLayShifSide = paradef->GetScntLayShifSide();
  pAnalysis = MultiSimAnalysis::AnPointer;
  // CmvHit_pointer = new CmvHit_Manager();
  // SipmHit_pointer = SipmHit_Manager::APointer; //new InoRPCStrip_Manager();


  
  // CmvStrip_pointer = CmvStrip_Manager::APointer; // added to get atime information

  //	}
	
}

void micalEventAction::CreateCmvHit() {

  
  CmvHit_pointer->CmvHit_list.clear(); 
  
  // Convert SiPM hit to Cmv Hit
  // SipmHit* foursipm[4]={0}; //GMA memory leakage ?
  cout <<"micalEventAction::CreateCmvHit() "<<endl;
  if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput ==3 || pAnalysis->InputOutput==5) {//0:GEN->RECO, 1:GEN->DIGI, 2:GEN->SIM, 3: SIM -> RECO, 4: SIM -> DIFI, 5 : DIGI -> RECO

    cout<<" SipmHit_pointer->SipmHit_list.size() "<<SipmHit_pointer->SipmHit_list.size()<<endl;
    for (unsigned int ijj=0; ijj<SipmHit_pointer->SipmHit_list.size(); ijj++) {
      SipmHit_pointer->SipmHit_list[ijj]->Print();  
    }
   
    
    for (unsigned int ij=0; ij<SipmHit_pointer->SipmHit_list.size(); ij++) {
      //  if(SipmHit_pointer->SipmHit_list[ij]->GetPulse()<160) continue;
      SipmHit* foursipm[4]={0}; //GMA memory leakage ?

      int isfoursipm[4]={0};
      int tmpstripid = -1;
      int tmpside = -1; //Used this to find global position of the layer
      int tmplayer=-1;
      if (!(SipmHit_pointer->SipmHit_list[ij]->GetUsed())) {
	tmpstripid = SipmHit_pointer->SipmHit_list[ij]->GetStripId();
	//	cout <<"tmp "<< tmpstripid<<" "<<int(SipmHit_pointer->SipmHit_list[ij]->GetUsed())<<" "<<SipmHit_pointer->SipmHit_list[ij]->GetSiPM()<<endl;
        int isipm = SipmHit_pointer->SipmHit_list[ij]->GetSiPM();

	isfoursipm[isipm]=1;
	foursipm[isipm] = SipmHit_pointer->SipmHit_list[ij];
	SipmHit_pointer->SipmHit_list[ij]->SetUsed(true);
				
	tmpside = SipmHit_pointer->SipmHit_list[ij]->GetPlane()-1; //We had added 1 while storing it.
	tmplayer = SipmHit_pointer->SipmHit_list[ij]->GetLayer();

	//	cout<<tmpside+1<<" "<<tmplayer<<" "<<isipm<<endl;
	
	//Look for all other SiPM of same strip
	for (unsigned int jk=ij+1; jk<SipmHit_pointer->SipmHit_list.size(); jk++) {
	  // if(SipmHit_pointer->SipmHit_list[jk]->GetPulse()<160) continue;
	  if (!(SipmHit_pointer->SipmHit_list[jk]->GetUsed())) {	
	    int tmpstripid2 = SipmHit_pointer->SipmHit_list[jk]->GetStripId();
	    //  cout <<"tmp2 "<< tmpstripid2 <<endl;
	    int isipm2 = SipmHit_pointer->SipmHit_list[jk]->GetSiPM();
	    int tmpside2 = SipmHit_pointer->SipmHit_list[jk]->GetPlane()-1;
	    int tmplayer2 = SipmHit_pointer->SipmHit_list[jk]->GetLayer();

	    //	    cout<<tmpside2+1<<" "<<tmplayer2<<" "<<isipm2<<endl;
	    if (tmpstripid !=tmpstripid2) continue; //this continues to next iteration
	
	    ////	if(tmpside !=tmpside2 || tmplayer !=tmplayer2) continue;
	    int isipm = SipmHit_pointer->SipmHit_list[jk]->GetSiPM();
	    //	    cout<<"isipm "<<isipm<<endl;
            isfoursipm[isipm]=1;
	    foursipm[isipm] = SipmHit_pointer->SipmHit_list[jk];
	    SipmHit_pointer->SipmHit_list[jk]->SetUsed(true);
	  }
	}
      }
      //      cout<<foursipm[0]<<" "<< foursipm[1]<< " "<< foursipm[2]<<" "<< foursipm[3]<<endl;
      //      cout<<"isfoursipm: "<< isfoursipm[0]<<" "<<isfoursipm[1]<<" "<<isfoursipm[2]<<" "<<isfoursipm[3]<<endl;
      if (tmpside>=0 && isfoursipm[0]+isfoursipm[1]+isfoursipm[2]+isfoursipm[3]>1) { //# atleast 2 sipms must have hit
	double pos[3];
	//	cout<<"posvec: ";
	for (int ix=0; ix<3; ix++) {pos[ix] = PhyVolGlPos[tmpside][tmplayer][ix];cout<<pos[ix]<<" ";}
	//	cout<<endl;

	//	cout<<"creating cmv hit "<<endl;
	if(foursipm[0]){	cout<<foursipm[0]->GetStripId()<<" "<<foursipm[0]->GetSiPM()<< endl;}
	if(foursipm[1]){	cout<<foursipm[1]->GetStripId()<<" "<<foursipm[1]->GetSiPM()<< endl;}
	if(foursipm[2]){	cout<<foursipm[2]->GetStripId()<<" "<<foursipm[2]->GetSiPM()<< endl;}
	if(foursipm[3]){	cout<<foursipm[3]->GetStripId()<<" "<<foursipm[3]->GetSiPM()<< endl;}
	CmvHit* tmpcmvHit = new CmvHit(foursipm[0], foursipm[1], foursipm[2], foursipm[3], pos);
	//	tmpcmvHit->Print();		     
	CmvHit_pointer->CmvHit_list.push_back(tmpcmvHit);
      }
    }
  }







  cout<<endl;

  //uncom
  //     //debug

  //     for (unsigned int ijj=0; ijj<CmvHit_pointer->CmvHit_list.size(); ijj++) {
				
  // int	tmpside0 = CmvHit_pointer->CmvHit_list[ijj]->GetPlane()-1; //We had added 1 while storing it.
  // int	tmplayer0 = CmvHit_pointer->CmvHit_list[ijj]->GetLayer();

  //    if(tmpside0 == 0 && tmplayer0 == 0){
  //    pAnalysis->hist44->Fill(    CmvHit_pointer->CmvHit_list[ijj]->GetTruePosZ()      );
  //    //  pAnalysis->hist55->Fill(    CmvHit_pointer->CmvHit_list[ijj]->GetTruePosZ()       );
    
  //    }

  //   }
  //   //debug



	
}

// cmv hit



void micalEventAction::FormCmvCluster() {
  cout<<"void micalEventAction::FormCmvCluster() {"<<endl;
  CmvCluster_pointer->CmvCluster_list.clear(); 

  // Form cmv cluster

  if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput ==3 || pAnalysis->InputOutput==5) {//0:GEN->RECO, 1:GEN->DIGI, 2:GEN->SIM, 3: SIM -> RECO, 4: SIM -> DIFI, 5 : DIGI -> RECO

    cout<<"CmvHit_pointer->CmvHit_list.size() "<<CmvHit_pointer->CmvHit_list.size()<<endl;
    for (unsigned int ijj=0; ijj<CmvHit_pointer->CmvHit_list.size(); ijj++) {
      CmvHit_pointer->CmvHit_list[ijj]->Print();
    }

  cout<<"...........First Sorting hits acc to stripno in each layer..................."<<endl;

  
    for (unsigned int ix=0; ix<CmvHit_pointer->CmvHit_list.size(); ix++) {
    
      int side1 = CmvHit_pointer->CmvHit_list[ix]->GetPlane()-1;
      int lay1 = CmvHit_pointer->CmvHit_list[ix]->GetLayer();
      int strp1  = CmvHit_pointer->CmvHit_list[ix]->GetStrip();

      
      for (unsigned int ixi=ix+1; ixi<CmvHit_pointer->CmvHit_list.size(); ixi++) {
      int side2 = CmvHit_pointer->CmvHit_list[ixi]->GetPlane()-1;
      int lay2 = CmvHit_pointer->CmvHit_list[ixi]->GetLayer();
      int strp2  = CmvHit_pointer->CmvHit_list[ixi]->GetStrip();

        
	if(side1==side2 && lay1==lay2){

	  
	  if(strp1<strp2){
	    swap(CmvHit_pointer->CmvHit_list[ix],CmvHit_pointer->CmvHit_list[ixi]);	  
	  }//  if(strp1<strp2){
	
	}// 	if(side1==side2 && lay1==lay2){

    	
      }// for (unsigned ixi=ix+1; ixi<CmvHit_pointer->CmvHit_list.size(); ixi++) {
    
    }//  for (unsigned ix=0; ix<CmvHit_pointer->CmvHit_list.size(); ix++) {


    cout<<"................sorted array: ..............."<<endl;
    for (unsigned int ix=0; ix<CmvHit_pointer->CmvHit_list.size(); ix++) {
      CmvHit* cmvhit = CmvHit_pointer->CmvHit_list[ix];
      cmvhit->SetUsed(false);
      cmvhit->Print();
      CmvHitBank[cmvhit->GetPlane()-1][cmvhit->GetLayer()].push_back(cmvhit);
    }//  for (unsigned int ix=0; ix<CmvHit_pointer->CmvHit_list.size(); ix++) {
  

    for(unsigned int tmpside=0;tmpside<6;tmpside++){
      
      for(unsigned int tmplay=0;tmplay<4;tmplay++){
	if(tmpside>0 && tmplay==3) continue;
	//	cout<<"tmpside tmplay hits "<<tmpside<<" "<<tmplay <<" "<<CmvHitBank[tmpside][tmplay].size()<< endl;   
    for(unsigned int jk=0; jk<CmvHitBank[tmpside][tmplay].size();jk++) {
      // cout<<"jk "<<jk<<endl;

      	if(!(CmvHitBank[tmpside][tmplay][jk]->GetUsed())) {

	  CmvCluster* tmpclust = new CmvCluster(CmvHitBank[tmpside][tmplay][jk]);

	  //Change Used from 0 to 1 to show that hit has been added
	  CmvHitBank[tmpside][tmplay][jk]->SetUsed(true);

	  //loop over other entries in same side and same lay
	  for(unsigned int kl=jk+1; kl<CmvHitBank[tmpside][tmplay].size(); kl++) {
	    //	    cout<<"kl "<<kl<<endl;

	    if( abs(CmvHitBank[tmpside][tmplay][kl-1]->GetStrip() - CmvHitBank[tmpside][tmplay][kl]->GetStrip()) >=2 ) continue;//Geometrically nearby  xx x xx 

	    
	    if(  !( CmvHitBank[tmpside][tmplay][kl]->GetUsed())  ){

	      tmpclust->AddHits(CmvHitBank[tmpside][tmplay][kl]);

	      
	      CmvHitBank[tmpside][tmplay][kl]->SetUsed(true);

	      
	    }//  if(  !( CmvHitBank[tmpside][tmplay][kl]->GetUsed())  ){
	  
	    
	    
	  }//	  for(unsigned int kl=jk+1; kl<CmvHitBank[tmpside][tmplay].size(); kl++) {


	  CmvCluster_pointer->CmvCluster_list.push_back(tmpclust);
	  
	}//  	if(!(CmvHitBank[tmpside][tmplay][jk]->GetUsed())) {




    }//  for(unsigned int jk=0; jk<CmvHitBank[tmpside][tmplay].size();jk++) {
      }// for(unsigned int tmplay=0;tmplay<4;tmplay++){
      
    }// for(unsigned int tmpside=0;tmpside<6;tmpside++){
    

    cout<<" CmvCluster_list.size() "<<CmvCluster_pointer->CmvCluster_list.size()<<endl;
    for (unsigned int ix=0; ix<CmvCluster_pointer->CmvCluster_list.size(); ix++) {
      CmvCluster* cmvcluster = CmvCluster_pointer->CmvCluster_list[ix];  
      cmvcluster->Print();
    }//  for (unsigned int ix=0; ix<CmvCluster_pointer->CmvCluster_list.size(); ix++) {
    




    
    

    /*  
    //old
    for (unsigned int ij=0; ij<CmvHit_pointer->CmvHit_list.size(); ij++) {
      CmvHit* clust[2]={0}; //GMA memory leakage ?
      int tmpstrip = -1;
      int tmpside = -1; //Used this to find global position of the layer
      int tmplayer=-1;
      int isclust[2]={0};

       
      if (!(CmvHit_pointer->CmvHit_list[ij]->GetUsed())) {
	//   CmvHit* clust[2]={0}; //GMA memory leakage ?
  
	tmpside = CmvHit_pointer->CmvHit_list[ij]->GetPlane()-1; //Used this to find global position of the layer
	tmplayer=CmvHit_pointer->CmvHit_list[ij]->GetLayer();
	tmpstrip = CmvHit_pointer->CmvHit_list[ij]->GetStrip();
	cout<<ij<<" GetUsed  "<< CmvHit_pointer->CmvHit_list[ij]->GetUsed()<<endl;
	//   if (!(CmvHit_pointer->CmvHit_list[ij]->GetUsed())) {

	//	tmpside = CmvHit_pointer->CmvHit_list[ij]->GetPlane()-1;
	//	tmplayer = CmvHit_pointer->CmvHit_list[ij]->GetLayer();
	//	tmpstrip = CmvHit_pointer->CmvHit_list[ij]->GetStrip();
	cout<<"tmpside: "<<tmpside<<" tmplayer: "<<tmplayer<<" tmpstrip: "<<tmpstrip<<endl;
	clust[0]  = CmvHit_pointer->CmvHit_list[ij];
	isclust[0] = 1;
	CmvHit_pointer->CmvHit_list[ij]->SetUsed(true);
	cout<<endl;
	//Look for neighboring strips in same layer and in same plane
	for (unsigned int jk=ij+1; jk<CmvHit_pointer->CmvHit_list.size(); jk++) {
	  if (!(CmvHit_pointer->CmvHit_list[jk]->GetUsed()) ) {
	    cout<<jk<<" GetUsed2  "<< CmvHit_pointer->CmvHit_list[jk]->GetUsed()<<endl;
	    int	tmpside2 = CmvHit_pointer->CmvHit_list[jk]->GetPlane()-1;
	    int    tmplayer2 = CmvHit_pointer->CmvHit_list[jk]->GetLayer();
	    int tmpstrip2 = CmvHit_pointer->CmvHit_list[jk]->GetStrip();

	    cout<<" tmpside2: "<<tmpside2<<" tmplayer2: "<<tmplayer2<<" tmpstrip2: "<<tmpstrip2<<endl;
	    
	 
	    
	    if (tmpside != tmpside2 || tmplayer!=tmplayer2 || abs(tmpstrip-tmpstrip2)>1 )continue;


	    clust[1] =  CmvHit_pointer->CmvHit_list[jk];
	    isclust[1] = 1;
	    CmvHit_pointer->CmvHit_list[jk]->SetUsed(true);

	  }
	}
      }
      
      if (tmpside>=0 && isclust[0]+isclust[1]>=1) { //#
    
	// cout<<"nearest neighbour found "<<tmpstrip<<" "<<tmpstrip2<<endl;
     
	CmvCluster* tmpcmvCluster = new CmvCluster(clust[0],clust[1]);
	tmpcmvCluster->Print();
	cout<<endl;
	//   if(tmplayer==0){
	//  tmpcmvCluster->SetClustersize(2);//}
	//   else if(tmplayer==1){  tmpcmvCluster->SetClustersizeL1(2);
	//  }
	//    else if(tmplayer==2){  tmpcmvCluster->SetClustersizeL2(2);  }
	CmvCluster_pointer->CmvCluster_list.push_back(tmpcmvCluster);


      }//  if (tmpside>=0 && isclust[0]+isclust[1]>=1) { //#
  

    } //  for (unsigned int ij=0; ij<CmvHit_pointer->CmvHit_list.size(); ij++) {
    
    cout<<endl;

    //... 

    
  
    cout<<"CmvCluster_list size"<<CmvCluster_pointer->CmvCluster_list.size()<<endl;
    for (unsigned int ab=0; ab<CmvCluster_pointer->CmvCluster_list.size(); ab++) {
      CmvCluster_pointer->CmvCluster_list[ab]->SetUsed(false);
      CmvCluster_pointer->CmvCluster_list[ab]->Print();
      cout<<endl;
    }





	
       
    
    //looping over clusters to combine 3 conscutive strip hits
    cout<<"........looping over cluster to combine 3 conscutive strip hits........"<<endl;
    for (unsigned int cd=0; cd<CmvCluster_pointer->CmvCluster_list.size(); cd++) {
      if (!(CmvCluster_pointer->CmvCluster_list[cd]->GetUsed())) {
	CmvCluster* clust2[2]={0}; //GMA memory leakage ?
	int isclust2[2]={0};
 
	cout<<cd<<" GetUsed  "<< CmvCluster_pointer->CmvCluster_list[cd]->GetUsed()<<endl;
     

	int tmpside = CmvCluster_pointer->CmvCluster_list[cd]->GetPlane()-1;
	int tmplayer = CmvCluster_pointer->CmvCluster_list[cd]->GetLayer();
	int tmpstrip = CmvCluster_pointer->CmvCluster_list[cd]->GetStrip();
	cout<<"tmpside: "<<tmpside<<" tmplayer: "<<tmplayer<<" tmpstrip: "<<tmpstrip<<endl;
	clust2[0]  = CmvCluster_pointer->CmvCluster_list[cd];
	isclust2[0] = 1;
	CmvCluster_pointer->CmvCluster_list[cd]->SetUsed(true);
	
	cout<<endl;

	//Look for neighboring strips in same layer and in same plane
	for (unsigned int dc=cd+1; dc<CmvCluster_pointer->CmvCluster_list.size(); dc++) {
	  if (!(CmvCluster_pointer->CmvCluster_list[dc]->GetUsed()) ) {
	    cout<<dc<<" GetUsed2  "<< CmvCluster_pointer->CmvCluster_list[dc]->GetUsed()<<endl;
	    int	tmpside2 = CmvCluster_pointer->CmvCluster_list[dc]->GetPlane()-1;
	    int    tmplayer2 = CmvCluster_pointer->CmvCluster_list[dc]->GetLayer();
	    int tmpstrip2 = CmvCluster_pointer->CmvCluster_list[dc]->GetStrip();

	    cout<<" tmpside2: "<<tmpside2<<" tmplayer2: "<<tmplayer2<<" tmpstrip2: "<<tmpstrip2<<endl;
	    
	 
	    
	    if (tmpside == tmpside2 && tmplayer==tmplayer2 && abs(tmpstrip-tmpstrip2)==2 ){

	      cout<<"Found 3 consecutive hits.."<<endl;
	      clust2[1] =  CmvCluster_pointer->CmvCluster_list[dc];

	      CmvCluster_pointer->CmvCluster_list[cd]->CombineClusts(clust2[1]);

	  
	      //clust2[1] =  CmvCluster_pointer->CmvCluster_list[dc];
	      //  if(tmplayer==tmplayer2==0){
	      //   CmvCluster_pointer->CmvCluster_list[cd]->SetClustersize(3);
	      //	}
	      //   else if(tmplayer==tmplayer2==1){  	  CmvCluster_pointer->CmvCluster_list[cd]->SetClustersizeL1(3);  }
	      //    else if(tmplayer==tmplayer2==2){  	  CmvCluster_pointer->CmvCluster_list[cd]->SetClustersizeL2(3);  }
    
 
     
	      //  CmvCluster* tmpcmvCluster = new CmvCluster(clust2[0],clust2[1]);
	      //  tmpcmvCluster->Print();
	      //  cout<<endl;		     
	      //  CmvCluster_pointer->CmvCluster_list.push_back(tmpcmvCluster);

	      CmvCluster_pointer->CmvCluster_list.erase(CmvCluster_pointer->CmvCluster_list.begin()+dc);
	      //	CmvCluster_pointer->CmvCluster_list.erase(CmvCluster_pointer->CmvCluster_list.begin()+cd);
	      dc--;
	      //	cd--;

	
	    }
	  }
	}
      }

    } //  for (unsigned int cd=0; cd<CmvCluster_pointer->CmvCluster_list.size(); cd++) {
    
    cout<<endl;

    

    cout<<"CmvCluster_list size"<<CmvCluster_pointer->CmvCluster_list.size()<<endl;
    for (unsigned int ab=0; ab<CmvCluster_pointer->CmvCluster_list.size(); ab++) {
      CmvCluster_pointer->CmvCluster_list[ab]->SetUsed(false);
      CmvCluster_pointer->CmvCluster_list[ab]->Print();
      cout<<endl;
    }




    //debug

    for (unsigned int ijj=0; ijj<CmvCluster_pointer->CmvCluster_list.size(); ijj++) {
				
      int	tmpside0 = CmvCluster_pointer->CmvCluster_list[ijj]->GetPlane()-1; //We had added 1 while storing it.
      int	tmplayer0 = CmvCluster_pointer->CmvCluster_list[ijj]->GetLayer();
      int clustersize0 =   CmvCluster_pointer->CmvCluster_list[ijj]->GetClustersize();
      if(tmpside0 == 1 && tmplayer0 == 0 && clustersize0==1){
	//  pAnalysis->hist55->Fill(    CmvCluster_pointer->CmvCluster_list[ijj]->GetTruePosX()      );
	// pAnalysis->hist55->Fill(    CmvCluster_pointer->CmvCluster_list[ijj]->GetTruePosX()       );
    

      }
    }

    //debug








	

    */
  

  } //  if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput ==3 || pAnalysis->InputOutput==5) {	
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalEventAction::~micalEventAction() {
  delete eventMessenger;




  // for (unsigned ij=0; ij<CmvHit_pointer->CmvHit_list.size(); ij++) {
  //   if (CmvHit_pointer->CmvHit_list[ij]) {
  //     cout <<"ij "<< ij<<" "<<CmvHit_pointer->CmvHit_list.size()<<endl;
  //     delete CmvHit_pointer->CmvHit_list[ij]; CmvHit_pointer->CmvHit_list[ij]=0;
  //   }
  // }

  // CmvHit_pointer->CmvHit_list.clear();
  // if (CmvHit_pointer) {delete CmvHit_pointer; CmvHit_pointer=0;}








  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Orighit algo has been included as  a function so that we can use more than once
//only by calling this function
//added by SSE 09/15
int micalEventAction::orighit_calc(vector<InoHit*> tmphitcluster){
  vector<iXYZ> xShwstrip;
  vector<iXYZ> yShwstrip;
	
  xShwstrip.clear();yShwstrip.clear(); 
  // cout<< " xShwstrip.size "<< xShwstrip.size()<< " yShwstrip.size "<< yShwstrip.size() <<endl;
  for( unsigned int mn=0; mn<tmphitcluster.size(); mn++) {
    x_stripno = tmphitcluster[mn]->GetXStripNum();
    y_stripno = tmphitcluster[mn]->GetYStripNum();
    z_plane=tmphitcluster[mn]->GetZPlane();
		
    if (x_stripno && y_stripno) {
      //cout<<" 09/15 orig func chk: "<< " "<<z_plane << " "<<x_stripno<< " "<<y_stripno<< " UID "<<tmphitcluster[mn]->GetUID()<<endl;
      for (unsigned nn=0; nn< xShwstrip.size(); nn++) {
	if (x_stripno == xShwstrip[nn].second &&
	    z_plane == xShwstrip[nn].first) {
	  // cout<< "09/15: ix working"<<endl;
	  x_stripno = -1; break;
	}
      }
      for (unsigned nn=0; nn< yShwstrip.size(); nn++) {
	if (y_stripno == yShwstrip[nn].second &&
	    z_plane == yShwstrip[nn].first) {
	  //cout<< "09/15: iy working"<<endl;
	  y_stripno = -1; break;
	}
      }
      if ((x_stripno>=0 || y_stripno>=0) ) {
	if(x_stripno>=0){
	  //cout<< " 09/15: ixstripno "<< ixstripno<<endl;
	  iXYZ ZXstrp(z_plane,x_stripno);xShwstrip.push_back(ZXstrp);
	}
	if(y_stripno>=0){
	  //cout<< " 09/15: iystripno "<< iystripno<<endl;
	  iXYZ ZYstrp(z_plane,y_stripno);yShwstrip.push_back(ZYstrp);
	}
      }//end of if ((ixstripno>=0 || iystripno>=0) )
    }
  }//end of for( int mn=0; mn<tcluster.size(); mn++) loop on total hits
  // cout<< " xShwstrip.size "<< xShwstrip.size()<< " yShwstrip.size "<< yShwstrip.size() <<endl;
  int x_nhits =0;
  int y_nhits =0;
  int Orighits_all=0;
  for(int kl=0; kl<nLayer; ++kl){

    if(xShwstrip.size() && yShwstrip.size()){
      x_nhits =0;
      y_nhits =0;
      for(unsigned nn=0; nn< xShwstrip.size(); nn++) {
	if(xShwstrip[nn].first==int(kl)){
	  x_nhits++;
	}
      }
      for (unsigned nn=0; nn< yShwstrip.size(); nn++){
	if(yShwstrip[nn].first==int(kl)){
	  y_nhits++;
	}
      }
      Orighits_all +=max(x_nhits,y_nhits);
    }
  }
  return Orighits_all;
}//SSE 09/15

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int Evt_Counts=0;
void micalEventAction::BeginOfEventAction(const G4Event* evt) {
  cout<<"......Begin of event action..........."<<endl;

  cout<<"Event_No"<<Evt_Counts<<endl;

  cout<<"check1"<<endl;
  inoStripX_pointer =  InoStripX_Manager::APointer;
  inoStripY_pointer = InoStripY_Manager::APointer;
  cout<<"check2"<<endl;


  CmvHit_pointer = new CmvHit_Manager();
  CmvHit_pointer->CmvHit_list.clear();
  cout<<"check3"<<endl;
  SipmHit_pointer = SipmHit_Manager::APointer; //new InoRPCStrip_Manager();
  cout<<"check4"<<endl;

  
  CmvStrip_pointer = CmvStrip_Manager::APointer; // added to get atime information

  for(int ijj=0;ijj<7;ijj++){
    for(int ikk=0;ikk<4;ikk++){
      CmvHitBank[ijj][ikk].clear();
    }
    
  }
  CmvCluster_pointer = new CmvCluster_Manager();
  CmvCluster_pointer->CmvCluster_list.clear();

  cout<<"check5"<<endl;

	
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    //HepRandom::showEngineStatus();
  }
  
  rang = 0.0;
  //initialisation per event
  Energycal0 = EnergyGap = 0.;
  TrackLcal0 = TrackLGap = 0.;
  Energycal1 = Energycal2 = 0.;
  TrackLcal1 = TrackLcal2 = 0.;
  
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  //  if(cal0CollID<0||cal1CollID<0||cal2CollID<0)
  
  if(cal0CollID<0||cal1CollID<0) {
    G4String colNam;
    cal0CollID = SDman->GetCollectionID(colNam="cal0Collect");
    cal1CollID = SDman->GetCollectionID(colNam="cal1Collect");//cmv
    //    cal2CollID = SDman->GetCollectionID(colNam="cal2Collect");
  }
  Evt_Counts++;


  for(int jki=0;jki<10;jki++){
	  pAnalysis->XPosdev_exclu[jki] = 100000; 
	  pAnalysis->YPosdev_exclu[jki] = 100000;
  }







 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void micalEventAction::EndOfEventAction(const G4Event* evt) {

  cout<<"......end of event action..........."<<endl;
  //  paradef = micalDetectorParameterDef::AnPointer;
  //StripXWidth = (1/m)*paradef->GetXStrwd();
  //StripYWidth = (1/m)*paradef->GetYStrwd();
  nLayer      = paradef->GetnLayer();
  LayerThickness = (1/m)*2*(paradef->GetParlay(2)+paradef->GetParirlay(2));
  
  typedef pair<int,int> ixyz;
  nevent++;
  /*MultiSimAnalysis */ pAnalysis = MultiSimAnalysis::AnPointer;
  InoRPCStrip_Manager* inoRPC_pointer = new InoRPCStrip_Manager(); //InoRPCStrip_Manager::APointer;

  pAnalysis->range = rang; //meghna
  pAnalysis->nhtcal0 = 0;
  
  G4int evtNb = evt->GetEventID();
  //  cout<<"# "<<evtNb<<endl;
  // pAnalysis->timeAsciiOutput<<"# "<<evtNb<<endl;
  if (evtNb%printModulo == 0) {
    if(pAnalysis->isXtermOut==1) {
      //    G4cout << "---> End of event: " << evtNb << G4endl;
      
      G4cout
	<< G4endl
	<< "    Absrober: total energy: " << std::setw(7)
	<< G4BestUnit(EnergyGap,"Energy")
	<< "        total track length: " << std::setw(7)
	<< G4BestUnit(TrackLGap,"Length")
	<< G4endl
	<< "RPC_gas: cal0 total energy: " << std::setw(7)
	<< G4BestUnit(Energycal0,"Energy")
	<< "        total track length: " << std::setw(7)
	<< G4BestUnit(TrackLcal0,"Length")
	<< G4endl;
    }
    pAnalysis->simtotabenr = EnergyGap;
    pAnalysis->simtotrpcenr= Energycal0;
    pAnalysis->simtotablen = TrackLGap;
    pAnalysis->simtotrpclen= TrackLcal0;
  }
  
  if (pAnalysis->InputOutput <=2) {
    G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
    if (pVisManager) {
      G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
      
      for (G4int ij=0; ij<n_trajectories; ij++) {
        G4VTrajectory* trj = ((*(evt->GetTrajectoryContainer()))[ij]);
        if (drawFlag == "all") { pVisManager->Draw(*trj); // GMA14 ,100);
        } else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.)) {
          pVisManager->Draw(*trj); // GMA14 ,100);
        } else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.)) {
          pVisManager->Draw(*trj); //GMA14 ,100);
        }
      }
      
      if(cal0CollID<0 || cal1CollID<0) return;
      //      if(cal0CollID<0) return;

      G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
      micalcal0HitsCollection* EHC0 = 0;
      micalcal1HitsCollection* EHC1 = 0;//cmv
      //  micalcal2HitsCollection* EHC2 = 0;
      if(HCE) {
        EHC0 = (micalcal0HitsCollection*)(HCE->GetHC(cal0CollID));
        EHC1 = (micalcal1HitsCollection*)(HCE->GetHC(cal1CollID));//cmv
        //    EHC2 = (micalcal2HitsCollection*)(HCE->GetHC(cal2CollID));
      }
      if(EHC0) {
        int n_hit = EHC0->entries();
        //    pAnalysis->ascii_output <<n_hit<<" ";
        
        G4double totE = 0;
        for (int ij=0; ij<n_hit; ij++) {
          totE +=(*EHC0)[ij]->GetEdep();
          //G4cout <<"Energy deposited in hit "<<i<<" "<<(*EHC0)[ij]->GetEdep()<< " at "<<(*EHC0)[ij]->GetPos()<<" at "<<(*EHC0)[ij]->GetTime()<<" "<<(*EHC0)[ij]->GetHitId() <<G4endl;
          if (pAnalysis->nhtcal0 <pAnalysis->nmxhit) {
            pAnalysis->calid0[pAnalysis->nhtcal0] = (*EHC0)[ij]->GetHitId();
            pAnalysis->calen0[pAnalysis->nhtcal0] = (*EHC0)[ij]->GetEdep()/keV;
            //G4ThreeVector MCpos = (*EHC0)[i]->GetPos();
            //cout<<"x "<<1.e-3*MCpos.x()<<"     y "<<1.e-3*MCpos.x()<<"     z "<<1.e-3*MCpos.z()<<endl;
            pAnalysis->nhtcal0++;
            //pAnalysis->ascii_output <<std::setw(12)<<(*EHC0)[i]->GetHitId()<<" "<<std::setw(9)<<(*EHC0)[i]->GetEdep()/keV;
          }
        }
        
        pAnalysis->caltot0 = totE/keV;
        //    pAnalysis->ascii_output <<"  "<<totE/keV<<G4endl;
        if(pAnalysis->isXtermOut==1) {
          G4cout << "     " << n_hit
                 << " hits are stored in micalEcalHitsCollection with total Energy  "<<totE<< " KeV"<<endl; // G4BestUnit(totE,"Energy") << G4endl;
        }
      }
    } // if (pVisManager)
  } // if (pAnalysis->InputOutput <=2)
  
  if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput==3 || pAnalysis->InputOutput==5) {
    //    InoPatternRecognition pattern; //GMA14
    //    pattern.RunPattRecog();
    
    InoTrackFinder trackfinder;
    // cout<<"Running the finder now ...."<<endl;
    trackfinder.RunTheFinder();	//VALGRIND
    // cout<<"Completed..."<<endl;
    //Now add other cluster close to vertex;
    vector <InoCluster*> totcluster = trackfinder.inoCluster_pointer->InoCluster_list;
    int totclustersize = totcluster.size();




    cout<<"totalclustersize = "<<totclustersize<<endl;


    InoTrack_Manager *pinotrack = InoTrack_Manager::APointer;

    inoTrackCand_pointer = new InoTrackCand_Manager();//earlier defined in InoTrackFitAlg constructor

    if (pinotrack) {

       CorrTimeError = pAnalysis->GetCorrTimeError();
        UnCorrTimeError = pAnalysis->GetUnCorrTimeError();
        timeerr=pow((pow(CorrTimeError,2.) + pow(UnCorrTimeError,2.)),0.5);
      bool ZIncreasesWithTime = true;
	
      double szxy=0, sz=0, sxy=0, sn=0, sz2=0;
      double errsq=timeerr*timeerr;
      //  double dZdT = 0.0; declared in hh

      // double ztinter = 0;
      for (unsigned iji=0; iji<pinotrack->InoTrack_list.size() ; iji++) {
	int cluster_size = pinotrack->InoTrack_list[iji]->ClustsInTrack.size();
	
	 

 
	for (int jk =0; jk<cluster_size;jk++) {
	  cout<<pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetZPos()<<" "<<pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetTime()<<endl;
	  szxy +=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetZPos()*pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetTime()/errsq;
	  sz +=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetZPos()/errsq;
	  sz2 +=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetZPos()*pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetZPos()/errsq;
	  sxy +=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetTime()/errsq;
	  sn +=1/errsq;
    
	}// for ( unsigned int jk =0; jk<cluster_size;jk++) {

   
	if ((sz2*sn - sz*sz) !=0) { 
	  dZdT = (szxy*sn - sz*sxy)/(sz2*sn - sz*sz);

	  ztinter  = sxy/sn - dZdT*sz/sn;


	}
	cout<<"dZdT:ztinter "<<dZdT<<" "<<ztinter<<endl;
  
	if(dZdT>0){
	  ZIncreasesWithTime= true;
	}

	else{
	  ZIncreasesWithTime=false;

	}
  
	cout<<"Zincreaseswith time: "<<ZIncreasesWithTime<<endl;
	//   //..
      }
      //

      //   int Nt;
      //   const int  layfirst =0; //Used first layer in track fitting
      //   const int  laylast =9; //Used last layer in track fitting
      //   int occulyr=-1;
      //   int nlayer=11;
      //   const float xyPosDev=3.0; // seven sigma 2.0; //maximum deviation of points from fit line
   
   
      //   double Tpos[nlayer]={0}; 
      //   bool Tusedpos[nlayer]={0};
      //   double Tdev[nlayer]; for (int ij=0; ij<nlayer; ij++) { Tdev[ij] = 100; Tpos[ij]=-5000;}

      //   double CorrTimeError = pAnalysis->GetCorrTimeError();
      //   double  UnCorrTimeError = pAnalysis->GetUnCorrTimeError();
      //   double  timeerr=pow((pow(CorrTimeError,2.) + pow(UnCorrTimeError,2.)),0.5);
 
      //   float errtco[nlayer];
      //   for(int ab=0;ab<nlayer;ab++){errtco[ab]=timeerr;  } 


      //   double xxerr[nlayer];
   
      //   for ( int ix=0; ix<nlayer; ix++) {
      // 	xxerr[ix] = errtco[ix]*errtco[ix];
      //   }


      //   double tchi2;
      //   int nmnhits =5;
      //   int mtchisq =2;
      //   double terrcst, terrcov, terrlin;


      //   Nt=0;
      //   int ntfail = 0;
      //   tchi2 = 0;
      //   double tresol = 0;
      //   double zval[nlayer], text[nlayer], texter[nlayer], tposinstr[nlayer], valz[nlayer];

      //   for (int ix=0; ix<nlayer; ix++) { text[ix]= texter[ix] =tposinstr[ix] =  100;}

   
      //   //..
      //   for (unsigned ij=0; ij<pinotrack->InoTrack_list.size() ; ij++) {
	    
      // 	int cluster_size=pinotrack->InoTrack_list[ij]->ClustsInTrack.size();
	    
      // 	for ( unsigned int jk =0; jk<cluster_size;jk++) {
      // 	  //	 cout<<"jk: "<<jk<<endl;
      // 	  //	 cout<<"pointer cluster"<<pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]<<endl;
 
		
      // 	  zval[jk]=pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetZPos();
      // 	  Tpos[jk]=pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetTime();
        
        
      // 	  Tusedpos[jk]=true;
 
      // 	  cout<<" "<<Tpos[jk]<<" "<<zval[jk]<<" "<<endl;
        
 
	      
      // 	}// for ( unsigned int jk =0; jk<nlayer;jk++) {
	    
   
	    
 
	    
      // 	//............RPC Straight line fit in T-Z plane.............
       
      // 	StraightLineFit tposfit(0, zval, Tpos,  xxerr, Tusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
    
	    
	    
      // 	tposfit.GetParameters(ntfail, tinters, tslope);
      // 	cout<<"Slope and intercept TX  "<<tslope<<" "<<tinters<<" "<<Tpos[9]<<" "<<zval[9]<<" "<<zval[9]*tslope+tinters<<"  "<<(-1.081)*tslope+tinters<<endl;

      // 	//	pAnalysis->hist11->Fill(Tpos[9]);
      // 	//	pAnalysis->hist22->Fill(zval[9]*tslope+tinters);     
      // 	tposfit.GetError(terrcst, terrlin, terrcov);
      // 	// cout<<"Error in slope and intercept TX "<<terrcst<<" "<<terrlin<<endl;
      // 	tposfit.GetChisqure(Nt, tchi2);
      // 	// cout<<"NDOF and chi-square T "<<Nt<<" "<<tchi2<<endl;
       
      // 	tposfit.GetFitValues(text,valz, Tdev, texter);
      // 	// for(int jk = 0;jk<nlayer;jk++){
      // 	// cout<<"get t fit values "<<text[jk]<<" "<<Tpos[jk]<<" "<<valz[jk]<<" "<<endl;
      // 	// }

      // 	// store events only when they satisfy:
       
      // 	if (Nt<=nmnhits/*-ntcor*/ && tchi2/(Nt-2)>mtchisq && ntfail!=0) {

      // 	  tslope=-1000000;
      // 	  tinters = -1000000;


      // 	}

      //   }//for (unsigned ij=0; ij<pinotrack->InoTrack_list.size() ; ij++) {








 
      //% 
      // cout<<"pAnalysis->InoTrack_listsize->Fill(pinotrack->InoTrack_list.size()); "<<pinotrack->InoTrack_list.size()<<endl;
      pAnalysis->InoTrack_listsize->Fill(pinotrack->InoTrack_list.size());        //asm : from here
      cout<<"Hola..."<<endl;
      //default version has isVisOut= 1, but we change that to 0 to obtain hit information in output without generating .inh file
      if (pAnalysis->isVisOut==1)  {
        pAnalysis->H->NFinders=pinotrack->InoTrack_list.size(); // number of finder sets
        for (unsigned ij=0; ij<pinotrack->InoTrack_list.size() ; ij++) {
          for ( unsigned int jk =0; jk<pinotrack->InoTrack_list[ij]->ClustsInTrack.size();jk++) {
            pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object //VALGRIND
            pAnalysis->Hp->TrackType=-4;// Track Type: -1: hits, -2: clulster, -3: triplet, -4: track
            pAnalysis->Hp->FindNum=ij;// track Number
            pAnalysis->Hp->ZZ=pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetZPlane();
            pAnalysis->Hp->XX=pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetXPos();
            pAnalysis->Hp->YY=pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetYPos();
	    cout<<".........HP........."<<endl;
	    cout<<pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetXPos()<<" "<<pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetZPlane()<<endl;

	  }
        } //asm : upto here
        
        if (pinotrack->InoTrack_list.size()>0) {
          for (unsigned int jk =0; jk<pinotrack->InoTrack_list[0]->ClustsInTrack.size();jk++) {
            //GMA where is it ptinting ?
            //cout<<"pAnalysis->TrkDist->Fill(pinotrack->InoTrack_list[0]->ClustsInTrack[jk]->GetZPlane());"<<endl;
            pAnalysis->TrkDist->Fill(pinotrack->InoTrack_list[0]->ClustsInTrack[jk]->GetZPlane());
            //cout<<"pAnalysis->EffDist->Fill(pinotrack->InoTrack_list[0]->ClustsInTrack[jk]->GetZPlane());"<<endl;
            pAnalysis->EffDist->Fill(pinotrack->InoTrack_list[0]->ClustsInTrack[jk]->GetZPlane());
          }
        }
      }

      //      cout<<"Cluster info"<<endl;
      // for (unsigned ij=0; ij<pinotrack->InoTrack_list.size() ; ij++) {
      //          for ( unsigned int jk =0; jk<pinotrack->InoTrack_list[ij]->ClustsInTrack.size();jk++) {
	    
	    
      // 	    cout<<pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetZPlane()<<" "<<pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetZPos()<<" "<<pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetXPos()<<" "<<pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetYPos()<<" "<<  pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetTime()    <<endl;
      // 	  }

      // }

 



      micalElectroMagneticField* inoical0Field=new micalElectroMagneticField();
 
      magfield = inoical0Field->GetConstantFieldvalue();
      magfield=0;
      cout<<" Mag field from event action "<<magfield<<endl;
    
      InoTrackFitAlg trackfitter;
      if(magfield==0){

	cout<<".......Entered RPC  straight line fit.........."<<endl;


	inoTrackCand_pointer->InoTrackCand_list.clear();
	double xvtx_parameter[6]={0};
	InoTrackCand* fTrackCand;
	cout<<"fTackCand pointer "<<fTrackCand<<endl;
     
	cout<<"inotrack list size"<<pinotrack->InoTrack_list.size()<<endl;

	unsigned int itrk=0;
	 occulyr=-1;
	for (unsigned int iji=0; iji<pinotrack->InoTrack_list.size() ; iji++) {
	  if(iji<pAnalysis->ntrkmx){ //ntrkmx is 20
	  double zcor;
	  int cluster_size=pinotrack->InoTrack_list[iji]->ClustsInTrack.size();


	  int nhits1 = pinotrack->InoTrack_list[iji]->GetEntries();
	  cout<<"nhits1: "<<nhits1<<endl;

	   Xpos[nlayer]={0.0}; 
	   Xusedpos[nlayer]={0};
	   for (int ij1=0; ij1<nlayer; ij1++) { Xdev[ij1] = 100; Xpos[ij1]=-5000;}
	   

	   Ypos[nlayer]={0.0};
	   Yusedpos[nlayer]={0};//=new float[nlayer];
	   for (int ij2=0; ij2<nlayer; ij2++) { Ydev[ij2] = 100; Ypos[ij2]=-5000.0;}
	    
	   
	   CorrTimeError = pAnalysis->GetCorrTimeError();
	   UnCorrTimeError = pAnalysis->GetUnCorrTimeError();
	   timeerr=pow((pow(CorrTimeError,2.) + pow(UnCorrTimeError,2.)),0.5);
 
   


	  Nx=0;
	  int nxfail = 0;
	  xchi2 = 1000;
	  double xresol = 0;
	  double zval[nlayer], xext[nlayer], xexter[nlayer], xposinstr[nlayer], valz[nlayer];
	  for (int ij3=0; ij3<nlayer; ij3++) { zval[ij3] = 0;}

	  Ny=0;
	  int nyfail = 0;
	  ychi2 = 1000;
	  double yresol = 0;
	  double yext[nlayer], yexter[nlayer], yposinstr[nlayer];



   
	  for (int ixx=0; ixx<nlayer; ixx++) { xext[ixx]= xexter[ixx] =xposinstr[ixx] =  100;}
	  for (int iyy=0; iyy<nlayer; iyy++) { yext[iyy]= yexter[iyy] =yposinstr[iyy] =  100;}


	  for (int nlay =0; nlay<nlayer;nlay++) {
	    ClustsInTrackBank[iji][nlay].clear();
	      
	  }


	  cout<<"# of Clusters: "<<pinotrack->InoTrack_list[iji]->ClustsInTrack.size()<<endl;

	  for (unsigned int ix=0; ix<pinotrack->InoTrack_list[iji]->ClustsInTrack.size(); ix++) {
	    InoCluster* clust = pinotrack->InoTrack_list[iji]->ClustsInTrack[ix];
	    clust->Print();
	    cout<<"clust->GetZPlane() "<<clust->GetZPlane()<<" "<<ClustsInTrackBank[iji][clust->GetZPlane()].size()<<endl;
	    ClustsInTrackBank[iji][clust->GetZPlane()].push_back(clust);
	 
	  
	  }// for (unsigned int ix=0; ix<Cluster_pointer->InoStripX_list.size(); ix++) {
	  cout<<"raj "<<nlayer<<endl;
	  for (int nlay =0; nlay<nlayer;nlay++) {
	    cout<<"ClustsInTrackBank[nlay].size() "<<ClustsInTrackBank[iji][nlay].size()<<endl;
	      
	  }
  float errxco[12]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};
  float erryco[12]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};
     
	  for (int nlay =0; nlay<nlayer;nlay++) {
        
	    //Only one clust per layer and that cluster should have less than 5 multiplicity:
	    //	    cout<<"nxstrip nystrip  "<<ClustsInTrackBank[nlay][0]->GetNXStripsInClust() << " "<< ClustsInTrackBank[nlay][0]->GetNYStripsInClust() <<endl;

	    if(ClustsInTrackBank[iji][nlay].size()>=1){
	      cout<<"ClustsInTrackBank[iji][nlay].size() "<<ClustsInTrackBank[iji][nlay].size()<<endl;
		
	      if (ClustsInTrackBank[iji][nlay].size()==1 &&  ClustsInTrackBank[iji][nlay][0]->GetNXStripsInClust()<5 && ClustsInTrackBank[iji][nlay][0]->GetNYStripsInClust()<5 ){
		cout<<"nxstrip nystrip  "<<ClustsInTrackBank[iji][nlay][0]->GetNXStripsInClust() << " "<< ClustsInTrackBank[iji][nlay][0]->GetNYStripsInClust() <<endl;

	      
		zval[nlay]=ClustsInTrackBank[iji][nlay][0]->GetZPos();
		Xpos[nlay]=ClustsInTrackBank[iji][nlay][0]->GetXPos();
		Ypos[nlay]=ClustsInTrackBank[iji][nlay][0]->GetYPos();
		 errxsq[nlay]=ClustsInTrackBank[iji][nlay][0]->GetXPosErr() * ClustsInTrackBank[iji][nlay][0]->GetXPosErr() ;
		 errysq[nlay]=ClustsInTrackBank[iji][nlay][0]->GetYPosErr() * ClustsInTrackBank[iji][nlay][0]->GetYPosErr() ;

        

		
		Xusedpos[nlay]=true;
		Yusedpos[nlay]=true;


	      }
	      else{
	      
		zval[nlay]=-5000.0;
		Xpos[nlay]=-5000.0;
		Ypos[nlay]=-5000.0;
		// errxsq[nlay]=5000.0;
		// errysq[nlay]=5000.0;
		errxsq[nlay]=errxco[nlay]*errxco[nlay]*0.03*0.03;
		errysq[nlay]=erryco[nlay]*erryco[nlay]*0.03*0.03;

		Xusedpos[nlay]=false;
		Yusedpos[nlay]=false;
	      }  
	    }
	      cout<<nlay<<" "<<Xpos[nlay]<<" "<<Ypos[nlay]<<" "<<zval[nlay]<<" "<<Xusedpos[nlay]<<" "<<Yusedpos[nlay]<<" "<<errxsq[nlay]<<" "<<errysq[nlay]<<endl;
	    
	 
	  }//  for (int nlay =0; nlay<10;nlay++) {

	  cout<<"cluster size: "<<cluster_size<<endl;
	  // for (int jk =0; jk<cluster_size;jk++) {
	  	
	  //   zval[jk]=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetZPos();
	  //   Xpos[jk]=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetXPos();
	  //   Ypos[jk]=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetYPos();
	  //   errxsq[jk]=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetXPosErr() * pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetXPosErr() ;
	  //   errysq[jk]=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetYPosErr() * pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetYPosErr() ;

	  //   Xusedpos[jk]=true;
	  //   Yusedpos[jk]=true;
	  //   cout<<jk<<" "<<Xpos[jk]<<" "<<Ypos[jk]<<" "<<zval[jk]<<" "<<Xusedpos[jk]<<" "<<Yusedpos[jk]<<" "<<errxsq[jk]<<" "<<errysq[jk]<<endl;
        
 
	      
	  // }// for ( unsigned int jk =0; jk<nlayer;jk++) {




	  //Finding Topmost Layer hit:

	  double zposmx = 10000;
	  double topmostlay;
	  for (int jk =0; jk<nlayer;jk++) {
	    if (Xusedpos[jk]==false || Yusedpos[jk]==false) continue;	  
	    topmostlay=ClustsInTrackBank[iji][jk][0]->GetZPos();
	    if(topmostlay<1000){
	      zposmx = topmostlay;
	    }
	  }

	  cout<<"TopMost layer having hit is: "<<zposmx<<endl;



	  
   
	  double posresol =true;
	  if(posresol){

	    itrk++;
	    cout<<"...Position resolution..."<<endl;
	    for(int jki=0;jki<nlayer;jki++){
	      occulyr = jki;

	      //  if (Xusedpos[jki]==false || Yusedpos[jki]==false) continue;
	      cout<<"Fitting Excluding Layer:"<<occulyr<<endl;

	      StraightLineFit xposresolfit(1, zval, Xpos,  errxsq, Xusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
	      xposresolfit.GetParameters(nxfail, xinters, xslope);
	      cout<<"Slope and intercept X  "<<xslope<<" "<<xinters<<endl;
	      xposresolfit.GetError(xerrcst, xerrlin, xerrcov);
	      cout<<"Error in slope and intercept X "<<xerrcst<<" "<<xerrlin<<endl;
	      xposresolfit.GetChisqure(Nx, xchi2);
	      cout<<"NDOF and chi-square X "<<Nx<<" "<<xchi2<<endl;
	      xposresolfit.GetFitValues(xext,valz, Xdev, xexter);
	      for(int jk = 0;jk<nlayer;jk++){
		cout<<"get xz fit values "<<xext[jk]<<" "<<Xpos[jk]<<" "<<valz[jk]<<" "<<Xdev[jk]<<endl;
	      }
	      
	      StraightLineFit yposresolfit(1, zval, Ypos,  errysq, Yusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);

	      yposresolfit.GetParameters(nyfail, yinters, yslope);
	      cout<<"Slope and intercept Y  "<<yslope<<" "<<yinters<<endl;
	      yposresolfit.GetError(yerrcst, yerrlin, yerrcov);
	      cout<<"Error in slope and intercept Y "<<yerrcst<<" "<<yerrlin<<endl;
	      yposresolfit.GetChisqure(Ny, ychi2);
	      cout<<"NDOF and chi-square Y "<<Ny<<" "<<ychi2<<endl;
	      yposresolfit.GetFitValues(yext,valz, Ydev, yexter);
	      for(int jk = 0;jk<nlayer;jk++){
		cout<<"get yz fit values "<<yext[jk]<<" "<<Ypos[jk]<<" "<<valz[jk]<<" "<<Ydev[jk]<<endl;
	      }

	      if (ClustsInTrackBank[iji][jki].size()>=1){
		pAnalysis->XPosdev_exclu[jki] = Xdev[jki]; 
		pAnalysis->YPosdev_exclu[jki] =  Ydev[jki];
		pAnalysis->inPosX[jki] = Xpos[jki];
		pAnalysis->inPosY[jki] = Ypos[jki];
		pAnalysis->extPosX[jki] = xext[jki];
		pAnalysis->extPosY[jki] = yext[jki];
		pAnalysis->nXStrips[jki] = ClustsInTrackBank[iji][jki][0]->GetNXStripsInClust();
		pAnalysis->nYStrips[jki] = ClustsInTrackBank[iji][jki][0]->GetNYStripsInClust();
		
		pAnalysis->Xchisq[jki] = xchi2;
		pAnalysis->Ychisq[jki] = ychi2;
		pAnalysis->Xndof[jki] = Nx-2;
		pAnalysis->Yndof[jki] = Ny-2;
		pAnalysis->nXfail[jki] = nxfail;
		pAnalysis->nYfail[jki] = nyfail;
	      }
	  
	      cout<<"pAnalysis->XPosdev_exclu[jki] "<<pAnalysis->XPosdev_exclu[jki]<<" pAnalysis->YPosdev_xclu[jki] "<<pAnalysis->YPosdev_exclu[jki]<<endl;
	    }//  for(int jki=0;jki<nLayer;jki++){








	    
	  } // if(posresol){

	  
	  occulyr = -1;
	  cout<<"occulyr: "<<occulyr<<endl;
     
	    
	  //............RPC Straight line fit in X-Z plane.............
       
	  StraightLineFit xposfit(1, zval, Xpos,  errxsq, Xusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
	  // cout<<"check 20"<<endl;
	    
	    
	  xposfit.GetParameters(nxfail, xinters, xslope);
	  cout<<"Slope and intercept X  "<<xslope<<" "<<xinters<<endl;
	  xposfit.GetError(xerrcst, xerrlin, xerrcov);
	  cout<<"Error in slope and intercept X "<<xerrcst<<" "<<xerrlin<<endl;
	  xposfit.GetChisqure(Nx, xchi2);
	  cout<<"NDOF and chi-square X "<<Nx<<" "<<xchi2<<endl;
       
	  xposfit.GetFitValues(xext,valz, Xdev, xexter);
	  for(int jk = 0;jk<nlayer;jk++){
	    cout<<"get xz fit values "<<xext[jk]<<" "<<Xpos[jk]<<" "<<valz[jk]<<" "<<Xdev[jk]<<endl;
	  }


	
	  //	.............RPC Straight line fit in y-zplane....................

	    
	  StraightLineFit yposfit(1, zval, Ypos,  errysq, Yusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
	  yposfit.GetParameters(nyfail, yinters, yslope);
	  cout<<"Slope and intercept Y  "<<yslope<<" "<<yinters<<endl;
	  yposfit.GetError(yerrcst, yerrlin, yerrcov);
	  cout<<"Error in slope and intercept Y "<<yerrcst<<" "<<yerrlin<<endl;
	  yposfit.GetChisqure(Ny, ychi2);
	  cout<<"NDOF and chi-square Y "<<Ny<<" "<<ychi2<<endl;
	  yposfit.GetFitValues(yext,valz, Ydev, yexter);
	  for(int jk = 0;jk<nlayer;jk++){
	    cout<<"get yz fit values "<<yext[jk]<<" "<<Ypos[jk]<<" "<<valz[jk]<<" "<<Ydev[jk]<<endl;
	  }
       
	  // // dz by dt calculation:
       
	  // double szxy=0, sz=0, sxy=0, sn=0, sz2=0;
	  // double errsq=timeerr*timeerr;
	  // double dZdT = 0.0;

          // double ztinter = 0;
 
	  // for (int jk =0; jk<cluster_size;jk++) {

	  //   szxy +=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetZPos()*pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetTime()/errsq;
	  //   sz +=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetZPos()/errsq;
	  //   sz2 +=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetZPos()*pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetZPos()/errsq;
	  //   sxy +=pinotrack->InoTrack_list[iji]->ClustsInTrack[jk]->GetTime()/errsq;
	  //   sn +=1/errsq;
    
	  // }// for ( unsigned int jk =0; jk<cluster_size;jk++) {

   
	  // if ((sz2*sn - sz*sz) !=0) { 
	  //   dZdT = (szxy*sn - sz*sxy)/(sz2*sn - sz*sz);

	  //   ztinter  = sxy/sn - dZdT*sz/sn;

	 
	  //   cout<<"dZdT: "<<dZdT<<endl;
	  // }
	  // cout<<"dZdT: "<<dZdT<<endl;
  
	  // if(dZdT>0){
	  //   ZIncreasesWithTime= true;
	  // }

	  // else{
	  //   ZIncreasesWithTime=false;

	  // }
  
	  cout<<"Zincreaseswith time: "<<ZIncreasesWithTime<<endl;
	  //

	  // store events only when they satisfy:
       
	  // if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {
	  //   if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) {

	  if (Ny>=nmnhits/*-ntcor*/  && nyfail==0) {
	    if (Nx>=nmnhits/*-ntcor*/  && nxfail==0) {

	      
	      cout<<"Slope and intercept X  "<<xslope<<" "<<xinters<<endl;
	      cout<<"Error in slope and intercept X "<<xerrcst<<" "<<xerrlin<<endl;
	      cout<<"NDOF and chi-square X "<<Nx<<" "<<xchi2<<endl;
	      for(int jk = 0;jk<nlayer;jk++){
		cout<<"get x/y fit values "<<xext[jk]<<" "<<Xpos[jk]<<" "<<valz[jk]<<" "<<yext[jk]<<" "<<Ypos[jk]<<" "<<valz[jk]<<" "<<endl;
	      }


	      cout<<"Slope and intercept Y  "<<yslope<<" "<<yinters<<endl;
	      cout<<"Error in slope and intercept Y "<<yerrcst<<" "<<yerrlin<<endl;
	      cout<<"NDOF and chi-square Y "<<Ny<<" "<<ychi2<<endl;


	      //storing the data :


	      double theta, phi;
       
	      // theta and phi from slopes in x-z and y-z plane
       
	      phi = atan2(yslope,xslope);
     
	      //since z is decreasing with time we have minus sign
	      ZIncreasesWithTime=false;
	      if(ZIncreasesWithTime==false){
	  
		theta = acos(-1.0/sqrt(1+pow(xslope,2)+pow(yslope,2)));
	   
		double PI = acos(-1.0);
		phi +=PI;
		if (phi > PI) { phi -=2*PI;}
        
	      }

	      else{
		theta= acos(1.0/sqrt(1+pow(xslope,2)+pow(yslope,2)));
	      }
     
	      cout<<iji<<" theta: "<<theta<<" phi: "<<phi<<" xslope "<<xslope<<" yslope  "<<yslope<<endl;//:
	 

	      G4ThreeVector dirVector0(0,0,1);       
	      dirVector0.setTheta(theta);
	      dirVector0.setPhi(phi);
	      cout<<"dirVector0 "<<dirVector0<<endl;
	      
	      //:

	  
	  
	      fTrackCand = new InoTrackCand(pinotrack->InoTrack_list[iji], ZIncreasesWithTime); 
	  
	      fTrackCand->SetTheta(theta);
	      fTrackCand->SetPhi(phi);

     
	      double therr =sqrt(pow((xslope*xerrcst),2)+pow((yslope*yerrcst),2))/(1+pow(xslope,2)+pow(yslope,2))*(pow(xslope,2)+pow(yslope,2)) ;
	      double pherr = sqrt(pow(yslope*xerrcst,2)+pow(xslope*yerrcst,2))/(pow(xslope,2)+pow(yslope,2));
	      cout<<" therr: "<<therr<<" pherr: "<<pherr<<endl;
	   
	      fTrackCand->SetThErr(therr);
	      fTrackCand->SetPhErr(pherr);

	      xvtx_parameter[0] =0;//momen
	      xvtx_parameter[1] =theta;
	      xvtx_parameter[2] =phi;
	      // 11th layer extrapolation
	      // xvtx_parameter[3] =xslope*(-1.66062)+xinters;//in meters
	      // xvtx_parameter[4] =yslope*(-1.66062)+yinters;//in meters
	      // //
	      // xvtx_parameter[5]=-1.66062*1000;//top layer z pos in mm



	      // xvtx_parameter[3] =xslope*(-1.65512)+xinters;//in meters
	      // xvtx_parameter[4] =yslope*(-1.65512)+yinters;//in meters
	      // //
	      // xvtx_parameter[5]=-1.65512*1000;//top layer z pos in mm
	      
	      //	      zposmx = paradef->GetParroom(2); //roof coordinates
	      xvtx_parameter[3] =xslope*(zposmx)+xinters;//in meters
	      xvtx_parameter[4] =yslope*(zposmx)+yinters;//in meters
	      //
	      xvtx_parameter[5]=zposmx*1000;//top layer z pos in mm



	      
	      fTrackCand->SetExtPara(xvtx_parameter);
	  
	      for(int zz=0;zz<6;zz++){
		cout<<"xvtx_parameter "<<xvtx_parameter[zz]<<endl;	    
	      }

	      //   double xexterr = xerrcst + 2*xerrcov*(-1.66062) + xerrlin*(-1.66062)*(-1.66062);
	      //   double yexterr = yerrcst + 2*yerrcov*(-1.66062) + yerrlin*(-1.66062)*(-1.66062);


	      double xexterr = xerrcst + 2*xerrcov*(zposmx) + xerrlin*(zposmx)*(zposmx);
	      double yexterr = yerrcst + 2*yerrcov*(zposmx) + yerrlin*(zposmx)*(zposmx);


	      
     
	      cout<<"error in extrapolated x and y point: "<<xexterr<<" "<<yexterr<<endl;
     
	      fTrackCand->SetVtxUError(pow(xexterr,0.5));//err in xpos
	      fTrackCand->SetVtxVError(pow(yexterr,0.5));//err in y pos	   
	      fTrackCand->SetVtxdUError(pow(xerrlin,0.5));// error in  dxdz
	      fTrackCand->SetVtxdVError(pow(yerrlin,0.5));// error in dydz
	      fTrackCand->SetVtxQPError(0);
	      fTrackCand->SetVtxUdUError(xerrcov);//store cov
	      fTrackCand->SetVtxUdVError(0);
	      fTrackCand->SetVtxVdVError(yerrcov);//
	      fTrackCand->SetVtxVdUError(0);
	      fTrackCand->SetVtxdUdVError(0);
	      cout<<"Pointer exists "<<fTrackCand<<endl;
	      inoTrackCand_pointer->InoTrackCand_list.push_back(fTrackCand);

	      fTrackCand->SetChi2(xchi2);
	      fTrackCand->SetChi22(ychi2);
	      fTrackCand->SetNDOF(Nx-2);
	      fTrackCand->SetNDOF2(Ny-2);
	      
	      for (int jk =0; jk<nlayer;jk++) {
   
		pAnalysis->XPosdev[jk] = Xdev[jk];
		pAnalysis->YPosdev[jk] = Ydev[jk];
	      }
	      
 
	      
 	    }//    if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) {
	      
 	  }//  if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {

	  else{

	    fTrackCand=0;
	    cout<<"Pointer is Null "<<fTrackCand<<endl;
	  }
	  }//if(iji<pAnalysis->ntrkmx){       
	}//for (unsigned iji=0; iji<pinotrack->InoTrack_list.size() ; iji++) {

	pAnalysis->ntrack = itrk;
	//
 




     

	
      }//if(magfield)
      else 
	{

       
	  cout<<"Going to Track Fitter..."<<endl;
	  //  InoTrackFitAlg trackfitter;
	  cout<<"Running Track Fit Alg ..."<<endl;
	  trackfitter.RunAlg(); //VALGRIND
   
	  cout<<"Track Fit Alg Completed "<<endl;




   
	}



  
      vector<ixyz> xtrkstrip;
      vector<ixyz> ytrkstrip;
      vector<ixyz> xshwstrip;
      vector<ixyz> yshwstrip;
      
      int tottrkXstrp = 0;
      int tottrkYstrp = 0;
      
      InoTrackCand_Manager *pfitTrack =  InoTrackCand_Manager::APointer;
      cout<<"check 56  "<<pfitTrack<<endl;
       
      if (pfitTrack) {
	// cout<<"check 57"<<endl;
        //G4cout <<"tmphitlist in eventaction "<< pinohit->InoHit_list.size()<<G4endl;
	cout <<"tmptracklist in eventaction "<< nevent<<" "<<pfitTrack->InoTrackCand_list.size()<<" "<<pAnalysis->ihist<<G4endl;
        {
          unsigned ij=0;
          
        
          for (unsigned jk=0; jk<pfitTrack->InoTrackCand_list.size() ; jk++) {

            // cout <<"Identical "<< int(pfitTrack->InoTrackCand_list[jk]->ClustsInTrack[0]->isIdentical(trackfinder.inoCluster_pointer->InoCluster_list[0]));
            G4HCofThisEvent * hce = evt->GetHCofThisEvent();
            //	    micalcal0HitsCollection* ehco = 0;
            G4ThreeVector MCpos;
            if(hce) {
              // SSE : Error while running code sepataterly (other than 0 option)
              //  ehco = (micalcal0HitsCollection*)(hce->GetHC(cal0CollID));
              //MCpos = (*ehco)[jk]->GetPos();
            }
            
            if (ij <pAnalysis->ntrkmx) {
	      
	   
	      // cout<<"VtxPlane: "<<pfitTrack->InoTrackCand_list[jk]->GetVtxPlane()<<endl;
              pAnalysis->itype[ij] =  pfitTrack->InoTrackCand_list[jk]->GetFitType();
              pAnalysis->nhits[ij] = (pfitTrack->InoTrackCand_list[jk]->GetNDOF()+5)/2;
              pAnalysis->chisq[ij] =  pfitTrack->InoTrackCand_list[jk]->GetChi2();
              pAnalysis->cvalue[ij] = pfitTrack->InoTrackCand_list[jk]->Getcval();
              pAnalysis->chisq2[ij] =  pfitTrack->InoTrackCand_list[jk]->GetChi22(); // for straightline fit
	      pAnalysis->ndof[ij] =  pfitTrack->InoTrackCand_list[jk]->GetNDOF();
	      pAnalysis->ndof2[ij] =  pfitTrack->InoTrackCand_list[jk]->GetNDOF2();

	      
              pAnalysis->fc_or_pc[ij] = pfitTrack->InoTrackCand_list[jk]->GetFCPC();
              pAnalysis->trkmm[ij] = pfitTrack->InoTrackCand_list[jk]->GetMomentum();
              pAnalysis->trkth[ij] = pfitTrack->InoTrackCand_list[jk]->GetTheta();
              pAnalysis->trkph[ij] = pfitTrack->InoTrackCand_list[jk]->GetPhi();
              // cout<<"-------------------------------------------------------------------------"<<endl;
              // cout<<"Reconstructed P = "<<pAnalysis->trkmm[ij]<<"  |  theta "<<pAnalysis->trkth[ij]*180/3.1415<<"  |  phi "<<pAnalysis->trkph[ij]*180/3.1415<<endl;
              // cout<<"-------------------------------------------------------------------------"<<endl;
              pAnalysis->therr[ij] = pfitTrack->InoTrackCand_list[jk]->GetThErr();
              pAnalysis->pherr[ij] = pfitTrack->InoTrackCand_list[jk]->GetPhErr();
	    
	      cout<<"atimeslope "<<dZdT<<endl;
	      pAnalysis->atimslope[ij] =dZdT;
	      pAnalysis->atiminter[ij] = ztinter;
	 
	      //Changed for miniICAL
	      //////////////////////////////////////////////////////////////////////////
	      // pAnalysis->momvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetMomentumCurve();
              // pAnalysis->thevx[ij] = acos(pfitTrack->InoTrackCand_list[jk]->GetDirCosZ());
              // pAnalysis->phivx[ij] = atan2(pfitTrack->InoTrackCand_list[jk]->GetDirCosV(),
              // 	   pfitTrack->InoTrackCand_list[jk]->GetDirCosU());
              // pAnalysis->posxvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxU();
              // pAnalysis->posyvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxV();
              // pAnalysis->poszvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxZ();
              
              pAnalysis->momvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetExtPara(0);
              pAnalysis->thevx[ij] = pfitTrack->InoTrackCand_list[jk]->GetExtPara(1);
              pAnalysis->phivx[ij] = pfitTrack->InoTrackCand_list[jk]->GetExtPara(2);
              
              pAnalysis->posxvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetExtPara(3);
              pAnalysis->posyvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetExtPara(4);
              pAnalysis->poszvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetExtPara(5);
              
	      cout<<"..........."<<      pfitTrack->InoTrackCand_list[jk]->GetExtPara(3)      <<"  "<<  pfitTrack->InoTrackCand_list[jk]->GetExtPara(4)<< "  "<<  pfitTrack->InoTrackCand_list[jk]->GetExtPara(5)<<endl;

	      pAnalysis->momrf[ij] = pfitTrack->InoTrackCand_list[jk]->GetRoofPara(0);
	      pAnalysis->therf[ij] = pfitTrack->InoTrackCand_list[jk]->GetRoofPara(1);
	      pAnalysis->phirf[ij] = pfitTrack->InoTrackCand_list[jk]->GetRoofPara(2);
		 
	      pAnalysis->posxrf[ij] = pfitTrack->InoTrackCand_list[jk]->GetRoofPara(3);
	      pAnalysis->posyrf[ij] = pfitTrack->InoTrackCand_list[jk]->GetRoofPara(4);
	      pAnalysis->poszrf[ij] = pfitTrack->InoTrackCand_list[jk]->GetRoofPara(5);
              
	      cout<<"....ROOF......."<<      pfitTrack->InoTrackCand_list[jk]->GetRoofPara(3)      <<"  "<<  pfitTrack->InoTrackCand_list[jk]->GetRoofPara(4)<< "  "<<  pfitTrack->InoTrackCand_list[jk]->GetRoofPara(5)<<endl;

             
	      
              pAnalysis->momend[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndMomentumCurve();
              pAnalysis->theend[ij] = acos(pfitTrack->InoTrackCand_list[jk]->GetEndDirCosZ());
              pAnalysis->phiend[ij] = atan2(pfitTrack->InoTrackCand_list[jk]->GetEndDirCosV(),
                                            pfitTrack->InoTrackCand_list[jk]->GetEndDirCosU());
              pAnalysis->tx_end[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndDirCosU();
              pAnalysis->ty_end[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndDirCosV();
              
              pAnalysis->posxend[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndU();
              pAnalysis->posyend[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndV();
              pAnalysis->poszend[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndZ();

	      
              /*
                pAnalysis->momds[ij] = pfitTrack->InoTrackCand_list[jk]->GetMomentumdS()
                + pfitTrack->InoTrackCand_list[ij]->GetdSExtra(); //GMA -ve value for lower momentum !!!!!!
                pAnalysis->momrg[jk] = pfitTrack->InoTrackCand_list[ij]->GetMomentumRange()
                +  pfitTrack->InoTrackCand_list[ij]->GetRangeExtra();
              */

	      
              pAnalysis->momds[ij] = pfitTrack->InoTrackCand_list[jk]->GetMomentumdS();  //GMA14 was bug
	      cout<<"Check 41"<<endl;
              pAnalysis->momrg[ij] = pfitTrack->InoTrackCand_list[jk]->GetMomentumRange();
	      cout<<"Check 45"<<endl;
              // cout<<"micalEventAction ij "<<ij<<endl;
              // cout<<"pAnalysis->trkmm "<<pAnalysis->trkmm[ij]<<endl;
              // cout<<"pAnalysis->momds "<<pAnalysis->momds[ij]<<endl;
              // cout<<"pAnalysis->momrg "<<pAnalysis->momrg[ij]<<endl;

	      cout<<"VtxPlane: "<<pfitTrack->InoTrackCand_list[jk]->GetVtxPlane()<<endl;
	      pAnalysis->vtxzplane[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxPlane();
	     
              pAnalysis->endzplane[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndPlane();

	     
              
              pAnalysis->xxerr[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxUError(); //GMA14 all these 15 tesrms


              pAnalysis->yyerr[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxVError();
              pAnalysis->txerr[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxdUError();
              pAnalysis->tyerr[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxdVError();
              pAnalysis->qperr[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxQPError();

	      pAnalysis->xxtxerr[ij]=pfitTrack->InoTrackCand_list[jk]->GetVtxUdUError();
	      pAnalysis->xxtyerr[ij]=pfitTrack->InoTrackCand_list[jk]->GetVtxUdVError();
	      pAnalysis->yytyerr[ij]=pfitTrack->InoTrackCand_list[jk]->GetVtxVdVError();
	      pAnalysis->yytxerr[ij]=pfitTrack->InoTrackCand_list[jk]->GetVtxVdUError();
	      pAnalysis->txtyerr[ij]=pfitTrack->InoTrackCand_list[jk]->GetVtxdUdVError();
	      
              pAnalysis->xxenderr[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndUError();
              pAnalysis->yyenderr[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndVError();
              pAnalysis->txenderr[ij] = pfitTrack->InoTrackCand_list[jk]->GetEnddUError();
              pAnalysis->tyenderr[ij] = pfitTrack->InoTrackCand_list[jk]->GetEnddVError();
              pAnalysis->qpenderr[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndQPError();
              
              pAnalysis->xxin[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxXX();
              pAnalysis->yyin[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxYY();
              pAnalysis->txin[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxTX();
              pAnalysis->tyin[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxTY();
              pAnalysis->tx[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxdU();
              pAnalysis->ty[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxdV();
              
              pAnalysis->mcxgnvx[ij]=1.e-3*MCpos.x();
              pAnalysis->mcygnvx[ij]=1.e-3*MCpos.y();
              //Initialise with this mainly for noise track
              
              // cout<<"---------------------------------------------------------------"<<endl;
              // cout<<"p = "<<pAnalysis->trkmm[ij]<<", momds = "<<pAnalysis->momds[ij]<<", E_mu = "<<pAnalysis->momin[0]<<endl;
              // cout<<"---------------------------------------------------------------"<<endl;
              
              pAnalysis->momgnvx[ij]  = pAnalysis->momgnend[ij] = 0.0;
              pAnalysis->thegnvx[ij] = pAnalysis->thegnend[ij] = 10.;
              pAnalysis->phignvx[ij] = pAnalysis->phignend[ij] = 10.;
              

 
              vector<InoCluster*> tmpclusts = pfitTrack->InoTrackCand_list[jk]->ClustsInTrack; // fTrack->ClustsInTrack;
              
              int plane = pfitTrack->InoTrackCand_list[jk]->GetVtxPlane();




	      
              // cout <<"genplane "<<plane<<endl;
              for (unsigned kl = 0; kl <tmpclusts.size(); kl++)	{
                if (tmpclusts[kl]->GetZPlane()==plane) {
                  pAnalysis->momgnvx[ij] = 0.001*tmpclusts[kl]->HitsInCluster[0]->GetMomentum();
		  pAnalysis->thegnvx[ij] = tmpclusts[kl]->HitsInCluster[0]->GetTheta();
		  pAnalysis->phignvx[ij] = tmpclusts[kl]->HitsInCluster[0]->GetPhi();
		  break;
		} // if (tmpclusts[kl]->GetZPlane()==plane)
	      } //  for (unsigned kl = 0; kl <tmpclusts.size(); kl++)
	      int trkclustersize = tmpclusts.size();
	      int tottrkhit=0;
	      int trkxstrp = 0;
	      int trkystrp = 0;
	      for (int kl=0; kl<trkclustersize; kl++) {
		vector<InoHit*>  tmphit = tmpclusts[kl]->HitsInCluster;
		int nhits = tmphit.size();
		int zplane =tmphit[0]->GetZPlane();
					
		for (int lm=0; lm<nhits; lm++) {
		  int ixstripno = tmphit[lm]->GetXStripNum();
		  int iystripno = tmphit[lm]->GetYStripNum();
		  if (tmphit[lm]->GetXStrip()) {
		    for (unsigned nn=0; nn< xtrkstrip.size(); nn++) {
		      //if (istripno == xtrkstrip[nn]) {istripno = -1; break;}
		      if (ixstripno == xtrkstrip[nn].second &&
			  zplane == xtrkstrip[nn].first) {ixstripno = -1; break;}
		    }
		  }
		  if (tmphit[lm]->GetYStrip()) {
		    for (unsigned nn=0; nn< ytrkstrip.size(); nn++) {
		      //if (istripno == ytrkstrip[nn]) {istripno = -1; break;}
		      if (iystripno == xtrkstrip[nn].second &&
			  zplane == xtrkstrip[nn].first) {iystripno = -1; break;}
		    }
		  }
		  if (ixstripno>=0 && iystripno>=0) {
		    trkxstrp++;  trkystrp++;
		    if (pAnalysis->cvalue[ij]>0){ixyz Zxstrip(zplane,ixstripno) ; xtrkstrip.push_back(Zxstrip);}
		    if (pAnalysis->cvalue[ij]>0){ixyz Zystrip(zplane,iystripno) ; ytrkstrip.push_back(Zystrip);}
		  }
		} // for (int lm=0; lm<nhits; lm++)
		//tmphit.clear();
		tottrkhit +=nhits;
	      } // for (int kl=0; kl<trkclustersize; kl++)
				
	      pAnalysis->ntrkcl[ij] = 1000*trkclustersize + tottrkhit;
	      pAnalysis->ntrkst[ij] = 1000*trkxstrp + trkystrp;
	      tottrkXstrp +=trkxstrp;
	      tottrkYstrp +=trkystrp;
	    
	      plane = pfitTrack->InoTrackCand_list[jk]->GetEndPlane();
	    
	      for (unsigned kl = 0; kl <tmpclusts.size(); kl++) {
		if (tmpclusts[kl]->GetZPlane()==plane) {
		  pAnalysis->momgnend[ij] = 0.001*tmpclusts[kl]->HitsInCluster[0]->GetMomentum();
		  pAnalysis->thegnend[ij] = tmpclusts[kl]->HitsInCluster[0]->GetTheta();
		  pAnalysis->phignend[ij] = tmpclusts[kl]->HitsInCluster[0]->GetPhi();
		  break;
		}
	      }
	      //tmpclusts.clear();
	    } // if (ij <pAnalysis->ntrkmx)
	    ij++;
	  } //for (unsigned jk=0; jk<pfitTrack->InoTrackCand_list.size() ; jk++)
      
	  cout <<"ntrack "<< ij<<endl;
	  pAnalysis->ntrkt = ij; //(pfitTrack->InoTrackCand_list.size() <=pAnalysis->ntrkmx) ? pfitTrack->InoTrackCand_list.size() : pAnalysis->ntrkmx;
	
	  //Now add other cluster close to vertex;
	  //	  vector <InoCluster*> totcluster = trackfinder.inoCluster_pointer->InoCluster_list;
	  //	  int totclustersize = totcluster.size();
	  //cout<<"totclustersize "<<totclustersize<<endl;
	  //int trackhits=0; //asm  temp variable
	  //Tag all clusters which belongs to any leading track

	  for (int jk=0; jk<totclustersize; jk++) {
	  
	    totcluster[jk]->SetInTrack(0);
	    totcluster[jk]->SetInShower(0);
	  
	    for (unsigned kl=0; kl<pfitTrack->InoTrackCand_list.size() ; kl++) {
	      if(((pfitTrack->InoTrackCand_list[0]->GetChi2()/pfitTrack->InoTrackCand_list[0]->GetNDOF())<30 ||
		  (abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())>0. &&
		   abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())<1000 ))) { // ||
		//		(pfitTrack->InoTrackCand_list[0]->Getcval()>0))) {

		for (unsigned lm=0; lm< pfitTrack->InoTrackCand_list[kl]->GetClusterEntries(); lm++) {
		  if (totcluster[jk]->isIdentical(pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm])) {
		  
		    if((kl>0&&abs(pfitTrack->InoTrackCand_list[kl]->GetVtxZ()-pfitTrack->InoTrackCand_list[0]->GetVtxZ())<.1)) {
		      totcluster[jk]->SetInShower(1);
		    
		      vector<InoHit*> tmphit = totcluster[jk]->HitsInCluster;
		      int nhits = tmphit.size();
		      int zplane =tmphit[0]->GetZPlane();
		    
		      for (int mn=0; mn<nhits; mn++) {
			int ixstripno = tmphit[mn]->GetXStripNum();
			int iystripno = tmphit[mn]->GetYStripNum();
			//		      if ((ixstripno==0 && iystripno==0) ||( ixstripno>=0 && iystripno>=0) ) {
			if (ixstripno>=0 && iystripno>=0) {
			  ixyz ZXstrp(zplane,ixstripno); ixyz ZYstrp(zplane,iystripno);
			  xshwstrip.push_back(ZXstrp); yshwstrip.push_back(ZYstrp);
			}
			if (tmphit[mn]->GetXStrip() && tmphit[mn]->GetYStrip() ) {
			  //tmphit[mn]->SetUID(-442);// SSE 08/15
			  if (pAnalysis->isVisOut) {
			    pAnalysis->H->NShowerHits++;
			    pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object  //VALGRIND
			    pAnalysis->Hp->ShowerHitNum++;
			    pAnalysis->Hp->ZZ=  zplane;//pfitTrack->InoTrackCand_list[i]->ClustsInTrack[j]->GetZPlane(); //os();
			    pAnalysis->Hp->XX=tmphit[mn]->GetXPos();
			    pAnalysis->Hp->YY=tmphit[mn]->GetYPos();
			    pAnalysis->Hp->TrackType=-442; // track type up -44 + 0 FIT UP
			    pAnalysis->Hp->Fup=2;// Shower hits from an event
			  } //isVis
			}
		      }
		    }  else {
		      totcluster[jk]->SetInTrack(1);
		    }
		  }
		}
	      }
	    }
	  } // for (int jk=0; jk<totclustersize; jk++)
	
	  {int nc_cl=0; int nc_ht=0; int ccc_cl=0; int ccc_ht=0 ; int t_ht=0;
	    for (int jk=0; jk<totclustersize; ++jk) {
	      if(totcluster[jk]->GetInTrack()==1&&totcluster[jk]->GetInShower()!=1) {
		nc_cl++; nc_ht += totcluster[jk]->GetHitEntries();
	      } else {
		ccc_cl++; ccc_ht+=totcluster[jk]->GetHitEntries();
	      }
	      t_ht+=totcluster[jk]->GetHitEntries();
	    }
	    //	    cout <<"before implemeting : " <<"hits not in track:"<< ccc_ht
	    //	 << " + track hits: "<<  nc_ht <<" total " << t_ht <<endl;//SSE 29 Oct 2015
	    if(pAnalysis->isVisOut==1){
	      pAnalysis->ascii_output <<endl;
	      pAnalysis->ascii_output <<"before implemeting : " <<"hits not in track:"<< ccc_ht << " + track hits: "<<  nc_ht <<" total " << t_ht <<endl;
	    }
	  }

	  //	  cout<< "totclustersize  " <<totclustersize<<  endl;
	  //====================================================================
	  //  calculation of hadron energy //SSE mofied and added here 291015
	  //Calculation of hadron momentum and direction  from simulation input Oct8, 15
	  //====================================================================
	
	  double E_had;
	  float Px_nu=0.0, Py_nu=0.0, Pz_nu=0.0; //SSE 081015
	  float Px_muon=0.0, Py_muon=0.0, Pz_muon=0.0; //SSE 081015
	  float  thetain_had=0.0, phiin_had=0.0;//,pmagin_had=0.0; //SSE 081015
	  //float  costhetain_had=0.0; //SSE 081015
	  //float costheta_had=0.0;
	  float dotAngle=0.0;
	  G4ThreeVector tmp3pin_nu(0.,0.,0.);
	  G4ThreeVector tmp3pin_muon(0.,0.,0.);
	  G4ThreeVector tmp3pin_had(0.,0.,0.);
	  if (abs(pAnalysis->pidin[1])==13||abs(pAnalysis->pidin[1])==14){
	    Px_muon=(pAnalysis->momin[1])*(TMath::Sin(pAnalysis->thein[1]))*(TMath::Cos(pAnalysis->phiin[1]));
	    Py_muon=(pAnalysis->momin[1])*(TMath::Sin(pAnalysis->thein[1]))*(TMath::Sin(pAnalysis->phiin[1]));
	    Pz_muon=(pAnalysis->momin[1])*(TMath::Cos(pAnalysis->thein[1]));
	    E_had=abs(pAnalysis->momin[0])- abs(pAnalysis->momin[1]);
	  } else {
	    E_had=pAnalysis->momin[0];
	  }

	  //=========================cal. of  hadron momentum and direction=================================
	
	  Px_nu=(pAnalysis->momin[0])*(TMath::Sin(pAnalysis->thein[0]))*(TMath::Cos(pAnalysis->phiin[0]));
	  Py_nu=(pAnalysis->momin[0])*(TMath::Sin(pAnalysis->thein[0]))*(TMath::Sin(pAnalysis->phiin[0]));
	  Pz_nu=(pAnalysis->momin[0])*(TMath::Cos(pAnalysis->thein[0]));
	
	  //G4ThreeVector 
	  tmp3pin_nu=G4ThreeVector (Px_nu,Py_nu, Pz_nu);
	  //G4ThreeVector 
	  tmp3pin_muon= G4ThreeVector (Px_muon,Py_muon, Pz_muon);
	  //G4ThreeVector 
	  tmp3pin_had=tmp3pin_nu-tmp3pin_muon;
	
	  //pmagin_had=tmp3pin_had.mag();
	  thetain_had=tmp3pin_had.theta();
	  //costhetain_had=tmp3pin_had.cosTheta();
	  phiin_had=tmp3pin_had.phi();
	
	  pAnalysis->theta_hadron_in=0;//SSE 081015
	  //pAnalysis->costheta_hadron_in=0;//SSE 081015
	  pAnalysis->phi_hadron_in=0;//SSE 081015
	  pAnalysis->theta_hadron_in=thetain_had ;//SSE 081015
	  //pAnalysis->costheta_hadron_in=costhetain_had ;//SSE 081015
	  pAnalysis->phi_hadron_in=phiin_had ;//SSE 081015
	
	  // cout<< " nu p in: x, y, z "<< tmp3pin_nu[0]<< " "<<tmp3pin_nu[1]<<" "<< tmp3pin_nu[2]<<endl;
	  // cout<< " mu p in: x, y, z "<< tmp3pin_muon[0]<< " "<<tmp3pin_muon[1]<<" "<< tmp3pin_muon[2]<<endl;
	  // cout<< " had p in: x, y, z "<< tmp3pin_had[0]<< " "<<tmp3pin_had[1]<<" "<< tmp3pin_had[2]<<endl;
	  // cout<< " had theta phi  input: x, y, z "<< thetain_had<< " "<<phiin_had<<endl;
	  //============= end of cal. of  hadron momentum and direction ==================================
	
	  // cout<<"E_had"<<" "<<E_had<<" "<<pAnalysis->momin[0]<<" "<<pAnalysis->momin[1]<<endl;//SSE 291015
	
	  //======================================================================
	  //Add flag through SetUID (Aug 2015; SSE) added here 291015
	  //
	  //
	  //=======================================================================
	  // cout<<"pinotrack->InoTrack_list.size() = "<<pinotrack->InoTrack_list.size()<<endl;
	  for (unsigned ixj=0; ixj<pinotrack->InoTrack_list.size() ; ixj++) {
	    // cout<<"Inside Loop 1, pinotrack->InoTrack_list[ixj]->ClustsInTrack.size() = "<<pinotrack->InoTrack_list[ixj]->ClustsInTrack.size()<<endl;
	    for ( unsigned int jxk =0; jxk<pinotrack->InoTrack_list[ixj]->ClustsInTrack.size();jxk++) {	
	      // cout<<"Inside Loop 2, pinotrack->InoTrack_list[ixj]->ClustsInTrack[jxk]->GetHitEntries() = "<<pinotrack->InoTrack_list[ixj]->ClustsInTrack[jxk]->GetHitEntries()<<endl;
	      for(unsigned int mxn=0;mxn<pinotrack->InoTrack_list[ixj]->ClustsInTrack[jxk]->GetHitEntries();mxn++){
		// cout<<"Inside Loop 3"<<endl;
		pinotrack->InoTrack_list[ixj]->ClustsInTrack[jxk]->HitsInCluster[mxn]->SetUID(-4);//SSE 08/15
	      }
	    }
	  }
	
	  // cout<<"pfitTrack->InoTrackCand_list.size() = "<<pfitTrack->InoTrackCand_list.size()<<endl;
	  for (unsigned kxl=0; kxl< pfitTrack->InoTrackCand_list.size() ; kxl++) {
	    // cout<<"Inside Loop 1, pfitTrack->InoTrackCand_list[kxl]->GetClusterEntries() = "<<pfitTrack->InoTrackCand_list[kxl]->GetClusterEntries()<<endl;
	    for ( unsigned int lxm =0; lxm<pfitTrack->InoTrackCand_list[kxl]->GetClusterEntries(); lxm++){
	      // cout<<"Inside Loop 2, pfitTrack->InoTrackCand_list[kxl]->ClustsInTrack[lxm]->GetHitEntries() = "<<pfitTrack->InoTrackCand_list[kxl]->ClustsInTrack[lxm]->GetHitEntries()<<endl;
	      for(unsigned int mxn=0; mxn<pfitTrack->InoTrackCand_list[kxl]->ClustsInTrack[lxm]->GetHitEntries();mxn++){
		// cout<<"Inside Loop 3"<<endl;
		if(pfitTrack->InoTrackCand_list[kxl]->GetFitType()==1){
		  // cout<<"UID == 10"<<endl;
		  pfitTrack->InoTrackCand_list[kxl]->ClustsInTrack[lxm]->HitsInCluster[mxn]->SetUID(10);	//SSE 08/15
		} else {
		  // cout<<"UID == 11"<<endl;
		  pfitTrack->InoTrackCand_list[kxl]->ClustsInTrack[lxm]->HitsInCluster[mxn]->SetUID(11);
		}
	      }
	    }
	  }
	  // cout<<"Outside Loop"<<endl;
	  //=========================================================


	  int Orighits_cluster=0;//SSE Hits in cluster + strip finding algo
	  vector<InoHit*> tmphit_trape;////SSE 301015  
	  int Hit_wo_ghst=0;
	    
	  int Orighits_trape=0;//hits in trapezoid + strip finding algo 
	  int wogh_Orighits=0;
	  int num_hadron=0;//SSE 250915
	  
	  float nhits_large_clust_selected=0;// SSE 091015
	  
	  double matrix_XX=0.0, matrix_YY=0.0, matrix_ZZ=0.0;//SSE 031115
	  double Matrix_xx=0.0, Matrix_yy=0.0, Matrix_zz=0.0;//SSE 031115
	  double Matrix_xy=0.0, Matrix_yz=0.0, Matrix_zx=0.0;//SSE 031115
	  
	  
	  int minplane_cluster=0;//SSE 091015
	  
	  
    
	  int allremhit = 0;
	  int allremcls = 0;
	  double x0=0,y0=0, z0=0;
	  if(pfitTrack->InoTrackCand_list.size() >0) {
	    x0 = pfitTrack->InoTrackCand_list[0]->GetVtxU();  // asm:  vertex of longest track
	    y0 = pfitTrack->InoTrackCand_list[0]->GetVtxV();   //asm
	    z0 = pfitTrack->InoTrackCand_list[0]->GetVtxZ();  
	    // cout<<" chi2 "<< (pfitTrack->InoTrackCand_list[0]->GetChi2()/pfitTrack->InoTrackCand_list[0]->GetNDOF())<<endl;
	    // cout<<	" abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())  "<< abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())<<endl;	
	    // cout<<	" pfitTrack->InoTrackCand_list[0]->Getcval()  "<< pfitTrack->InoTrackCand_list[0]->Getcval()<<endl;	
	  }
	  // cout<<" x0 "<< x0<< " y0 "<<y0<<" "<<pfitTrack->InoTrackCand_list.size()<<endl;//SSE 281015
	  float mXX=0; float  mYY=0; float mZZ=0;
	  float Mxx =0; float Myy=0; float Mzz =0;
	  float Mxy =0; float Mxz=0; float Myz =0;
	  //	commented by SSE because cval= 0 or -0

	  //   &&(pfitTrack->InoTrackCand_list[0]->Getcval()>0)){
	  if( (x0!=0 || y0!=0 ) && 
	      (pfitTrack->InoTrackCand_list.size()>0) &&
	      (pfitTrack->InoTrackCand_list[0]->GetChi2()/pfitTrack->InoTrackCand_list[0]->GetNDOF())<30 &&  //add a function for nhits
	      abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())>0. && 
	      abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())<1000000) {
	    
	    // cout<<"in loop:  x0 "<< x0<< " y0 "<<y0<<endl;//SSE 281015
	    for (unsigned kl=0; kl<1 ; kl++) {
	      int plane = pfitTrack->InoTrackCand_list[kl]->fVertex->GetPlane();
	      //TVector3 direction = pfitTrack->InoTrackCand_list[kl]->fVertex->GetDirCosine();
	    
	      double theta = pfitTrack->InoTrackCand_list[0]->GetTheta(); if (abs(theta<1.e-4)) {theta=0;};
	      double phi = pfitTrack->InoTrackCand_list[0]->GetPhi();   if (abs(phi<1.e-4)) {phi=0;};
	      double dxdz = tan(theta)*cos(phi);     if( abs(dxdz)<1.e-4) {dxdz=0;}
	      double dydz = tan(theta)*sin(phi);     if( abs(dydz)<1.e-4) {dydz=0;}
	      // cout<< " theta "<< theta<<" phi  "<< phi<< " dxdz "<< dxdz<<" dydz "<<dydz<<endl;//SSE 28Oct 
	    
	      int dir = (cos(theta)>0) ? 1 : -1;  //asm: check this
	      int plusplane=5;
	      int minusplane=-1;
	      double basewindow1 = 0.05, basewindow=0.05;// 10 cm
	    
	      double slopex = 1.0  ; //dz/dx(dy) = 1
	      double slopey = 1.0  ;
	    
	      if (abs(dxdz)>1) {slopex = 1.0;} else {slopex = 2.0;}
	      if (abs(dydz)>1) {slopey = 1.0;} else {slopey = 2.0;}
	      if(pAnalysis->isXtermOut==2) {
		pAnalysis->ascii_output <<"-------------------------------------------------------------"<<endl;
		pAnalysis->ascii_output <<" x0:y0 "<< x0 <<":"<<y0<< "  dxdz:dydz "<<  dxdz <<":" <<dydz<<endl;
	      }
	      int collected_hits=0;
	      do {
		collected_hits=0;//Trk->SetStraight();
		//inoTrack_pointer->InoTrack_list.push_back(Trk);
		basewindow1 +=0.05;
		if (basewindow1<0.2){
		  plusplane=5; minusplane=-2;
		} else if (basewindow1>=0.2){
		  plusplane+=1; minusplane = -3;
		}
	      
		for (int jxk=minusplane; jxk<plusplane; jxk++) {
		  int newplane = plane + jxk*dir;
		  if (newplane<0 && newplane>=nLayer) continue;
		
		  if (jxk<0) {
		    basewindow =basewindow1+0.05;
		  } else if(jxk>=0) {
		    basewindow=basewindow1;
		  }
		
		  for (int lm=0; lm<totclustersize; lm++) {
		    // cout<<"Loop lm at 739 = "<<lm<<endl;
		    if (totcluster[lm]->GetZPlane() !=newplane) continue;
		    if (totcluster[lm]->GetInTrack()) continue; //remove clusters in tracks
		  
		    double xmn = x0 +  LayerThickness*(jxk+1)*dir*(dxdz - dir*slopex) - basewindow;
		    double xmx = x0 +  LayerThickness*(jxk+1)*dir*(dxdz + dir*slopex) + basewindow;
		  
		    double ymn = y0 + LayerThickness*(jxk+1)*dir*(dydz - dir*slopey) - basewindow;
		    double ymx = y0 + LayerThickness*(jxk+1)*dir*(dydz + dir*slopey) + basewindow;
		  
		    if (totcluster[lm]->GetBegXPos() > xmn &&
			totcluster[lm]->GetEndXPos() < xmx &&
			totcluster[lm]->GetBegYPos() > ymn &&
			totcluster[lm]->GetEndYPos() < ymx &&
			!totcluster[lm]->GetInShower() && 
			!totcluster[lm]->GetInTrack()) {
		    
		      allremhit +=totcluster[lm]->GetHitEntries();
		      allremcls +=1;   totcluster[lm]->SetInShower(1);
		    
		      vector<InoHit*> tmphit = totcluster[lm]->HitsInCluster;
		      int nhits = tmphit.size();
		      int zplane =tmphit[0]->GetZPlane();
		    
		    
		      for (int mn=0; mn<nhits; mn++) {
			// cout<<" UID " <<tmphit[mn]->GetUID()<<endl;
			int ixstripno = tmphit[mn]->GetXStripNum();
			int iystripno = tmphit[mn]->GetYStripNum();
		      
			if (tmphit[mn]->GetXStrip() && tmphit[mn]->GetYStrip() ) {
			  if (tmphit[mn]->GetXPos() < xmn || tmphit[mn]->GetXPos() > xmx) continue;
			  if (tmphit[mn]->GetYPos() < ymn || tmphit[mn]->GetYPos() > ymx) continue;

			  for (unsigned nn=0; nn< xtrkstrip.size(); nn++) {
			    if (ixstripno == xtrkstrip[nn].second &&
				zplane == xtrkstrip[nn].first) {ixstripno = -1; break;}
			  }
			
			  for (unsigned nn=0; nn< ytrkstrip.size(); nn++) {
			    if (iystripno == ytrkstrip[nn].second &&
				zplane == ytrkstrip[nn].first) {
			      iystripno = -1; break;
			    }
			  }
			  for (unsigned nn=0; nn< xshwstrip.size(); nn++) {
			    if (ixstripno == xshwstrip[nn].second &&
				zplane == xshwstrip[nn].first) {
			      ixstripno = -1; break;
			    }
			  }
			 
			  for (unsigned nn=0; nn< yshwstrip.size(); nn++) {
			    if (iystripno == yshwstrip[nn].second &&
				zplane == yshwstrip[nn].first) {
			      iystripno = -1; break;
			    }
			  }
			
			  if((ixstripno>=0 || iystripno>=0) ) {
			    tmphit_trape.push_back(tmphit[mn]);//SSE 291015

			    if(ixstripno>=0) {
			      ixyz ZXstrp(zplane,ixstripno);xshwstrip.push_back(ZXstrp);
			    }
			    if(iystripno>=0){
			      ixyz ZYstrp(zplane,iystripno);yshwstrip.push_back(ZYstrp);
			    }
			    collected_hits++;
			    
			  }
		       
			  //tmphit[mn]->SetUID(-442); //SSE 08/15	
			  if (pAnalysis->isVisOut) {
			    pAnalysis->H->NShowerHits++;
			    pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object  //VALGRIND
			    pAnalysis->Hp->ShowerHitNum++;
			    pAnalysis->Hp->ZZ=  zplane;//pfitTrack->InoTrackCand_list[i]->ClustsInTrack[j]->GetZPlane(); //os();
			    pAnalysis->Hp->XX=tmphit[mn]->GetXPos();
			    pAnalysis->Hp->YY=tmphit[mn]->GetYPos();
			    pAnalysis->Hp->TrackType=-442; // track type up -44 + 0 FIT UP
			    pAnalysis->Hp->Fup=2;// Shower hits from an event
			    mXX=tmphit[mn]->GetXPos();
			    mYY=tmphit[mn]->GetYPos();
			    mZZ=tmphit[mn]->GetZPos();
			    Mxx += mXX*mXX;
			    Mxy += mXX*mYY;
			    Mxz += mXX*mZZ;
			    Myy += mYY*mYY;
			    Myz += mYY*mZZ;
			    Mzz += mZZ*mZZ;
			  
			  } //isVis
			} //if (tmphit[mn]->GetXStrip() && tmphit[mn]->GetYStrip() )
		      } //for (int mn=0; mn<nhits; mn++)
		      //	      tmphit.clear();
		    } //if (totcluster[lm]->GetBegXPos() > xmn .....)
		    // cout<< " trapezoidal cluster size Oct 27, 2015  size x & y"<< xshwstrip.size ()<< " "<<yshwstrip.size () <<endl;
		  } //for (int lm=0; lm<totclustersize; lm++)//asm: Dec7 : need to cleat x/ytrkstrip  for the next plane
		} //for (int jxk=minusplane; jxk<plusplane; jxk++)
		
	      } while(collected_hits>4||basewindow<=0.35);
	    
	    } //for (unsigned kl=0; kl<1 ; kl++) c: //sse 
	    // WE need to call orighit_calc function here
	    
	    Orighits_trape = orighit_calc(tmphit_trape);//SSE 30102015
	    
	    //  histogram filling for Orighits
	    // pAnalysis->hist_orighits_mod_E->Fill(Orighits_mod, E_had);
	    
	    
	    
	   
	    //=================================================
	    //Added on 031115 
	    //to Find hadron shower direction
	    //=================================================
	    
	    for( unsigned int ijk=0; ijk<tmphit_trape.size(); ijk++) {
	      //need min two planes including muon vetrex plane
	      if (ijk>0 && pfitTrack->InoTrackCand_list[0]->fVertex->GetPlane()!=tmphit_trape[ijk]->GetZPlane()){
		minplane_cluster++;
		break;
	      }
	    } 
	    
	    if (minplane_cluster>0){// lagest cluster must contains min two planes.
	      for( unsigned int mn=0; mn<tmphit_trape.size(); mn++) {
		
		matrix_XX=(tmphit_trape[mn]->GetXPos() - x0);
		matrix_YY=(tmphit_trape[mn]->GetYPos() - y0);
		matrix_ZZ=(tmphit_trape[mn]->GetZPos() - z0);

		// cout<< " xx  "<< tmphit_trape[mn]->GetXPos()
		// << " yy  "<< tmphit_trape[mn]->GetYPos()
		// << " zz  "<< " "<<tmphit_trape[mn]->GetZPos()<<endl<<endl;
		
		// cout<< " matrix_XX "<< matrix_XX<< " matrix_YY "<< matrix_YY<< "  matrix_ZZ "<<  matrix_ZZ<<endl;			 
        
		Matrix_xx +=matrix_XX*matrix_XX;
		Matrix_xy +=matrix_XX*matrix_YY;
		Matrix_yy +=matrix_YY*matrix_YY;
		Matrix_yz +=matrix_YY*matrix_ZZ;
		Matrix_zz +=matrix_ZZ*matrix_ZZ;
		Matrix_zx +=matrix_ZZ*matrix_XX;
		
		//cout<< " Matrix_xx  " << Matrix_xx<< " Matrix_xy "<<Matrix_xy <<
		//	" Matrix_yy " << Matrix_yy<< " Matrix_yz "<<Matrix_yz <<
		//	" Matrix_zz " << Matrix_zz<< " Matrix_zx "<<Matrix_zx << endl;
		
	      }


	      // cout<< " Matrix_xx  " << Matrix_xx<< " Matrix_xy "<<Matrix_xy <<
	      // 	" Matrix_yy " << Matrix_yy<< " Matrix_yz "<<Matrix_yz <<
	      // 	" Matrix_zz " << Matrix_zz<< " Matrix_zx "<<Matrix_zx << endl;
	

		
	    }// end of if (minplane_cluster>0
	    
	  } else { // of (x0!=0||y0!=0)&&(pfitTrack->InoTrackCand_list.size()!=0) // finding a good recosntructed muon 
	    
	    //========================================================================
	    //  To remove event if a layer has more than 100 hits //SSE 09/15
	    //========================================================================
	    
	  

	    InoHit_Manager* tmphitlist = InoHit_Manager::APointer;
	    // cout<<"2E_had"<<" "<<E_had<<endl;//SSE
	    int totalhits = 0; // tmphitlist->InoHit_list.size();
	    vector<InoHit*> tmphadcluster1;
	    
	    for (int ixj =0; ixj<nLayer;ixj++){
	      HitBank_All[ixj].clear();
	    }
	    // cout <<"xx1 "<<endl;
	    // cout<<"tmphitlist->InoHit_list.size() = "<<tmphitlist->InoHit_list.size()<<endl;
	    for( unsigned int iyj=0; iyj<tmphitlist->InoHit_list.size(); iyj++) {
	      // cout<<"Inside Loop "<<endl;
	      // cout <<"iyj "<< iyj<<endl;
	      int nxstrip=0, nystrip=0;
	      
	      int irpcmod = tmphitlist->InoHit_list[iyj]->GetRPCmod();
	      // cout <<"ixxx "<<irpcmod <<endl;
	      // cout<<"<inoRPC_pointer- "<<(inoRPC_pointer)<<endl;
	      // cout<<inoRPC_pointer->InoRPCStrip.size()<<endl;

	      for (unsigned int ix=0; ix<inoRPC_pointer->InoRPCStrip.size(); ix++) {
		// cout <<"ix "<< ix<<" "<<irpcmod <<endl;
		if (inoRPC_pointer->InoRPCStrip[ix].first==irpcmod) {
		  nxstrip = inoRPC_pointer->InoRPCStrip[ix].second/100;
		  nystrip = inoRPC_pointer->InoRPCStrip[ix].second%100;
		  break;
		}
	      }
	      // cout <<"x2x2 "<<endl;
	      if (nxstrip*nystrip <100) {
		int plane = tmphitlist->InoHit_list[iyj]->GetZPlane();
		HitBank_All[plane].push_back(tmphitlist->InoHit_list[iyj]);
		tmphadcluster1.push_back(tmphitlist->InoHit_list[iyj]);
	      }
	    }
	    // cout <<"xx2 "<<endl;
	    //  pAnalysis->total_inohits      =0;
	    // pAnalysis->orighits_new	=0;//SSE 250915
	    //pAnalysis->orighits_mod	=0;//SSE 250915
	    //pAnalysis->hit_wo_ghst	=0;//SSE 250915
	    //pAnalysis->hit_wogh_orighits       =0;//SSE 250915
	    //pAnalysis->e_hadron=0;//SSE 250915
	    //pAnalysis->nhits_largest_cluster=0;//SSE 250915
	    //pAnalysis->theta_hadron_shw=0;//SSE 071015
	    // pAnalysis->costheta_hadron_shw=0;//SSE 081015
	    //pAnalysis->phi_hadron_shw=0;//SSE 071015
	    //pAnalysis-> dot_angle_had_shw=0.0;
	    // pAnalysis-> nhits_largest_cluster_selected=0; //SSE 091015	
	    vector<InoHit*> tmphadcluster;
	    vector < vector <InoHit* > > allhadcluster;
	    
	    //pAnalysis->total_inohits=tmphadcluster1.size();
	    // cout<< " before removal :  "<<tmphitlist->InoHit_list.size()<<endl;
	    

	    // cout<< " tmphadcluster1.size() "<<tmphadcluster1.size()<<endl;
	    for (unsigned int ixx=0; ixx<tmphadcluster1.size(); ixx++) {
	      // cout<<" tmphadcluster1[ixx]->GetUID()  "<< tmphadcluster1[ixx]->GetUID()<<endl;
	      if (tmphadcluster1[ixx]->GetUID()==10 || tmphadcluster1[ixx]->GetUID()==11) {
		tmphadcluster1.erase(tmphadcluster1.begin()+ixx);
		ixx--;
	      }
	    }
	    // cout<< " tmphadcluster1.size() "<<tmphadcluster1.size()<<endl;
	    //-------------------------------------------------
	    //orighit calculation 
	    //
	    //-------------------------------------------------
	    //Orighits_new = orighit_calc(tmphadcluster1);
	    // cout<<" orighit function: new orighit "<<Orighits_new<<endl;

	    //cout<< " after removal :  "<<tmphitlist->InoHit_list.size()<<" totalhits before removal "<< totalhits<<endl;
	    //totalhits = tmphitlist->InoHit_list.size();	
	    totalhits=tmphadcluster1.size();
	    // cout<<"tmphadcluster1.size() "<<tmphadcluster1.size()<<endl;
	    //----------------------------------------------
	    //Begin cluster algo
	    //
	    //---------------------------------------------
	    
	    int tmphitcnt = 0;
	    // cout <<"tmphitcnt "<< tmphitcnt<<" "<<totalhits<<endl;
	    while (totalhits>tmphitcnt) {
	      //cout<<"total hits before for (int ix=0; ix<totalhits; ix++) loop "<<totalhits<<endl;
	      for (int ix=0; ix<totalhits; ix++) {
		tmphadcluster.clear();
		InoHit* tmphit = tmphadcluster1[ix];//SSE
		//		cout <<"ixxxxx "<<ix<<" "<< tmphit->GetUID()<<endl;
		
		if (tmphit->GetUID()>1) continue;
		tmphadcluster.push_back(tmphit);
		tmphit->SetUID(2);
		tmphitcnt++;
		// cout<<"tmphitcint "<<tmphitcnt<<endl;
		// cout <<"ix "<<ix<<"  tmphit->GetZPlane() "<< tmphit->GetZPlane()<<" "<<tmphit->GetUID()<<" "<<tmphit->GetXStripNum()<<" "<<tmphit->GetYStripNum()<<endl;
	    
		int initialclustersize = 0;
		int finalclustersize = tmphadcluster.size();
		while (initialclustersize != finalclustersize) {
		  initialclustersize = finalclustersize;
		  // cout<< " initial cluster size "<<initialclustersize<<endl;//SSE
		  for (int iy=0; iy<initialclustersize; iy++) {
		    InoHit* tmphit2 = tmphadcluster[iy];
		    //cout <<"iy "<< iy<<" "<<tmphit2->GetZPlane()<<" "<<tmphit2->GetUID()<<" X strp "<<tmphit2->GetXStripNum()<<" YStrp "<<tmphit2->GetYStripNum()<<" XPOS " <<tmphit2->GetXPos()<<" YPOS " <<tmphit2->GetYPos()<<endl;
		    if (tmphit2->GetUID() >2) continue;
		    tmphit2->SetUID(3);
		    
		    for (int iz=0; iz<totalhits; iz++) {
		      //InoHit* tmphit3 = tmphitlist->InoHit_list[iz];////this line is replaced by next line
		      InoHit* tmphit3 = tmphadcluster1[iz];//SSE
		      // cout <<"iz "<< iz<<" "<<tmphit3->GetUID()<<endl;
		      if (tmphit3->GetUID() >1) continue;
		      //cout <<"iz "<< iz<<" "<<tmphit3->GetUID()<<" "<<tmphit2->GetUID()<<" "<<tmphit3->GetXStripNum()<<" "<<tmphit3->GetYStripNum()<<endl;
		      //  cout<<"dif "<< tmphit2->GetZPlane()<<" "<<tmphit3->GetZPlane()<<" "
		      //   <<tmphit2->GetXPos()<<" "<<tmphit3->GetXPos()<<" "
		      //  <<tmphit2->GetYPos()<<" "<<tmphit3->GetYPos()<<endl;
		      if (abs(tmphit2->GetZPlane()- tmphit3->GetZPlane()) < 5 &&
			  abs(tmphit2->GetXPos()- tmphit3->GetXPos())<0.15 &&
			  abs(tmphit2->GetYPos()- tmphit3->GetYPos())<0.15) {// strip distance 15 cm and 5 layers
			tmphadcluster.push_back(tmphit3);
			tmphit3->SetUID(2);
			tmphitcnt++;
			finalclustersize++;
			// cout <<"2: iz "<< iz<<" "<<tmphit3->GetUID()<<" "<<tmphit2->GetUID()<<" "<<tmphit3->GetXStripNum()<<" "<<tmphit3->GetYStripNum()<<endl;
			
			//   cout <<" tmphadcluster size tmphitcnt final ini size "<< tmphadcluster.size()<<" "<< tmphitcnt<<" "<<finalclustersize<<" "<< initialclustersize<<endl;
		      }
		    }
		  }
		  // cout<< finalclustersize <<" "<<initialclustersize<<endl;
		}
		allhadcluster.push_back(tmphadcluster);
		// cout<<"allhadcluster "<<allhadcluster.size()<<endl;
	      }
	      // cout<<"for (int ix=0; ix<totalhits; ix++) out of loop"<<endl;
	    }
	    // cout <<"tmphitcnt2 "<< tmphitcnt<<" "<<totalhits<<endl;
	    //According to the number of hit, sort them out
	    int nclhit=allhadcluster.size();//no. of clusters 
	    // cout<<" total no. of cluster"<<nclhit<<endl;
	    if (nclhit>1000) nclhit=1000;
	    int in_r[1000];
	    int clhits[1000];
	    for (int ix=0; ix<1000;ix++){
	      clhits[ix]=0;//added by SSE, initialisation needed, otherwise gives problem during storing num_hadron		
	    }
	    
	    for (int ix=0; ix<nclhit; ix++) {
	      in_r[ix] = ix;//store cluster number
	      clhits[ix] = allhadcluster[ix].size();//no of hits in "ix th" cluster
	      // cout<< "ix :"<< ix<< " hit in ix th cluster  "<< clhits[ix]<<endl;//SSE
	    }

	    for (int ix=0; ix<nclhit; ix++){
	      for (int iy=ix+1; iy<nclhit; iy++) {
		if (clhits[iy]>clhits[ix]) {
		  int tmpid = in_r[ix];
		  in_r[ix] = in_r[iy];// 0 th element of in_r stores the largest cluster's serial no
		  in_r[iy] = tmpid;
		  
		  tmpid =clhits[ix];
		  clhits[ix] = clhits[iy];
		  clhits[iy] = tmpid;
		  //	cout<< "clhits[ix] "<<clhits[ix]<<" clhits[iy] "<<clhits[iy]<<endl;	
		}
	      }
	    }

	    num_hadron=clhits[0]; 
	    vector<InoHit*> tmphadcluster2;//vector to store hits from largest cluster
	    float cluster_plane[1000];
	    int ii=0;
	    for (unsigned int ixx=0; ixx<allhadcluster.size(); ixx++) {
	      int ix= in_r[ixx];
	      //cout<<"ixxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx "<< ix<<" size  "<<allhadcluster[ix].size()<<endl;
	      for (unsigned int iy=0; iy<allhadcluster[ix].size(); iy++) {
		
		if (ixx==0) {
		  allhadcluster[ix][iy]->SetUID(3333);
		  ii++;
		  tmphadcluster2.push_back(allhadcluster[ix][iy]);
		  cluster_plane[ii]=allhadcluster[ix][iy]->GetZPlane(); 
		}
		//cout<<"ix "<<ix<<" "<< iy<<" "<<allhadcluster[ix][iy]->GetZPlane()<<" "
		//  <<allhadcluster[ix][iy]->GetXPos()<<" "
		// <<allhadcluster[ix][iy]->GetYPos()<<" "
		
		
		// <<allhadcluster[ix][iy]->GetXStripNum()<< " "
		//<<allhadcluster[ix][iy]->GetYStripNum()<<endl;
	      }
	    }

	    //pAnalysis->hist_nhits_LargestCluster_E->Fill(num_hadron, E_had);//SSE

	    //==================================================================
	    // orighit calculation taking the output of cluster algorithm as 
	    // input of orighit algo  
	    //SSE 09/15
	    //==================================================================
	    Orighits_cluster = orighit_calc(tmphadcluster2);
	    
	    //histogram filling for Orighits
	    //pAnalysis->hist_orighits_mod_E->Fill(Orighits_cluster, E_had);
	    //
	    //==================================================================
	    //choosing the largest cluster contains atleast two plane
	    //=================================================================
	    
	    for( unsigned int ijk=0; ijk<tmphadcluster2.size(); ijk++) {
	      if (ijk>0 && cluster_plane[0]!=cluster_plane[ijk]) {
		minplane_cluster++;
		break;
	      }
	    } 
	    //cout<< " mincluster_plane "<<minplane_cluster<<endl;
	    
	    //===================================================================
	    //Algorithm to find direction of hadronic shower
	    //=================================================================== 
	    
	    if (minplane_cluster>0){// lagest cluster must contains min two planes.
	      nhits_large_clust_selected=num_hadron;
	      
	      //Calculation of center of largsest cluster	
	      
	      float center_x=0.0, center_y=0.0, center_z=0.0;
	      
	      for( unsigned int mn=0; mn<tmphadcluster2.size(); mn++) {
		center_x +=tmphadcluster2[mn]->GetXPos();
		center_y +=tmphadcluster2[mn]->GetYPos();
		center_z +=tmphadcluster2[mn]->GetZPos();
	      }
	      center_x= center_x/(tmphadcluster2.size()); 
	      center_y= center_y/(tmphadcluster2.size()); 
	      center_z= center_z/(tmphadcluster2.size()); 
	      //cout<< " center_x "<< center_x<< " center_y "<< center_y<<" center_z "<< center_z<<endl;
	      // and Calculating matrix element to find hadron direction.
	      
	      for( unsigned int mn=0; mn<tmphadcluster2.size(); mn++) {
		
		matrix_XX=(tmphadcluster2[mn]->GetXPos() - center_x);
		matrix_YY=(tmphadcluster2[mn]->GetYPos() - center_y);
		matrix_ZZ=(tmphadcluster2[mn]->GetZPos() - center_z);
		
		//	cout<< " xx  "<< tmphadcluster2[mn]->GetXPos()
		//	<< " yy  "<< tmphadcluster2[mn]->GetYPos()
		//	<< " zz  "<< " "<<tmphadcluster2[mn]->GetZPos()<<endl<<endl;
		
		//	cout<< " matrix_XX "<< matrix_XX<< " matrix_YY "<< matrix_YY<< "  matrix_ZZ "<<  matrix_ZZ<<endl;			 
        
		Matrix_xx +=matrix_XX*matrix_XX;
		Matrix_xy +=matrix_XX*matrix_YY;
		Matrix_yy +=matrix_YY*matrix_YY;
		Matrix_yz +=matrix_YY*matrix_ZZ;
		Matrix_zz +=matrix_ZZ*matrix_ZZ;
		Matrix_zx +=matrix_ZZ*matrix_XX;
		
		//cout<< " Matrix_xx  " << Matrix_xx<< " Matrix_xy "<<Matrix_xy <<
		//	" Matrix_yy " << Matrix_yy<< " Matrix_yz "<<Matrix_yz <<
		//	" Matrix_zz " << Matrix_zz<< " Matrix_zx "<<Matrix_zx << endl;
		
	      }
	      
	      // cout<< " Matrix_xx  " << Matrix_xx<< " Matrix_xy "<<Matrix_xy <<
	      // 	" Matrix_yy " << Matrix_yy<< " Matrix_yz "<<Matrix_yz <<
	      // 	" Matrix_zz " << Matrix_zz<< " Matrix_zx "<<Matrix_zx << endl;
	      
	    } //if (minplane_cluster>0)// SSE 03 Nov 2015
	    
	    //==================================================================================
	    //				Algorithm of ghost hits removal
	    //=================================================================================
	    //create hitbank using hits from largest cluster from cluster algo.
	    //initialise the HitBank_largestClust otherwise some info stored init
	    for (int ixj =0; ixj<nLayer;ixj++){
	      HitBank_largestClust[ixj].clear();
	    }
	    
	    // Filling HIt Bank (needed in GHR algo)
	    for( unsigned int mn=0; mn<tmphadcluster2.size(); mn++) {
	      HitBank_largestClust[tmphadcluster2[mn]->GetZPlane()].push_back(tmphadcluster2[mn]); //filling HitBank
	    }
	    
	    int iFlag_had[500][300];
	    //Initialisation
	    int iXHFill[500][500];// 1st array is for layer, 2nd on n_had (n_had no. of layer have one hit) 
	    int iYHFill[500][500];
	    
	    int diff_old=1000;
	    int diff=1000;
	    int n_had=0,n_had1=0;
	    int layer_prev=-10000;
	    int layer_post=-10000;
	    
	    int count_hit_layer[250];
	    
	    //initilisation of iFlag_had
	    for(int ixj=0; ixj<500; ++ixj) {
	      for(unsigned int jxk=0; jxk<HitBank_largestClust[ixj].size(); ++jxk) {
		iFlag_had[ixj][jxk]=0;				
	      }
	    }
	    //initialisation of iXHFill[ixj][jxk] and iYHFill[ixj][jxk], array of layer and n_had		
	    for(unsigned int ixj=0; ixj<500; ++ixj) {
	      for(int jxk=0; jxk<500; ++jxk) {
		iXHFill[ixj][jxk]=-100; 	 iYHFill[ixj][jxk]=-100; 				
	      }
	    }
	    for (int iter=0;iter<1;iter++){
	      layer_prev=-10000;//initialise previous layer
	      layer_post=-10000;//initialise post layer
	      
	      //To find which layer has only one hit
	      for(int ixj=0; ixj<nLayer; ++ixj) {//loop on layer number 
		count_hit_layer[ixj]=0;
		for(unsigned int jxk=0;jxk<HitBank_largestClust[ixj].size(); ++jxk){
		  if (iFlag_had[ixj][jxk]==0 && HitBank_largestClust[ixj][jxk]->GetUID()==3333){
		    count_hit_layer[ixj]=count_hit_layer[ixj]+1;//count number of hit in a layer
		  }
		}
		if (count_hit_layer[ixj]==1) {
		  // cout<<ij<< " layer has  count "<< count_hit_layer[ij]<<endl;
		}
	      }
	      
	      //Start finding ghost hit 
	      //	Start loop in upward direction 
	      
	      for(int iyj=0; iyj<nLayer; ++iyj) {//loop on layer number
		if((count_hit_layer[iyj]==1 && iyj>=0) || (count_hit_layer[iyj]>1 && iyj==layer_post)) {
		  //if a layer having one hit or connected to the hit of post layer
		  //start the loop on hits of that layer
		  //start connecting to hit in next layer and also identify ghost hits
		  //
		  for(unsigned int jyk=0; jyk<HitBank_largestClust[iyj].size(); ++jyk) {
		    if (count_hit_layer[iyj]==1 && iFlag_had[iyj][jyk]==0) {
		      n_had=n_had+1; 
		      iFlag_had[iyj][jyk]=(1+10*iter);//Flag changing
		      iXHFill[iyj][n_had]= HitBank_largestClust[iyj][jyk]->GetXStripNum() ;
		      iYHFill[iyj][n_had]= HitBank_largestClust[iyj][jyk]->GetYStripNum() ;
		      //cout<< "layer "<<iyj << " has 1 hit " <<iXHFill[iyj][n_had] << " "<<iYHFill[iyj][n_had] << " n_had  "<< n_had<<endl; 
		    }

		    diff_old=1000;		
		    //cout	<< " start connecting hit (of next layer with)  iXHFill[iyj][n_had]" 
		    //	<<iXHFill[iyj][n_had]<<" "<<iYHFill[iyj][n_had]<<endl;	
	     		
		    //Start the if on connected hit
		    if(HitBank_largestClust[iyj][jyk]->GetXStripNum()==iXHFill[iyj][n_had] &&
		       HitBank_largestClust[iyj][jyk]->GetYStripNum()==iYHFill[iyj][n_had] &&
		       (iFlag_had[iyj][jyk]==0 ||iFlag_had[iyj][jyk]==(1+10*iter))) {

		      if (iFlag_had[iyj][jyk]==0) {iFlag_had[iyj][jyk]=(2+10*iter);} //Flag changing
		      for(unsigned int kl=0; iyj<nLayer-1 && kl<HitBank_largestClust[iyj+1].size(); ++kl) {
			//start loop on hit points of consecutive post layer
			//    cout	<< "next layer: "<< (iyj+1)
			//		<< " hit no. "<< kl
			//		<<" X: "<<HitBank_largestClust[iyj+1][kl]->GetXStripNum() << "Y: "<< HitBank_largestClust[iyj+1][kl]->GetYStripNum()<<endl;
			//if next layer's hit is non-connected then find nearest hit
			if(iFlag_had[iyj+1][kl]==0) {           
			  diff=abs(HitBank_largestClust[iyj][jyk]->GetXStripNum()-HitBank_largestClust[iyj+1][kl]->GetXStripNum())
						 +abs(HitBank_largestClust[iyj][jyk]->GetYStripNum()-HitBank_largestClust[iyj+1][kl]->GetYStripNum());
			  // 		cout<<" post: diff "<< diff<< " diff _old "<< diff_old<<endl;
			  if(diff<diff_old){
			    diff_old=diff;
			    layer_post=iyj+1;		
			    iXHFill[iyj][n_had]= HitBank_largestClust[iyj][jyk]->GetXStripNum();
			    iYHFill[iyj][n_had]= HitBank_largestClust[iyj][jyk]->GetYStripNum(); 
			    iXHFill[iyj+1][n_had]= HitBank_largestClust[iyj+1][kl]->GetXStripNum() ;
			    iYHFill[iyj+1][n_had]= HitBank_largestClust[iyj+1][kl]->GetYStripNum();
			  }
			}//end of if(iFlag_had[iyj-1][kl]==0)
		      }//end loop on hit points of  post layer
	
		      //   cout	<<" *layer "<<iyj<<" n_had "<< n_had<<" "<<iXHFill[iyj][n_had]<<" "<<iYHFill[iyj][n_had]
		      //		<<" post layer "<< iyj+1<< " "<<iXHFill[iyj+1][n_had]<< " "<<iYHFill[iyj+1][n_had]<<endl;
		    } else { //end if on connected hit
		      int irpcmod =HitBank_largestClust[iyj][jyk]->GetRPCmod();
		      int nxstrip=0;
		      int nystrip=0;
		      for (unsigned int ix=0; ix<inoRPC_pointer->InoRPCStrip.size(); ix++) {
			if (inoRPC_pointer->InoRPCStrip[ix].first==irpcmod) {
			  nxstrip = inoRPC_pointer->InoRPCStrip[ix].second/100;
			  nystrip = inoRPC_pointer->InoRPCStrip[ix].second%100;
			  break;
			}
		      }
			
		      if (nxstrip>1 && nystrip>1 &&
			  ((HitBank_largestClust[iyj][jyk]->GetXStripNum()==iXHFill[iyj][n_had]&&
			    HitBank_largestClust[iyj][jyk]->GetYStripNum()!=iYHFill[iyj][n_had])||
			   (HitBank_largestClust[iyj][jyk]->GetXStripNum()!=iXHFill[iyj][n_had]&& 
			    HitBank_largestClust[iyj][jyk]->GetYStripNum()==iYHFill[iyj][n_had]))&& 
			  iFlag_had[iyj][jyk]==0){//if on unconnected hit
			iFlag_had[iyj][jyk]=-1-10*iter;//flag realted to ghost hits
			// cout<<"layer : "<<iyj<< " "<< HitBank_largestClust[iyj][jyk]->GetXStripNum() <<" "<<HitBank_largestClust[iyj][jyk]->GetYStripNum()<< "is ghost hit"<<endl;
		      }//endi of else if
		      //		}//if(HitBank_largestClust[iyj][jyk]->GetUID()==3333)
		    }
		  }//end of loop on hit bank on the layer which has one hit or connected hit
		}//end of if(HitBank_largestClust[iyj].size()==1 && iyj>0)
	      }//end ....loop on layer no.
	      //cout<<"End of loop on post-layer, no. of layer having one hit " <<n_had<<endl<<endl;


	      //-------------------------------------------------
	      //start the loop on the downwards direction
	      
	      for(int iyj=nLayer-1; iyj>=0; --iyj) {//loop on layer number
		//if((HitBank[iyj].size()==1 && iyj<250)||(HitBank[iyj].size()>1 && iyj==layer_prev)) {//if on layer having one hit
		if((count_hit_layer[iyj]==1 && iyj<nLayer)||(count_hit_layer[iyj]>1 && iyj==layer_prev)) {//if on layer having one hit
		  // cout<<" Enter SEcond loop for dowards layer"<<endl;
		  for(unsigned int jyk=0; jyk<HitBank_largestClust[iyj].size(); ++jyk) {
		    //	if(HitBank_largestClust[iyj][jyk]->GetUID()==3333){		
		    //if (HitBank_largestClust[iyj].size()==1){
		    if (count_hit_layer[iyj]==1 && iFlag_had[iyj][jyk]==1+10*iter) {
		      n_had1=n_had1+1; iFlag_had[iyj][jyk]=1+10*iter;
		      iXHFill[iyj][n_had1]= HitBank_largestClust[iyj][jyk]->GetXStripNum() ;	
		      iYHFill[iyj][n_had1]= HitBank_largestClust[iyj][jyk]->GetYStripNum() ;
		      //  cout<< "layer "<<iyj << " have one hit "<<iXHFill[iyj][n_had1] << " "<<iYHFill[iyj][n_had1] <<endl; 
		    }	
		    //checking previous layer 
		    diff_old=1000;
		    // cout<< "connect hit of prev layer with  iXHFill[iyj][n_had]" <<iXHFill[iyj][n_had1]<<" "<<iYHFill[iyj][n_had1]<<endl;	
		    if (HitBank_largestClust[iyj][jyk]->GetXStripNum()==iXHFill[iyj][n_had1]&& HitBank_largestClust[iyj][jyk]->GetYStripNum()==iYHFill[iyj][n_had1] 
			&& (iFlag_had[iyj][jyk]==0|| iFlag_had[iyj][jyk]==(1+10*iter)||iFlag_had[iyj][jyk]==(2+10*iter))){//if on connected hit
		      if (iFlag_had[iyj][jyk]==0)iFlag_had[iyj][jyk]=3+10*iter;
		      for(unsigned int kl=0; iyj>0 && kl<HitBank_largestClust[iyj-1].size(); ++kl) {
			//  cout<< "prev. layer: "<< iyj-1<< " hit no. "<< kl
			//	  <<" "<<HitBank_largestClust[iyj-1][kl]->GetXStripNum() 
			//	  << " "<< HitBank_largestClust[iyj-1][kl]->GetYStripNum()<<"iFlag"<< iFlag_had[iyj-1][kl]<<endl;
			
			if(iFlag_had[iyj-1][kl]==0){    
			  diff=abs(HitBank_largestClust[iyj][jyk]->GetXStripNum()-HitBank_largestClust[iyj-1][kl]->GetXStripNum())
			    +abs(HitBank_largestClust[iyj][jyk]->GetYStripNum()-HitBank_largestClust[iyj-1][kl]->GetYStripNum());
			  //cout<<" calculate difference "<<endl;
			  //cout<<"pre : diff "<< diff<< "diff _old "<< diff_old<<endl;
			  if(diff<diff_old){
			    diff_old=diff;
			    layer_prev=iyj-1;		
			    iXHFill[iyj][n_had1]=HitBank_largestClust[iyj][jyk]->GetXStripNum() ; 
			    iYHFill[iyj][n_had1]= HitBank_largestClust[iyj][jyk]->GetYStripNum(); 
			    iXHFill[iyj-1][n_had1]=HitBank_largestClust[iyj-1][kl]->GetXStripNum();
			    iYHFill[iyj-1][n_had1]=HitBank_largestClust[iyj-1][kl]->GetYStripNum();
			  }
			} else if (iFlag_had[iyj-1][kl]==1+10*iter || iFlag_had[iyj-1][kl]==2+10*iter) {
			  //if two hits are already connected by upwards case
			  //then take those are nearest
			  //which help to avoid the case same hit of a particular layer
			  //connect two diiferent hits in previous(or succesive layer)
			  //during upwards & downwards direction checking
			  
			  iXHFill[iyj][n_had1]=HitBank_largestClust[iyj][jyk]->GetXStripNum() ; 
			  iYHFill[iyj][n_had1]= HitBank_largestClust[iyj][jyk]->GetYStripNum(); 
			  iXHFill[iyj-1][n_had1]=HitBank_largestClust[iyj-1][kl]->GetXStripNum();
			  iYHFill[iyj-1][n_had1]=HitBank_largestClust[iyj-1][kl]->GetYStripNum();
			  diff_old=-0.1;//arbitrary -ve number so that diff from prev. if never less than that.
			  layer_prev=iyj-1;
			  //cout<<"forced to particular hit"<<endl;
			  //continue;
			  //break;
			}
			
		      }//loop on hit points of pre layer
		      //cout	<<" ***layer "<<iyj<<" n_had1 "<<n_had1 <<" "<<iXHFill[iyj][n_had1]<<" "<<iYHFill[iyj][n_had1]
		      //		<<" prev layer "<< iyj-1<< " "<<iXHFill[iyj-1][n_had1]<< " "<<iYHFill[iyj-1][n_had1]<<endl;
		    } else { //end if on connected hit
		      int irpcmod =HitBank_largestClust[iyj][jyk]->GetRPCmod();
		      int nxstrip=0;
		      int nystrip=0;
		      for (unsigned int ix=0; ix<inoRPC_pointer->InoRPCStrip.size(); ix++) {
			if (inoRPC_pointer->InoRPCStrip[ix].first==irpcmod) {
			  nxstrip = inoRPC_pointer->InoRPCStrip[ix].second/100;
			  nystrip = inoRPC_pointer->InoRPCStrip[ix].second%100;
			  break;
			}
		      }
		      if (nxstrip>1 && nystrip>1 &&
			  ((HitBank_largestClust[iyj][jyk]->GetXStripNum()==iXHFill[iyj][n_had1]&& 
			    HitBank_largestClust[iyj][jyk]->GetYStripNum()!=iYHFill[iyj][n_had1])||
			   (HitBank_largestClust[iyj][jyk]->GetXStripNum()!=iXHFill[iyj][n_had1]&& 
			    HitBank_largestClust[iyj][jyk]->GetYStripNum()==iYHFill[iyj][n_had1])) && 
			  iFlag_had[iyj][jyk]==0){//if on unconnected hit
			iFlag_had[iyj][jyk]=-1-10*iter;
			// cout<< "layer: "<< iyj<< " "<<HitBank_largestClust[iyj][jyk]->GetXStripNum() 
			//     <<" "<<HitBank_largestClust[iyj][jyk]->GetYStripNum()<< "is ghost hit"<<endl;
		      }//end of else if
		    }
		    //	}//if(HitBank_largestClust[iyj][jyk]->GetUID()==3333)
		  }//end of loop on hit bank	
		}//end of if(HitBank[iyj].size()==1 && iyj>0)
		//cout<<"layer "<<iyj<<" "<<iXHFill[iyj][n_had]<<" "<<iYHFill[iyj][n_had]<<" prev "<< iXHFill[iyj-1][n_had]<< " "<<iYHFill[iyj-1][n_had] <<endl;
	      }//end ....loop on layer no.
	      // cout<<" one  in a layer " <<n_had1<<endl;
	    } //for loop on iter //End of ghost removal algorithm
	    
	    int count_h=0;//counter to calculate remaining hit after ghost hit removal algo
	    vector<InoHit*> tmphadcluster3;//Create vector of InoHit with hits remained after removing ghost hits
	    for(int iyj=0; iyj<nLayer; ++iyj) {//loop on layer number 
	      count_hit_layer[iyj]=0;
	      for(unsigned int jyk=0;jyk<HitBank_largestClust[iyj].size(); ++jyk){
		// cout<< "layer "<<iyj<< "   x: "<<HitBank_largestClust[iyj][jyk]->GetXStripNum() <<" y:"<<HitBank_largestClust[iyj][jyk]->GetYStripNum()
		//	  << "    iFlag     "<< iFlag_had[iyj][jyk] <<" GetUID() "<<HitBank_largestClust[iyj][jyk]->GetUID() <<endl;
		
		if (iFlag_had[iyj][jyk]>=0) {
		  HitBank_largestClust[iyj][jyk]->SetUID(4444);
		  tmphadcluster3.push_back(HitBank_largestClust[iyj][jyk]);
		  count_h=count_h+1;
		}
		//cout<< " GetUID() "<<HitBank_largestClust[iyj][jyk]->GetUID() <<endl;
	      }
	    }

	    Hit_wo_ghst=count_h;
	    //cout<< " Hit_wo_ghst "<< count_h<< " "<<tmphadcluster3.size()<<endl;
	    //==========================================================================
	    //			APPLY orighit algorithm after removing ghost hits
	    //==========================================================================
	    wogh_Orighits = orighit_calc(tmphadcluster3);
	    // cout<<" chk: orighit function after ghost hit removal"<<wogh_Orighits<<endl;
	    
	    //cout<< " hit_wo_ghst  "<< " "<<Hit_wo_ghst<<" pAnalysis->momin[0] "<<pAnalysis->momin[0]<< " pAnalysis->momin[1] "<<pAnalysis->momin[1] <<" E_had"<<E_had<<endl;
	    //pAnalysis->hh_woghst_E->Fill(Hit_wo_ghst, E_had);//SSE
	    
	    //histogram filling for Orighits 
	    // pAnalysis->hist_wogh_orighits_E->Fill(wogh_Orighits , E_had);
	    cout<<" Orighits_cluster "<<Orighits_cluster<<"  "<<Hit_wo_ghst<<" wogh_Orighits  "<<wogh_Orighits<<endl;
	    
	  } //else x0= // loop jyk  // else of finding a good recosntructed muon 
	  
	  
	  //=========================================================================
	  //hadron shower direction
	  //calculation of eigen value and eigen vector
	  //==============================================================================
	  double Rx_had=0.0,  Ry_had=0.0,  Rz_had=0.0;
	  float R_had=0.0, theta_had=0.0, phi_had=0.0; //SSE 071015
	  
	  double matrix_elemnt[3][3]={{Matrix_xx,Matrix_xy,Matrix_zx}, {Matrix_xy, Matrix_yy,  Matrix_yz}, {Matrix_zx, Matrix_yz,  Matrix_zz}};
	  //double matrix_elemnt[3][3]={{1,0,0}, {0,0,1}, {0,1,0}}; //for testing
	  TMatrixD mat_dir(3,3);
	  for (int row=0; row<3; row++) { 
	    for (int column=0;column<3;column++){ 
	      mat_dir[row][column] = matrix_elemnt[row][column];
	    } 
	  }
	  
	  // mat_dir.Print();
	  const TMatrixDEigen eigen(mat_dir);
	  //const TVectorD eigenVal = eigen.GetEigenValues();
	  const TMatrixD eigenVal = eigen.GetEigenValues();
	  const TMatrixD eigenVect=eigen.GetEigenVectors();
	  // eigenVal.Print();
	  
	  //cout<< " eigen values " <<eigenVal[0]<< " "<<eigenVal[1] << " " << eigenVal[2] << endl;
	  // eigenVect.Print();
	  
	  //have to check whether sorting required
	  for (int row=0; row<3; row++) {
	    for (int column=0;column<1;column++){
	      if (row==0)Rx_had=eigenVal[0][0]*eigenVect[row][column];
	      if (row==1)Ry_had=eigenVal[0][0]*eigenVect[row][column];
	      if (row==2)Rz_had=eigenVal[0][0]*eigenVect[row][column];
	      //cout<< eigenVect[row][column]<<endl;
	    }
	  }
	  G4ThreeVector tmp3v(Rx_had,Ry_had,Rz_had);
	  R_had = tmp3v.mag();
	  theta_had = tmp3v.theta();
	  phi_had = tmp3v.phi();
	  //costheta_had=tmp3v.cosTheta();
	  dotAngle= tmp3v.cosTheta(tmp3pin_had);
	  //cout<< " Rx_had "<< Rx_had << " Ry_had "<< Ry_had<< " Rz_had "<<Rz_had<<endl;
	  //cout<< " R "<<R_had  <<" had shower theta "<< theta_had<< " had shower phi  "<<phi_had <<endl;
	  //cout<< " had shower costheta "<< costheta_had<< " dotAngle "<< dotAngle<<endl;
	  
	  pAnalysis->e_hadron=E_had;//SSE 250915
	  pAnalysis->orighits_trape=Orighits_trape; //SSE 250915
	  pAnalysis->nhits_largest_cluster=num_hadron;//SSE 250915
	  pAnalysis->orighits_cluster=Orighits_cluster;////SSE 250915
	  pAnalysis->hit_wo_ghst=Hit_wo_ghst;//SSE 250915
	  pAnalysis->hit_wogh_orighits=wogh_Orighits ;//SSE 250915
	  if(theta_had!=0)pAnalysis->theta_hadron_shw=theta_had ;//SSE 071015
	  //pAnalysis->costheta_hadron_shw=costheta_had ;//SSE 071015
	  if(phi_had!=0)pAnalysis->phi_hadron_shw=phi_had ;//SSE 071015
	  if( dotAngle!=0)pAnalysis-> dot_angle_had_shw=dotAngle ;//SSE 081015
	  pAnalysis->nhits_largest_cluster_selected=nhits_large_clust_selected;//SSE 
	  if(theta_had!=0) {
	    for(int eval=0; eval<3; eval++) {
	      pAnalysis->had_eigen_val[eval] = eigenVal[eval][eval];
	    }
	  }
	  pAnalysis->ntotcl = 1000*allremcls + allremhit;
	  pAnalysis->ntotst = 1000*(xtrkstrip.size()-tottrkXstrp) + (ytrkstrip.size()-tottrkYstrp);
	  //pAnalysis->inohits =0;
	  pAnalysis->orighits =0;
	  pAnalysis->x_hits =0;
	  pAnalysis->y_hits =0;
	
	  pAnalysis->inohits =0;
	  int nc_cl=0; int nc_ht=0; int ccc_ht=0 ; int t_ht=0;
	  for (int kl=0; kl<totclustersize; ++kl){
	    
	    if(totcluster[kl]->GetInShower()==0&&totcluster[kl]->GetInTrack()==0){
	      nc_cl++ ; nc_ht += totcluster[kl]->GetHitEntries();
	    } else if(totcluster[kl]->GetInShower()==1){
	      ccc_ht+=totcluster[kl]->GetHitEntries();
	    } else if(totcluster[kl]->GetInTrack()==1) {
	      t_ht+=totcluster[kl]->GetHitEntries();
	    }
	  }
	
	  pAnalysis->inohits =ccc_ht;
	  if(pAnalysis->isXtermOut==2){
	    pAnalysis->ascii_output <<"hadron shower collected: " <<"shower:"<< ccc_ht << " + : track "<<  t_ht <<" + not selected " <<nc_ht <<endl;
	  }
	}

	//===========================================
	//previous orighit calc
	//=====================================

	int Orighits=0;
	int xxnhits =0;
	int yynhits =0;
	// cout<< " xshwstrip.size "<< xshwstrip.size()<< " yshwstrip.size "<< yshwstrip.size() <<endl;	  
	for(int kl=0; kl<nLayer; ++kl){
	  
	  if(xshwstrip.size() && yshwstrip.size()){
	    xxnhits =0;
	    yynhits =0;
	    for(unsigned nn=0; nn< xshwstrip.size(); nn++) {
	      if(xshwstrip[nn].first==int(kl)){
		xxnhits++;
	      }
	    }
	    for (unsigned nn=0; nn< yshwstrip.size(); nn++){
	      if(yshwstrip[nn].first==int(kl)){
		yynhits++;
	      }
	    }
	    //if(showerXHit[kl].size()||showerYHit[kl].size()){
	    //xxnhits += showerXHit[kl].size();
	    //yynhits += showerYHit[kl].size();
	    //pAnalysis->x_hits += xxnhits;
	    //pAnalysis->y_hits += yynhits;
	    Orighits +=max(xxnhits,yynhits);
	  }
	}
	
	// cout<<" orighits "<<Orighits<<endl;
	//histogram filling for Orighits //SSE
	//pAnalysis->hh_E->Fill(Orighits, E_had);//SSE
	//	TProfile* profx =pAnalysis->hh_E->ProfileX("Profile_orighits",-10, 90, "o");//firstybin, lastybin, "o" : original range//SSE
	//		profx->Draw();//SSE
	
	if( abs(pAnalysis->x_hits) <100000) pAnalysis->x_hits=xshwstrip.size();else{pAnalysis->x_hits=-1;}//asm//temporarily commented out
	if( abs(pAnalysis->y_hits) <100000) pAnalysis->y_hits=yshwstrip.size();else{pAnalysis->y_hits=-1;}
	if( abs(pAnalysis->x_hits)<100000 && abs(pAnalysis->y_hits) <100000)pAnalysis->orighits=Orighits;else{pAnalysis->orighits=-1;}

	if(pAnalysis->isXtermOut==2){
	  pAnalysis->ascii_output <<  "(x_hits y_hit orighits inohits)  (" << pAnalysis->x_hits << " " <<pAnalysis->y_hits << " " << pAnalysis->orighits<< " " << pAnalysis->inohits<< ")"<<endl;
	}
      } else { //if (pfitTrack)
        cout <<"XXXXXXXXXXXXXXXXXXXXXXX InoTrackCand_Manager::APointer is not found "<<endl; 
        
      } //loop jk
      //      pAnalysis->pEventTree->Fill(); //VALGRIND
      int nUp=0;
      int nDown=0;
      
      if (pAnalysis->isVisOut) {
        for (unsigned kl=0; kl< pfitTrack->InoTrackCand_list.size() ; kl++) {
          pAnalysis->H->NRecTracks++;
          pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object/
          pAnalysis->Hp->TrackType=-140;// Track Type: -1: hits, -2: clulster, -3: triplet, -4: track -14: particle info
          pAnalysis->Hp->ZZ =((((0.5*(paradef->GetnIRLayer()+1))*LayerThickness)+pAnalysis->poszvx[kl]*cm)/(LayerThickness*m));// vertex z incase of particle ifo
          
          pAnalysis->Hp->XX=pAnalysis->posxvx[kl]/100; // vertex x incase of particle info
          pAnalysis->Hp->YY=pAnalysis->posyvx[kl]/100; // vertex y incase of particle info
          pAnalysis->Hp->pmag=pAnalysis->trkmm[kl]; // vertex y incase of particle info
          pAnalysis->Hp->pt=pAnalysis->trkth[kl]; // vertex y incase of particle info
          pAnalysis->Hp->pp=pAnalysis->trkph[kl]; // vertex y incase of particle info
          
          for ( unsigned int lm =0; lm<pfitTrack->InoTrackCand_list[kl]->GetClusterEntries(); lm++){
            
            pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object  //VALGRIND
            for(unsigned int mn=0;mn<pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm]->GetHitEntries();mn++){
              
              int zplane = pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm]->HitsInCluster[mn]->GetZPlane();
              if (zplane >=250) {zplane -=250;}
              pAnalysis->Hp->ZZ=  zplane;//pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm]->GetZPlane(); //os();
              pAnalysis->Hp->XX=pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm]->HitsInCluster[mn]->GetXPos();
              pAnalysis->Hp->YY=pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm]->HitsInCluster[mn]->GetYPos();
              if(pfitTrack->InoTrackCand_list[kl]->GetFitType()==1){
                
                pAnalysis->Hp->TrackType=-440; // track type up -44 + 0 FIT UP
                pAnalysis->Hp->Fup=0;// Finder set up
                pAnalysis->Hp->FitUpNum = nUp;
                pAnalysis->H->NFitUp=nUp+1;
              } else {
                pAnalysis->Hp->TrackType=-441; // track type down -44 + 1 FIT DOWN
                pAnalysis->Hp->Fup=1;//Finder set down
                pAnalysis->Hp->FitDownNum = nDown;
                pAnalysis->H->NFitDown=nDown+1;
              }
            }
            //--------------------------
          }
	  if (pfitTrack->InoTrackCand_list[kl]->GetFitType()==1) {
	    nUp++;
	  } else {
	    nDown++;
	  }
	}
      } // if (pAnalysis->isVisOut



   
    } //if (pinotrack)
  } //if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput==3 || pAnalysis->InputOutput==5) 
  
 
  //   if (inoRPC_pointer) { delete inoRPC_pointer; inoRPC_pointer=0;}
  //Create Cmvhits from SiPM 
  CreateCmvHit();
  
  //Fill cmv hit data:
  //22032022

  pAnalysis->pRootFile->cd();
  pAnalysis->cmv_nhit = CmvHit_pointer->CmvHit_list.size();
  cout<<"  CmvHit_pointer->CmvHit_list.size():  "<<  CmvHit_pointer->CmvHit_list.size()<<endl;
  cout <<"  pAnalysis->cmv_nclusthit "<<  pAnalysis->cmv_nhit<<endl;

  unsigned int hitId =0;
  if (pAnalysis->cmv_nhit >pAnalysis->cmv_nhtmx) pAnalysis->cmv_nhit =pAnalysis->cmv_nhtmx;
  for (unsigned ij=0; ij<CmvHit_pointer->CmvHit_list.size() && ij<pAnalysis->cmv_nhit; ij++) {

    hitId = CmvHit_pointer->CmvHit_list[ij]->GetPlane();
    hitId<<=2;
    hitId+=CmvHit_pointer->CmvHit_list[ij]->GetLayer();
    hitId<<=7;
    hitId+=CmvHit_pointer->CmvHit_list[ij]->GetStrip();
    hitId<<=2; //just shifted by 2 as we did this while forming strip	
    //pAnalysis->cmv_clustpdgid[ij] =CmvHit_pointer->CmvHit_list[ij]->GetpdgId();

    cout<<"clustid "<<hitId<<" "<<endl;


	
    pAnalysis->cmv_hitid[ij] =hitId;

    cout<< pAnalysis->cmv_hitid[ij]<<endl;

    pAnalysis->cmv_hitLeTim[ij] =CmvHit_pointer->CmvHit_list[ij]->GetLeTime();
    pAnalysis->cmv_hitRiTim[ij] =CmvHit_pointer->CmvHit_list[ij]->GetRiTime();

    pAnalysis->cmv_hitLePul[ij] =CmvHit_pointer->CmvHit_list[ij]->GetLePulse();
    pAnalysis->cmv_hitRiPul[ij] =CmvHit_pointer->CmvHit_list[ij]->GetRiPulse();

	
    pAnalysis->cmv_hitTrueposx[ij] =CmvHit_pointer->CmvHit_list[ij]->GetTruePosX();
    pAnalysis->cmv_hitTrueposy[ij] =CmvHit_pointer->CmvHit_list[ij]->GetTruePosY();
    pAnalysis->cmv_hitTrueposz[ij] =CmvHit_pointer->CmvHit_list[ij]->GetTruePosZ();

    pAnalysis->cmv_hitRecoposx[ij] =CmvHit_pointer->CmvHit_list[ij]->GetRecoPosX();
    pAnalysis->cmv_hitRecoposy[ij] =CmvHit_pointer->CmvHit_list[ij]->GetRecoPosY();
    pAnalysis->cmv_hitRecoposz[ij] =CmvHit_pointer->CmvHit_list[ij]->GetRecoPosZ();

    //	pAnalysis->cmv_clustsiz[ij] =CmvHit_pointer->CmvHit_list[ij]->GetHitsize();
    // pAnalysis->cmv_clustmom[ij] =CmvHit_pointer->CmvHit_list[ij]->GetSimMom();
    // pAnalysis->cmv_clustthe[ij] =CmvHit_pointer->CmvHit_list[ij]->GetSimThe();
    // pAnalysis->cmv_clustphi[ij] =CmvHit_pointer->CmvHit_list[ij]->GetSimPhi();

	
    if (ij >=pAnalysis->cmv_nhtmx) break; //redundant
  }
  //			pAnalysis->pEventTree->Fill();



  




  //22032022
  FormCmvCluster();


  //store clust data

  pAnalysis->pRootFile->cd();
  pAnalysis->cmv_nclust = CmvCluster_pointer->CmvCluster_list.size();
  cout<<"  CmvCluster_pointer->CmvCluster_list.size():  "<<  CmvCluster_pointer->CmvCluster_list.size()<<endl;
  cout <<"  pAnalysis->cmv_nclustclust "<<  pAnalysis->cmv_nclust<<endl;

  unsigned int clustId =0;
  if (pAnalysis->cmv_nclust >pAnalysis->cmv_nclustmx) pAnalysis->cmv_nclust =pAnalysis->cmv_nclustmx;
  for (unsigned ij=0; ij<CmvCluster_pointer->CmvCluster_list.size() && ij<pAnalysis->cmv_nclust; ij++) {

    clustId = CmvCluster_pointer->CmvCluster_list[ij]->GetPlane();
    clustId<<=2;
    clustId+=CmvCluster_pointer->CmvCluster_list[ij]->GetLayer();
    clustId<<=7;
    clustId+=CmvCluster_pointer->CmvCluster_list[ij]->GetStrip();
    clustId<<=2; //just shifted by 2 as we did this while forming strip	
    //pAnalysis->cmv_clustpdgid[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetpdgId();

    cout<<"clustid "<<clustId<<" "<<endl;


	
    pAnalysis->cmv_clustid[ij] =clustId;

    cout<< pAnalysis->cmv_clustid[ij]<<endl;

    // pAnalysis->cmv_clustLeTim[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetLeTime();
    // pAnalysis->cmv_clustRiTim[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetRiTime();

    // pAnalysis->cmv_clustLePul[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetLePulse();
    // pAnalysis->cmv_clustRiPul[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetRiPulse();

	
    pAnalysis->cmv_clustTrueposx[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetTruePosX();
    pAnalysis->cmv_clustTrueposy[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetTruePosY();
    pAnalysis->cmv_clustTrueposz[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetTruePosZ();

    pAnalysis->cmv_clustRecoposx[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetRecoPosX();
    pAnalysis->cmv_clustRecoposy[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetRecoPosY();
    pAnalysis->cmv_clustRecoposz[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetRecoPosZ();

    pAnalysis->cmv_clustsiz[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetClusterSize();
    // pAnalysis->cmv_clustmom[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetSimMom();
    // pAnalysis->cmv_clustthe[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetSimThe();
    // pAnalysis->cmv_clustphi[ij] =CmvCluster_pointer->CmvCluster_list[ij]->GetSimPhi();

	
    if (ij >=pAnalysis->cmv_nclustmx) break; //redundant
  }
  //			pAnalysis->pEventTree->Fill();











      
  CMVD_Extrapolation();
  cout<<"cmvd extrapolation ended "<<endl;

  //22032022
  //Store Extrapol positions:


  pAnalysis->pRootFile->cd();

  // pAnalysis->cmv_nexphit = CmvLayExtra_pointer->CmvLayExtra_list.size();
  // cout<<"  CmvLayExtra_pointer->CmvLayExtra_list.size():  "<<  CmvLayExtra_pointer->CmvLayExtra_list.size()<<endl;
  // cout <<"  pAnalysis->cmv_nexpthit "<<  pAnalysis->cmv_nexphit<<endl;


  // if (pAnalysis->cmv_nexphit >pAnalysis->cmv_nexphtmx) pAnalysis->cmv_nexphit =pAnalysis->cmv_nexphtmx;
  // for (unsigned ijj=0; ijj<CmvLayExtra_pointer->CmvLayExtra_list.size()/* && ij<pAnalysis->cmv_nclusthit*/; ijj++) {

  //   CmvLayExtra_pointer->CmvLayExtra_list[ijj]->Print();
		
  //   //pAnalysis->cmv_clustpdgid[ijj] =CmvLayExtra_pointer->CmvLayExtra_list[ijj]->GetpdgId();


  //   pAnalysis->cmv_expid[ijj] =CmvLayExtra_pointer->CmvLayExtra_list[ijj]->GetId();

  //   cout<< pAnalysis->cmv_hitid[ijj]<<endl;


	
  //   pAnalysis->cmv_Expposx[ijj] =CmvLayExtra_pointer->CmvLayExtra_list[ijj]->GetExtXPos();
  //   pAnalysis->cmv_Expposy[ijj] =CmvLayExtra_pointer->CmvLayExtra_list[ijj]->GetExtYPos();
  //   pAnalysis->cmv_Expposz[ijj] =CmvLayExtra_pointer->CmvLayExtra_list[ijj]->GetExtZPos();
  //   pAnalysis->distofclosapp[ijj]=CmvLayExtra_pointer->CmvLayExtra_list[ijj]->GetClosDist();
  //   pAnalysis->planeedge[ijj]=CmvLayExtra_pointer->CmvLayExtra_list[ijj]->GetEdge();
  //   pAnalysis->cmv_DCAposx[ijj]=CmvLayExtra_pointer->CmvLayExtra_list[ijj]->GetDCAXPos(); 
  //   pAnalysis->cmv_DCAposy[ijj]=CmvLayExtra_pointer->CmvLayExtra_list[ijj]->GetDCAYPos();
  //   pAnalysis->cmv_DCAposz[ijj]=CmvLayExtra_pointer->CmvLayExtra_list[ijj]->GetDCAZPos();
  //   if (ijj >=pAnalysis->cmv_nexphtmx) break; //redundant
  // }
  







  //22032022








  pAnalysis->pEventTree->Fill(); //VALGRIND
  cout<<"checkq"<<endl;
  //.............................End/.............................


  if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput==3 || pAnalysis->InputOutput==5) {
    
    //  if(magfield==0){
   
    for (unsigned int ij=0; ij< inoTrackCand_pointer->InoTrackCand_list.size(); ij++) {
      if (inoTrackCand_pointer->InoTrackCand_list[ij]) {
  	delete inoTrackCand_pointer->InoTrackCand_list[ij];
  	inoTrackCand_pointer->InoTrackCand_list[ij]=0;
      }
    }
    inoTrackCand_pointer->InoTrackCand_list.clear();
    if (inoTrackCand_pointer) {
      delete inoTrackCand_pointer;
      inoTrackCand_pointer=0;
    }



    //  }//if(magfield==0)


  }//if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput==3 || pAnalysis->InputOutput==5) 
  cout<<"checkw"<<endl;
  //
  


  for (unsigned ij=0; ij<CmvHit_pointer->CmvHit_list.size(); ij++) {
    if (CmvHit_pointer->CmvHit_list[ij]) {
      cout <<"ij "<< ij<<" "<<CmvHit_pointer->CmvHit_list.size()<<endl;
      delete CmvHit_pointer->CmvHit_list[ij]; CmvHit_pointer->CmvHit_list[ij]=0;
    }
  }

  CmvHit_pointer->CmvHit_list.clear();
  if (CmvHit_pointer) {delete CmvHit_pointer; CmvHit_pointer=0;}




  for (unsigned ij=0; ij<CmvCluster_pointer->CmvCluster_list.size(); ij++) {
    if (CmvCluster_pointer->CmvCluster_list[ij]) {
      cout <<"ij "<< ij<<" "<<CmvCluster_pointer->CmvCluster_list.size()<<endl;
      delete CmvCluster_pointer->CmvCluster_list[ij]; CmvCluster_pointer->CmvCluster_list[ij]=0;
    }
  }

  CmvCluster_pointer->CmvCluster_list.clear();
  if (CmvCluster_pointer) {delete CmvCluster_pointer; CmvCluster_pointer=0;}






  //


  // for (unsigned ij=0; ij<CmvLayExtra_pointer->CmvLayExtra_list.size(); ij++) {
  //   if (CmvLayExtra_pointer->CmvLayExtra_list[ij]) {
  //     //  cout <<"ij "<< ij<<" "<<CmvLayExtra_pointer->CmvLayExtra_list.size()<<endl;
  //     delete CmvLayExtra_pointer->CmvLayExtra_list[ij]; CmvLayExtra_pointer->CmvLayExtra_list[ij]=0;
  //   }
  // }

  // CmvLayExtra_pointer->CmvLayExtra_list.clear();
  // if (CmvLayExtra_pointer) {delete CmvLayExtra_pointer; CmvLayExtra_pointer=0;}






  //





  

  for (unsigned ij=0; ij<CmvStrip_pointer->CmvStrip_list.size(); ij++) {
    if (CmvStrip_pointer->CmvStrip_list[ij]) {
      //  cout <<"ij "<< ij<<" "<<CmvStrip_pointer->CmvStrip_list.size()<<endl;
      delete CmvStrip_pointer->CmvStrip_list[ij]; CmvStrip_pointer->CmvStrip_list[ij]=0;
    }
  }

  CmvStrip_pointer->CmvStrip_list.clear();
  if (CmvStrip_pointer) {delete CmvStrip_pointer; CmvStrip_pointer=0;}

  for (unsigned ij=0; ij<SipmHit_pointer->SipmHit_list.size(); ij++) {
    if (SipmHit_pointer->SipmHit_list[ij]) {
      //  cout <<"ij "<< ij<<" "<<SipmHit_pointer->SipmHit_list.size()<<endl;
      delete SipmHit_pointer->SipmHit_list[ij]; SipmHit_pointer->SipmHit_list[ij]=0;
    }
  }

  SipmHit_pointer->SipmHit_list.clear();
  if (SipmHit_pointer) {delete SipmHit_pointer; SipmHit_pointer=0;}


 
  for (unsigned ij=0; ij<inoStripX_pointer->InoStripX_list.size(); ij++) {
    if (inoStripX_pointer->InoStripX_list[ij]) {
      cout <<"ijkl "<< ij<<" "<<inoStripX_pointer->InoStripX_list.size()<<endl;
      delete inoStripX_pointer->InoStripX_list[ij]; inoStripX_pointer->InoStripX_list[ij]=0;
    }
  }

  inoStripX_pointer->InoStripX_list.clear();
  if (inoStripX_pointer) {delete inoStripX_pointer; inoStripX_pointer=0;}
  
  for (unsigned ij=0; ij<inoStripY_pointer->InoStripY_list.size(); ij++) {
    if (inoStripY_pointer->InoStripY_list[ij]) {
      delete inoStripY_pointer->InoStripY_list[ij]; inoStripY_pointer->InoStripY_list[ij]=0;
    }
  }
  inoStripY_pointer->InoStripY_list.clear();
  if (inoStripY_pointer) {delete inoStripY_pointer; inoStripY_pointer=0;}
  // these pointers u delete only at end of run . it is deleted in destructor cal0sd

  cout<<"checkw"<<endl;


  if (inoRPC_pointer) { delete inoRPC_pointer; inoRPC_pointer=0;}



 for (unsigned int ij=0; ij<7; ++ij) {
   for (unsigned int jk=0; jk<4; ++jk) {
    CmvClusterBank[ij][jk].clear();
   }
 }






  cout<<"hey"<<endl;
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
