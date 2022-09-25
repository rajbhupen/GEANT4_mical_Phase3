#ifndef CMVCLUSTER_H
#define CMVCLUSTER_H
//ClusterCam
//GMA use proper fView value to calcualte X/Y-direction(length) peoperly

//#include "MessageService/MsgService.h"
#include "CmvHit.h"
#include "micalDetectorParameterDef.hh"
#include "MultiSimAnalysis.hh"
#include <vector>
using std::vector;

class CmvHit;
class CmvCluster
{

public:

  CmvCluster(CmvHit* hit);
  virtual ~CmvCluster();
   void Print();
  void AddHits(CmvHit* hit);
  bool ContainsHit(CmvHit* hit);

   CmvHit* GetHit(unsigned int i) const;
   

  int    GetPlane() const { return fPlane;}; 
  int    GetLayer() const { return fLayer;}; 
  int    GetStrip() const { return fStrip;};
  int    GetBegStrip() const { return fBegStrip;};
  int    GetEndStrip()const  { return fEndStrip;};
  int    GetClusterSize(){return fClusterSize;}
  bool   GetUsed() const {return isUsed;};		



  double GetTruePosX() const { return fTruePosX;}
  double GetTruePosY() const { return fTruePosY;}
  double GetTruePosZ() const { return fTruePosZ;}



  double GetBegTruePosX() const  { return fBegTruePosX;}
  double GetBegTruePosY() const { return fBegTruePosY;}
  double GetBegTruePosZ() const { return fBegTruePosZ;}
  
  double GetEndTruePosX() const { return fEndTruePosX;}
  double GetEndTruePosY() const { return fEndTruePosY;}
  double GetEndTruePosZ() const { return fEndTruePosZ;}



  double GetRecoPosX() const { return fRecoPosX;}
  double GetRecoPosY() const { return fRecoPosY;}
  double GetRecoPosZ() const { return fRecoPosZ;}



  double GetBegRecoPosX()const  { return fBegRecoPosX;}
  double GetBegRecoPosY()const  { return fBegRecoPosY;}
  double GetBegRecoPosZ() const { return fBegRecoPosZ;}
  
  double GetEndRecoPosX() const { return fEndRecoPosX;}
  double GetEndRecoPosY() const { return fEndRecoPosY;}
  double GetEndRecoPosZ() const { return fEndRecoPosZ;}



  double GetWRecoPosX() const { return fWRecoPosX;}
  double GetWRecoPosY() const { return fWRecoPosY;}
  double GetWRecoPosZ() const { return fWRecoPosZ;}





  
  double GetPosXErr() const { return fXPosErr;}
  double GetPosYErr() const { return fYPosErr;}
  double GetPosZErr() const { return fZPosErr;}
	
  /* double GetLeTime() { return fLeTime;}	 */
  /* double GetRiTime() { return fRiTime;}	 */

  /* double GetLePulse() { return fLePulse;}	 */
  /* double GetRiPulse() { return fRiPulse;}	 */
  double GetPulse() const {return fPulse;}

  double GetBegPulse() const {return fBegPulse;}
  double GetEndPulse() const {return fEndPulse;}

  
  double GetGenMom() const { return fMomentum;}
  double GetGenThe() const { return fTheta;}
  double GetGenPhi() const { return fPhi;}	


  void SetUsed(bool fd) {isUsed=fd;};


  vector <CmvHit*> HitsInCmvCluster;
private:
  micalDetectorParameterDef* paradef;
  MultiSimAnalysis* pAnalysis;

  int fPlane;
  int fLayer;
  int fStrip;  int fBegStrip;  int fEndStrip;  

  
  int fClusterSize;
  bool isUsed;
  
  //Global X-position from Simulation
  double fTruePosX;   double fBegTruePosX;   double fEndTruePosX;
  
  //Global Y-position from Simulation
  double fTruePosY;   double fBegTruePosY;   double fEndTruePosY;

  //Global Z-position from Simulation
  double fTruePosZ;   double fBegTruePosZ;   double fEndTruePosZ;


  //Global X-position from Reconstruction
  double fRecoPosX;   double fBegRecoPosX;   double fEndRecoPosX;
  
  //Global Y-position from Reconstruction
  double fRecoPosY;   double fBegRecoPosY;   double fEndRecoPosY;

  //Global Z-position from Reconstruction
  double fRecoPosZ;   double fBegRecoPosZ;   double fEndRecoPosZ;


    //energy weighted reco X-position
  double fWRecoPosX;
  
  //energy weighted reco Y-position
  double fWRecoPosY;  

  //Energy weighted reco  Z-position 
  double fWRecoPosZ;   


  

  // Timing 
  double fLeTime;    double fBegLeTime;  double fEndLeTime;
  
  
 
  double fPulse;  double fBegPulse;  double fEndPulse;
  
 
  
  
  double fXPosErr; //Error in X-position (m)
  double fYPosErr; //Error in Y-position (m)
  double fZPosErr; //Error in Y-position (m)

  double fMomentum; //Momentum of track which is behind this cluster (use only one SiPM)
  double fTheta;  //Theta
  double fPhi;    // Phi

  double PhotonSpeed;
  double layergap;


};

#endif
