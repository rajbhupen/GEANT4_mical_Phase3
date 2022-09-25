#ifndef CMVHIT_H
#define CMVHIT_H
//HitCam
//GMA use proper fView value to calcualte X/Y-direction(length) peoperly

//#include "MessageService/MsgService.h"
#include "SipmHit.h"
#include "micalDetectorParameterDef.hh"
#include "MultiSimAnalysis.hh"

class CmvHit
{

 public:
  CmvHit();
  CmvHit(SipmHit* L0,  SipmHit* L1, SipmHit* R0, SipmHit* R1, double* xx);
  CmvHit(CmvHit* hit);
  ~CmvHit();

  SipmHit* GetL0SiPM() const {return fLeMiSipm;};
  SipmHit* GetL1SiPM() const {return fLePlSipm;};
  SipmHit* GetR0SiPM() const {return fRiMiSipm;};
  SipmHit* GetR1SiPM() const {return fRiPlSipm;};
	
  int    GetStripId() const { return fStripId;}
  int    GetPlane() const { return (fStripId>>9)&0x07;}; 
  int    GetLayer() const { return (fStripId>>7)&0x03;}; 
  int    GetStrip() const { return fStripId&0x7f;};
	
  void Print();
  bool GetUsed() const {return isUsed;};	
  double GetPulse() const {
    double pulse =0;
    int count =0;
    if (fLeMiSipm) {count++; pulse+= fLeMiSipm->GetPulse();}
    if (fLePlSipm) {count++; pulse+= fLePlSipm->GetPulse();}
    if (fRiMiSipm) {count++; pulse+= fRiMiSipm->GetPulse();}
    if (fRiPlSipm) {count++; pulse+= fRiPlSipm->GetPulse();}
    if (count ==0) count = 1;
    return pulse;
    //return pulse/count; //why average ?
  };

  void SetTruePosX(double fd) { fXTruePosX=fd;}
  void SetTruePosY(double fd) { fYTruePosY=fd;}
  void SetTruePosZ(double fd) { fZTruePosZ=fd;}

  void SetRecoPosX(double fd) { fXRecoPosX=fd;}
  void SetRecoPosY(double fd) { fYRecoPosY=fd;}
  void SetRecoPosZ(double fd) { fZRecoPosZ=fd;}

  void SetPosXErr(double fd) { fXPosErr=fd;}
  void SetPosYErr(double fd) { fYPosErr=fd;}
  void SetPosZErr(double fd) { fZPosErr=fd;}
	
  void SetLeTime(double fd) { fLeTime=fd;}	
  void SetRiTime(double fd) { fRiTime=fd;}	

  void SetLePulse(double fd) { fLePulse=fd;}	
  void SetRiPulse(double fd) { fRiPulse=fd;}	
	
  void SetGenMom(double fd) { fMomentum=fd;}
  void SetGenThe(double fd) { fTheta=fd;}
  void SetGenPhi(double fd) { fPhi=fd;}	
	
  double GetTruePosX() { return fXTruePosX;}
  double GetTruePosY() { return fYTruePosY;}
  double GetTruePosZ() { return fZTruePosZ;}

  double GetRecoPosX() { return fXRecoPosX;}
  double GetRecoPosY() { return fYRecoPosY;}
  double GetRecoPosZ() { return fZRecoPosZ;}

  double GetPosXErr() { return fXPosErr;}
  double GetPosYErr() { return fYPosErr;}
  double GetPosZErr() { return fZPosErr;}
	
  double GetLeTime() { return fLeTime;}	
  double GetRiTime() { return fRiTime;}	

  double GetLePulse() { return fLePulse;}	
  double GetRiPulse() { return fRiPulse;}	

	
  double GetGenMom() { return fMomentum;}
  double GetGenThe() { return fTheta;}
  double GetGenPhi() { return fPhi;}	

  bool isIdentical(CmvHit* hit);
  void SetUsed(bool fd) {isUsed=fd;};
  //private:
  micalDetectorParameterDef* paradef;
  MultiSimAnalysis* pAnalysis;

  SipmHit* fLeMiSipm;   //Strip for L-Y-SiPM 
  SipmHit* fLePlSipm;   //Strip for L+Y-SiPM
  SipmHit* fRiMiSipm;   //Strip for R-Y-SiPM 
  SipmHit* fRiPlSipm;   //Strip for R+Y-SiPM

  int fStripId; // All four SiPM should have come fro mthe same StripID

  double fXTruePosX;      //Global X-position from Simulation
  double fYTruePosY;      //Global Y-position from Simulation
  double fZTruePosZ;      //Global Z-position from Simulation

  double fXRecoPosX;      //Global X-position from Reconstruction
  double fYRecoPosY;      //Global Y-position from Reconstruction
  double fZRecoPosZ;      //Global Z-position from Reconstruction

  double fLeTime;      // Timing of Left-SiPMs (ns)
  double fRiTime;      // Timing of Right-SiPMs (ns)

  double fLePulse;      // Timing of Left-SiPMs (ns)
  double fRiPulse;      // Timing of Right-SiPMs (ns)
	
  double fXPosErr; //Error in X-position (m)
  double fYPosErr; //Error in Y-position (m)
  double fZPosErr; //Error in Y-position (m)

  double fMomentum; //Momentum of track which is behind this hit (use only one SiPM)
  double fTheta;  //Theta
  double fPhi;    // Phi

  double PhotonSpeed;
  double layergap;
  bool isUsed;
  bool debug=false;	
};

#endif
