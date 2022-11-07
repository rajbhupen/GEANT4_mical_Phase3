#ifndef SIPMHIT_H
#define SIPMHIT_H
//CandStripHandle
#include "CmvStrip.h"
//#include "TRandom3.h"
#include "TRandom.h"
class SipmHit
{

 public:
  SipmHit();
  SipmHit(SipmHit* cd);
  SipmHit(CmvStrip* str, int Sipm); // SiPM : 0:L-x, 1:L+x, 2:R-x, 3:R+x
  ~SipmHit();
  SipmHit *DupHandle() const;
  void Print();
  void Trace(const char *c = "") const;

  //  CmvStrip* getCmvStrip() { return fCmvStrip;}
  int    GetpdgId() const { return fpdgSipm;}	
  int    GetId() const { return fSipmId;}
  int    GetStripId() const { return (fSipmId>>2)&0xfff;};
  int    GetPlane() const { return (fSipmId>>11)&0x07;};
  int    GetLayer() const { return (fSipmId>>9)&0x03;};
  int    GetStrip() const { return (fSipmId>>2)&0x7f;};
  int    GetSiPM() const { return  fSipmId&0x03;};
	
  double GetXPos() const {return fXPos;};
  double GetYPos() const {return fYPos;};
  double GetZPos() const {return fZPos;};
	
  double GetXLocPos() const {return fXLocPos;};
  double GetYLocPos() const {return fYLocPos;};
  double GetZLocPos() const {return fZLocPos;};
  int GetTimePulse() const {return  iTimePulse;};
  int GetTime() const {return  iTimePulse&0x0fffff;};
  int GetPulse() const {
    //	  cout<<"check: "<<iTimePulse<<"hey"<<endl;
    return (iTimePulse>>20)&0x0fff;};
   
	
	
  double GetSimMom() const {return fSimMom;};
  double GetSimThe() const {return fSimThe;};
  double GetSimPhi() const {return fSimPhi;};
	
  bool GetUsed() const {return isUsed;};
	
  void SetpdgId(int ipdg) { fpdgSipm=ipdg;}
  void SetId(int id) { fSipmId=id;}
	
  void SetXPos(double fd) {fXPos=fd;}
  void SetYPos(double fd) {fYPos=fd;}
  void SetZPos(double fd) {fZPos=fd;}
	
  void SetXLocPos(double fd) {fXLocPos=fd;}
  void SetYLocPos(double fd) {fYLocPos=fd;}

  void SetZLocPos(double fd) {fZLocPos=fd;}
	
  void SetTimePulse(int fd) {iTimePulse = fd;};	
  //	void SetTime(int fd)  { iTimePulse+=fd;}
  //	void SetPulse(int fd) {iTimePulse += (fd<<20);}
	
  void SetSimMom(double fd) {fSimMom=fd;};
  void SetSimThe(double fd) {fSimThe=fd;};
  void SetSimPhi(double fd) {fSimPhi=fd;};

  void Update(int edep, int time);
  void SetUsed(bool fd) {isUsed=fd;};
	
  //	CmvStrip* fCmvStrip;
  int    fSipmId;   //StripID
  int    fpdgSipm;
  double fXPos;
  double fYPos;
  double fZPos;
  double fXLocPos;
  double fYLocPos;
  double fZLocPos;
  int iTimePulse;
 	
  double fSimMom;
  double fSimThe;
  double fSimPhi;

  bool isUsed;
  //	TRandom3* gRandom3;

  bool debug = true;
};


#endif                                              // INOSTRIP_H
