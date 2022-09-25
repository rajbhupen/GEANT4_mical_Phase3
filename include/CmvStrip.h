#ifndef CMVSTRIP_H
#define CMVSTRIP_H
//CandStripHandle
#include <iostream>
#include "MultiSimAnalysis.hh"
using namespace std;

class CmvStrip
{

 public:
  CmvStrip();
  CmvStrip(CmvStrip* cd);
  ~CmvStrip();
  CmvStrip *DupHandle() const;
  void Print();
  void Trace(const char *c = "") const;

  int    GetpdgId() const { return fpdgStrip;}	
  int    GetId() const { return fStrip;}
  int    GetPlane() const { return (fStrip>>11)&0x07;};
  int    GetLayer() const { return (fStrip>>9)&0x03;};
  int    GetStrip() const { return (fStrip>>2)&0x7f;};
	
  double GetXPos() const {return fXPos;};
  double GetYPos() const {return fYPos;};
  double GetZPos() const {return fZPos;};
	
  double GetXLocPos() const {return fXLocPos;};
  double GetYLocPos() const {return fYLocPos;};
  double GetZLocPos() const {return fZLocPos;};
  
  double GetTime() const {return  fTime;};
  double GetPulse() const {return fPulse;}; 

  double GetSimMom() const {return fSimMom;};
  double GetSimThe() const {return fSimThe;};
  double GetSimPhi() const {return fSimPhi;};


  void SetpdgId(int ipdg) { fpdgStrip=ipdg;}
  void SetId(int id) { fStrip=id;}
	
  void SetXPos(double fd) {fXPos=fd;}
  void SetYPos(double fd) {fYPos=fd;}
  void SetZPos(double fd) {fZPos=fd;}
	
  void SetXLocPos(double fd) {fXLocPos=fd;}
  void SetYLocPos(double fd) {fYLocPos=fd;}
  void SetZLocPos(double fd) {fZLocPos=fd;}
  
  void SetTime(double fd) {fTime=fd;}
  void SetPulse(double fd) {fPulse=fd;}; 

  void SetSimMom(double fd) {fSimMom=fd;};
  void SetSimThe(double fd) {fSimThe=fd;};
  void SetSimPhi(double fd) {fSimPhi=fd;};

  int    fStrip;   //StripID
  int    fpdgStrip;
  double fXPos;
  double fYPos;
  double fZPos;
  double fXLocPos;
  double fYLocPos;
  double fZLocPos;
  double fTime;
  double fPulse;

  double fSimMom;
  double fSimThe;
  double fSimPhi;
};
#endif                                              // INOSTRIP_H
