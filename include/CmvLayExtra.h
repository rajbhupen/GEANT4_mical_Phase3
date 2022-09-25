#ifndef CMVLAYEXTRA_H
#define CMVLAYEXTRA_H
//CandStripHandle
#include <iostream>
#include "MultiSimAnalysis.hh"
using namespace std;

class CmvLayExtra
{

public:
  CmvLayExtra();
  CmvLayExtra(CmvLayExtra* cd);
  ~CmvLayExtra();
  CmvLayExtra *DupHandle() const;
  void Print();
  void Trace(const char *c = "") const;
	
  int    GetId() const { return flay;}
  int    GetPlane() const { return (flay>>2)&0x07;};
  int    GetLayer() const { return flay&0x03;};
  
  bool GetUsed() const {return isUsed;};	
  double GetExtXPos() const {return fExtXPos;};
  double GetExtYPos() const {return fExtYPos;};
  double GetExtZPos() const {return fExtZPos;};
  double GetDCAXPos() const {return fDCAXPos;};
  double GetDCAYPos() const {return fDCAYPos;};
  double GetDCAZPos() const {return fDCAZPos;};


  double GetClosDist() const {return fClosDist;};
  double GetEdge() const {return fEdge;};
  void SetId(int id) { flay=id;}
  void SetUsed(bool fd) {isUsed=fd;};
  void SetExtXPos(double fd) {fExtXPos=fd;}
  void SetExtYPos(double fd) {fExtYPos=fd;}
  void SetExtZPos(double fd) {fExtZPos=fd;}
  void SetDCAXPos(double fd) {fDCAXPos=fd;}
  void SetDCAYPos(double fd) {fDCAYPos=fd;}
  void SetDCAZPos(double fd) {fDCAZPos=fd;}


  void SetClosDist (double fd){fClosDist = fd;};
 void SetEdge(int fd){fEdge = fd;};


  
  int    flay ;   //StripID
 
  double fExtXPos;
  double fExtYPos;
  double fExtZPos;

  double fDCAXPos;
  double fDCAYPos;
  double fDCAZPos;

  double fClosDist;
  bool isUsed;
  int fEdge;
 
};
#endif                                              // INOSTRIP_H
