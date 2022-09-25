#include <cassert>
#include <iostream>
#include "CmvStrip.h"
#include <cmath>
#include "TMath.h"
#include <iostream>
using namespace std;
//______________________________________________________________________

//______________________________________________________________________
CmvStrip::CmvStrip() {
  fpdgStrip = -100;
  fXPos=-1000.0;
  fYPos=-1000.0;
  fZPos=-1000.0;
	
  fXLocPos=-1000.0;
  fYLocPos=-1000.0;
  fZLocPos=-1000.0;
  
  fTime=-1000.0;
  fPulse=-1000.0;

  fSimMom=-1000.0;
  fSimThe=-1000.0;
  fSimPhi=-1000.0;
}

//______________________________________________________________________
CmvStrip::CmvStrip(CmvStrip* cd) {
  fStrip    = cd->fStrip;    
  fpdgStrip	= cd->fpdgStrip; 
  fXPos			= cd->fXPos;     
  fYPos			= cd->fYPos;     
  fZPos			= cd->fZPos;     
  fXLocPos	= cd->fXLocPos;  
  fYLocPos	= cd->fYLocPos;
  fZLocPos	= cd->fZLocPos;
						                 
  fTime			= cd->fTime;     
  fPulse		= cd->fPulse;    
						                 
  fSimMom		= cd->fSimMom;   
  fSimThe		= cd->fSimThe;   
  fSimPhi		= cd->fSimPhi;
}

//______________________________________________________________________
CmvStrip::~CmvStrip() {
}

//______________________________________________________________________
CmvStrip *CmvStrip::DupHandle() const {
  return (new CmvStrip(*this));
}


//______________________________________________________________________
void CmvStrip::Trace(const char *c) const {
  std::cout<<"**********Begin CmvStrip::Trace(\"" << c << "\")" << std::endl
	   << "Information from CmvStrip's CandHandle: " << std::endl;
  //  CandHandle::Trace(c);
  std::cout
    << "**********End CmvStrip::Trace(\"" << c << "\")" << std::endl;
}

void CmvStrip::Print() {
  cout<< "CmvStrip():" 
    // <<std::setw(4) <<jk <<" "
      << " Detid "<< GetId()
      << " Plane "<< std::setw(2)<<  GetPlane()
      << " Layer "<<std::setw(2)<<  GetLayer()
      << " Strip "<<std::setw(3)<<  GetStrip()
      << " pdgId "<<std::setw(3)<<  GetpdgId()
      << " X_Loc="   <<std::setw(8)<<  GetXLocPos() 
      << " Y_Loc=" <<std::setw(8)<<  GetYLocPos()
      << " Z_Loc=" <<std::setw(8)<<  GetZLocPos()  
      << " X_Pos=" <<std::setw(8)<<  GetXPos()
      << " Y_Pos=" <<std::setw(8)<<  GetYPos()
      << " Z_Pos=" <<std::setw(8)<<  GetZPos()
      << " chg="   <<std::setw(8)<<  GetPulse()
      << " time="  <<std::setw(8)<<  GetTime()
      << endl;
  //  cout<<"......................................................................"<<endl;	

}
