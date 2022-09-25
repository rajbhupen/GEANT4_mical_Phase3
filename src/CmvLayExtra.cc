#include <cassert>
#include <iostream>
#include "CmvLayExtra.h"
#include <cmath>
#include "TMath.h"
#include <iostream>
using namespace std;
//______________________________________________________________________

//______________________________________________________________________
CmvLayExtra::CmvLayExtra() {
        
  fExtXPos=-1000000.0;
  fExtYPos=-1000000.0;
  fExtZPos=-1000000.0;

  fDCAXPos=-1000000.0;
  fDCAYPos=-1000000.0;
  fDCAZPos=-1000000.0;
	
  fClosDist = -1000000.0;
  isUsed = 0;
  fEdge=-1000000;
}

//______________________________________________________________________
CmvLayExtra::CmvLayExtra(CmvLayExtra* cd) {
  flay    = cd->flay;    

  fExtXPos			= cd->fExtXPos;     
  fExtYPos			= cd->fExtYPos;     
  fExtZPos			= cd->fExtZPos;

  fDCAXPos			= cd->fDCAXPos;     
  fDCAYPos			= cd->fDCAYPos;     
  fDCAZPos			= cd->fDCAZPos;

  fClosDist = cd->fClosDist;
  isUsed = cd->isUsed;
  fEdge = cd->fEdge;
}

//______________________________________________________________________
CmvLayExtra::~CmvLayExtra() {
}

//______________________________________________________________________
CmvLayExtra *CmvLayExtra::DupHandle() const {
  return (new CmvLayExtra(*this));
}


//______________________________________________________________________
void CmvLayExtra::Trace(const char *c) const {
  std::cout<<"**********Begin CmvLayExtra::Trace(\"" << c << "\")" << std::endl
	   << "Information from CmvLayExtra's CandHandle: " << std::endl;
  //  CandHandle::Trace(c);
  std::cout
    << "**********End CmvLayExtra::Trace(\"" << c << "\")" << std::endl;
}

void CmvLayExtra::Print() {
  cout<< "CmvLayExtra():" 
    // <<std::setw(4) <<jk <<" "
      << " Detid "<< GetId()
      << " Plane "<< std::setw(2)<<  GetPlane()
      << " Layer "<<std::setw(2)<<  GetLayer()
      << " Ext X_Pos=" <<std::setw(8)<<  GetExtXPos()
      << " Ext Y_Pos=" <<std::setw(8)<<  GetExtYPos()
      << " Ext Z_Pos=" <<std::setw(8)<<  GetExtZPos()
      << " Dist of Closest Approach=" <<std::setw(8)<<  GetClosDist()
      << " min Edge "<< GetEdge() 
      << " DCA X_Pos=" <<std::setw(8)<<  GetDCAXPos()
      << " DCA Y_Pos=" <<std::setw(8)<<  GetDCAYPos()
      << " DCA Z_Pos=" <<std::setw(8)<<  GetDCAZPos()

      << endl;
  //  cout<<"......................................................................"<<endl;	

}
