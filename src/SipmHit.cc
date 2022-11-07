//Assume: nearer sipm has 33 npho for a muon passage.
//1pe = 0.27 pC @2.5 ov
//peaks in Led spectrum were fitted with gaussian whose sigma sq is s1^2+npho*s2^2 (s1 = 0.015pC,s2=0.006pC )
#include <cassert>
#include <iostream>
#include "SipmHit.h"

using namespace std;
//______________________________________________________________________

//______________________________________________________________________
SipmHit::SipmHit() {
  //	fCmvStrip = 0;
  fpdgSipm = -100;
  fXPos=-1000.0;
  fYPos=-1000.0;
  fZPos=-1000.0;
  fXLocPos=-1000.0;
  fYLocPos=-1000.0;
  fZLocPos=-1000.0;
  iTimePulse=0;
  
  fSimMom=-1000.0;
  fSimThe=-1000.0;
  fSimPhi=-1000.0;
  isUsed = 0;
  //	gRandom3 = new TRandom3();
}

//______________________________________________________________________
SipmHit::SipmHit(SipmHit* cd) {
  //	fCmvStrip = cd->fCmvStrip;
  fpdgSipm	= cd->fpdgSipm; 
  fSipmId    = cd->fSipmId;    
  fXPos			= cd->fXPos;     
  fYPos			= cd->fYPos;     
  fZPos			= cd->fZPos;     
  fXLocPos	= cd->fXLocPos;  
  fYLocPos	= cd->fYLocPos;  
  fZLocPos	= cd->fZLocPos;  						                 
  iTimePulse	= cd->iTimePulse;     
  
  fSimMom		= cd->fSimMom;   
  fSimThe		= cd->fSimThe;   
  fSimPhi		= cd->fSimPhi;
  isUsed    = cd->isUsed;
}

SipmHit::SipmHit(CmvStrip* str, int Sipm) {
  //	fCmvStrip = str->fCmvStrip;
	
  fpdgSipm	= str->GetpdgId(); 
  fSipmId = str->GetId() + Sipm;
  fXPos			= str->GetXPos();     
  fYPos			= str->GetYPos();     
  fZPos			= str->GetZPos();     
  fXLocPos	= str->GetXLocPos();  
  fYLocPos	= str->GetYLocPos();  
  fZLocPos	= str->GetZLocPos(); 						                 					                 
  fSimMom		= str->GetSimMom();   
  fSimThe		= str->GetSimThe();   
  fSimPhi		= str->GetSimPhi();

  double time = str->GetTime();
  double pulse = str->GetPulse();
  if(debug) cout<<"fSipmId "<<fSipmId<<endl;

  MultiSimAnalysis* pAnalysis = MultiSimAnalysis::AnPointer;
  micalDetectorParameterDef* paradef = micalDetectorParameterDef::AnPointer;
  

  int  Sipm_Pedestal = paradef->Sipm_Pedestal;
  int iplane = (fSipmId>>11)&0x07;
  
  if(debug)  cout<<"iplane:"<<iplane<<endl;
  double layhalflength = paradef->partopscint[1];//2300;
  double dist;
  //  double  ScntLayShifSide = paradef->GetScntLayShifSide();
  //  double SidePlaneHalfLength = paradef->GetSidePlaneHalfLength();
  if (iplane==4){
    layhalflength = paradef->partopscint[1]+50;//2350 i.e 4.7m;
    //   dist= (Sipm<=1) ? layhalflength+fYLocPos : layhalflength-fYLocPos; //Assume in cm
    if(debug)  cout<<"dist fYLocPos layhalflength "<<iplane<<" "<<dist<<" "<<fYLocPos<<" "<<layhalflength<<endl;
  }
  if (iplane==1){
    layhalflength = paradef->partopscint[1]-50;//2250 i.e 4.5m;
    //   dist= (Sipm<=1) ? layhalflength+fYLocPos : layhalflength-fYLocPos; //Assume in cm
    if(debug) cout<<"dist fYLocPos layhalflength "<<iplane<<" "<<dist<<" "<<fYLocPos<<" "<<layhalflength<<endl;
  }

  if (iplane==6 || iplane==7){
    layhalflength = 1000;//
  } 
  //  if(debug) cout<<"dist fYLocPos layhalflength "<<iplane<<" "<<dist<<" "<<fYLocPos<<" "<<layhalflength<<endl;
  

      
  dist= (Sipm<=1) ? layhalflength+fYLocPos : layhalflength-fYLocPos; //Assume in cm
    
  if(debug) cout<<"dist fXLocPos layhalflength "<<iplane<<" "<<dist<<" "<<fYLocPos<<" "<<layhalflength<<endl;
     
  if(debug) cout <<"SipmHit::SipmHit time edep dist "<<iplane<<"  " <<time<<" "<<pulse<<" "<<dist<<endl; //pulse is just edep in keV
  
  //	if (Sipm>1 ) { dist = 457.8 - dist;} // may be different for different CMV, change accordingly
  //Double Exponential loss of signal in fibre
  double lam1 = 883.3;
  double lam2 = 10121.0;
  double frac = 0.96;
  
  double signal =pulse*(lam1*exp(-dist/lam1) + frac*lam2*exp(-dist/lam2)) /(lam1+lam2*frac); // in KeV
  
  //Average number of ohotons will be vary for different scintillator, SiPM
  double lambda = signal* 0.0165; // Assume 1cm scint (2MeV), observed signal in 33 photon in nearer SiPM//33/2000
  //lambda *=0.1; // GMA Manually reduced it (Too much energy deposit in the scintillator veto detector : Probably multiple insertion of the chamber
  // lambda *=1;//material changed from lead to hydrocarbon
  //Poission fluctuations of observed SiPM
  int npho = gRandom->Poisson(lambda);
  
  //Charge in pC
  double picoC = npho*0.27;//0.24 
  //Smeaing of charge
  if(debug)	cout <<"picoC "<<picoC<<endl;
  double sigma1sq = pow(0.015,2);
  double sigma2sq = pow(0.006,2);
  double gausigmasq = sigma1sq + (npho*sigma2sq);
  picoC += gRandom->Gaus(0,pow(gausigmasq,0.5)); // Random noise in the signal
  //sqrt(a1^2+a2^2*npho)
  if (picoC<0) { picoC = 0;}
  //Digitisation of signal, assuming least count 0.01pC
  int iPul = picoC*100+ Sipm_Pedestal;//100= pedestal
  // int iPul = pulse;
  //Assuming only 12bit adc
    
  if (iPul > 4095) iPul = 4095;
  
  if(debug) cout<<"pulse: "<<pulse<<" signal: "<<signal<<" lambda: "<<lambda<<" npho: "<<npho<<" picoC: "<<picoC<<" ipul: "<<iPul<<" "<<endl;
 
  //Time at SiPM with the assumption that velocity inside fibre is 16.2cm/ns
  int sipmtime = time + (dist/pAnalysis->GetPhotonSpeedVal());
  if(debug) cout <<"SipmHit::SipmHit "<<iPul<<" "<<picoC<<" n"<<npho<<" "<<lambda<<" "<<signal<<" p"<<pulse<<" "<<time<<" "<<sipmtime<<endl;
  //Time resolution, 3ns for to photon and 1ns for 50 photon
  double timereso =0.707* exp(1.151-0.0031*npho);
  //   double  timereso =1;
 
  if(debug) cout<<"timeresol  "<<timereso<<endl;
  sipmtime += gRandom->Gaus(0, timereso);
  //// cout<<"sipmtime "<<sipmtime<<endl;
  //TDC least count 100ps
  int iTime = sipmtime/pAnalysis->GetCMVadctons();
  //       int iTime = sipmtime;
  if(debug) cout<<"itime "<<iTime<<endl;
  if (iTime <0) { iTime = 0;}
  if (iTime >1048575) { iTime =1048575;} // 20bits for TDC 2^20-1 ==1048575
       

  iTimePulse	= iPul;
  if(debug)  cout<<"iTimePulse "<<iTimePulse<<" "<<iPul<<" "<<iTime<<endl;
  
  iTimePulse<<=20;
  if(debug) cout<<"iTimePulse "<<iTimePulse<<endl;
  iTimePulse +=iTime;
 if(debug) cout<<"iTimePulse "<<iTimePulse<<endl;
  


  
  isUsed = false;

  if(debug) cout <<"iTimePulse iTime iPul  "<<iTimePulse<<"  " <<iTime<<" "<<iPul<<" "<<endl; //pulse is just edep in keV
  if(debug) cout<<"iTimePulse to time and pulse: "<< iTimePulse<<" "<<((iTimePulse)&0x0fffff) <<" "<< ((iTimePulse>>20)&0x0fff) <<endl;

      
}

void SipmHit::Update(int edep, int time) {
  int oldenr = (iTimePulse>>20)&0x0fff;
  int oldtime = iTimePulse&0xfffff;
  micalDetectorParameterDef* paradef = micalDetectorParameterDef::AnPointer;
  int  Sipm_Pedestal = paradef->Sipm_Pedestal;
  
  int enr = oldenr+edep-Sipm_Pedestal;//both signal and noise contains pedestal shift 
  if(debug) cout<<"enr: "<<enr<<" "<<oldenr<<" "<<edep<<endl;
  if(debug) cout<<"time:"<<time<<" "<<oldtime<<endl;
  if (enr>4097) enr = 4097;
  if (time >oldtime) time=oldtime;
  SetTimePulse((enr<<20)+time);
}

//______________________________________________________________________
SipmHit::~SipmHit() {
  //	delete gRandom3;
}

//______________________________________________________________________
SipmHit *SipmHit::DupHandle() const {
  return (new SipmHit(*this));
}


//______________________________________________________________________
void SipmHit::Trace(const char *c) const {
  std::cout<<"**********Begin SipmHit::Trace(\"" << c << "\")" << std::endl
	   << "Information from SipmHit's CandHandle: " << std::endl;
  //  CandHandle::Trace(c);
  std::cout
    << "**********End SipmHit::Trace(\"" << c << "\")" << std::endl;
}

void SipmHit::Print() {
  cout<< "SipmHit():" 
      << " Plane "<< std::setw(3)<<  GetPlane()
      << " Layer "<< std::setw(3)<<  GetLayer()
      << " Strip No "<< std::setw(3)<<  GetStrip()
      << " StripId  "<< std::setw(3)<<  GetStripId()
      << " Id "<< std::setw(3)<<  GetId()
      << " SiPM No "<< std::setw(3)<<  GetSiPM()
      << " pdgId No "<< std::setw(3)<<  GetpdgId()
      << " X_Loc=" <<std::setw(8)<<  GetXLocPos() 
      << " Y_Loc=" <<std::setw(8)<<  GetYLocPos()
      << " Z_Loc=" <<std::setw(8)<<  GetZLocPos()	  
      << " X_Pos=" <<std::setw(8)<<  GetXPos()
      << " Y_Pos=" <<std::setw(8)<<  GetYPos()
      << " Z_Pos=" <<std::setw(8)<<  GetZPos()
      << " chg="   <<std::setw(8)<<  GetPulse()
      << " time="  <<std::setw(8)<<  GetTime()
      << " timepulse="  <<std::setw(8)<<  GetTimePulse()
    
      << endl;
  //  cout<<"......................................................................"<<endl;
}
