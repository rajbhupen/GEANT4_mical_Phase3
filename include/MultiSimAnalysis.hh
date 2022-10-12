////////////////////////////////////////////////////////////////////////////////
//
//  Event class for event information and store in root file as 
//  tuple and/or histograms
//
////////////////////////////////////////////////////////////////////////////////


#ifndef MULTISIM_H
#define MULTISIM_H
#include <vector>
using std::vector;
#include "micalDetectorParameterDef.hh"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TTree.h"
#include "TFile.h"
#include "globals.hh"
#include "Hits.h"
#include "HitPos.h"
#include "TProfile.h"
#include "G4ThreeVector.hh"

//#include <iostream.h>
//#include <fstream.h>

#include <iostream>
#include <fstream>
using namespace std;

struct vectGr{
  float x;
  float y;
  float z;
  float dx;
  float dy;
  float dz;
};

class MultiSimAnalysis
{
  public:
    MultiSimAnalysis(G4String sprefix);
    void OpenRootfiles(G4String inf, G4String out, G4String colfile, G4String corrfile);

  void CloseRootfiles();
   ~MultiSimAnalysis();
 void SetCorrTimeError(G4double val);
  void SetUnCorrTimeError(G4double val);
  void SetTimeToDigiConvVal(G4double val);
  void SetSignalSpeedVal(G4double val);
	void SetPhotonSpeedVal(G4double val);
	void SetCMVadctons(G4double val);
	
  double GetCorrTimeError() {return CorrTimeError;}
  double GetUnCorrTimeError() {return UnCorrTimeError;}
  double GetTimeToDigiConvVal() {return TimeToDigiConv;}
  double GetSignalSpeedVal() {return SignalSpeed;}
	double GetPhotonSpeedVal() {return PhotonSpeed;}
	double GetCMVadctons() {return CMVadctons;}

	
  public:
    static MultiSimAnalysis* AnPointer;
    ofstream ascii_output;
  // ofstream timeAsciiOutput;
    ofstream B_ascii_output;
    G4int       InputOutput;
    G4String    text_inputFile; // input file for ascii output or .inh file
  G4int collatedIn;
    G4int corr_pos_timeIn;
    TFile *pVisFile;
    TFile *pRootFile;
    TFile *inputRootFile;  
  TFile* collatedRootFile;


  //Raj: Added to obtain multiplicity dependent error 
  TFile* InCorrFile;
  TH2D*	h_xposerrsq;
  TH2D*	h_yposerrsq;
  double xposerrsq[4][12];
  double yposerrsq[4][12];

  
    TTree *pEventTree;
    TTree *inputEventTree;

    TTree *rtree; 

    TH1F  *pDepthDose;
    TH2F  *pRadialMomentum;
    TH3F  *pEnergyDeposit;

  TH1F * hist11;
  TH1F * hist22;
  TH1F * hist33;
  TH1F * hist55;
  TH1F * hist44;

   TH1F * hist_ellip0;
 TH1F * hist_ellip1;
   TH1F * hist_ellip2;
   TH1F * hist_ellip3;
  
    TH1F  *pPosX;
    TH1F  *pPosY;  
    TH1F  *pPosZ;

    TH2F  *pPosXX;
    TH2F  *pPosYY;  
    TH2F  *pPosZZ;

    TH1F  *pdedz[20];  
    TH1F  *hitDist;    //asm
    TH1F  *TrkDist;    //asm
    TH1F  *EffDist;    //asm
    TH1F  *InoTrack_listsize;   //asm
   
  TH1F *strpXtime;
  TH1F *strpYtime;
  TH1F *strpXtimeCorr;
  TH1F *strpYtimeCorr;
  TH1F *hitXtime;
  TH1F *hitYtime;
  TH1D* smagFieldX;
  TH1D* smagFieldY;
  TH2D* smag2dX;
  TH2D* smag2dY;
  TH1D* rmagFieldX;
  TH1D* rmagFieldY;
  TH2D* rmag2dX;
  TH2D* rmag2dY;
  TH2D* smag2dXYpixel_iron;
  TH2D* smag2dXYpixel_air;
  TH2D* rmag2dXYpixel_iron;
  TH2D* rmag2dXYpixel_air;
  TH2D* xyvsbxin;
  TH2D* xyvsbyin;
  TH2D* xyvsbxdiff;
  TH2D* xyvsbydiff;
  TH2D* xyvsbxindiff;
  TH2D* xyvsbyindiff;
  TH2D* xyvsbxout;
  TH2D* xyvsbyout;

  // Collated Histograms
  TH2D* inefficiency_corx[20];
  TH2D* inefficiency_uncx[20];
  TH2D* inefficiency_uncy[20];
  TH2D* triggereffi_xevt[20];
  TH2D* triggereffi_yevt[20];
  TH2D* strp_xmulsim_cor[20];
  TH2D* strp_ymulsim_cor[20];

  TH2D* block_xmulsim[20][16][16];
  TH2D* block_ymulsim[20][16][16];


  TH1D* hdifftime1[20];
  TH1D* hdifftime2[20];
  TH1D* hxtime_ext[20];
  TH1D* hytime_ext[20];
 TH1D* hxpos_ext[20];
  TH1D* hypos_ext[20];
  TH1D* hxpos_ext_kalman[20];
  TH1D* hypos_ext_kalman[20];
  TH1D* h_hit_time_ext[20];



  TH1D* xtdc_minus_ref[20][8];
  TH1D* ytdc_minus_ref[20][8];

  
  
  int  isVizualisation;
  int  isVisOut;
  int  isXtermOut;
  Hits *H;
  HitPos *Hp;
  int EveCnt;
  int nloops;
  
  TH1F   *ShwXw;
  TH1F   *ShwYw;
  TH2F   *RC;
  TH1F   *DGap;
  
  TH1F *DeadStripX;
  TH1F *NoisyStripX;
  TH1F *DeadStripY;
  TH1F *NoisyStripY;
  TH1F *DiffTime;
  
//    TProfile *hprof;//SSE	
//    TH2D   *hh_E;//SSE
//    TH2D   *hh_woghst_E;//SSE
//    TH2D   *hist_orighits_new_E;//SSE
//    TH2D   *hist_orighits_mod_E;//SSE 09/15
//    TH2D   *hist_nhits_LargestCluster_E;//SSE
//    TH2D   *hist_wogh_orighits_E;//SSE

	int  ievent; //Event counter
	int  FirstEvt; //First to read from root file
	//20/02/2009 plot histogrammes, while run the code
	static const int nhistmx=1000;
	int   ihist;
  
	TH3F* gens_list[6][nhistmx];
	
	vector<vectGr> gens_vect[6];
	
  vector<double> ScntStrpXPos[4];
	vector<int>ScntStrpNo[4];
  
  vector<G4ThreeVector> ScntPos[4][4];
  
  //20/02/2009    Ranges for histogrammes
	int histxmn, histxmx, histymn, histymx, histzmn, histzmx;
	
	//12/07/09 for reco output
  
	// Common for all types of output
	
	unsigned   irun;                // Run number of these events
	unsigned   ievt;                //Event number
	unsigned   ngent;
  
	G4float	ievt_wt;		//*GMa
	int		intxn_id;		//*GMa
  
	static const unsigned int ngenmx=50;
	G4int   pidin[ngenmx]; 	  //PID of incident particle
	G4float momin[ngenmx]; 	  //Energy of incident particle
	G4float thein[ngenmx];	  //Initial polar angle of incident particle
	G4float phiin[ngenmx];     //Initial azimuthal angle of incident particle 
	G4float posxin[ngenmx];	  //Initial X-position
	G4float posyin[ngenmx];     //Initial Y-position
	G4float poszin[ngenmx];     //Initial Z-position
  G4int ngenerated;
  G4int naperture;
  // Reconstructed tracks
	
	unsigned   ntrkt;
	static  const unsigned int ntrkmx=20;
	
	G4int   itype[ntrkmx]; 	  // Fitting type, now only forward and backward
	G4int   nLayer;           //No. of layers having hit from simulation
	G4int   nhits[ntrkmx]; 	  // NUmber of hits in the track from  trackfitter
	G4int   nhits_finder[ntrkmx]; // NUmber of hits in the track from trackfinder
	G4int   vtxzplane[ntrkmx];
	G4int   endzplane[ntrkmx];
	G4float chisq[ntrkmx];        //chisquare 
	G4float cvalue[ntrkmx];       //velocity of muon 
	G4float chisq2[ntrkmx];
    G4int ndof[ntrkmx];
  G4int ndof2[ntrkmx];
  G4int  fc_or_pc[ntrkmx];
  //Backward direction
  //bit 5 :LocalPos[1];
  //bit 6 : CheckMat[1]
  //bit 7 : LocalPos[0]
  //bit 8 : CheckMat[0]
  //bit 9-14 : hadron cluster (only for backward)
  //bit 15-16 : Fiducial volume (top and bottom layer)
  //bit 17-20 : Fiducial volume (X/Y side)
  //Forward direction
  //bit 1 :LocalPos[1];
  //bit 2 : CheckMat[1]
  //bit 3 : LocalPos[0]
  //bit 4 : CheckMat[0]

	G4float trkmm[ntrkmx]; 	  //Extrapolated (half gap) Measured momentum of reconstrued track
	G4float trkth[ntrkmx];	  //Extra ...Measured polar angle of track
	G4float trkph[ntrkmx];        //Extra ....Measured azimuthal angle of track
	
	G4float momvx[ntrkmx]; 	  //Measured momentum of reconstrued track
	G4float thevx[ntrkmx];	  //Measured polar angle of track
	G4float phivx[ntrkmx];        //Measured azimuthal angle of track
	G4float posxvx[ntrkmx];	  //Starting X-position
	G4float posyvx[ntrkmx];       //Starting Y-position
	G4float poszvx[ntrkmx];       //Starting Z-position


  //ROOF RSA

	G4float momrf[ntrkmx]; 	  //Measured momentum of reconstrued track
	G4float therf[ntrkmx];	  //Measured polar angle of track
	G4float phirf[ntrkmx];        //Measured azimuthal angle of track
	G4float posxrf[ntrkmx];	  //Starting X-position
	G4float posyrf[ntrkmx];       //Starting Y-position
	G4float poszrf[ntrkmx];       //Starting Z-position


	static const unsigned int nlayermx=10;
  G4float XPosdev_exclu[nlayermx];
 G4float YPosdev_exclu[nlayermx];


  G4float XPosdev[nlayermx];
  G4float YPosdev[nlayermx];


  //for pos resolution: 
  unsigned   ntrack;
  
  G4float Xchisq[10]; 
  G4float Ychisq[10];
  G4int Xndof[10]; 
  G4int Yndof[10]; 

   G4int nXfail[10]; 
  G4int nYfail[10];
  
  Float_t inPosX[10];
  Float_t extPosX[10];
  Float_t inPosY[10];
  Float_t extPosY[10];
  Int_t nXStrips[10];
  Int_t nYStrips[10];

 

	static const unsigned int cmv_nhtmx=10000;

  unsigned int cmv_nhit;
  //cmvd cluster hit:
  unsigned int cmv_hitid[ cmv_nhtmx];
  G4int cmv_hitpdgid[ cmv_nhtmx];
  G4float cmv_hitLeTim[ cmv_nhtmx];
  G4float cmv_hitRiTim[ cmv_nhtmx];
  G4float cmv_hitLePul[ cmv_nhtmx];
  G4float cmv_hitRiPul[ cmv_nhtmx];
  G4float cmv_hitTrueposx[ cmv_nhtmx];
  G4float cmv_hitTrueposy[ cmv_nhtmx];
  G4float cmv_hitTrueposz[ cmv_nhtmx];
  G4float cmv_hitRecoposx[ cmv_nhtmx];
  G4float cmv_hitRecoposy[ cmv_nhtmx];
  G4float cmv_hitRecoposz[ cmv_nhtmx];
  

  //  G4int cmv_hitsiz[ cmv_nhtmx];

  //cmv_ cluster:


  static const unsigned int cmv_nclustmx=10000;
  unsigned int cmv_nclust;
  //cmvd cluster hit:
  unsigned int cmv_clustid[ cmv_nclustmx];
  G4int cmv_clustpdgid[ cmv_nclustmx];
  G4float cmv_clustLeTim[ cmv_nclustmx];
  G4float cmv_clustRiTim[ cmv_nclustmx];
  G4float cmv_clustLePul[ cmv_nclustmx];
  G4float cmv_clustRiPul[ cmv_nclustmx];
  G4float cmv_clustTrueposx[ cmv_nclustmx];
  G4float cmv_clustTrueposy[ cmv_nclustmx];
  G4float cmv_clustTrueposz[ cmv_nclustmx];
  G4float cmv_clustRecoposx[ cmv_nclustmx];
  G4float cmv_clustRecoposy[ cmv_nclustmx];
  G4float cmv_clustRecoposz[ cmv_nclustmx];
  G4int cmv_clustsiz[ cmv_nclustmx];









  
  //cmvd extrapol :

  static const unsigned int cmvdnlayermx=32;

  unsigned int CMVDnLayer;
  
  
  G4float CMVDTruePosX[32];
  G4float CMVDTruePosY[32];
  G4float CMVDTruePosZ[32];

  G4float CMVDRecoPosX[32];
  G4float CMVDRecoPosY[32];
  G4float CMVDRecoPosZ[32];

  G4float CMVDWRecoPosX[32];
  G4float CMVDWRecoPosY[32];
  G4float CMVDWRecoPosZ[32];
  
  G4float CMVDExpPosX[32];
  G4float CMVDExpPosY[32];
  G4float CMVDExpPosZ[32];

  G4int CMVDExpnHit[32];
  G4int CMVDisExpIn[32];
  G4int CMVDisHit[32];


  G4int CMVDpl2[32];
  G4int CMVDClustsInLay[32];
  G4int CMVDHitsInLay[32];
  G4int CMVDClusterSize[32];
  
  
  //

 
  
	
	G4float momend[ntrkmx]; 	  //Measured momentum of reconstrued track at end point
	G4float theend[ntrkmx];	  //Measured polar angle of track  at end point
	G4float phiend[ntrkmx];       //Measured azimuthal angle of track  at end point
	G4float posxend[ntrkmx];	  //End X-position
	G4float posyend[ntrkmx];      //End Y-position
	G4float poszend[ntrkmx];      //End Z-position
	G4float tx_end[ntrkmx];      //End tx
	G4float ty_end[ntrkmx];      //End ty
	
	G4float momds[ntrkmx];	  //path length
	G4float momrg[ntrkmx];	  //path length multiplied with density
	
	G4float mcxgnvx[ntrkmx];
	G4float mcygnvx[ntrkmx];
	G4float momgnvx[ntrkmx]; 	  //Generated track momentum at reconstructed vtx
	G4float thegnvx[ntrkmx];	  //Generated polar angle at reconstructed vtx
	G4float phignvx[ntrkmx];      //Generated azimuthal angle at reconstructed vtx
	
	G4float momgnend[ntrkmx]; 	  //Generated track momentum at reconstructed term
	G4float thegnend[ntrkmx];	  //Generated polar angle at reconstructed term
	G4float phignend[ntrkmx];      //Generated azimuthal angle at reconstructed term
  
  G4float tx[ntrkmx];
  G4float ty[ntrkmx];
  G4float xxin[ntrkmx];
  G4float yyin[ntrkmx];    
  G4float txin[ntrkmx];
  G4float tyin[ntrkmx];
  
  G4float therr[ntrkmx];
  G4float pherr[ntrkmx];
  G4float atimslope[ntrkmx];
  G4float atiminter[ntrkmx];
  
  G4float xxerr[ntrkmx];
  G4float yyerr[ntrkmx];
  G4float txerr[ntrkmx];
  G4float tyerr[ntrkmx];
  G4float qperr[ntrkmx];

 G4float xxtxerr[ntrkmx];
  G4float xxtyerr[ntrkmx];
  G4float yytyerr[ntrkmx];
  G4float yytxerr[ntrkmx];
    G4float txtyerr[ntrkmx];
  
  G4float xxenderr[ntrkmx];
  G4float yyenderr[ntrkmx];
  G4float txenderr[ntrkmx];
  G4float tyenderr[ntrkmx];
  G4float qpenderr[ntrkmx];
  
	G4int   ntrkcl[ntrkmx]; 	  //1000*cluster + hits associated with this track
	G4float ntrkst[ntrkmx];	  //1000*Xstrip + YStrip associated with this track
	
	G4int   ntotcl; 	  //1000*cluster + hits close to any track vertex
	G4float ntotst;	  //1000*Xstrip + YStrip close to any track vertex
	
	G4int inoclust;        //Total number of cluster reconstructed in TrackFinder FormCluster method.  //asmita_h 
	G4int origclust;       //Total number of strips taking max among xstrips or ystrips from each plane. //asmita_h
	G4int inohits;        //Total number of hits reconstructed in TrackFinder FormHits method.
	G4int orighits;       //Total number of strips taking max among xstrips or ystrips from each plane.  
	G4int x_hits;         //Sum of x strip hit from all layers
	G4int y_hits;         //Sum of y strip hit from all layers  
	G4float hPathlength;  //Pathlength //used only for pion track rigth now. 
	
	G4int inohits_old;        //Total number of hits reconstructed in TrackFinder FormHits method.
	G4int orighits_old;       //Total number of strips taking max among xstrips or ystrips from each plane.  
	G4int x_hits_old;         //Sum of x strip hit from all layers
	G4int y_hits_old;         //Sum of y strip hit from all layers  
	
  G4int orighits_trape;       // SSE 09/15 apply orighit algo on hits in trapezoid  
  G4int orighits_cluster;       //Apply orighit algo on hits  from the largest cluster SSE 09/15  
  //G4int total_inohits;//SSE// Total no of hits.	 
  G4int hit_wo_ghst;//SSE//Apply ghost hit removal (GHR) algo on the hits from the largest cluster
  G4int hit_wogh_orighits;//SSE//Hits obtained from GHR algo + orighits algo
  G4int nhits_largest_cluster;//SSE //Total no. of hits from the largest cluster
  G4float e_hadron;//SSE
  G4float theta_hadron_shw;//SSE 07Oct25
 // G4float costheta_hadron_shw;//SSE 08Oct25
  G4float phi_hadron_shw;//SSE 07Oct15
  G4float theta_hadron_in;//SSE 08Oct25 //from intital hadron direction
 // G4float costheta_hadron_in;//SSE 08Oct25 //from intital hadron direction
  G4float had_eigen_val[3]; // ADB 300316
  G4float phi_hadron_in;//SSE 08Oct15
  G4float dot_angle_had_shw;//SSE 08Oct15
	G4int nhits_largest_cluster_selected;//SSE 091015
	G4float caltot0;
  //    G4float caltot1;
  //    G4float caltot2;
	
	G4int   nmxhit;
	G4int   nhtcal0;
	
  //    G4int   nhtcal1;
  //    G4int   nhtcal2;
	
	G4float calen0[10000];
	unsigned int calid0[10000]; //This 10000 are same as 10000 in CellDetID[10000] in micalcal0SD.hh & numberInCell(10000) in micalcal0SD.cc
	
  //calid0 : 0-1 : 2bit : timing (will not use in future)
  //       2-8   7bit: Y-strip  0-127 (max)
  //       9-15  7bit: X-strip  0-127 (max)
  //       16-18  3bit: Chamber # 0-7
  //       19-21  3bit: Module #  0-7
  //       22-29  8bit: Layer $ 0-255(max)
  //       30-31  2bit : nInT (INO detector)   

  //Details coding and endoding example are in micalcal0SD.cc

  //    G4int   calid1[1000];
  //    G4float calen1[1000];

  //    G4int   calid2[1000];
  //    G4float calen2[1000];


  // For Simulation output
  unsigned int nsimht;
  static const unsigned int nsimhtmx=4000;
  unsigned int detid[nsimhtmx];
  int   simpdgid[nsimhtmx];
  float simtime[nsimhtmx];
  float simenr[nsimhtmx];
  float simvx[nsimhtmx]; 
  float simvy[nsimhtmx];
  float simvz[nsimhtmx];
  float simpx[nsimhtmx];
  float simpy[nsimhtmx];
  float simpz[nsimhtmx];

  float simlocvx[nsimhtmx];
  float simlocvy[nsimhtmx];
  float simlocvz[nsimhtmx];
  float simtotabenr;
  float simtotrpcenr;
  float simtotablen;
  float simtotrpclen;
  float range;

  //For Digitisation
  
  unsigned int ndigiht;
  static const unsigned int ndigihtmx=5000;
  int trigx;
  int trigy;
  unsigned int stripid[ndigihtmx];
  int   digipdgid[ndigihtmx];
  int   digitime[ndigihtmx];
  int   digitruetime[ndigihtmx];
  float digienr[ndigihtmx];
  float digivx[ndigihtmx]; 
  float digivy[ndigihtmx];
  float digivz[ndigihtmx];
  float digipx[ndigihtmx];
  float digipy[ndigihtmx];
  float digipz[ndigihtmx];
  
  //cmv
  
  //simulation
	static const unsigned int cmv_nsimhtmx=1000;
  unsigned int cmv_nsimhit;
  int cmv_detid[cmv_nsimhtmx];
  int  cmv_simpdgid[cmv_nsimhtmx];
  float cmv_simtime[cmv_nsimhtmx];
  float cmv_simenr[cmv_nsimhtmx];
  float cmv_simposx[cmv_nsimhtmx]; 
  float cmv_simposy[cmv_nsimhtmx];
  float cmv_simposz[cmv_nsimhtmx];
  
	Float_t cmv_simpx[cmv_nsimhtmx];
	Float_t cmv_simpy[cmv_nsimhtmx];
	Float_t cmv_simpz[cmv_nsimhtmx];
  
  float cmv_simlocx[cmv_nsimhtmx];
  float cmv_simlocy[cmv_nsimhtmx];
	 float cmv_simlocz[cmv_nsimhtmx];
  //digitization
	// for data only one 16 bit integer + one 32bit integer arrays will be there
	static const unsigned int cmv_ndigihtmx=1000;

  unsigned int cmv_ndigihit;
  int  cmv_digipdgid[cmv_ndigihtmx];
	unsigned int cmv_sipmid[cmv_ndigihtmx];
	// 3 bit for CMVD
	// 2 bit for Layer
	// 7 bit for Strip
	// 2 bit for SiPM
  unsigned int cmv_digitimpul[cmv_ndigihtmx]; // 12bit for energy and 20 bit for time

   unsigned int cmv_digitime[cmv_ndigihtmx];
   unsigned int cmv_digipul[cmv_ndigihtmx];
  float cmv_digiposx[cmv_ndigihtmx]; 
  float cmv_digiposy[cmv_ndigihtmx];
  float cmv_digiposz[cmv_ndigihtmx];
  
	float cmv_digimom[cmv_ndigihtmx];
	float cmv_digithe[cmv_ndigihtmx];
	float cmv_digiphi[cmv_ndigihtmx];
  
  float cmv_digilocx[cmv_ndigihtmx];
  float cmv_digilocy[cmv_ndigihtmx];
  float cmv_digilocz[cmv_ndigihtmx];
  //cmv sipm noise
  TFile* sipmnoise;
  TH1F* noise_hist[10][4];

  
  //cmv
  
  //For Nuance input

  /*
  struct cal0_info {
    G4int  calid0;
    G4float cal_en0;
  };

  struct cal1_info {
    G4int  calid1;
    G4float cal_en1;
  };

  struct cal2_info {
    G4int  calid2;
    G4float cal_en2;
  };
  */
private:
  micalDetectorParameterDef *paradef;
  double CorrTimeError;
  double UnCorrTimeError;
  double TimeToDigiConv;
  double SignalSpeed;
	double PhotonSpeed;
	double CMVadctons;
	
  double h2dDimX;
  double h2dDimY;
  int nbinxMag2d;
  int nbinyMag2d;

  double magZmax;
  int nbinxMagZ;
  int numberInLA;

  // micalElectroMagneticField* inoical0Field;


  
};

#endif
