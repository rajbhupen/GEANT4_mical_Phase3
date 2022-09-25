
// $Id: micalEventAction.hh,v 1.8 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef micalEventAction_h
#define micalEventAction_h 1

#include "MultiSimAnalysis.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "TCanvas.h"
#include "micalDetectorParameterDef.hh"
#include "InoHit.h"//SSE 09/15
#include "vect_manager.h"
#include "CmvHit.h"
class micalEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalEventAction : public G4UserEventAction
{
public:
  micalEventAction();
  ~micalEventAction();

public:
  void  BeginOfEventAction(const G4Event*);
  void    EndOfEventAction(const G4Event*);
    
  void Addcal0(G4double de, G4double dl) {Energycal0 += de; TrackLcal0 += dl;};
  void AddGap(G4double de, G4double dl) {EnergyGap += de; TrackLGap += dl;};
      
    void Addcal1(G4double de, G4double dl) {Energycal1 += de; TrackLcal1 += dl;};
   void Addcal2(G4double de, G4double dl) {Energycal2 += de; TrackLcal2 += dl;};
               
  void SetDrawFlag   (G4String val)  {drawFlag = val;};
  void SetPrintModulo(G4int    val)  {printModulo = val;};
  double rang; 
   
  //void shower_construct(G4int x_stripno, G4int y_stripno,  G4int z_plane, G4int counter ){
  int orighit_calc(vector<InoHit*> tmphitcluster);
        

  //InoHit_Manager* inoHit_pointer;//SSE 
  vector<InoHit*> HitBank_largestClust[500];//SSE 
  vector<InoHit*> HitBank_All[500];//SSE 


bool LinePlaneInt(double* Line, double* Plane, double* Point);
  bool ClosDistbtwLineEdge(double* Line, double* Plane, std::vector<double>& edge, double* GOnLine, double* GOnedge);
  
  void CMVD_Extrapolation();
  void CreateCmvHit();
  void FormCmvCluster();
  InoTrackCand_Manager* inoTrackCand_pointer;

  vector <CmvHit*> CmvHitBank[7][4]; //7:sides-top,left,right,back, 4:layers
private:
  micalDetectorParameterDef* paradef;
  double StripXWidth;
  double StripYWidth;
  double LayerThickness;
  int    nLayer;

  G4double  Energycal0, EnergyGap;
  G4double  TrackLcal0, TrackLGap;
     
  typedef pair<int,int> iXYZ;//SSE
  //  vector<iXYZ> xShwstrip;//SSE
  //  vector<iXYZ> yShwstrip;	//SSE
  //  int Orighits_all;//SSE
  int  x_stripno,  y_stripno , z_plane;//SSE
  // int TotalHits;//SSE
                   
     G4double  Energycal1, Energycal2;
     G4double  TrackLcal1, TrackLcal2;

  G4String  drawFlag;
  G4int     printModulo;
                             
  micalEventActionMessenger*  eventMessenger;

  G4int cal0CollID;
    G4int cal1CollID; //cmv
  //  G4int cal2CollID;

  G4int nevent;
  //  const TCanvas* c0;

	double partopscint[3];     //halflength of scintillators of thickness 1cm
	float AirGapScintTop;
double PhyVolGlPos[7][4][3];
  //	double PhyVolGlPos[4][4][3];
  int NoScntStrpTop;
  int NoScntStrpSide;

  int NoScintStrip;

  int NoScntStrpSideSmallay;
  double SideSmallPlaneHalfLength;
  double ScntLayShifSide;
   double SidePlaneHalfLength;
    double TopPlaneHalfLength;
  double magfield;
  

	int iflag = -1;

	CmvHit_Manager* CmvHit_pointer;
  SipmHit_Manager* SipmHit_pointer;
	MultiSimAnalysis *pAnalysis;
 	CmvStrip_Manager* CmvStrip_pointer;
	CmvLayExtra_Manager* CmvLayExtra_pointer;
     InoStripX_Manager* inoStripX_pointer;
     InoStripY_Manager* inoStripY_pointer;


  //     InoTrackCand_Manager* inoTrackCand_pointer;//defined for straightline fit

 	CmvCluster_Manager* CmvCluster_pointer;

	  double dZdT = 0.0;

          double ztinter = 0;


  
	
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
