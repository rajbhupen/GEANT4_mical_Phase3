/*

tar -cvf allnew.tar anal_basic_c217_newformat.cc EveTree.C EveTree.h effic_self_iter1x_150126.txt effic_self_iter1y_150126.txt effic_trig_iter2x_150126.txt effic_trig_iter2y_150126.txt time_corr_2479_150228.txt pos_time_inuse_150228.txt test.log StraightLineFit.cc StraightLineFit.h

proj-clhep.web.cern.ch/proj-clhep/DISTRIBUTION/
tar -zxvf clhep-2.1.4.2.tgz

./configure
gmake
gmake install

//g++ -g -pthread -m64 -Wno-deprecated -I${ROOTSYS}/include -I${CLHEP_BASE_DIR}/include -o anal_basic_c217_newformat EveTree.C StraightLineFit.cc anal_basic_c217_newformat.cc  `root-config --cflags --libs` -lMinuit -L${CLHEP_BASE_DIR}/lib -lCLHEP

g++ -g -pthread -m64 -Wno-deprecated -I${ROOTSYS}/include -o anal_basic_c217_newformat EveTree.C StraightLineFit.cc anal_basic_c217_newformat.cc  `root-config --cflags --libs` -lMinuit



./anal_basic_c217_newformat
testnew.log

*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <new>
#include<climits>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "EveTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"

#include "TStyle.h"
#include "TAttFill.h"
#include "TPaveStats.h"
#include "TMinuit.h"
#include "TPostScript.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom.h"
#include "TPaletteAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"

#include "StraightLineFit.h"
using namespace std;

const double cval = 29.979; // velocity of light in cm/ns

const int  nlayer =12; //Maximum RPC layers
const int  nstrip =32; // Strips in RPC chamber

bool  isonlyEff =true; //Only efficiency, no iteration
bool  firstiter = true;
const int  firstlayer =0; //Used first layer for efficiency calculation
const int  lastlayer =11; //Used last layer for efficiency calculation

const int  layfirst =0; //Used first layer in track fitting
const int  laylast =11; //Used last layer in track fitting

const int  firstXstrip =0;  //1st strip in X
const int  lastXstrip =31; // last strip in X

const int  firstYstrip =0;  //1st strip in Y
const int  lastYstrip =31; // last strip in Y

const double  mxchisq =2.0;   //Maximum value of Normalised chi^2 (position fit);
const double  mxtimechisq=2.0; // Maximum value of Normalised chi^2 (time fit);
const double  accprng =1.0; //Acceptance region with respect tot strip centre
//const double  effirng =0.25; //Additional region for efficeincy calculation
                           //Extrapolation error < 0.2 of strip width
const int xtcorstr=0;    //Starting layer for X-time correction
const int xtcorend=11;   // End layer for X-time correction
const int ytcorstr=0;    //Starting layer for Y-time correction
const int ytcorend=11;   //End layer for Y-time correction

double find_max_gap(int ix) {
  if (ix >=layfirst && ix <=laylast) return 1.5;
  if (ix <layfirst) return 1.5+0.3*(layfirst-ix);
  if (ix >layfirst) return 1.5+0.3*(ix-laylast);
  return 1.5;
}

const double pival=acos(-1);
const int nmxhits=4;
const int nmxusedhits=3;

double xposerrsq[nmxhits][nlayer] = {
  {0.0745597, 0.0601179, 0.0501015, 0.0394726, 0.0350947, 0.055039, 0.0390593, 0.0372429, 0.0398186, 0.0568854, 0.04426, 0.0647753},
  {0.0459906, 0.0543877, 0.0371751, 0.0559855, 0.0453876, 0.0360321, 0.0656442, 0.0641278, 0.043064, 0.053713, 0.0444662, 0.058869},
  {0.162488, 0.109485, 0.0785132, 0.0824725, 0.0803276, 0.0997201, 0.0863911, 0.0762118, 0.0899414, 0.0973464, 0.0703962, 0.0933884},
  {0.957458, 0.751783, 0.666736, 0.637636, 0.647876, 0.757628, 0.472752, 0.623146, 0.766528, 0.807635, 0.673675, 0.884292}
};

double yposerrsq[nmxhits][nlayer] = {
  {0.0459599, 0.0382816, 0.0381181, 0.0500019, 0.0393566, 0.0400298, 0.0522574, 0.0332281, 0.0377729, 0.0391096, 0.0295895, 0.0408621},
  {0.033369, 0.0253477, 0.0479287, 0.042458, 0.0278329, 0.0277059, 0.0528327, 0.038016, 0.0188255, 0.0222188, 0.0273548, 0.0275486},
  {0.092609, 0.112, 0.0768894, 0.0688024, 0.118428, 0.108849, 0.0775983, 0.0854793, 0.0995298, 0.110006, 0.071602, 0.0937076},
  {0.689821, 0.759425, 0.598038, 0.607062, 0.748, 0.710522, 0.621418, 0.675569, 0.720422, 0.719581, 0.702963, 0.877633}
};

Double_t gausX(Double_t* x, Double_t* par){
  return par[0]*(TMath::Gaus(x[0], par[1], par[2], kFALSE)); //kTRUE));
}

const int nmxchn = 200;
int nchannel =0;
double m_data[nmxchn];
double m_xpos[nmxchn];

void fcnsg(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t flag) {
  
  double fval=0;
  double x[2];
  for (int ij=0; ij<nchannel; ij++) {
    x[0] = m_xpos[ij];
    fval += pow( (m_data[ij] - gausX(x, par)), 2.) / TMath::Max(1., m_data[ij]);
  }
  
  f = fval;
}


Double_t fitspec(Double_t* x, Double_t* par) {
  int nx=int (x[0]/nstrip);
  double yy = x[0] - nx*nstrip;
  
  double yval = par[0] + par[1] + par[2]*sin(pival*yy/32);
  //  double yval = par[2]*pow(sin(pival*yy/32), par[3]);
  //  double yval = par[2]*sin(pival*yy/32);
  return yval;
}

//bias in time resolution due to uncertainties in other layers
double bias_intime_xreso2[12]={0.404097, 0.294606, 0.250296, 0.168577, 0.128642, 0.114688, 0.120965, 0.138543, 0.158922, 0.251805, 0.316103, 0.46363};
double bias_intime_yreso2[12]={0.508908, 0.390572, 0.342497, 0.227587, 0.167038, 0.142859, 0.140968, 0.173794, 0.230109, 0.359566, 0.454661, 0.635356};

double bias_inpos_xreso2[12]={0.0240879, 0.0176354, 0.0150597, 0.0098424, 0.00719886, 0.00583462, 0.00570846, 0.00690567, 0.00941585, 0.0144282, 0.0174868, 0.0240973};
double bias_inpos_yreso2[12]={0.0242652, 0.0176931, 0.0146815, 0.0100788, 0.00727195, 0.00586189, 0.00567934, 0.00681836, 0.00907941, 0.0137591, 0.0167822, 0.0235308};

//Time shift correction based on postion of position of track in a strip
double parabolic(double x, double* par) {
  double yy = par[0]+par[1]*abs(x)+par[2]*x*x;
  //  cout<<" yyy "<< par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<x <<" "<<fabs(x)<<" "<<abs(x)<<" "<<yy<<endl;
  return yy;
}

const int nmxtimehit=4;
const int nmxusedtimehit=2;

double strpos_vs_time[02*nmxtimehit][nlayer][3] = { // 0:1x, 1:2x, 3:3x, 4:4x, 5:1y, 6:2y, 7:3y & 8:4y
{{-0.108298, 0.234881, 0.673403}, 
{-0.0428534, 0.178566, 1.14984}, 
{-0.00179806, -0.058392, 1.48385}, 
{-0.0537578, -0.236439, 2.14217}, 
{-0.0253757, -0.187405, 3.09114}, 
{-0.0696471, -0.183559, 3.78667}, 
{0.00250928, -0.276274, 2.84852}, 
{0.00835515, -0.302841, 2.87873}, 
{-0.167504, 0.238886, 1.81847}, 
{-0.0788643, -0.159428, 1.91216}, 
{-0.0894312, 0.0856115, 1.11812}, 
{-0.0779776, 0.280962, 0.816331}}, 

{{-0.214927, 1.42247, -1.16535}, 
{-0.529811, 2.10891, -1.24669}, 
{-0.334925, 1.51751, -0.815884}, 
{-0.23738, 0.973092, 0.138824}, 
{-0.531924, 2.05702, -0.448773}, 
{-0.856707, 2.48331, 0.00234152}, 
{-0.336572, 0.837857, 0.942916}, 
{-0.34803, 1.08262, 0.572813}, 
{-0.489068, 1.92972, -0.656595}, 
{-0.30959, 1.66386, -0.946176}, 
{-0.20136, 1.233, -0.594369}, 
{-0.320393, 1.41313, -0.852566}}, 

 {{-1.13823, 0.415615, -0.0956477}, 
{-0.938263, 0.444087, -0.333743}, 
{-0.691977, 0.240515, -0.2864}, 
{-0.543465, 0.348418, -0.170337}, 
{-0.97317, 0.539603, -0.902608}, 
{-1.39851, 0.387069, -0.574954}, 
{-0.521833, 0.18826, 0.461717}, 
{-0.582583, 0.382221, -0.114176}, 
{-1.00909, -0.0518222, -0.00510963}, 
{-0.798994, 0.286262, -0.0887585}, 
{-0.551053, 0.278699, -0.448895}, 
{-0.607885, 0.123587, -0.00279212}}, 

{{-1.00401, 1.8883, -2.63737}, 
{-0.758471, 0.544713, -0.793395}, 
{-0.685073, -0.094587, -0.216575}, 
{-0.617107, 0.49322, -1.20662}, 
{-0.946502, 0.632696, -0.874812}, 
{-1.19493, 1.00998, -1.33925}, 
{-0.69642, 0.0394355, -0.33321}, 
{-0.640847, 0.168112, -0.245247}, 
{-0.670089, 0.702369, -1.27111}, 
{-0.736155, 0.147291, 0.0974525}, 
{-0.692997, 0.0315732, -0.6247}, 
{-0.607172, -0.669918, 0.388367}},

{{-0.00329939, 0.242462, 1.06257}, 
{-0.0902285, 0.0624302, 2.03717}, 
{-0.0578305, 0.0933354, 1.96958}, 
{-0.0585525, 0.401106, 1.45776}, 
{-0.0667166, -0.319948, 3.5735}, 
{-0.120356, -0.369183, 4.42829}, 
{-0.298342, 0.70353, 1.81377}, 
{-0.150394, -0.433648, 4.82061}, 
{-0.102069, 0.0859978, 2.08871}, 
{-0.124022, 0.0217614, 1.94001}, 
{-0.0913154, -0.114765, 2.79464}, 
{-0.114681, 0.102416, 1.63289}}, 

{{-0.272125, 1.68014, -1.17568}, 
{-0.487944, 2.29154, -1.52053}, 
{-0.323852, 1.38964, -0.240395}, 
{-0.34994, 1.6006, -0.671717}, 
{-0.718347, 2.41741, -0.471764}, 
{-0.721208, 2.50796, -0.40962}, 
{-0.267072, 1.43869, -0.586951}, 
{-0.649484, 2.28938, -0.168097}, 
{-0.95087, 3.55793, -1.35827}, 
{-0.881629, 3.59384, -2.48311}, 
{-0.435896, 1.95779, -1.00669}, 
 {-0.390051, 1.95935, -1.36556}},

{{-1.18847, 0.0757272, 0.117748}, 
{-1.19774, 0.267215, 0.0940863}, 
{-0.844947, 0.544917, -0.412739}, 
{-0.731456, 0.455506, -0.742424}, 
{-1.49598, 0.688469, 0.132483}, 
{-1.38607, 0.569129, -0.348203}, 
{-0.514269, 0.568307, -1.28993}, 
{-1.04816, 0.668811, -0.60223}, 
{-1.29873, 0.645937, -0.852456}, 
{-1.51009, 0.171537, 0.353457}, 
{-0.923922, 0.40255, -0.626108}, 
{-1.04374, 0.076331, -0.0552816}},

{{-0.92741, -0.989598, 1.67875}, 
{-0.670879, -0.485287, 1.01524}, 
{-0.791233, 0.0619398, -0.343801}, 
{-0.781125, -0.409228, -0.243095}, 
{-0.975324, 0.466529, -0.269322}, 
{-1.04694, 0.801549, -0.997753}, 
{-0.556878, 0.961198, -2.14108}, 
{-0.928105, 1.69794, -2.72732}, 
{-1.01487, 0.592797, 0.340055}, 
{-1.00587, -0.0115691, 1.0955}, 
{-0.855476, 0.00859748, -0.137312},
{-0.803551, 0.247394, -1.63056}}
};

double biasinxtime[nlayer]={0};

double biasinytime[nlayer]={0};

const double stripwidth = 3.0; // cm
const double stripgap = 0.2; // cm
const double layergap = 16.0; // cm
const double fact = stripwidth/layergap;
const double sigmaz = 0.1; //in cm as the smallest division in scale
const double sigmarpc = 1.5; //in ns RPC time resolution
const int nmnhits =6; //Minimum number of layers required for position and time measurements

const int nmnentry = 10; // Statistics requred in a strip for time correction and use in next iteration

bool isTiming =true; // true //Don't do timing fit etc
bool useallpos=false; //true; //While use all strip irrespecive of aligned or not
bool usealltime =false; //While use timing of all strip irrespecive of aligned or not
bool timefit;

const float xyPosDev=7.0; // seven sigma 2.0; //maximum deviation of points from fit line

const int trigly1 =2; //0;extpo
const int trigly2 =4;
const int trigly3 =7;
const int trigly4 =9; //11;

void GetXposInStrip(double* ext, double* off, double* pos) {
  for (int ij=0; ij<nlayer; ij++) {
    int istr = int(ext[ij]+off[ij]);
    if (ext[ij]+off[ij]<0.0) {
      pos[ij] = ext[ij]+off[ij] - istr + 0.5;
    } else {
      pos[ij] = ext[ij]+off[ij] - istr - 0.5;
    }
    //    cout <<"ij "<<ij<<" "<< ext[ij]<<" "<<off[ij]<<" "<<pos[ij]<<endl;
  }
}

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -1;
  for (int ix=1; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return 1000;
}

int main() {
  // ----------------- Don't Change Anything Starting from Here ------------       ********************   ---------------------------

  static unsigned int mypow_2[32];
  for (int ij=0; ij<32; ij++) {
    mypow_2[ij] = pow(2, ij);
  }

  gStyle->SetPalette(1,0);
  gStyle->SetFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatStyle(1001);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);

  gStyle->SetStatFont(22);        // Times New Roman
  gStyle->SetTextFont(22);        // Times New Roman
  gStyle->SetTitleFont(22,"XYZ"); // Times New Roman
  gStyle->SetLabelFont(22,"XYZ"); // Times New Roman
  gStyle->SetLabelSize(0.06, "XYZ"); // Times New Roman  
  gStyle->SetNdivisions(606, "XYZ");

  gStyle->SetOptTitle(0);
  gStyle->SetFuncWidth(1);
  gStyle->SetFuncColor(2);
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(101);
  gStyle->SetOptLogy(0);
  gStyle->SetStatW(.18);
  gStyle->SetStatH(.08);
  gStyle->SetPadTopMargin(.02); //0.09
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadLeftMargin(0.02);
  gStyle->SetPadRightMargin(0.02);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.12);
  latex.SetTextFont(42);
  latex.SetTextAlign(1); //(31); // align right

  int EvtIndx,yyindx,xxindx,BID;

  bool posinuse[2][nlayer][nstrip]; //Reject strips in timing, which are not aligned
  bool timeinuse[2][nlayer][nstrip]; //Reject strips in timing, which are not aligned
  
  for ( int ij=0; ij<2; ij++) {
    for ( int jk=0; jk<nlayer; jk++) {
      for ( int kl=0; kl<nstrip; kl++) {
	posinuse[ij][jk][kl] = timeinuse[ij][jk][kl]=true;
      }
    }
  }
  // To calcualte average X-time and Y-time to synchronise them.
  int ntotxtime=0;
  double totxtime=0;
  int ntotytime=0;
  double totytime=0;

  int isequence[nlayer]={0,1,2,3,4,5,6,7,8,9,10,11};
  int nentrymx=-1;
  int isalign=0;

  char outfile[100];
  char outfilx[100];

  char name[100];
  char title[100];
  char infile[200];
  char datafile[100];
  char rootfiles[100];
  cout <<"Give the input file name"<<endl;;
  cin>> rootfiles;
  //  sprintf(rootfiles, "test_2479.log");
  
  if (isonlyEff) { 
    isalign = 0;
  } else {
    //    cout <<"Do you want to align ? yes/no 1/0" <<endl;;
    //    cout <<"for alignment, it can be any number greater than zero"<<endl;
    //    cin>> isalign;
    isalign =1;
  }
  
  //# of iteration where all layers are included
  // For the time being it is not implemented, but can be used
  // Can be use in other way, first few iteration with combined+individual, then
  // Only combined
  //Total # of iteration all layers + individual layers
  const int nmxiter = (isalign>0) ? 1 : 1; //Less than 12, otherwise change the pad in canvas
  const int nlayerit = nlayer; // (isalign >0) ? nlayer : 1;
  
  int ievt, nhits;
  double xslope, xinters, yslope, yinters, timexslope, timeyslope, timex2slope, timey2slope;
  float txslop, tyslop, timexinters, timeyinters, errtimexslope,errtimeyslope,errtimexinters, errtimeyinters;
  double xchi2, ychi2, xt0chi2, yt0chi2;
  int nxstrip, nystrip, Nx, Ny, nxtime, nytime, ntxyla;
  double zen;

  int len = strlen(rootfiles);
  strncpy(outfilx, rootfiles, len-4);
  outfilx[len-4]='\0';
  sprintf(outfilx, "%s_basic_", outfilx);
  len = strlen(outfilx);
  outfilx[len]='\0';

  sprintf(outfile, "%s%i.root", outfilx, isalign);
  TFile* fileOut = new TFile(outfile, "recreate");

  sprintf(outfile, "%s%i.ps", outfilx, isalign);
  TPostScript ps(outfile,111);  
  ps.Range(20,30); //ps.Range(10,20);
  
  sprintf(outfile, "%s%i.txt", outfilx, isalign);
  ofstream file_out(outfile);

  sprintf(outfile, "%s%i_str.txt", outfilx, isalign);
  ofstream file_outstr(outfile);

  TTree* T1 = new TTree("T1", "store"); //("Tree Name","Tree Title")

  T1->Branch("Evt",&ievt,"ievt/I");
  T1->Branch("nhits", &nhits, "nhits/I");
  T1->Branch("xslope", &xslope, "xslope/D");
  T1->Branch("xinters", &xinters, "xinters/D");
  T1->Branch("yslope", &yslope, "yslope/D");
  T1->Branch("yinters", &yinters, "yinters/D");

  T1->Branch("xchi2", &xchi2, "xchi2/D"); 
  T1->Branch("ychi2", &ychi2, "ychi2/D"); 
  T1->Branch("zen",&zen,"zen/D");
  
  if (isTiming) {
    
    T1->Branch("txslop", &txslop, "txslop/F");
    T1->Branch("tyslop", &tyslop, "tyslop/F");
    T1->Branch("xt0chi2", &xt0chi2, "xt0chi2/D");
    T1->Branch("yt0chi2", &yt0chi2, "yt0chi2/D");
    T1->Branch("nxtime", &nxtime, "nxtime/I");
    T1->Branch("nytime", &nytime, "nytime/I");
    T1->Branch("ntxyla", &ntxyla, "ntxyla/I");
    T1->Branch("timexslope",&timexslope,"timexslope/D");
    T1->Branch("timeyslope",&timeyslope,"timeyslope/D");
    T1->Branch("errtimexslope",&errtimexslope,"errtimexslope/F");
  }

  int narray =60;

  TH1F* costhe[6];
  costhe[0] = new TH1F("costhe_all", "costhe_all", narray, 0., narray);//filling as dn/d(theta)
  costhe[1] = new TH1F("costhe_trig", "costhe_trig", narray, 0., narray);
  costhe[2] = new TH1F("costhe_selecx", "costhe_selecx", narray, 0., narray);
  costhe[3] = new TH1F("costhe_selecy", "costhe_selecy", narray, 0., narray);
  costhe[4] = new TH1F("costhe_accep", "costhe_accep", narray, 0., narray);
  costhe[5] = new TH1F("costhe_accep_H", "costhe_accep_H", 10, 0.5, 1.0);//10Nov Honda also give 0 to 60 deg in same bin filling as dn/d(costheta)

  TH1F* phiang[6];
  
  phiang[0] = new TH1F("phiang_all", "phiang_all", 36, -pival, pival);
  phiang[1] = new TH1F("phiang_trig", "phiang_trig", 36, -pival, pival);
  phiang[2] = new TH1F("phiang_selecx", "phiang_selecx", 36, -pival, pival);
  phiang[3] = new TH1F("phiang_selecy", "phiang_selecy", 36, -pival, pival);
  phiang[4] = new TH1F("phiang_accep", "phiang_accep", 36, -pival, pival);

  TH2F* dir_cxchi = new TH2F("dir_cxchi", "dir_cxchi", 100, 0., 100., 200, -20., 20.);
  TH2F* dir_cychi = new TH2F("dir_cychi", "dir_cychi", 100, 0., 100., 200, -20., 20.);
  
  TH2F* dir_cxy = new TH2F("dir_cxy", "dir_cxy", 100, -2., 2., 100, -2., 2.);
  
  TH1F* xlayer_reso[nlayer][2*nmxiter];
  TH1F* ylayer_reso[nlayer][2*nmxiter];

  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<2*nmxiter; jk++) {
      sprintf(title, "xlayer_reso_l%i_i%i", ij, jk);
      xlayer_reso[ij][jk]=new TH1F(title, title, 150, -2.0, 2.0);
      sprintf(title, "ylayer_reso_l%i_i%i", ij, jk);
      ylayer_reso[ij][jk]=new TH1F(title, title, 150, -2.0, 2.0);
    }
  }

  TH1F* time_xreso[nlayer][2*nmxiter];
  TH1F* time_yreso[nlayer][2*nmxiter];
  
  TH1F* time_xstrreso[nlayer][nstrip][2*nmxiter];
  TH1F* time_ystrreso[nlayer][nstrip][2*nmxiter];

  if (isTiming) {
    for (int kl=0; kl<2*nmxiter; kl++ ){
      for (int ij=0; ij<nlayer; ij++) {
	sprintf(title, "time_xreso_l%i_i%i", ij, kl);
	time_xreso[ij][kl] = new TH1F(title, title, 150, -10.0, 10.0); //(-7.5, 7.5);
	
	sprintf(title, "time_yreso_l%i_i%i", ij, kl);
	time_yreso[ij][kl] = new TH1F(title, title, 150, -10.0, 10.0); //(-7.5, 7.5);   

	for (int jk=0; jk<nstrip; jk++) {
	  sprintf(title, "time_xstrreso_l%i_s%i_i%i", ij, jk, kl);
	  time_xstrreso[ij][jk][kl] = new TH1F(title, title, 120, -7.5, 7.5);
	  
	  sprintf(title, "time_ystrreso_l%i_s%i_i%i", ij, jk, kl);
	  time_ystrreso[ij][jk][kl] = new TH1F(title, title, 120, -7.5, 7.5);
	  
	}
      }
    }
  }

  TH1F* occu_x[nlayer];
  TH1F* occu_y[nlayer];
  TH2F* raw_occu[nlayer];

  TH1F* xlayer_mult[nlayer];
  TH1F* ylayer_mult[nlayer];

  for (int ij=0; ij<nlayer; ij++) {
    sprintf(title, "occu_x%i", ij);
    occu_x[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "occu_y%i", ij);
    occu_y[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
    
    sprintf(title, "raw_occu_l%i", ij);
    raw_occu[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "xlayer_mult_l%i", ij);
    xlayer_mult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);
    
    sprintf(title, "ylayer_mult_l%i", ij);
    ylayer_mult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

  }

  TH1F* h_chisqx = new TH1F("chisqx", "chisqx", 120, 0.0, 150.0);
  TH1F* h_reduchisqx = new TH1F("reduchisqx", "reduced chisqx", 90, 0.0, 30.0/*45.0*/);
  TH1F* h_chisqy = new TH1F("chi2y", "chisqy", 120, 0.0, 150.0);
  TH1F* h_reduchisqy = new TH1F("reduchisqy", "reduced chisqy", 90, 0.0, 30.0/*45.0*/);
  TH1F* h_xndf = new TH1F("xndf", "xndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_yndf = new TH1F("yndf", "yndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_xprob = new TH1F("xprob", "xprob", 120, 0.0, 1.0);
  TH1F* h_yprob = new TH1F("yprob", "yprob", 120, 0.0, 1.0);
  
  TH1F* h_tchisqx = new TH1F("tchisqx", "tchisqx", 120, 0.0, 90.0);
  TH1F* h_treduchisqx = new TH1F("treduchisqx", "reduced tchisqx", 90, 0.0, 30.0);
  TH1F* h_tchisqy = new TH1F("tchi2y", "tchisqy", 120, 0.0, 90.0);
  TH1F* h_treduchisqy = new TH1F("treduchisqy", "reduced tchisqy", 90, 0.0, 30.0);
  TH1F* h_txndf = new TH1F("txndf", "txndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_tyndf = new TH1F("tyndf", "tyndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_xtprob = new TH1F("xtprob", "xtprob", 120, 0.0, 1.0);
  TH1F* h_ytprob = new TH1F("ytprob", "ytprob", 120, 0.0, 1.0);


  TH1F* pos_xslope = new TH1F("pos_xslope", "pos_xslope", 120, -3.6, 3.6);
  TH1F* pos_yslope = new TH1F("pos_yslope", "pos_yslope", 120, -3.6, 3.6);

  TH1F* pos_theta = new TH1F("pos_theta", "pos_theta", 120, 0.0, 80.0);
  TH1F* pos_phi = new TH1F("pos_phi", "pos_phi", 120, -180.0, 180.0);

  TH1F* dir_cx = new TH1F("dir_cx", "dir_cx", 180, -1.5, 3.5);
  TH1F* dir_cy = new TH1F("dir_cy", "dir_cy", 180, -1.5, 3.5);
  
  int xhits[nlayer],yhits[nlayer];   //number of hits after noise rejection for position fit
  int xallhits[nlayer],yallhits[nlayer]; // raw : for one hit, return strip number otherwise -10-multiplicity
  
  bool passxtime[nlayer], passytime[nlayer]; // passed through timing criteria in all hitted strips in x/ystr_x/ytdev
  bool passxtmy[nlayer], passytmx[nlayer]; //for X/Y strips, use boundary of Y/X extrapolation
  int istrxtime[nlayer], istrytime[nlayer]; //strip with early timing

  float errxco[nlayer]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};
  float erryco[nlayer]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};
  
  double xxerr[nlayer], yyerr[nlayer];
  for ( int ix=0; ix<nlayer; ix++) {xxerr[ix] = yyerr[ix] = errxco[ix]*errxco[ix];}

  double xrms[nlayer]={0};
  double yrms[nlayer]={0};

#include "effic_self_iter1x_150126.txt"
#include "effic_self_iter1y_150126.txt"

#include "effic_trig_iter2x_150126.txt"
#include "effic_trig_iter2y_150126.txt"

#include "time_corr_2479_150228.txt"
#include "pos_time_inuse_150228.txt"

//-------------------------------------------------------------------------

  double timesx[nlayer]={0};
  double timesy[nlayer]={0};

  double widthx[nlayer]={0};
  double widthy[nlayer]={0};

  //  bool filloccu = true; // For first entries fill occuplancy plot and remove noisy channels
  
  double errcst, errcov, errlin;
  double dist[nlayer], xtime[nlayer],ytime[nlayer],  xtdev[nlayer],ytdev[nlayer];
  double rawxtime[nlayer], rawytime[nlayer], rawxtime1[nlayer], rawytime1[nlayer], rawxtime2[nlayer], rawytime2[nlayer];

  bool xusedtime[nlayer], yusedtime[nlayer];

  double xval, yval; 
  int nxtfail, nytfail, nentry, nTotalp, nTotalt;
  double xc0inters,yc0inters,DDx,errcst_tx,errcov_tx,errlin_tx,DDy,errcst_ty,errcov_ty,errlin_ty;
  int ntcormx, iiterrs, lstr, lend,occulyr ,isfill,ntrigX,ntrigY;

  double xtext[nlayer], xtexter[nlayer];
  double ytext[nlayer], ytexter[nlayer];

  //*********************************************************************************************
  //               Iteration Starts
  //*********************************************************************************************
      
  nTotalp = nTotalt = 0;
  for (int iiter=0; iiter<nmxiter; iiter++) {
    ntotxtime = totxtime = ntotytime = totytime = 0;

    //    firstiter=0; // before fit, calculate time shift in each layer with all hists
    
    //    int ntcormx = 2;
    //Check this condition properly, why should we use this
    ntcormx =(layfirst-firstlayer>0 || lastlayer-laylast>0) ? 1 : 2;

    if (isonlyEff) ntcormx=1;

    for (int ntcor=0; ntcor< ntcormx; ntcor++) {
      iiterrs=nmxiter*ntcor+iiter;
      
      lstr = max(firstlayer,0);  //Always 0
      //      lstr = max(firstlayer,2);  //GMA140130

      lend = min(lastlayer+1,nlayerit); //nlayerit=nlayer=12     //Always 12
      
      if (ntcor==0) { lend = lstr+1;}
      
      file_out <<"lay1 "<< iiter<<" "<<ntcor<<" "<<iiterrs<<" "<<lstr<<" "<<lend<<endl;
      
      for (int laye= lstr; laye<lend; laye++) {
	occulyr = (ntcor==1) ? isequence[abs(laye)] : nlayer;
	
        isfill = (iiter==nmxiter-1 && ntcor==0) ? true : false; //Fill rootuple and other histogramme for the final iteration with all layers
	
	ifstream file_db;
	file_db.open(rootfiles);  
	
	while(!(file_db.eof())) {
	  file_db >> datafile>>nentrymx;
	  cout<<"datafile "<<datafile<<endl;
	  if (strstr(datafile,"#")) continue;
	  if(file_db.eof()) break;
	  //	  sprintf(infile, "/data/gobinda/ino/rpcdata/daq/%s", datafile);
	  sprintf(infile, "daq/%s", datafile);

          TFile *fileIn = new TFile(infile, "read");
	  TTree *event_tree= (TTree*)fileIn->Get("evetree");

	  EveTree *event = new EveTree(event_tree);
          event->Loop();
	  
          nentry = event_tree->GetEntries();
	  //	  nentry = 1000;
          nentry = min(nentry,nentrymx);
	  cout <<infile<<" has "<< nentry<<" events "<<endl;
	  int ntiming=0;
          for(int ievt=0;ievt<nentry;ievt++) {    //ij is event loop variable. while checking put ij<3 or a small number.

	    xslope= xinters= yslope= yinters= timexslope= timeyslope= timex2slope= timey2slope=-100.;
	    xchi2= ychi2= xt0chi2= yt0chi2=-100.;
	    nxstrip= nystrip=Nx=Ny=nxtime=nytime=ntxyla=0;
	    
	    fileIn->cd();
	    event_tree->GetEntry(ievt);    // while running for entire event file put "ij<numentries".
	    
	    vector<int> xpts[nlayer];
	    vector<int> ypts[nlayer];

	    vector<int> xptsall[nlayer]; //For trigger criteria
	    vector<int> yptsall[nlayer];	    
	    
	    fileOut->cd();
	    
	    for(int jk=0;jk<nlayer;jk++) {
	      for(int kl=0; kl<nstrip; kl++) {
		if(event->xLayer[jk]->TestBitNumber(kl)) {
		  xptsall[jk].push_back(kl);
		}
		
		if(event->yLayer[jk]->TestBitNumber(kl)) {
		  yptsall[jk].push_back(kl);
		}
	      }
	    }
	    if (firstiter) { 
	      nTotalp++;
	      //	      cout<<"XSIDE "<<endl;
	      for(int jk=0;jk<nlayer;jk++) {
		for (int kl=0; kl< xptsall[jk].size(); kl++) {
		  occu_x[jk]->Fill(xptsall[jk][kl]);
		  for (int lm=0; lm< yptsall[jk].size(); lm++) {
		    raw_occu[jk]->Fill(xptsall[jk][kl], yptsall[jk][lm]);
		  }
		}
		//		cout<<endl;
	      }
	      
	      //	      cout<<"YSIDE "<<endl;
	      for(int jk=0;jk<nlayer;jk++) {
		//		cout<<"L "<<jk<<" ";
		for (int kl=0; kl< yptsall[jk].size(); kl++) {
		  //		  cout<<yptsall[jk][kl]<<", ";
		  occu_y[jk]->Fill(yptsall[jk][kl]);
		}
		//		cout<<endl;
	      }
	    }


	    // Store total number of hits in X/Y layer irrespective of noise etc.
	    for (int iz=0; iz<nlayer; iz++) {
	      xallhits[iz] = (xptsall[iz].size()==1) ? xptsall[iz][0] : -10-xptsall[iz].size(); //avoid confusion with no hits and only hits at strip# 1
	      yallhits[iz] = (yptsall[iz].size()==1) ? yptsall[iz][0] : -10-yptsall[iz].size();
	    }

	    ////////////////////////////////////////////
	    //
	    //  First clean up noisey layers and then
	    //  Strip efficiency and resolutions etc etc
	    //
	    ////////////////////////////////////////////
	    
	    xslope = 0;
	    xinters = 0;
	    yslope = 0;
	    yinters = 0;
	    ntrigX = ntrigY = 0;

	    for (int iz=0; iz<nlayer; iz++) {
	      xhits[iz] = yhits[iz] =0;
	      if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && xptsall[iz].size()>0) {ntrigX++;} //18 01 2012 
	      for (int ix=0; ix<xptsall[iz].size(); ix++) {
		bool failed=false;
		if (!posinuse[0][iz][xptsall[iz][ix]]) {failed=true;}
		if (!failed){
		  xpts[iz].push_back(xptsall[iz][ix]);
		}
	      }
	      xhits[iz] = xpts[iz].size();
	      
	      if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && yptsall[iz].size()>0) { ntrigY++;} //18 01 2012
	      for (int iy=0; iy<yptsall[iz].size(); iy++) {
		bool failed=false;
		if (!posinuse[1][iz][yptsall[iz][iy]]) {failed=true;}
		if (!failed){
		  ypts[iz].push_back(yptsall[iz][iy]);
		}
	      }
	      yhits[iz] = ypts[iz].size();
	    }
	    
	    for (int ij=0; ij<nlayer; ij++) {
	      xlayer_mult[ij]->Fill(xpts[ij].size());
	      ylayer_mult[ij]->Fill(ypts[ij].size());
	    }


	    //	    if(ntrigX<4) continue;  
	    
	    double Xpos[nlayer]; 
	    bool Xusedpos[nlayer];
            double Xdev[nlayer]; for (int ij=0; ij<nlayer; ij++) { Xdev[ij] = 100; Xpos[ij]=0.0;}
	    
            for (int ij=0;ij<nlayer;ij++) {
	      if (xhits[ij]<=0 || xhits[ij]>nmxhits) {
		Xpos[ij]= -100;
	      } else {
		for (int ix=0; ix<xhits[ij]; ix++) {
		  if (ix<xhits[ij]-1 && abs(xpts[ij][ix]-xpts[ij][ix+1])>1) { Xpos[ij]=-100; break;}
		  Xpos[ij] +=xpts[ij][ix];
		}
		if (Xpos[ij]>=0.0) {
		  Xpos[ij]  = Xpos[ij]/xhits[ij] + 0.5 - xoff[ij];
		  xxerr[ij] = xposerrsq[xhits[ij]-1][ij];
		}
	      }
	    }

	    //Sort out hits, which can be used for fit
	    for (int ij=0;ij<nlayer;ij++) {
	      Xusedpos[ij] = (Xpos[ij]>-100 && xhits[ij]<=nmxusedhits) ? true : false; //Xpos[ij] : -101;
	    }
	    



            Nx=0;
	    int nxfail = 0;
	    xchi2 = 0;
	    double xresol = 0;
	    double zval[nlayer], xext[nlayer], xexter[nlayer], xposinstr[nlayer];
	    for (int ix=0; ix<nlayer; ix++) { zval[ix]=ix;}
	    for (int ix=0; ix<nlayer; ix++) { xext[ix]= xexter[ix] =xposinstr[ix] =  100;}
	    StraightLineFit xposfit(1, zval, Xpos,  xxerr, Xusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
	    xposfit.GetParameters(nxfail, xinters, xslope);
	    //	    xposfit.GetError(errcst, errlin, errcov);
	    xposfit.GetChisqure(Nx, xchi2);
	    xposfit.GetFitValues(xext, Xdev, xexter);

	    GetXposInStrip(xext, xoff, xposinstr);
	    
	    if (nxfail==0 && isfill) {
              h_chisqx->Fill(xchi2);
              if (Nx-2>0) {
		h_reduchisqx->Fill(xchi2/(Nx-2));
		double probx = TMath::Prob(xchi2, Nx-2);
		h_xprob->Fill(probx);
	      }
              h_xndf->Fill(Nx);
            }

	    double Ypos[nlayer];
	    bool Yusedpos[nlayer];//=new float[nlayer];
            double Ydev[nlayer]; for (int ij=0; ij<nlayer; ij++) { Ydev[ij] = 100; Ypos[ij]=0.0;}
	    
            for (int ij=0;ij<nlayer;ij++) {
	      if (yhits[ij]<=0 || yhits[ij]>nmxhits) {
		Ypos[ij]=-100;
	      } else {
		for (int iy=0; iy<yhits[ij]; iy++) {
		  if (iy<yhits[ij]-1 && abs(ypts[ij][iy]-ypts[ij][iy+1])>1) { Ypos[ij]=-100; break;}
		  Ypos[ij] +=ypts[ij][iy];
		}
		if (Ypos[ij]>=0.0) {
		  Ypos[ij]  = Ypos[ij]/yhits[ij] + 0.5 - yoff[ij];
		  yyerr[ij] = yposerrsq[yhits[ij]-1][ij];
		}
	      }
	    }

	    //Sort out hits, which can be used for fit
	    for (int ij=0;ij<nlayer;ij++) {
	      Yusedpos[ij] = (Ypos[ij]>-99 && yhits[ij]<=nmxusedhits) ? true : false; //Ypos[ij] :  -101;
	    }

            Ny=0;
	    int nyfail = 0;
	    ychi2 = 0;
	    double yresol = 0;
	    double yext[nlayer], yexter[nlayer], yposinstr[nlayer];;
            for (int ix=0; ix<nlayer; ix++) { yext[ix]= yexter[ix] =yposinstr[ix] =  100;}
	    
	    StraightLineFit yposfit(1, zval, Ypos,  yyerr, Yusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
	    yposfit.GetParameters(nyfail, yinters, yslope);
	    //	    yposfit.GetError(errcst, errlin, errcov);
	    yposfit.GetChisqure(Ny, ychi2);
	    yposfit.GetFitValues(yext, Ydev, yexter);
	    GetXposInStrip(yext, yoff, yposinstr);

	    if (nyfail==0 && isfill) {
	      h_chisqy->Fill(ychi2);
	      if (Ny-2>0) {
		h_reduchisqy->Fill(ychi2/(Ny-2));
		double probx =TMath::Prob(ychi2, Ny-2);
		h_yprob->Fill(probx);
	      }
	      h_yndf->Fill(Ny);
	    }

	    nhits = 100*Nx + Ny;

            if (Nx>= nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) {//4Nov,2011
	      if (occulyr>=nlayer) {
		pos_xslope->Fill(stripwidth*xslope);
		for (int ij=0; ij<nlayer; ij++) {
		  if (abs(Xdev[ij])<6.0 && Xusedpos[ij]) { 
		    xlayer_reso[ij][iiterrs]->Fill(Xdev[ij]);
		  }
		}
	      } else {
		if (abs(Xdev[occulyr])<6.0 && Xusedpos[occulyr]) {
		  xlayer_reso[occulyr][iiterrs]->Fill(Xdev[occulyr]);
		} //if (iiter==nmxiter-1) 
	      }
	    }


            

            if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {//4Nov
              if (occulyr>=nlayer) {
		pos_yslope->Fill(stripwidth*yslope);
                for (int ij=0; ij<nlayer; ij++) {
                  if (abs(Ydev[ij])<6.0  && Yusedpos[ij]) {
                    ylayer_reso[ij][iiterrs]->Fill(Ydev[ij]);
                  }
                }
	      } else {
                if (abs(Ydev[occulyr])<6.0 && Yusedpos[occulyr]) {
                  ylayer_reso[occulyr][iiterrs]->Fill(Ydev[occulyr]);
                }
	      }
	    }

	        if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {//4Nov
	      if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) {

                double fitthe1 = acos(sqrt(1./(1+pow(fact*xslope,2.)+pow(fact*yslope,2.)))); //10Nov
		double fitthe = (180./pival)*fitthe1; //acos(sqrt(1./(1+pow(fact*xslope,2.)+pow(fact*yslope,2.)))); 
		
		//       x1 x2 x3 ...............................x30   x31
		//       ^  ^  ^  ^  ^  ^  ^  ^  W  ^  ^  ^  ^  ^  ^  ^  ^
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y1
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y2
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//  S<-- |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.--> N
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y30
		//       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y31
                //                               | 
		//                               E
		
    //		double fitphi = atan2(-xslope, yslope);  // With respect to South 

		double fitphi = atan2(xslope, yslope);  // With respect to South 


		if (isfill) { 
		  pos_theta->Fill(fitthe);
		  pos_phi->Fill(fitphi*180./pival);

		  costhe[0]->Fill(fitthe, 1.0);
		  phiang[0]->Fill(fitphi, 1.0);
		  
		  if (xpts[trigly1].size()>0 &&
		      xpts[trigly2].size()>0 &&
		      xpts[trigly3].size()>0 &&
		      xpts[trigly4].size()>0) {
		    costhe[1]->Fill(fitthe, 1.0);
		    phiang[1]->Fill(fitphi, 1.0);
		    
		    if (Xpos[trigly1] >-99 && Xpos[trigly2] >-99 &&
			Xpos[trigly3] >-99 && Xpos[trigly4] >-99) {
		      costhe[2]->Fill(fitthe, 1.0);
		      phiang[2]->Fill(fitphi, 1.0);
		      
		      if (Ypos[trigly1] >-99 && Ypos[trigly2] >-99 &&
			  Ypos[trigly3] >-99 && Ypos[trigly4] >-99) {
			costhe[3]->Fill(fitthe, 1.0);
			phiang[3]->Fill(fitphi, 1.0);
			
			if (abs(Xdev[trigly1]) < 1.5 && Xusedpos[trigly1] && abs(Ydev[trigly1]) < 1.5 &&
			    abs(Xdev[trigly2]) < 1.5 && Xusedpos[trigly2] && abs(Ydev[trigly2]) < 1.5 &&
			    abs(Xdev[trigly3]) < 1.5 && Xusedpos[trigly3] && abs(Ydev[trigly3]) < 1.5 &&
			    abs(Xdev[trigly4]) < 1.5 && Xusedpos[trigly4] && abs(Ydev[trigly4]) < 1.5) {
			  costhe[4]->Fill(fitthe, 1.0);
			  phiang[4]->Fill(fitphi, 1.0);
			  costhe[5]->Fill(cos(fitthe1),1.0);//10Nov
			  zen = fitthe;
			}
		      }
		    }
		  }
		} // if (isfill)
	      } // if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0)
            } // if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0)

            // in the root structure; otherwise it will fill after time fit
	    //GMA            if ( !isTiming && isfill &&  nxfail==0 &&  nyfail==0 && Nx>nmnhits && Ny>nmnhits && xchi2/(Ny-2) >mxchisq && ychi2/(Nx-2) >mxchisq) T1->Fill();
            // Now filling T1 contains only proper track fit data
	    nxtime=0;
	    nytime=0;
	    ntxyla=0;
	    
	    if (isTiming&& Nx>=nmnhits/*-ntcor*/ && Ny>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && ychi2/(Ny-2)<mxchisq && nxfail==0 && nyfail==0) {	      
              //////////////////////////////////////////////
              //                                          //
              //  Timing informations and directionality  //
              //                                          //
              //////////////////////////////////////////////
	      
	      //Store time informations
	      
	      for (int jk=0; jk<32; jk++) {
		if (jk<nlayer) {
		  timesx[jk] = (event->tdcdata[jk] < 50.) ? -100. : 0.1*event->tdcdata[jk]; //<500.
		  widthx[jk] = 0.0; //
		} else if (jk>=16 && jk<nlayer+16) {
		  int jkm=jk-16;
		  timesy[jkm] = (event->tdcdata[jk] < 50.) ? -100. : 0.1*event->tdcdata[jk]; //<500.
		  widthy[jkm] = 0.0; //
		}
	      }

	      //	      cout<<"timing "<<ievt<<" "<<ntiming++<<endl;
	      ntiming++;
	      
	      for (int ij=0; ij<nlayer; ij++) { 
		dist[ij]= -100.; istrxtime[ij] = istrytime[ij] = -1;
		passxtime[ij] = passytime[ij] = false; 
		passxtmy[ij] = passytmx[ij] = true;
		
		xtime[ij] = rawxtime[ij] = rawxtime1[ij] = rawxtime2[ij] = -100;
		ytime[ij] = rawytime[ij] = rawytime1[ij] = rawytime2[ij] = -100;
		xusedtime[ij] = yusedtime[ij] = false;
		xtdev[ij] = 100;
		ytdev[ij] = 100;
	      }
	      
              xval=-100; yval=-100;
              int init=-1;

              for( int ij=0;ij<nlayer;ij++) {
		if(abs(Xdev[ij]) < 2.0 && abs(Ydev[ij]) < 2.0) {
		  
		  if (xpts[ij].size()==0 || xpts[ij].size() >nmxhits || ypts[ij].size()==0 || ypts[ij].size() >nmxhits) continue;
		  if (init<0) {
		    xval = xext[ij];
		    yval = yext[ij];
		    dist[ij] = 0.0;
		    init = ij;
		  } else {
		    dist[ij] = sqrt( pow((xext[ij] - xval)*stripwidth, 2.) + 
				     pow((yext[ij] - yval)*stripwidth, 2.) +
                                     pow((ij - init)*layergap, 2.));
		  }
		  
		  //calcualte strips where signal come earlier and the value of offset
		  double tshft=1000.0;
		  bool tpass=true;
		  passxtime[ij] = true;
		  passxtmy[ij] = true; // (yext[ij] > firstYstrip && yext[ij] < lastYstrip) ? true : false;
		  int istr = istrxtime[ij] = xpts[ij][0];
		  for (int ix=0; ix<xpts[ij].size(); ix++) {
		    if ((!timeinuse[0][ij][xpts[ij][ix]]) || 
			xpts[ij][ix]< firstXstrip || 
			xpts[ij][ix]> lastXstrip) tpass= passxtime[ij] = false;
		    
		    if (xtoffset[ij][xpts[ij][ix]]<tshft) { 
		      tshft =xtoffset[ij][xpts[ij][ix]]; istr = istrxtime[ij] = xpts[ij][ix];
		    }
		  }

		  rawxtime[ij] = timesx[ij];

		  timesx[ij] -= slope_path*yext[ij]; //(5./32.)*yext[ij];
		  
		  rawxtime1[ij] = timesx[ij];
		  timesx[ij] -=tshft;
		  rawxtime2[ij] = timesx[ij];
		  
		  istr = int(yext[ij]+0.5);
		  if (istr<0) istr=0;
		  if (istr>=nstrip) istr = nstrip-1;
		  double dx = yext[ij]-istr;

		  // Linear extrapolation using only two points
		  if ((istr==0 && dx<=0.0) || (istr==nstrip-1 && dx>=0.0)) { 
		    timesx[ij] -=xt_slope_cor[ij][istrxtime[ij]][istr];
		  } else if (dx>0) {
		    timesx[ij] -=(1-dx)*xt_slope_cor[ij][istrxtime[ij]][istr]+dx*xt_slope_cor[ij][istrxtime[ij]][istr+1]; 
		  } else {
		    timesx[ij] -=abs(dx)*xt_slope_cor[ij][istrxtime[ij]][istr-1]+(1-abs(dx))*xt_slope_cor[ij][istrxtime[ij]][istr]; 
		  }
		  
		  if ((xpts[ij].size() >0 && xpts[ij].size()<=nmxtimehit) && (usealltime || (ij==occulyr) || (tpass && passxtmy[ij]))) {
		    int isiz = max(0, min(nmxtimehit,int(xpts[ij].size()))-1);
		    xtime[ij] = timesx[ij] - parabolic(xposinstr[ij], strpos_vs_time[isiz][ij]);
		    xusedtime[ij] = (xpts[ij].size()<=nmxusedtimehit) ? true : false; //  xtime[ij] : -101;
		  } else {rawxtime[ij] = rawxtime1[ij] = rawxtime2[ij] = -90;}
		  
		  tshft=1000.0;
		  tpass=true;
		  passytime[ij] = true;
		  passytmx[ij] = true; //(xext[ij] > firstXstrip && xext[ij] < lastXstrip) ? true : false;
		  istr = istrytime[ij] = ypts[ij][0];
		  for (int iy=0; iy<ypts[ij].size(); iy++) {
		    if ((!timeinuse[1][ij][ypts[ij][iy]]) || 
			ypts[ij][iy]< firstYstrip || 
			ypts[ij][iy]> lastYstrip) tpass = passytime[ij] = false;
		    
		    if (ytoffset[ij][ypts[ij][iy]]<tshft) { 
		      tshft =ytoffset[ij][ypts[ij][iy]]; istr = istrytime[ij] = ypts[ij][iy];
		    }
		  }
		  
		  rawytime[ij] = timesy[ij];
		  timesy[ij] -=ytimeshift; // shift y-time to match with x-time

		  timesy[ij] -= 5. - slope_path*xext[ij]; //(5./32.)*xext[ij];
		  rawytime1[ij] = timesy[ij];
		  
		  timesy[ij] -=tshft;
		  rawytime2[ij] = timesy[ij];

		  istr = int(xext[ij]+0.5);
		  if (istr<0) istr=0;
		  if (istr>=nstrip) istr = nstrip-1;

		  double dy = xext[ij]-istr;

		  // Linear extrapolation using only two points
		  if ((istr==0 && dy<=0.0) || (istr==nstrip-1 && dy>=0.0)) { 
		    timesy[ij] -=yt_slope_cor[ij][istrytime[ij]][istr];
		  } else if (dy>0) {
		    timesy[ij] -=(1-dy)*yt_slope_cor[ij][istrytime[ij]][istr]+dy*yt_slope_cor[ij][istrytime[ij]][istr+1]; 
		  } else {
		    timesy[ij] -=abs(dy)*yt_slope_cor[ij][istrytime[ij]][istr-1]+(1-abs(dy))*yt_slope_cor[ij][istrytime[ij]][istr]; 
		  }
		  
		  if ((ypts[ij].size() >0 && ypts[ij].size()<=nmxtimehit) && (usealltime || (ij==occulyr) || (tpass && passytmx[ij]))) { 
		    int isiz = max(0, min(nmxtimehit,int(ypts[ij].size()))-1)+nmxtimehit;
		    ytime[ij] = timesy[ij] - parabolic(yposinstr[ij], strpos_vs_time[isiz][ij]);
		    yusedtime[ij] = (ypts[ij].size()<=nmxusedtimehit) ?  true : false; //ytime[ij] : -101;
		    
		    if (ypts[ij].size()>nmxusedtimehit && yusedtime[ij]) cout <<"ypts "<<ij<<" "<< ypts[ij].size()<<" "<<nmxusedtimehit<<" "<<int(yusedtime[ij])<<endl;
 		  } else { rawytime[ij] = rawytime1[ij] = rawytime2[ij] = -90;}
		} //if(abs(Xdev[ij]) < xyPosDev && abs(Ydev[ij]) < xyPosDev)
	      } //for(int ij=0;ij<nlayer;ij++)
	      
	      int tmpxtent[nlayer];
	      int tmpytent[nlayer];
	      
	      timexslope = 0;
              xc0inters = 0;
	      xt0chi2 = 1.e+6;
              nxtfail = 0;
	      
	      timeyslope = 0;
              yc0inters = 0;
	      yt0chi2 = 1.e+6;
              nytfail = 0;

	      for (int ix=0; ix<nlayer; ix++) { xtext[ix]= xtexter[ix] = 1000;}

	      // Xtime fit
	      int iTimeSlopeConst = (iiter<nmxiter-1 && ntcor==1) ? 0 : -1;
	      StraightLineFit xtimefit(iTimeSlopeConst, dist, xtime,  timeserrx2, xusedtime, occulyr, occulyr, xtcorstr, xtcorend, float(7.0));
	      xtimefit.GetParameters(nxtfail, xc0inters, timexslope);
	      //	    xtimefit.GetError(errcst, errlin, errcov);
	      xtimefit.GetChisqure(nxtime, xt0chi2);
	      xtimefit.GetFitValues(xtext, xtdev, xtexter);

              if (nxtime>=nmnhits && nxtfail==0) {//4Nov For resolution 8layers are considered
		if (occulyr >=nlayer) {
		  if (xt0chi2/(nxtime-2)<mxtimechisq) { 
		    dir_cx->Fill(-cval*timexslope);
		  }

                  for (int ij=0; ij<nlayer; ij++) {
		    if (dist[ij]<0 || xtime[ij] <-50) continue;
		    if (xusedtime[ij]) { 
		      if (istrxtime[ij]>=0 && istrxtime[ij] <nstrip && passxtmy[ij] ) { 
			time_xstrreso[ij][istrxtime[ij]][iiterrs]->Fill(xtdev[ij]); 
		      }
		      if (passxtime[ij] && passxtmy[ij]) {
			if (xusedtime[ij]) {time_xreso[ij][iiterrs]->Fill(xtdev[ij]);}  // Filling and saving. Not fitted.
		      }
		    }
		  }
		} else {
		  if (dist[occulyr] >=0 && xtime[occulyr] >-50.0) {
		    if (xusedtime[occulyr]) {
		      if (istrxtime[occulyr]>=0 && istrxtime[occulyr] <nstrip) { 
			if (passxtmy[occulyr]) {time_xstrreso[occulyr][istrxtime[occulyr]][iiterrs]->Fill(xtdev[occulyr]);}
		      }
		    }
		    if (passxtime[occulyr] && passxtmy[occulyr]) {
		      if (xusedtime[occulyr]) {time_xreso[occulyr][iiterrs]->Fill(xtdev[occulyr]);}
		    }
		  } //if (dist[occulyr] >=0 && xtime[occulyr] >-50.0)
		} //else of if (occulyr >=nlayer) 
	      } // if (nxtime>=nmnhits/*-ntcor*/ && nxtfail==0)
	      
	      for (int ix=0; ix<nlayer; ix++) { ytext[ix]= ytexter[ix] = 100;}
	      if (nxtfail==0 && isfill) {
		h_tchisqx->Fill(xt0chi2);
		if (nxtime-2>0) {
		  h_treduchisqx->Fill(xt0chi2/(nxtime-2));
		  double probx =TMath::Prob(xt0chi2, nxtime-2);
		  h_xtprob->Fill(probx);
		}
		h_txndf->Fill(nxtime);
	      }

	      int tmpntxy = 0;
	      StraightLineFit ytimefit(iTimeSlopeConst, dist, ytime, timeserry2, yusedtime, occulyr, occulyr, ytcorstr, ytcorend, float(7.0));
	      ytimefit.GetParameters(nytfail, yc0inters, timeyslope);
	      //	    ytimefit.GetError(errcst, errlin, errcov);
	      ytimefit.GetChisqure(nytime, yt0chi2);
	      ytimefit.GetFitValues(ytext, ytdev, ytexter);

	      //150101		  tmpytent[ij]= 0;
	      
	      if (nytime>=nmnhits && nytfail==0) {
		if (occulyr >=nlayer) {
		  if (yt0chi2/(nytime-2)<mxtimechisq) { 
		    dir_cy->Fill(-cval*timeyslope);
		  }
                  for (int ij=0; ij<nlayer; ij++) {
		    if (dist[ij]<0 || ytime[ij] <-50) continue;
		    if (yusedtime[ij]) { 
		      if (istrytime[ij]>=0 && istrytime[ij] <nstrip && passytmx[ij]) {
			time_ystrreso[ij][istrytime[ij]][iiterrs]->Fill(ytdev[ij]);
		      }
		    }
		    if (passytime[ij] && passytmx[ij]) {
		      if (yusedtime[ij]) {time_yreso[ij][iiterrs]->Fill(ytdev[ij]);}
		    }
		  }
		} else {
		  if (yusedtime[occulyr]) { 
		      if (istrytime[occulyr]>=0 && istrytime[occulyr] <nstrip) {
			if (passytmx[occulyr]) {time_ystrreso[occulyr][istrytime[occulyr]][iiterrs]->Fill(ytdev[occulyr]);}
		      }
		  }
		  if (passytime[occulyr] && passytmx[occulyr]) {
		    if (yusedtime[occulyr]) {time_yreso[occulyr][iiterrs]->Fill(ytdev[occulyr]);}
		  }
		} // else of if (occulyr >=nlayer) 
              } //if (nytime>=nmnhits/*-ntcor*/ && nytfail==0)


	      
	      if (nytfail==0 && isfill) {
		h_tchisqy->Fill(yt0chi2);
		if (nytime-2>0) {
		  h_treduchisqy->Fill(yt0chi2/(nytime-2));
		  double probx =TMath::Prob(yt0chi2, nytime-2);
		  h_ytprob->Fill(probx);
		}
		h_tyndf->Fill(nytime);
	      }





	    } // if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<2.0 && ychi2/(Ny-2)<2.0 && nxfail==0 && nyfail==0); Position fit cut
	    if (isTiming && isfill &&  nxfail==0 &&  nyfail==0 && nxtime>2 && nytime>2) {
	      txslop = -1./timexslope/cval;
	      tyslop = -1./timeyslope/cval;
	      
	      if (nxtime>10 && nytime>10) T1->Fill(); //07/09/2011
            }
	    
	    if((ievt%10000)==0) cout<<"Processed " <<event->ENum<<" events so far for iteration# "<< iiter<<" occu "<<occulyr<<endl;
	  }  // EVENT loop
	  fileIn->cd();
	  delete event_tree;
	  delete fileIn;
	} // while(!(file_db.eof()))
	
	file_db.close();
	fileOut->cd();
	
	if (firstiter) { 
	  double scale  = 100./max(1, nTotalp);
	  
	  gStyle->SetOptLogy(0);
	  gStyle->SetPadLeftMargin(0.12);
	  gStyle->SetPadBottomMargin(0.14);
	  gStyle->SetOptStat(0);
	  ps.NewPage();
	  
	  TCanvas *c1=new TCanvas ("c1","Occupancy",500,700);
	  c1->Divide(3,4);
	  for (int ij=0; ij<nlayer; ij++) {
	    c1->cd(ij+1);
	    occu_x[ij]->Scale(scale);
	    occu_x[ij]->GetXaxis()->SetTitle(occu_x[ij]->GetTitle());
	    occu_x[ij]->GetXaxis()->SetTitleSize(0.07);
	    occu_x[ij]->GetXaxis()->CenterTitle();
	    occu_x[ij]->Draw();
	  }
	  c1->Update();

	  ps.NewPage();
	  for (int ij=0; ij<nlayer; ij++) {
	    c1->cd(ij+1);
	    occu_y[ij]->Scale(scale);
	    occu_y[ij]->GetXaxis()->SetTitle(occu_y[ij]->GetTitle());
	    occu_y[ij]->GetXaxis()->SetTitleSize(0.07);
	    occu_y[ij]->GetXaxis()->CenterTitle();
	    occu_y[ij]->Draw();
	  }
	  c1->Update();
	  if (c1) { delete c1; c1=0;}
	  
	  gStyle->SetOptLogz(1);
	  gStyle->SetPadLeftMargin(0.14);
	  gStyle->SetPadRightMargin(0.16);
	  
	  ps.NewPage();
	  TCanvas *c2=new TCanvas ("c2","Occupancy",500,700);
	  c2->Divide(3,4);
	  for (int ij=0; ij<nlayer; ij++) {
	    c2->cd(ij+1);
	    //  gPad->SetLogz(1);
	    raw_occu[ij]->Scale(scale);
	    raw_occu[ij]->GetXaxis()->SetTitle("XStrip #");
	    raw_occu[ij]->GetXaxis()->SetTitleSize(0.07);
	    raw_occu[ij]->GetXaxis()->CenterTitle();
	    raw_occu[ij]->GetYaxis()->SetTitle("YStrip #");
	    raw_occu[ij]->GetYaxis()->SetTitleSize(0.07);
	    raw_occu[ij]->GetYaxis()->CenterTitle();

	    raw_occu[ij]->Draw("colz");
	  }
	  c2->Update();
	  
	  if (c2) { delete c2; c2=0;}
	  
	  ps.NewPage();
	  gStyle->SetOptStat(1100);
	  gStyle->SetOptTitle(1);
	  gStyle->SetOptLogy(1);
	  gStyle->SetStatW(.40); //40);
	  gStyle->SetStatH(.20); //30);
	  gStyle->SetTitleFontSize(0.07);
	  gStyle->SetPadLeftMargin(0.09);
	  gStyle->SetPadBottomMargin(0.06);
	  gStyle->SetPadTopMargin(0.09); //(0.03);
	  gStyle->SetPadRightMargin(0.01);
	  
	  TCanvas* c4c = new TCanvas("c4c", "c4c", 700, 900);
	  c4c->Divide(3,4);
	  gStyle->SetStatY(.99); gStyle->SetStatTextColor(3);
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    xlayer_mult[ij]->Scale(1./nTotalp);
	    xlayer_mult[ij]->SetLineColor(3); xlayer_mult[ij]->Draw();
	  }
	  c4c->Update();
	  gStyle->SetStatY(.79); gStyle->SetStatTextColor(4);
	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    ylayer_mult[ij]->Scale(1./nTotalp);
	    ylayer_mult[ij]->SetLineColor(4); ylayer_mult[ij]->Draw("sames");
	  }
	  
	  c4c->Update();
	  

	  ps.NewPage(); 
	  gStyle->SetOptStat(1100);
	  gStyle->SetStatY(.99);
	  gStyle->SetStatH(.20);
	  gStyle->SetPadLeftMargin(0.15);
	  ps.NewPage(); 
	  // gStyle->SetPadLeftMargin(0.15);
	  c4c->cd(1);gPad->SetLeftMargin(0.11); h_xprob->Draw();
	  c4c->cd(2); gPad->SetLeftMargin(0.11);h_reduchisqx->Draw();
	  c4c->cd(3); gPad->SetLeftMargin(0.11);h_xndf->Draw();
	  
	  c4c->cd(4); gPad->SetLeftMargin(0.11);h_yprob->Draw();
	  c4c->cd(5); gPad->SetLeftMargin(0.11);h_reduchisqy->Draw();
	  c4c->cd(6); gPad->SetLeftMargin(0.11);h_yndf->Draw();
	  
	  if (isTiming) {
	    c4c->cd(7); gPad->SetLeftMargin(0.11);h_xtprob->Draw();
	    c4c->cd(8); gPad->SetLeftMargin(0.11); h_treduchisqx->Draw();
	    c4c->cd(9); gPad->SetLeftMargin(0.11); h_txndf->Draw();
	    
	    c4c->cd(10); gPad->SetLeftMargin(0.11); h_ytprob->Draw();
	    c4c->cd(11);gPad->SetLeftMargin(0.11); h_treduchisqy->Draw();
	    c4c->cd(12); gPad->SetLeftMargin(0.11);  h_tyndf->Draw();
	  }
	  c4c->Update();

	  ps.NewPage();
	  gStyle->SetOptStat(111110);
	  gStyle->SetOptFit(101);
	  gStyle->SetOptLogy(0);
	  
	  gStyle->SetPadBottomMargin(0.11);
	  gStyle->SetPadTopMargin(0.09);
	  gStyle->SetPadLeftMargin(0.07);
	  gStyle->SetPadRightMargin(0.02);
	  gStyle->SetOptTitle(1);
	  gStyle->SetTitleFontSize(0.07);
	  gStyle->SetStatW(.30);
	  gStyle->SetStatH(.14);
	  gStyle->SetStatY(.99);
	  gStyle->SetStatX(.99);
	  
	  TCanvas *c2a=new TCanvas ("c2a","Slope",500,700);
	  c2a->Divide(2,3);
	  const int nxplot=6;
	  //	latex.SetTextSize(0.08);
	  TH1F* histax[nxplot]={0};
	  
	  double ncx = 0;
	  double ncy = 0;
	  for (int ij=0; ij<dir_cx->GetNbinsX(); ij++) {
	    if (dir_cx->GetBinCenter(ij)<0.0) {
	      ncx +=dir_cx->GetBinContent(ij);
	      ncy +=dir_cy->GetBinContent(ij);
	    } else { break;}
	  }

	  for (int ix=0; ix<nxplot; ix++) {
	    switch(ix) {
	    case 0 : histax[ix] = (TH1F*)pos_xslope->Clone(); break;
	    case 1 : histax[ix] = (TH1F*)pos_yslope->Clone(); break;
	    case 2 : histax[ix] = (TH1F*)pos_theta->Clone(); break;
	    case 3 : histax[ix] = (TH1F*)pos_phi->Clone(); break;
	    case 4 : histax[ix] = (TH1F*)dir_cx->Clone(); break;
	    case 5 : histax[ix] = (TH1F*)dir_cy->Clone(); break;
	    default : histax[ix] = (TH1F*)pos_xslope->Clone(); break;
	    }
	    c2a->cd(ix+1);
	    histax[ix]->GetXaxis()->CenterTitle();
	    histax[ix]->GetXaxis()->SetTitleOffset(0.7);
	    histax[ix]->GetXaxis()->SetTitleSize(0.06);
	    histax[ix]->GetXaxis()->SetLabelOffset(-0.01);
	    histax[ix]->GetXaxis()->SetLabelSize(0.06);
	    histax[ix]->Draw();
	    if (ix==4) {latex.DrawLatex(0.12, 0.66,Form("#scale[0.6]{%g\%}", int(1000000*ncx/max(1.,dir_cx->GetEntries()))/10000.));}
	    if (ix==5) {latex.DrawLatex(0.12, 0.66,Form("#scale[0.6]{%g\%}", int(1000000*ncy/max(1.,dir_cy->GetEntries()))/10000.));}
	    
	  }
	  c2a->Update();
	  if (c2a) { delete c2a; c2a=0;}
	  for (int ix=0; ix<nxplot; ix++) {
	    if (histax[ix]) { delete histax[ix]; histax[ix]=0;}
	  }

	  firstiter = false;
	}
	
	int istr = (occulyr >= nlayer) ? 0 : occulyr;
	int iend = (occulyr >= nlayer) ? nlayer : occulyr+1;

        for (int iocc = istr; iocc<iend; iocc++) {

	  double absWidth=min(2.5,max(1.9, 0.15*(nmxiter -iiter)));
	  double widthScale=min(0.70,max(0.40, 0.05*(nmxiter -iiter))); //Initially one can allow very large asymmetric distribution and and gradually make it smaller
	  
	  ps.NewPage();
	  gStyle->SetOptStat(1110);
	  gStyle->SetOptFit(101);
	  gStyle->SetOptLogy(1);
	  
	  gStyle->SetPadBottomMargin(0.11);
	  gStyle->SetPadTopMargin(0.08);
	  gStyle->SetPadLeftMargin(0.07);
	  gStyle->SetPadRightMargin(0.02);
	  gStyle->SetOptTitle(1);
	  gStyle->SetTitleFontSize(0.07);
	  gStyle->SetStatW(.20);
	  gStyle->SetStatH(.12);
	  gStyle->SetStatY(.99);
	  gStyle->SetStatX(.99);
	  

	  ps.NewPage();
	  TCanvas *c2=new TCanvas ("c2","Residue",500,700);
	  c2->Divide(2,2);
	  
	  c2->cd(1);
	  
	  double alw= (ntcor==0) ? -2.05 : -2.8; //xlayer_reso[iocc][iiterrs]->GetXaxis()->GetXmin();
	  double ahg= (ntcor==0) ?  2.05 :  2.8; //xlayer_reso[iocc][iiterrs]->GetXaxis()->GetXmax();
	  TF1* fitfx = new TF1("fitfx", gausX, alw, ahg, 3);

	  double  parx[3]={xlayer_reso[iocc][iiterrs]->GetMaximum(), xlayer_reso[iocc][iiterrs]->GetMean(), max(0.15,0.7*xlayer_reso[iocc][iiterrs]->GetRMS())};
	  fitfx->SetParameters(parx);
	  //	    fitfx->SetParLimits(0, 0.17*parx[0], 1.5*parx[0]);
	  fitfx->SetParLimits(1, parx[1]-1., parx[1]+1.);
	  fitfx->SetParLimits(2, 0.12, 1.5*parx[2]);
	  
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitle("X-residues (pitch)");
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->CenterTitle();
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
	  
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);

          xlayer_reso[iocc][iiterrs]->Fit(fitfx, "WW:BMRQ");
	  
	  c2->cd(2);
	  
	  TF1* fitfy = new TF1("fitfy", gausX, alw, ahg, 3);
	  double  pary[3]={ylayer_reso[iocc][iiterrs]->GetMaximum(), ylayer_reso[iocc][iiterrs]->GetMean(), max(0.15,0.7*ylayer_reso[iocc][iiterrs]->GetRMS())};
	  fitfy->SetParameters(pary);
	  //	    fitfy->SetParLimits(0, 0.15*pary[0], 1.5*pary[0]);
	  fitfy->SetParLimits(1, pary[1]-1., pary[1]+1.);
	  fitfy->SetParLimits(2, 0.60, 1.5*pary[2]);
	  
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitle("Y-residues (pitch)");
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->CenterTitle();
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
	  
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
          ylayer_reso[iocc][iiterrs]->Fit(fitfy, "WW:BMRQ");
	  
	  //time resolution only for plot, but not for any calculation
	  alw= (ntcor==0) ? -6.0 :-9.0;
	  ahg= (ntcor==0) ?  6.0 : 9.0;
	  c2->cd(3);
	  TF1* fittfx = new TF1("fittfx", gausX, alw, ahg, 3);

	  double  partx[3]={time_xreso[iocc][iiterrs]->GetMaximum(), time_xreso[iocc][iiterrs]->GetMean(), max(0.15,0.7*time_xreso[iocc][iiterrs]->GetRMS())};
	  fittfx->SetParameters(partx);
	  //	    fittfx->SetParLimits(0, 0.17*partx[0], 1.5*partx[0]);
	  fittfx->SetParLimits(1, partx[1]-1., partx[1]+1.);
	  fittfx->SetParLimits(2, 0.6, 1.5*partx[2]);
	  
	  time_xreso[iocc][iiterrs]->GetXaxis()->SetTitle("X - #Deltat (ns)");
	  time_xreso[iocc][iiterrs]->GetXaxis()->CenterTitle();

	  time_xreso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	  time_xreso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
	  time_xreso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	  time_xreso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);

          time_xreso[iocc][iiterrs]->Fit(fittfx, "WW:BMRQ");
	  
	  c2->cd(4);
	  TF1* fittfy = new TF1("fittfy", gausX, alw, ahg, 3);

	  double  party[3]={time_yreso[iocc][iiterrs]->GetMaximum(), time_yreso[iocc][iiterrs]->GetMean(), max(0.15,0.7*time_yreso[iocc][iiterrs]->GetRMS())};
	  fittfy->SetParameters(party);
	  //	    fittfy->SetParLimits(0, 0.17*party[0], 1.5*party[0]);
	  fittfy->SetParLimits(1, party[1]-1., party[1]+1.);
	  fittfy->SetParLimits(2, 0.6, 1.5*party[2]);
	  
	  time_yreso[iocc][iiterrs]->GetXaxis()->SetTitle("Y - #Deltat (ns)");
	  time_yreso[iocc][iiterrs]->GetXaxis()->CenterTitle();
	  time_yreso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	  time_yreso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
	  time_yreso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	  time_yreso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);

          time_yreso[iocc][iiterrs]->Fit(fittfy, "WW:BMRQ");
	  
	  c2->Update();
	  if (c2) { delete c2; c2=0;}
	  
	  file_out <<"lay "<< iocc<<" it "<< iiterrs <<" X-sh "<< fitfx->GetParameter(1)<<" "<<fitfx->GetParameter(2)<<" Y-sh "<< fitfy->GetParameter(1)<<" "<<fitfy->GetParameter(2)<<" X-sh "<< fittfx->GetParameter(1)<<" "<<fittfx->GetParameter(2)<<" Y-sh "<< fittfy->GetParameter(1)<<" "<<fittfy->GetParameter(2)<<endl; 
	  
	  if (iiter<nmxiter-1 && ntcor==1) {
	    if (fabs(xlayer_reso[iocc][iiterrs]->GetMean()-
		     fitfx->GetParameter(1)) 
		< xlayer_reso[iocc][iiterrs]->GetRMS()) {
	      xoff[iocc] += fitfx->GetParameter(1);
	    }
	    
	    if (fabs(ylayer_reso[iocc][iiterrs]->GetMean()-
		     fitfy->GetParameter(1)) 
		< ylayer_reso[iocc][iiterrs]->GetRMS()) {
	      yoff[iocc] += fitfy->GetParameter(1);
	    }
	  }
	  
	  xrms[iocc] = fitfx->GetParameter(2);
	  yrms[iocc] = fitfy->GetParameter(2);	
	  
	  delete fitfx; fitfx=0;
	  delete fitfy; fitfy=0;
	  
	  if (isTiming) { 
	    ps.NewPage();
	    gStyle->SetOptTitle(0);
	    gStyle->SetOptStat(0);
	    gStyle->SetOptFit(0);
	    gStyle->SetOptLogy(1);
	    gStyle->SetPadTopMargin(.001);
	    gStyle->SetPadBottomMargin(0.001); //0.07
	    gStyle->SetPadLeftMargin(0.001);
	    gStyle->SetPadRightMargin(0.001);
	    TCanvas* c1 = new TCanvas("c1", "c1", 700, 900);
	    c1->Divide(8,8,1.e-6, 1.e-6);

	    TH1F* time_shift[nstrip][2]={0};
	    double fitmean[nstrip][2]={0};
	    double fitrms[nstrip][2]={0};
	    double fitchi[nstrip][2]={0};
	    double statmean[nstrip][2]={0};
	  
	    TF1* fity[nstrip][2]={0}; 
	  
	    const int nsgpr=3;
	    double fitres[nsgpr];
	    double parerr[nsgpr];
	    double fchisq;
            for (int ij=0; ij<2; ij++) {
              for (int jk=0; jk<nstrip; jk++) {
	      
		c1->cd(nstrip*ij+jk+1);
		if (ij==0) {
		  time_shift[jk][ij] = (TH1F*)time_xstrreso[iocc][jk][iiterrs]->Clone();
		} else {
		  time_shift[jk][ij] = (TH1F*)time_ystrreso[iocc][jk][iiterrs]->Clone();
		} 
		
		if (time_shift[jk][ij]->Integral()>2) {
		  time_shift[jk][ij]->GetXaxis()->SetLabelSize(.07);
		  
		  double  par[3]={time_shift[jk][ij]->GetMaximum(),time_shift[jk][ij]->GetMean(), time_shift[jk][ij]->GetRMS()};
		  
		  int nbinx = time_shift[jk][ij]->GetNbinsX();
		  int nlow=0;
                  for (int kl=0; kl<nbinx; kl++) {
		    if (time_shift[jk][ij]->GetBinContent(kl+1) >0) {
		      nlow=kl; break;
		    }
		  }
		  int nhig = 0;
                  for (int kl=nbinx; kl>0; kl--) {
		    if (time_shift[jk][ij]->GetBinContent(kl+1) >0) {
		      nhig=kl; break;
		    }
		  }
		  float amean = 0.5*(nlow + nhig);
		  nlow = int(TMath::Max(0., nlow - 1.0*(amean - nlow)));
		  nhig = int(TMath::Min(nbinx-1., nhig + 1.0*(nhig - amean)));
		  
		  nchannel = 0;
                  for (int kl=nlow; kl<=nhig; kl++) {
		    if (nchannel <nmxchn) {
		      m_data[nchannel] = time_shift[jk][ij]->GetBinContent(kl+1);
		      m_xpos[nchannel] = time_shift[jk][ij]->GetBinCenter(kl+1);
		      nchannel++;
		    }
		  }
		  double alw= time_shift[jk][ij]->GetBinCenter(nlow+1);
		  double ahg = time_shift[jk][ij]->GetBinCenter(nhig+1);
		  
		  //		  time_shift[jk][ij]->GetXaxis()->SetRangeUser(alw, ahg);
		  //		  time_shift[jk][ij]->GetXaxis()->SetLabelSize(.1);
		  time_shift[jk][ij]->Draw();

		  TMinuit *gMinuit = new TMinuit(nsgpr);
		  gMinuit->SetPrintLevel(-1);

		  TString hname[nsgpr] = {"height", "mean", "rms"};
		  
		  int nmx = time_shift[jk][ij]->GetMaximumBin();
		  double hgh = 0.35*(time_shift[jk][ij]->GetBinContent(nmx-1) + 
				     time_shift[jk][ij]->GetBinContent(nmx) + 
				     time_shift[jk][ij]->GetBinContent(nmx+1));
		  double strt[nsgpr] = {hgh,time_shift[jk][ij]->GetMean(), max(0.6, min(3.0, 0.9*time_shift[jk][ij]->GetRMS()))};
		  
		  double alow[nsgpr] = {0.5*strt[0], strt[1]-2.0, max(0.7*strt[2],0.5)};
		  double ahig[nsgpr] = {2.0*strt[0], strt[1]+2.0, min(1.1*strt[2]+0.1,3.5)};
		  double step[nsgpr] = {0.5, 0.01, 0.01};
		  
		  gMinuit->SetFCN(fcnsg);
		  
		  double arglist[10];
		  int ierflg = 0;
		  arglist[0] =  1 ;
		  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
		  
                  for (int kl=0; kl<nsgpr; kl++) {
		    gMinuit->mnparm(kl, hname[kl], strt[kl], step[kl], alow[kl], ahig[kl],ierflg);
		  }
		  
		  arglist[0] = 0;
		  //	      gMinuit->mnexcm("MIGRAD", arglist, 0, ierflg);
		  gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);
		  
		  arglist[0] = 0;
		  gMinuit->mnexcm("IMPROVE", arglist, 0, ierflg);
		  
		  TString chnam;
		  double parv,err,xlo,xup, plerr, mierr, eparab, gcc;
		  int iuit;
		  
                  for (int kl=0; kl<nsgpr; kl++) {
		    gMinuit->mnpout(kl, chnam, parv, err, xlo, xup, iuit);
		    gMinuit->mnerrs(kl, plerr, mierr, eparab, gcc);
		    fitres[kl] = parv;
		    parerr[kl] = err;
		  }
		  double  fedm, errdef;
		  int  nparx, istat, fitndof;
		  gMinuit->mnstat(fchisq, fedm, errdef, fitndof, nparx, istat);
		  
		  //		  cout <<"chisq "<< fchisq<<" "<<fitndof<<" "<<nparx<<" "<<istat<<" "<<fitres[0]<<" "<<fitres[1]<<" "<<fitres[2]<<" "<<time_shift[jk][ij]->GetMean()<<" "<<time_shift[jk][ij]->GetRMS()<<" "<<time_shift[jk][ij]->GetEntries()<<endl;
		  if (istat==0) fchisq =100000000.0;
		  //		  time_shift[jk][ij]->Fit("gaus");

		  sprintf(name, "fity_%i_%i", ij, jk);
		  //	fity[jk][ij] = new TF1(name, gausX, alw, ahg,  3);
		  //	fity[jk][ij]->SetParameters(par);
		  //	time_shift[jk][ij]->Fit(name);
		  
		  //	time_shift[jk][ij]->Draw();
		  
		  fity[jk][ij] = new TF1(name, gausX, alw, ahg, 3);
		  fity[jk][ij]->SetParameters(fitres);
		  fity[jk][ij]->SetLineColor(2);
		  fity[jk][ij]->SetLineWidth(1);
		  fity[jk][ij]->Draw("same");

		  fitmean[jk][ij] = fitres[1]; // fity[jk][ij]->GetParameter(1);
		  fitrms[jk][ij] = fitres[2];
		  fitchi[jk][ij] = fchisq;

		  statmean[jk][ij] = time_shift[jk][ij]->GetMean();

		  latex.DrawLatex(0.32, 0.36,Form("%g", int(1000*fitres[1])/1000.));
		  latex.DrawLatex(0.32, 0.26,Form("%g", int(1000*statmean[jk][ij])/1000.));
		  latex.DrawLatex(0.32, 0.16,Form("%g", int(1000*fitres[2])/1000.));
		  latex.DrawLatex(0.32, 0.06,Form("%g", int(1000*time_shift[jk][ij]->GetRMS())/1000.));
		  delete gMinuit; gMinuit=0;
		} // if (time_shift[jk][ij]->GetEntries()>15)
	      } // for (int jk=0; jk<nstrip; jk++)
	    } //for (int ij=0; ij<2; ij++) 
	    c1->Update();
	    if (c1) { delete c1; c1=0;} 
	    
            for (int ij=0; ij<2; ij++) {
              for (int jk=0; jk<nstrip; jk++) {
		if (time_shift[jk][ij]) { delete time_shift[jk][ij]; time_shift[jk][ij]=0;}
		if (fity[jk][ij]) { delete fity[jk][ij]; fity[jk][ij]=0;}
	      }
	    }
	    file_out <<"lay3 "<< iiter<<" "<<ntcor<<" "<<iiterrs<<" "<<lstr<<" "<<lend<<" "<<laye<<" "<<occulyr<<" "<<iocc<<" "<<nlayer<<endl;
	  } //if (isTiming)
	} // for (int iocc = istr; iocc<iend; iocc++)
      } // for (int laye=0; laye<nlayerit; laye++) 



    } // or (int ntcor=0; ntcor< (iiter<nmxiter-1) ? 2 : 1; ntcor++)

  } // for (int iiter=0; iiter<nmxiter; iiter++)





  ps.Close();
  fileOut->cd();
  fileOut->Write();
  file_out.close();
  file_outstr.close();
  fileOut->Close();

  return 0;
  //--------------------------------------------------------------
}


/*

scp INORUN_20110514_0912.ine INORUN_20110514_0757.ine INORUN_20110514_0641.ine INORUN_20110514_0526.ine INORUN_20110514_0412.ine INORUN_20110514_0255.ine INORUN_20110514_0141.ine INORUN_20110514_0257.ine INORUN_20110514_0024.ine INORUN_20110514_0026.ine INORUN_20110513_2309.ine INORUN_20110513_2154.ine INORUN_20110513_2035.ine INORUN_20110513_2040.ine INORUN_20110513_1809.ine INORUN_20110513_1806.ine INORUN_20110513_1921.ine INORUN_20110514_0917.ine 


*/



