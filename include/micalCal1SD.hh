//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B2TrackerSD.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file B2TrackerSD.hh
/// \brief Definition of the B2TrackerSD class

#ifndef micalcal1SD_h
#define micalcal1SD_h 1
#include "G4Box.hh"
#include "G4VSolid.hh"
#include "G4GenericMessenger.hh"
#include "micalCal1SDMessenger.hh"
#include "G4VSensitiveDetector.hh"
#include "micalDetectorParameterDef.hh"
#include "MultiSimAnalysis.hh"
#include "G4SystemOfUnits.hh"
#include "micalCal1Hit.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include <vector>
#include<iostream>
#include "vect_manager.h"
using namespace std;



class G4Step;
class G4HCofThisEvent;
class G4Box;
class micalcal1SDMessenger;


class G4VSolid;

//class B1RunAction;//..
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// B1Tracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero 
/// energy deposit.

class micalcal1SD : public G4VSensitiveDetector
{
  public:
    micalcal1SD(const G4String name);
     ~micalcal1SD();
 
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

//unsigned long detid;
// unsigned int ScntStrpid;

//	void SetStripWidth(G4int val){SWidth = val;}
//	G4int GetStripWidth() const { return SWidth; }
	
	void SetPhotonSpeed(G4double val);
	void SetCMVadctons(G4double val);

        
  private:
  	micalcal1HitsCollection* cal1Collection;
	micalDetectorParameterDef* paradef;
	micalcal1SDMessenger* cal1SDMessenger;
	MultiSimAnalysis *pAnalysis;
	

	
	const G4int numberInCell;
	G4int Counter;
	G4int InCell;
     
 // G4Box* fEnvelopeBox;
	double partopscint[3];     //halflength of scintillators of thickness 1cm
	float AirGapScintTop;
	G4double Phys_TopScint_GPos[4][3];
	G4double Phys_SideScint_R_GPos[3][3];
	G4double Phys_SideScint_L_GPos[3][3];
	G4double Phys_SideScint_D_GPos[3][3];
	double  PhyVolGlPos[4][4][3];
  int NoScntStrpTop;
  int NoScntStrpSide;
  int Sipm_Pedestal;
  int Cmv_Threshold;
  
	unsigned int CellDetID[20000]; //16bit
	

	CmvStrip_Manager* CmvStrip_pointer;
	SipmHit_Manager* SipmHit_pointer;
 void DefineCommands();
// G4GenericMessenger* fMessenger;
 // static const B1RunAction* fRunAction;//..

  bool debug = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
