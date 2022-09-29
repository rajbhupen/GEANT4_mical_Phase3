#include "micalPrimaryGeneratorAction.hh"

#include "micalDetectorConstruction.hh"
#include "micalPrimaryGeneratorMessenger.hh"
#include "micalDetectorParameterDef.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4Box.hh"

#include "math.h"
#include "CLHEP/Random/RandGauss.h"

//using namespace std;
//#include "micalDetectorParameterDef.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
micalPrimaryGeneratorAction *micalPrimaryGeneratorAction::AnPointer;
micalPrimaryGeneratorAction::micalPrimaryGeneratorAction(
							 micalDetectorConstruction* micalDC, MultiSimAnalysis *panalysis)
  :micalDetector(micalDC), pAnalysis(panalysis) { 
  AnPointer =this;
  // G4cout<<" Initialized micalPrimaryGeneratorAction Constructor"<<endl;
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  //Default settings :
  SetRunNumber(0);
  SetInputFlag(0);
  SetFirstEvt(1);
  SetRndmFlag("on");
  SetPartId(13);
  SetIncEnergy(2.0*GeV);
  SetIncEnergySmr(0*MeV);
  SetIncDirection(G4ThreeVector(0.0,0.0,-1.0));
  SetIncThetaSmr(0*mrad);
  SetIncPhiSmr(0*mrad);
  SetIncPosition(G4ThreeVector(0.0*cm,0.0*cm,0*cm));
  SetIncVxSmr(0*cm);
  SetIncVySmr(0*cm);
  SetIncVzSmr(0*cm);
  
  enerin[0] = 0.3;
  enerin[1] = 0.4;
  enerin[2] = 0.5;
  enerin[3] = 0.6;
  enerin[4] = 0.7;
  enerin[5] = 0.8;
  enerin[6] = 0.9;
  enerin[7] = 1.0;
  enerin[8] = 1.1;
  enerin[9] = 1.2;
  enerin[10] = 1.3;
  enerin[11] = 1.4;
  enerin[12] = 1.5;
  enerin[13] = 1.8;
  enerin[14] = 2.0;
  enerin[15] = 2.5;
  enerin[16] = 3.0;
  enerin[17] = 4.0;
  enerin[18] = 5.0;
  enerin[19] = 6.0;
  
  pivalGA = acos(-1.0);
  initialise = 0;
  initialiseCor = 0;
  g_nevt=-1;
  
  //create a messenger for this class
  gunMessenger = new micalPrimaryGeneratorMessenger(this);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalPrimaryGeneratorAction::~micalPrimaryGeneratorAction() {
  // G4cout<<"micalPrimaryGeneratorAction Distructor"<<endl;
  if (particleGun)	{delete particleGun;}
  if (gunMessenger) {delete gunMessenger;}
}

void micalPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  // cout<<"micalPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {..."<<endl;
  // int nevt, npart, ipart;
  // G4int      pid;
  G4double   vx=0, vy=0, vz=0;
  // G4double   px=0, py=0, pz=0;
  // G4double   atime;
  
  if (initialise==0) {
    // gunMessenger = new micalPrimaryGeneratorMessenger(this);
    paradef = micalDetectorParameterDef::AnPointer;
  
  ShiftInX = paradef->GetShiftInX();
  ShiftInY = paradef->GetShiftInY();
ShiftInZ = paradef->GetShiftInZ(0);
  for(int ij=0; ij<3; ij++) {

      parino[ij] = paradef->GetParino(ij);
      pargas[ij] = paradef->GetPargas(ij);
           parchm[ij] = paradef->GetParchm(ij);
      StackPosInWorld[ij] = paradef->GetStackPosInRoom(ij) + paradef->GetINOroomPos(ij);
      cout<<"GetStackPosInRoom "<<ij<<" "<<paradef->GetStackPosInRoom(ij)<<" "<<paradef->GetINOroomPos(ij)<<endl;
    }
    cout<<"StackPosInWorld[3] "<<StackPosInWorld[0]<<" "<<StackPosInWorld[1]<<" "<<StackPosInWorld[2]<<endl;
    
    WorldXDim = paradef->GetParworld(0);
    WorldYDim = paradef->GetParworld(1);
    WorldZDim = paradef->GetParworld(2);
    nINODet = paradef->GetNumino();
    nLayer=paradef->GetnLayer();
    nIRLayer=paradef->GetnIRLayer();
    //cmv

       
  int jmax;
  for(int i =0;i<4;i++){
   
    jmax =  (i==0)? 4:3;
 
    for(int j=0;j<jmax;j++){
      for(int k=0;k<3;k++){
	ScintlayGlPos[i][j][k] = paradef->ScintLayGPos[i][j][k];

      }//k
      // cout<<endl;
    }//j
    // cout<<endl<<endl;
  }//i


  for(int op=0; op<3;op++) {
    partopscint[op] = paradef->partopscint[op];
  }
	
  AirGapScintTop= paradef->AirGapScintTop;
	







  
 NoScintStrip = 88;
 layhalflength = 0.5*((NoScintStrip*2*partopscint[0])+((NoScintStrip+1)*AirGapScintTop));//2289mm for top and 1041mm for sides

    //cmv
    for(int ij=0; ij<nLayer; ij++) {
      RPCLayerPosZ[ij] = paradef->GetRPCLayerPosZ(ij);
      cout<<"rpclayerposz: "<<ij<<" "<<StackPosInWorld[2]+ShiftInZ+RPCLayerPosZ[ij]<<" "<<RPCLayerPosZ[ij]<<endl;
      LayerZdim[ij] = paradef->GetLayerZdim(ij);
    }
    for(int ij=0; ij<nIRLayer; ij++) {
      IRONLayerPosZ[ij] = paradef->GetIRONLayerPosZ(ij);
      IronLayerZdim[ij] = paradef->GetIronLayerZdim(ij);
    }
    g_nevt=-1;
    initialise = 1;
    FirstEvt=pAnalysis->FirstEvt;
  }

  //this function is called at the begining of event
  // default particle kinematic

	pAnalysis->irun = RunNumber;  //Keep an option that in a file, one may have more than two run number
  g_nevt++;
  cout<<"g_nevt "<<g_nevt<<" "<<pAnalysis->InputOutput<<endl;
  if (pAnalysis->InputOutput<3) {

  //inputFlag = 0 for G4MC & =0 for Nuance 2 :GINIE
  if (InputFlag==0) { // G4MC case
    pAnalysis->ievt=g_nevt;
    
		cout<<"g_nevt"<<g_nevt<<endl;
    pAnalysis->ngent = ((unsigned)particleGun->GetNumberOfParticles() <=pAnalysis->ngenmx) ?
      particleGun->GetNumberOfParticles() : pAnalysis->ngenmx;
    
    for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++)	{
      //Option to have multiple particle
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      if (G4UniformRand()>0.5) partId *= -1;
      G4ParticleDefinition* particle = particleTable->FindParticle(partId);
      particleGun->SetParticleDefinition(particle);
      
       G4ThreeVector ini_Dir(incDirection);
      // G4ThreeVector ini_Dir(0,0,1);
      // ini_Dir.setTheta(2.792);
      // ini_Dir.setPhi(1.0472);
       // cout<<"incDirection "<<ini_Dir<<endl;
      double in_Energy = incEnergy*MeV;
			vx = incPosition.x()*mm;
			vy = incPosition.y()*mm;
			vz = incPosition.z()*mm;
	        
			//			cout<<"theta "<<ini_Dir.theta()<<" phi  "<<ini_Dir.phi()<<endl;
				
			
			
      if (rndmFlag=="on") {
				G4double theta=0;
				if(incThetaSmr>=0    && incDirection.theta() !=0) {
					theta = G4RandGauss::shoot(0,incThetaSmr);
				} else if(incThetaSmr<0 || incDirection.theta()==0) {
					theta = incThetaSmr*(2*G4UniformRand()-1);
				}
				
				G4double phi=0;
				if(incPhiSmr>=0    && incDirection.theta() !=0) {
					phi = G4RandGauss::shoot(0,incPhiSmr);
				} else if(incPhiSmr<0 || incDirection.theta()==0) {
					phi = incPhiSmr*(2*G4UniformRand()-1);
				}
	
				ini_Dir.setTheta(ini_Dir.theta()+theta);
				ini_Dir.setPhi(ini_Dir.phi()+phi);

				theta = ini_Dir.theta();
				//bool tmpDD = true;//false;
				//while(tmpDD) {
				//double xx11 = G4UniformRand();// - 1;
				//if(xx11<1.00 && xx11>0.5) {
				//theta = acos(-xx11);
				//break;
				//}
				//}
				//	theta = acos(-1.0);//-0.95);
				//phi = pivalGA*(2*G4UniformRand()-1);
				
				// theta = acos(-1);
				// phi = 0;
				
				ini_Dir.setTheta(theta);
				ini_Dir.setPhi(phi);
				cout<<"theta "<<theta<<" phi  "<<phi<<endl;
				
				if(incEnergySmr>=0) {
					in_Energy += G4RandGauss::shoot(0,incEnergySmr);
					if (in_Energy <1*MeV) in_Energy=1*MeV;
				} else if(incEnergySmr<0) {
					in_Energy += incEnergySmr*(2*G4UniformRand()-1);
					if (in_Energy <1*MeV) in_Energy=1*MeV;
				}
				
				//	in_Energy = 1*GeV;//(3.5 + 3.0*(2*G4UniformRand()-1))*GeV;
				// in_Energy = 0.7*GeV;
				// in_Energy = incEnergy*MeV;
				// in_Energy = enerin[g_nevt/5000]*GeV;
				// cout<<vx<<","<< vy<<","<< vz<<endl;
				if(incVxSmr>=0&&incVySmr>=0&&incVzSmr>=0) {
					vx += G4RandGauss::shoot(0,incVxSmr)*mm;
					vy += G4RandGauss::shoot(0,incVySmr)*mm;
					vz += G4RandGauss::shoot(0,incVzSmr)*mm;
				} else {
					vx += incVxSmr*(2*G4UniformRand()-1.)*mm;
					vy += incVySmr*(2*G4UniformRand()-1.)*mm;
					vz += incVzSmr*(2*G4UniformRand()-1.)*mm;
				}
      }
      
			cout<<"PGA after = "<<vx<<" "<<vy<<" "<<vz<<" "<<ini_Dir<<" "<<in_Energy<<endl;	
			//      cout <<"ini_Dir "<<ini_Dir<<endl;
      particleGun->SetParticleMomentumDirection(ini_Dir);
      particleGun->SetParticleMomentum(in_Energy);
      particleGun->SetParticlePosition(G4ThreeVector(vx, vy, vz));
      particleGun->GeneratePrimaryVertex(anEvent);
      if (ij < (int)pAnalysis->ngenmx) {
				pAnalysis->pidin[ij] = particleGun->GetParticleDefinition()->GetPDGEncoding();
				pAnalysis->posxin[ij]= particleGun->GetParticlePosition().x()/m;
				pAnalysis->posyin[ij]= particleGun->GetParticlePosition().y()/m;
				pAnalysis->poszin[ij]= particleGun->GetParticlePosition().z()/m;
				double momentum = particleGun->GetParticleMomentum()/GeV; 
				if(particle->GetPDGCharge()==0) {
					pAnalysis->momin[ij] = momentum;
				} else {
					pAnalysis->momin[ij] = momentum*(particle->GetPDGCharge());
				}
				pAnalysis->thein[ij] = particleGun->GetParticleMomentumDirection().theta();
				pAnalysis->phiin[ij] = particleGun->GetParticleMomentumDirection().phi();
				
				// cout<<"PGA = "<<ij<<" "<<particleGun->GetParticlePosition()/m<<" "<<pAnalysis->momin[ij]<<" "<<particleGun->GetParticleMomentum()/GeV<<" "<<pAnalysis->thein[ij]*180/pivalGA<<" "<<pAnalysis->phiin[ij]*180/pivalGA<<endl;
      } // if (ij < (int)pAnalysis->ngenmx)
    } // for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++)
  } else if (InputFlag==1) {
 
    // Cosmic Flux from data
    pAnalysis->ievt=g_nevt;
    pAnalysis->ngent = ((unsigned)particleGun->GetNumberOfParticles() <=pAnalysis->ngenmx) ?
      particleGun->GetNumberOfParticles() : pAnalysis->ngenmx;
    for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++) {
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* particle = particleTable->FindParticle(partId);
      particleGun->SetParticleDefinition(particle);
      
      G4ThreeVector ini_Dir(incDirection);
      double rand_theta;
      double theta;
			double theta_gen;
			double phi_gen;
      double phi;
      int brkpt = 1;
      double Point2[3];
      // double energy;
      double vertexX;
      double vertexY;
      double vertexZ;
      double Ini_Theta = 0;
      double Ini_Phi = 0;

      while(brkpt) {
	vx = 2.4*pargas[0]*(G4UniformRand()-0.5);
	vy = 2.4*pargas[1]*(G4UniformRand()-0.5);
	vz = RPCLayerPosZ[toptrgly];
	//   for Standalone Theta generation
	rand_theta = G4UniformRand();
	theta = acos(pow((1-rand_theta*norm1),(1./(PowCosTheta+1.0))))*rad;
	phi = pivalGA*(2*G4UniformRand()-1)*rad;

	double Line1[6];
	double Plane1[6];
	double Point1[3] = {-100000.,-100000.,-100000.};
	Line1[0] = vx;
	Line1[1] = vy;
	Line1[2] = vz;
	Line1[3] = -sin(theta)*cos(phi);
	Line1[4] = -sin(theta)*sin(phi);
	Line1[5] = -cos(theta);
	Plane1[0] =  0;
	Plane1[1] =  0;
	Plane1[2] =  RPCLayerPosZ[bottomtrgly];;
	Plane1[3] = 0;
	Plane1[4] = 0;
	Plane1[5] = 1;
	
	int trgCheck = LinePlaneInt(Line1,Plane1,Point1);
	if(trgCheck == 1) {
	  if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	    double Line2[6];
	    double Plane2[6];
	    for(int xxi=0;xxi<3;xxi++) {Point2[xxi] = -100000000.;}	    
	    for(int lmn=0;lmn<3;lmn++) {Point2[lmn] = Point1[lmn];}
	    Line2[0] = StackPosInWorld[0] + vx;
	    Line2[1] = StackPosInWorld[1] + vy;
	    Line2[2] = StackPosInWorld[2] + vz;
	    Line2[3] = -sin(theta)*cos(phi);
	    Line2[4] = -sin(theta)*sin(phi);
	    Line2[5] = -cos(theta);

	    Plane2[0] =  0;
	    Plane2[1] =  0;
	    Plane2[2] =  WorldZDim - 1*mm;
	    Plane2[3] = 0;
	    Plane2[4] = 0;
	    Plane2[5] = 1;

	    int TopPlane = LinePlaneInt(Line2,Plane2,Point2);
	    if(TopPlane ==1) {
	      if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
		vertexX = Point2[0];
		vertexY = Point2[1];
		vertexZ = Point2[2];
		theta_gen = theta;
		phi_gen = phi;
		Ini_Theta = theta + pivalGA*rad;
		Ini_Phi = phi - pivalGA*rad;
		if(Ini_Phi < -pivalGA) {
		  Ini_Phi = Ini_Phi + 2*pivalGA;
		} else if(Ini_Phi > pivalGA) {
		  Ini_Phi = Ini_Phi - 2*pivalGA;
		}
		brkpt = 0;
	      } // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
	    } // if(TopPlane ==1) {
	  } // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	} // if(trgCheck == 1) {
      } // while(brkpt) {

      double in_Energy = GetCosmicEnergy(ELowLim,EUpLim);
      
      ini_Dir.setTheta(Ini_Theta);
      ini_Dir.setPhi(Ini_Phi);
      particleGun->SetParticleMomentumDirection(ini_Dir);
      particleGun->SetParticleEnergy(in_Energy);
      particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
      particleGun->GeneratePrimaryVertex(anEvent);
      if (ij < (int)pAnalysis->ngenmx) {
	pAnalysis->pidin[ij] = particleGun->GetParticleDefinition()->GetPDGEncoding();
	pAnalysis->posxin[ij]= particleGun->GetParticlePosition().x();
	pAnalysis->posyin[ij]= particleGun->GetParticlePosition().y();
	pAnalysis->poszin[ij]= particleGun->GetParticlePosition().z();
	if(particle->GetPDGCharge()==0) {
	  pAnalysis->momin[ij] = particleGun->GetParticleMomentum()/GeV;
	} else {
	  pAnalysis->momin[ij] = (particleGun->GetParticleMomentum()/GeV)*(particle->GetPDGCharge());
	}
	pAnalysis->thein[ij] = particleGun->GetParticleMomentumDirection().theta();
	pAnalysis->phiin[ij] = particleGun->GetParticleMomentumDirection().phi();
      } // if (ij < (int)pAnalysis->ngenmx) {    
    } // for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++) {
  } else if (InputFlag==2) {
    // Cosmic Flux 3D hist Corsika
    if(initialiseCor==0) {
      FileCORSIKA->cd();
      corsikaFlux = (TH3D*)FileCORSIKA->Get("flux");
      initialiseCor = 1;
    }
    pAnalysis->ievt=g_nevt;
    pAnalysis->ngent = ((unsigned)particleGun->GetNumberOfParticles() <=pAnalysis->ngenmx) ? particleGun->GetNumberOfParticles() : pAnalysis->ngenmx;
    for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++) {
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* particle = particleTable->FindParticle(partId);
      particleGun->SetParticleDefinition(particle);
      G4ThreeVector ini_Dir(incDirection);
      double vertexX;
      double vertexY;
      double vertexZ;
      double Ini_Theta = 0;
      double Ini_Phi = 0;
      double Ini_Enrgy = 0;
      double costheta;
      double phi;
      double logEnrgy;
      double theta;
      double enrgy;
      int brkpt = 1;
      double Point2[3];
      while(brkpt) {
	vx = 2.4*pargas[0]*(G4UniformRand()-0.5);
	vy = 2.4*pargas[1]*(G4UniformRand()-0.5);
	vz = RPCLayerPosZ[toptrgly];
	
	corsikaFlux->GetRandom3(logEnrgy,costheta,phi);
	phi = phi*pivalGA/180;
	phi = phi*rad;
	theta = acos(costheta);
	enrgy = pow(10,logEnrgy);
	double Line1[6];
	double Plane1[6];
	double Point1[3] = {-100000.,-100000.,-100000.};
	Line1[0] = vx;
	Line1[1] = vy;
	Line1[2] = vz;
	Line1[3] = -sin(theta)*cos(phi);
	Line1[4] = -sin(theta)*sin(phi);
	Line1[5] = -cos(theta);
	Plane1[0] =  0;
	Plane1[1] =  0;
	Plane1[2] =  RPCLayerPosZ[bottomtrgly];;
	Plane1[3] = 0;
	Plane1[4] = 0;
	Plane1[5] = 1;
	
	int trgCheck = 1; //LinePlaneInt(Line1,Plane1,Point1);
	if(trgCheck == 1) {
	  // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	  {
	    double Line2[6];
	    double Plane2[6];
	    for(int xxi=0;xxi<3;xxi++) {Point2[xxi] = -100000000.;}	    
	    for(int lmn=0;lmn<3;lmn++) {Point2[lmn] = Point1[lmn];}
	    Line2[0] = StackPosInWorld[0] + vx;
	    Line2[1] = StackPosInWorld[1] + vy;
	    Line2[2] = StackPosInWorld[2] + vz;
	    Line2[3] = -sin(theta)*cos(phi);
	    Line2[4] = -sin(theta)*sin(phi);
	    Line2[5] = -cos(theta);

	    Plane2[0] =  0;
	    Plane2[1] =  0;
	    Plane2[2] =  WorldZDim - 1*mm;
	    Plane2[3] = 0;
	    Plane2[4] = 0;
	    Plane2[5] = 1;
	    int TopPlane = 1; //LinePlaneInt(Line2,Plane2,Point2);
	    if(TopPlane ==1) {
	      // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
	      {
		vertexX = 0;//Point2[0];
		vertexY = 0;//Point2[1];
		vertexZ = 0;//Point2[2];
		Ini_Theta = theta;
		Ini_Phi = phi;
		if(Ini_Phi < -pivalGA) {
		  Ini_Phi = Ini_Phi + 2*pivalGA;
		} else if(Ini_Phi > pivalGA) {
		  Ini_Phi = Ini_Phi - 2*pivalGA;
		}
		Ini_Enrgy = enrgy;
		brkpt = 0;
	      } // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
	    } // if(TopPlane ==1) {
	  } // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	} // // if(trgCheck == 1) {
      } // while(brkpt) {	
      
      ini_Dir.setTheta(Ini_Theta);
      ini_Dir.setPhi(Ini_Phi);
      particleGun->SetParticleMomentumDirection(ini_Dir);
      particleGun->SetParticleEnergy(Ini_Enrgy);
      particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
      particleGun->GeneratePrimaryVertex(anEvent);
      if (ij < (int)pAnalysis->ngenmx) {
	pAnalysis->pidin[ij] = particleGun->GetParticleDefinition()->GetPDGEncoding();
	pAnalysis->posxin[ij]= particleGun->GetParticlePosition().x();
	pAnalysis->posyin[ij]= particleGun->GetParticlePosition().y();
	pAnalysis->poszin[ij]= particleGun->GetParticlePosition().z();
	if(particle->GetPDGCharge()==0) {
	  pAnalysis->momin[ij] = particleGun->GetParticleMomentum()/GeV;
	} else {
	  pAnalysis->momin[ij] = (particleGun->GetParticleMomentum()/GeV)*(particle->GetPDGCharge());
	}
	pAnalysis->thein[ij] = particleGun->GetParticleMomentumDirection().theta();
	pAnalysis->phiin[ij] = particleGun->GetParticleMomentumDirection().phi();
      } // if (ij < (int)pAnalysis->ngenmx) {
    } // for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++) {
  } else if (InputFlag==3) {
    cout<<"(InputFlag==3) "<<endl;
    if (initialiseCor==0) {
      TreeCORSIKA = (TTree*)FileCORSIKA->Get("corsikaTreeAll");
      TreeCORSIKA->SetBranchAddress("iC_nevt",&iC_nevt);
      TreeCORSIKA->SetBranchAddress("iC_eventweight",&iC_eventweight);
      TreeCORSIKA->SetBranchAddress("iC_npart",&iC_npart);
      TreeCORSIKA->SetBranchAddress("iC_cpid",iC_cpid);
      TreeCORSIKA->SetBranchAddress("iC_cvx",iC_cvx);
      TreeCORSIKA->SetBranchAddress("iC_cvy",iC_cvy);
      TreeCORSIKA->SetBranchAddress("iC_cpx",iC_cpx);
      TreeCORSIKA->SetBranchAddress("iC_cpy",iC_cpy);
      TreeCORSIKA->SetBranchAddress("iC_cpz",iC_cpz);
      initialiseCor=1;
    } // if (initialiseCor==0) {
    if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
      TreeCORSIKA->GetEntry(FirstEvt+g_nevt);
      pAnalysis->ievt		= iC_nevt;
      pAnalysis->ievt_wt	= iC_eventweight;
      pAnalysis->ngent = (iC_npart < (int)pAnalysis->ngenmx) ? iC_npart : pAnalysis->ngenmx;
      for (unsigned int count=0; count<pAnalysis->ngent; count++) {
	cout<<"count "<<count<<endl;
	if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
	  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	  G4ParticleDefinition* particle = particleTable->FindParticle(iC_cpid[count]);
	  particleGun->SetParticleDefinition(particle);
	  particleGun->SetParticlePosition(G4ThreeVector(0,0,0));
	  G4ThreeVector tmp3v(iC_cpx[count]*GeV,iC_cpy[count]*GeV,-iC_cpz[count]*GeV);
	  cout<<tmp3v<<endl;
	  
	  double Point2[3];
	  // double energy;
	  double vertexX;
	  double vertexY;
	  double vertexZ;
	  double Ini_Theta = 0;
	  double Ini_Phi = 0;
	  G4ThreeVector ini_Dir(tmp3v.unit());
	  cout<<"Theta::"<<ini_Dir.theta()<<"  "<<ini_Dir<<endl;
	  
	
	  double theta, phi;
	  int brkpt = 1;
	  int brkcnt = 0;
	  // if (ini_Dir.phi()<0) continue;

	  // if( ini_Dir.theta()>2.3 || ini_Dir.theta()< 1.9) {
	  //   //  goto label;
	  //   cout<<"Theta outside range"<<endl;
	  // continue;
	  // }

	    
	  


  cout<<"Theta:"<<ini_Dir.theta()<<"  "<<ini_Dir<<endl;


	  
   cout<<"toptrgly "<<toptrgly<<endl;
	  pAnalysis->ngenerated = 0;
	  pAnalysis->naperture = 0;
	  //	  double stackposy = parchm[1];
	  // if (G4UniformRand()>0.5) stackposy*= -1;
	  // cout<<"stackposy: "<<stackposy<<endl;
	  while(brkpt) {
	    
	    //	    vx = 1.2*pargas[0]*(2*G4UniformRand()-1.0);
	    //	    vy = 1.2*pargas[1]*(2*G4UniformRand()-1.0);
	    //	    vz = RPCLayerPosZ[toptrgly];
	    vx =StackPosInWorld[0]+ 6000*(2*G4UniformRand()-1.0);
	      vy =StackPosInWorld[1]+ 6000*(2*G4UniformRand()-1.0);
	      //vy = -7500+1000*(2*G4UniformRand()-1.0); // from back side
	    //   	     vy =-12000+ 3000*(2*G4UniformRand()-1.0);
	      //  vz =   ScintlayGlPos[0][3][2] + 831 + 50*(2*G4UniformRand()-1.0)  ;
	      vz=paradef->GetINOroomPos(2)+1500;
	      
	    //	    vx = 0;
	    //	    vy = StackPosInWorld[1];
	    //	    vz = 0;
	    // vx =(-3800 + 100*(2*G4UniformRand()-1.0))*mm;
	    //  vy = (-11730 + 1000*(2*G4UniformRand()-1.0))*mm;
	    //  vz = (-250 +   100*(2*G4UniformRand()-1.0))*mm;
	      
	    cout<<vx<<" "<<vy<<" "<<vz<<endl;
	    // if(abs(vx)<pargas[0] && abs(vy)<pargas[1]) {
	    // cout<<"vertex inside"<<endl;
	    // pAnalysis->ngenerated++;
	      // }
	    phi = ini_Dir.phi();
	    theta = ini_Dir.theta();
	  
	    double Line1[6];
	    double Plane1[6];
	    double Point1[3] = {-100000.,-100000.,-100000.};
	    Line1[0] = vx;
	    Line1[1] = vy;
	    Line1[2] = vz;
	    Line1[3] = -sin(theta)*cos(phi);
	    Line1[4] = -sin(theta)*sin(phi);
	    Line1[5] = -cos(theta);
	    Plane1[0] = StackPosInWorld[0];//=  INOroomPos[0]+StackPosInRoom[0];
	    Plane1[1] = StackPosInWorld[1];//=  INOroomPos[1]+StackPosInRoom[1];
	    Plane1[2] = StackPosInWorld[2] + ShiftInZ +  RPCLayerPosZ[0];// INOroomPos[2]+StackPosInRoom[2]+ShiftInZ+  RPCLayerPosZ[bottomtrgly];
	    Plane1[3] = 0;
	    Plane1[4] = 0;
	    Plane1[5] = 1;
	    int trgCheck = LinePlaneInt(Line1,Plane1,Point1);
	    
	    cout<<"trgcheck"<<trgCheck<<endl;
	    cout<<"Point1:"<<Point1[0]<<" "<<Point1[1]<<" "<<Point1[2]<<endl;
	    if(trgCheck == 1) {
	      if( (Point1[0]<StackPosInWorld[0]+ShiftInX+pargas[0] &&
		   Point1[0]>StackPosInWorld[0]+ShiftInX-pargas[0] &&
		   Point1[1]<StackPosInWorld[1]+parchm[1]+ShiftInY+pargas[1] &&
		   Point1[1]>StackPosInWorld[1]+parchm[1]+ShiftInY-pargas[1]) ||
		  (Point1[0]<StackPosInWorld[0]+ShiftInX+pargas[0] &&
		   Point1[0]>StackPosInWorld[0]+ShiftInX-pargas[0] &&
		   Point1[1]<StackPosInWorld[1]-parchm[1]+ShiftInY+pargas[1] &&
		   Point1[1]>StackPosInWorld[1]-parchm[1]+ShiftInY-pargas[1] ))

		 
		{
		  cout<<" point1 inside"<<endl;
	        	pAnalysis->naperture++;
		double Line2[6];
		double Plane2[6];
		//		for(int lmn=0;lmn<3;lmn++) {Point2[lmn] = Point1[lmn];}
		//			Line2[0] = StackPosInWorld[0] + vx;
		//		Line2[1] = StackPosInWorld[1] + vy + stackposy;
		//		Line2[2] = StackPosInWorld[2] + vz;
		Line2[0] = Point1[0];
		Line2[1] = Point1[1];
		Line2[2] = Point1[2];
		
		Line2[3] = -sin(theta)*cos(phi);
		Line2[4] = -sin(theta)*sin(phi);
		Line2[5] = -cos(theta);
		
		//	Plane2[0] =  ScintlayGlPos[0][3][0] ;
		//		Plane2[1] =  ScintlayGlPos[0][3][1] ;
		//		Plane2[2] =   ScintlayGlPos[0][3][2] + 500*mm;// WorldZDim - 1*mm;
		Plane2[0] =  StackPosInWorld[0]; 
		Plane2[1] =  StackPosInWorld[1] ;
		Plane2[2] = StackPosInWorld[2]+ShiftInZ + RPCLayerPosZ[3];// WorldZDim - 1*mm;
		
		
		cout<<Plane2[0]<<" "<<Plane2[1]<<" "<<Plane2[2]<<endl;
		Plane2[3] = 0;
		Plane2[4] = 0;
		Plane2[5] = 1;
		int TopPlane = LinePlaneInt(Line2,Plane2,Point2);
		cout<<"topplane"<<TopPlane<<endl;
		if(TopPlane ==1) {
		  cout<<Point2[0]<<" "<<Point2[1]<<" "<<Point2[2]<<endl;
		  // if(Point2[0]<(ScintlayGlPos[0][3][0]+layhalflength) && Point2[0]>(ScintlayGlPos[0][3][0]-layhalflength)  /*WorldXDim*/ && Point2[1]<(ScintlayGlPos[0][3][1]+2300) && Point2[1]>(ScintlayGlPos[0][3][1]-2300) /*WorldYDim*/) {
		  if( (Point2[0]<StackPosInWorld[0]+ShiftInX+pargas[0] &&
		       Point2[0]>StackPosInWorld[0]+ShiftInX-pargas[0] &&
		       Point2[1]<StackPosInWorld[1]+parchm[1]+ShiftInY+pargas[1] &&
		       Point2[1]>StackPosInWorld[1]+parchm[1]+ShiftInY-pargas[1])
		      ||
		      (Point2[0]<StackPosInWorld[0]+ShiftInX+pargas[0] &&
		       Point2[0]>StackPosInWorld[0]+ShiftInX-pargas[0] &&
		       Point2[1]<StackPosInWorld[1]-parchm[1]+ShiftInY+pargas[1] &&
		       Point2[1]>StackPosInWorld[1]-parchm[1]+ShiftInY-pargas[1] ))

		 
		{
		     pAnalysis->ngenerated++;
		    cout<<"inside scintillator"<<endl;
		    //vertexX = Point2[0];
		    //vertexY = Point2[1];
		    //vertexZ = Point2[2];

		   vertexX = vx;
		    vertexY = vy;
		    vertexZ = vz;
		    Ini_Theta = theta;
		    Ini_Phi = phi;
		    if(Ini_Phi < -pivalGA) {
		      Ini_Phi = Ini_Phi + 2*pivalGA;
		    } else if(Ini_Phi > pivalGA) {
		      Ini_Phi = Ini_Phi - 2*pivalGA;
		    }
		    brkpt = 0;

		} // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
		} // if(TopPlane ==1) {
	      } // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	    }// if(trgCheck == 1) {
	    brkcnt++;
	    cout<<endl;
	    if(brkcnt>1000) {
	      cout<<"brkcnt "<<brkcnt<<endl;
	      vertexX = 0;
	      vertexY = 0; //Point2[1];
	      vertexZ = WorldZDim - 1*mm; //Point2[2];
	      Ini_Theta = theta;
	      Ini_Phi = phi;
	      if(Ini_Phi < -pivalGA) {
		Ini_Phi = Ini_Phi + 2*pivalGA;
	      } else if(Ini_Phi > pivalGA) {
		Ini_Phi = Ini_Phi - 2*pivalGA;
	      }
	      brkpt = 0;
	      cout<<"brkt "<<brkcnt<<" "<<Ini_Theta*180/pivalGA <<" "<<Ini_Phi*180/pivalGA<<endl;
	    }
	  } // while(brkpt) {	
	  
	  // Ini_Enrgy = enrgy;
	  cout<<"PGA after = "<<vertexX<<" "<<vertexY<<" "<<vertexZ<<" "<<ini_Dir<<" "<<" "<<ini_Dir.theta()<<" "<<ini_Dir.phi()<<" "<<tmp3v.mag()<<endl;
	  ini_Dir.setTheta(Ini_Theta);
	  ini_Dir.setPhi(Ini_Phi);
	  particleGun->SetParticleMomentumDirection(ini_Dir);
	  particleGun->SetParticleMomentum(tmp3v.mag());
	  particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
	  particleGun->GeneratePrimaryVertex(anEvent);
	  if (count < (int)pAnalysis->ngenmx) {
	    pAnalysis->pidin [count]= particleGun->GetParticleDefinition()->GetPDGEncoding();//_cpid[count];
	    pAnalysis->posxin[count]= particleGun->GetParticlePosition().x();
	    pAnalysis->posyin[count]= particleGun->GetParticlePosition().y();
	    pAnalysis->poszin[count]= particleGun->GetParticlePosition().z();
	    if(particle->GetPDGCharge()==0){
	      pAnalysis->momin[count] = particleGun->GetParticleMomentum()/GeV;
	    } else {
	      pAnalysis->momin[count] = (particleGun->GetParticleMomentum())*(particle->GetPDGCharge())/GeV;
	    }
	    pAnalysis->thein[count] = particleGun->GetParticleMomentumDirection().theta();
	    pAnalysis->phiin[count] = particleGun->GetParticleMomentumDirection().phi();
	  } // if (count < (int)pAnalysis->ngenmx) {
	} // if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {

	//  label:
	//	cout<<" "<<endl;
      } // for (int count=0; count<iC_npart; count++) {
    } // if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
    // Corsika Event By Event
  } else if (InputFlag==4) { //Modified for flux_Pethu_FLUKA_SIBYLL.root file.

   cout<<"toptrgly "<<toptrgly<<" bottomtrgly " <<bottomtrgly<<endl;
   // cout<<" InputFlag==4 "<<endl;
    // system("free");

    if (initialiseCor==0) {
      FileFLUX->cd();
      muFlux = (TH3F*)FileFLUX->Get("muFlux");
      mupFlux = (TH3F*)FileFLUX->Get("mupFlux");
      munFlux = (TH3F*)FileFLUX->Get("munFlux");
      initialiseCor=1;
    } // if (initialiseCor==0) {
    // cout<<" InputFlag==4 1 "<<endl;
    // system("free");

    //    if(1) {
      pAnalysis->ievt		= g_nevt;
      pAnalysis->ievt_wt	= 1.;
      pAnalysis->ngent = particleGun->GetNumberOfParticles();	// number of particles
      // cout << " npart " <<  pAnalysis->ngent << endl;
       if(pAnalysis->InputOutput==0 || pAnalysis->InputOutput==1 || pAnalysis->InputOutput==2) {
	 // if(1) {
	for (unsigned int count=0; count<pAnalysis->ngent; count++) {
	  // cout<<"seed "<<gRandom->GetSeed()<<endl;
	  
	  double Pxx,Pyy,Pzz;
	  muFlux->GetRandom3(Pxx,Pyy,Pzz);
	  int gBin = muFlux->FindBin(Pxx,Pyy,Pzz);
	  cout<<"gBin "<<gBin<<endl;
	  double mupCnt = mupFlux->GetBinContent(gBin);
	  double munCnt = munFlux->GetBinContent(gBin);
	  cout<<mupCnt<<" "<<munCnt<<endl;
	  int partID;
	  if(G4UniformRand()*(mupCnt+munCnt)>munCnt) {
	    partID = -13;
	  } else {
	    partID = 13;
	  }
	  cout << " partID " << partID << endl;
	  
	  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	  G4ParticleDefinition* particle = particleTable->FindParticle(partID);
	  particleGun->SetParticleDefinition(particle);
	  particleGun->SetParticlePosition(G4ThreeVector(0,0,0));

	  G4ThreeVector tmp3v(Pxx*GeV,Pyy*GeV,Pzz*GeV);
	  
	  double Point2[3];
	  // double energy;
	  double vertexX;
	  double vertexY;
	  double vertexZ;
	  double Ini_Theta = 0;
	  double Ini_Phi = 0;
	  G4ThreeVector ini_Dir(tmp3v.unit());
	  double theta, phi;
	  int brkpt = 1;
	  int brkcnt = 0;

	  pAnalysis->ngenerated = 0;
	  pAnalysis->naperture = 0;
	  
	  while(brkpt) {

	    vx = pargas[0]*(2*G4UniformRand()-1.0);
            vy =(2*(pargas[1]+parchm[1])*G4UniformRand())+(-parchm[1]-paradef->GetParFrpBox(1)+paradef->GetShiftInY() );
	    vz = RPCLayerPosZ[toptrgly];

	    cout<<"vy "<<vy<<endl;	    
	    if(vy<(parchm[1]+paradef->GetShiftInY()-paradef->GetParFrpBox(1)) && vy > (-parchm[1]+paradef->GetParFrpBox(1)+paradef->GetShiftInY()+2*pargas[1] ) ) continue;//Dead space between 2-stacks

	    cout<<"..inside top RPC.."<<endl;
	    if(fabs(vx)<pargas[0]) {
	      pAnalysis->ngenerated++;
	    }
	    phi = ini_Dir.phi();
	    theta = ini_Dir.theta();
	    double Line1[6];
	    double Plane1[6];
	    double Point1[3] = {-100000.,-100000.,-100000.};
	    Line1[0] = vx;
	    Line1[1] = vy;
	    Line1[2] = vz;
	    Line1[3] = -sin(theta)*cos(phi);
	    Line1[4] = -sin(theta)*sin(phi);
	    Line1[5] = -cos(theta);
	    Plane1[0] =  0;
	    Plane1[1] =  0;
	    Plane1[2] =  RPCLayerPosZ[bottomtrgly];;
	    Plane1[3] = 0;
	    Plane1[4] = 0;
	    Plane1[5] = 1;
	    int trgCheck = LinePlaneInt(Line1,Plane1,Point1);
	    if(trgCheck == 1) {
	   
	      if(fabs(Point1[0])<pargas[0]){

		if( (Point1[1]>(parchm[1]+paradef->GetShiftInY()-paradef->GetParFrpBox(1) ) && Point1[1] < (parchm[1]+paradef->GetShiftInY()+paradef->GetParFrpBox(1)+2*pargas[1] ) )    ||  (Point1[1] < (-parchm[1]+paradef->GetParFrpBox(1)+paradef->GetShiftInY()+2*pargas[1]) && Point1[1]>  (-parchm[1]-paradef->GetParFrpBox(1)+paradef->GetShiftInY())          ) ) { //Last condition for deadspace between 2 stacks
		 cout<<"..inside bottom RPC.."<<endl;
		pAnalysis->naperture++;
		double Line2[6];
		double Plane2[6];
		for(int lmn=0;lmn<3;lmn++) {Point2[lmn] = Point1[lmn];}
		Line2[0] = StackPosInWorld[0] + paradef->GetShiftInX() + vx;
		Line2[1] = StackPosInWorld[1] + vy;
		Line2[2] = StackPosInWorld[2]+ paradef->GetShiftInZ(0) + vz;
		Line2[3] = -sin(theta)*cos(phi);
		Line2[4] = -sin(theta)*sin(phi);
		Line2[5] = -cos(theta);
		
		Plane2[0] =  0;
		Plane2[1] =  0;
		Plane2[2] =paradef->GetINOroomPos(2)+paradef->GetParairroom(2)-1*mm;//room top // WorldZDim - 1*mm;
		Plane2[3] = 0;
		Plane2[4] = 0;
		Plane2[5] = 1;
		int TopPlane = LinePlaneInt(Line2,Plane2,Point2);
		if(TopPlane ==1) {
		  if(fabs(Point2[0])<WorldXDim && fabs(Point2[1])<WorldYDim) {
		    vertexX = Point2[0];
		    vertexY = Point2[1];
		    vertexZ = Point2[2];
		    Ini_Theta = theta;
		    Ini_Phi = phi;
		    if(Ini_Phi < -pivalGA) {
		      Ini_Phi = Ini_Phi + 2*pivalGA;
		    } else if(Ini_Phi > pivalGA) {
		      Ini_Phi = Ini_Phi - 2*pivalGA;
		    }
		    brkpt = 0;

		  } // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
		} // if(TopPlane ==1) {
	      }
	      } // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	    }// if(trgCheck == 1) {
	    brkcnt++;
	    if(brkcnt>1000) {
	      cout<<"brkcnt "<<brkcnt<<endl;
	      vertexX = 0;
	      vertexY = 0; //Point2[1];
	      vertexZ = WorldZDim - 1*mm; //Point2[2];
	      Ini_Theta = TMath::Pi()*0.5;
	      // Ini_Theta = theta;
	      Ini_Phi = phi;
	      if(Ini_Phi < -pivalGA) {
		Ini_Phi = Ini_Phi + 2*pivalGA;
	      } else if(Ini_Phi > pivalGA) {
		Ini_Phi = Ini_Phi - 2*pivalGA;
	      }
	      brkpt = 0;
	      cout<<"brkt "<<brkcnt<<" "<<Ini_Theta*180/pivalGA <<" "<<Ini_Phi*180/pivalGA<<endl;
	    }
	  } // while(brkpt) {	
	   cout << " count " << count << " vertexX " << vertexX << " vertexY " << vertexY << " vertexZ " << vertexZ << endl;
	   cout << " mom theta phi " << tmp3v.mag() << " " << Ini_Theta << " " << Ini_Phi << endl;
	  
	  // Ini_Enrgy = enrgy;
	  // double Ini_Enrgy = (G4UniformRand()*EUpLim+ELowLim)*MeV;
	  // cout << " Ini_Enrgy " << Ini_Enrgy << endl;
	  // cout << " ELowLim " << ELowLim << " EUpLim " << EUpLim << endl;	  
	  
	  ini_Dir.setTheta(Ini_Theta);
	  ini_Dir.setPhi(Ini_Phi);
	  particleGun->SetParticleMomentumDirection(ini_Dir);
	  particleGun->SetParticleMomentum(tmp3v.mag());
	  particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
	  particleGun->GeneratePrimaryVertex(anEvent);
	  if (count < (int)pAnalysis->ngenmx) {
	    pAnalysis->pidin [count]= particleGun->GetParticleDefinition()->GetPDGEncoding();
	    pAnalysis->posxin[count]= particleGun->GetParticlePosition().x();
	    pAnalysis->posyin[count]= particleGun->GetParticlePosition().y();
	    pAnalysis->poszin[count]= particleGun->GetParticlePosition().z();
	    if(particle->GetPDGCharge()==0){
	      pAnalysis->momin[count] = particleGun->GetParticleMomentum()/GeV;
	    } else {
	      pAnalysis->momin[count] = (particleGun->GetParticleMomentum())*(particle->GetPDGCharge())/GeV;
	    }
	    pAnalysis->thein[count] = particleGun->GetParticleMomentumDirection().theta();
	    pAnalysis->phiin[count] = particleGun->GetParticleMomentumDirection().phi();
	  } // if (count < (int)pAnalysis->ngenmx) {
	} // if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
       // } else {
       // 	pAnalysis->ngent = 0;}
       } // if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
    // cout<<" InputFlag==4 end "<<endl;
    // system("free");
    

  
  }else if (InputFlag==5) {
    if (initialiseCor==0) {
      TreeCORSIKA = (TTree*)FileCORSIKA->Get("corsikaTreeAll");
      TreeCORSIKA->SetBranchAddress("iC_nevt",&iC_nevt);
      TreeCORSIKA->SetBranchAddress("iC_eventweight",&iC_eventweight);
      TreeCORSIKA->SetBranchAddress("iC_npart",&iC_npart);
      TreeCORSIKA->SetBranchAddress("iC_cpid",iC_cpid);
      TreeCORSIKA->SetBranchAddress("iC_cvx",iC_cvx);
      TreeCORSIKA->SetBranchAddress("iC_cvy",iC_cvy);
      TreeCORSIKA->SetBranchAddress("iC_cpx",iC_cpx);
      TreeCORSIKA->SetBranchAddress("iC_cpy",iC_cpy);
      TreeCORSIKA->SetBranchAddress("iC_cpz",iC_cpz);
      initialiseCor=1;
    } // if (initialiseCor==0) {
    if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
      TreeCORSIKA->GetEntry(FirstEvt+g_nevt);
      pAnalysis->ievt		= iC_nevt;
      pAnalysis->ievt_wt	= iC_eventweight;
      pAnalysis->ngent = (iC_npart < (int)pAnalysis->ngenmx) ? iC_npart : pAnalysis->ngenmx;
      for (unsigned int count=0; count<pAnalysis->ngent; count++) {
	if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
	  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	  G4ParticleDefinition* particle = particleTable->FindParticle(iC_cpid[count]);
	  particleGun->SetParticleDefinition(particle);
	  particleGun->SetParticlePosition(G4ThreeVector(0,0,0));
	  G4ThreeVector tmp3v(iC_cpx[count]*GeV,iC_cpy[count]*GeV,iC_cpz[count]*GeV);
	  // cout<<"Momin "<<tmp3v.mag()<<endl;
	  double Point2[3];
	  // double energy;
	  double vertexX;
	  double vertexY;
	  double vertexZ;
	  double Ini_Theta = 0;
	  double Ini_Phi = 0;
	  G4ThreeVector ini_Dir(tmp3v.unit());
	  double theta, phi;
	  int brkpt = 1;
	  int brkcnt = 0;

	  pAnalysis->ngenerated = 0;
	  pAnalysis->naperture = 0;

	  vx = 0;//500*(2*G4UniformRand()-1.0);
	  vy = 0;//500*(2*G4UniformRand()-1.0);
	  vz = RPCLayerPosZ[toptrgly] + 101*mm;
	  
	  vertexX = vx + StackPosInWorld[0];
	  vertexY = vy + StackPosInWorld[1];
	  vertexZ = vz + StackPosInWorld[2];
	  
	  phi = ini_Dir.phi();
	  theta = ini_Dir.theta();

	  Ini_Theta = theta;
	  Ini_Phi = phi;
	  if(Ini_Phi < -pivalGA) {
	    Ini_Phi = Ini_Phi + 2*pivalGA;
	  } else if(Ini_Phi > pivalGA) {
	    Ini_Phi = Ini_Phi - 2*pivalGA;
	  }
	  	  
	  // Ini_Enrgy = enrgy;
	    
	  ini_Dir.setTheta(Ini_Theta);
	  ini_Dir.setPhi(Ini_Phi);
	  particleGun->SetParticleMomentumDirection(ini_Dir);
	  particleGun->SetParticleMomentum(tmp3v.mag());
	  particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
	  particleGun->GeneratePrimaryVertex(anEvent);
	  // cout<<"particleGun->PartPos "<<particleGun->GetParticlePosition()<<endl;
	  // cout<<"particleGun->Mom "<<particleGun->GetParticleMomentum()/GeV<<endl;
	  if (count < (int)pAnalysis->ngenmx) {
	    pAnalysis->pidin [count]= particleGun->GetParticleDefinition()->GetPDGEncoding();//_cpid[count];
	    pAnalysis->posxin[count]= particleGun->GetParticlePosition().x();
	    pAnalysis->posyin[count]= particleGun->GetParticlePosition().y();
	    pAnalysis->poszin[count]= particleGun->GetParticlePosition().z();
	    if(particle->GetPDGCharge()==0){
	      pAnalysis->momin[count] = particleGun->GetParticleMomentum()/GeV;
	    } else {
	      pAnalysis->momin[count] = (particleGun->GetParticleMomentum())*(particle->GetPDGCharge())/GeV;
	    }
	    pAnalysis->thein[count] = particleGun->GetParticleMomentumDirection().theta();
	    pAnalysis->phiin[count] = particleGun->GetParticleMomentumDirection().phi();
	  } // if (count < (int)pAnalysis->ngenmx) {
	} // if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
      } // for (int count=0; count<iC_npart; count++) {
    } // if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
    // Corsika Event By Event
  }
  else
    	// cmvd case
	if(InputFlag==6){
	  
 if (initialiseCor==0) {
      TreeCORSIKA = (TTree*)FileCORSIKA->Get("corsikaTreeAll");
      TreeCORSIKA->SetBranchAddress("iC_nevt",&iC_nevt);
      TreeCORSIKA->SetBranchAddress("iC_eventweight",&iC_eventweight);
      TreeCORSIKA->SetBranchAddress("iC_npart",&iC_npart);
      TreeCORSIKA->SetBranchAddress("iC_cpid",iC_cpid);
      TreeCORSIKA->SetBranchAddress("iC_cvx",iC_cvx);
      TreeCORSIKA->SetBranchAddress("iC_cvy",iC_cvy);
      TreeCORSIKA->SetBranchAddress("iC_cpx",iC_cpx);
      TreeCORSIKA->SetBranchAddress("iC_cpy",iC_cpy);
      TreeCORSIKA->SetBranchAddress("iC_cpz",iC_cpz);
      initialiseCor=1;
    } // if (initialiseCor==0) {
    if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
      TreeCORSIKA->GetEntry(FirstEvt+g_nevt);
      pAnalysis->ievt		= iC_nevt;
      pAnalysis->ievt_wt	= iC_eventweight;
      pAnalysis->ngent = (iC_npart < (int)pAnalysis->ngenmx) ? iC_npart : pAnalysis->ngenmx;
      for (unsigned int count=0; count<pAnalysis->ngent; count++) {
	if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {

	  
	  particleGun->SetParticlePosition(G4ThreeVector(0,0,0));
	  G4ThreeVector tmp3v(iC_cpx[count]*GeV,iC_cpy[count]*GeV,iC_cpz[count]*GeV);
	  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    if (G4UniformRand()>0.5) partId *= -1;
      cout<<"partId "<<partId<<endl;
      G4ParticleDefinition* particle = particleTable->FindParticle(partId);
      particleGun->SetParticleDefinition(particle);
        double Pie = acos(-1.0);
	 vx =  2000*(2*G4UniformRand()-1)*mm;
	 vy = 2300*(2*G4UniformRand()-1)*mm;
	vy +=-10730*mm;
	 vz = -500*mm;


	 //	G4ThreeVector ini_Dir();

        // cout<<"incDirection "<<ini_Dir<<endl;
     
	// 		vx = incPosition.x()*mm;
	// 		vy = incPosition.y()*mm;
	// 		vz = incPosition.z()*mm;
	        
	 //	cout<<"theta "<<ini_Dir.theta()<<" phi  "<<ini_Dir.phi()<<endl;
				
	 /*				
      if (rndmFlag=="on") {
			

	G4double theta=0;
				if(incThetaSmr>=0    && incDirection.theta() !=0) {
					theta = G4RandGauss::shoot(0,incThetaSmr);
				} else if(incThetaSmr<0 || incDirection.theta()==0) {
					theta = incThetaSmr*(2*G4UniformRand()-1);
				}
				
				G4double phi=0;
				if(incPhiSmr>=0    && incDirection.theta() !=0) {
					phi = G4RandGauss::shoot(0,incPhiSmr);
				} else if(incPhiSmr<0 || incDirection.theta()==0) {
					phi = incPhiSmr*(2*G4UniformRand()-1);
				}
	
				ini_Dir.setTheta(ini_Dir.theta()+theta);
				ini_Dir.setPhi(ini_Dir.phi()+phi);

		       
				cout<<"theta "<<theta<<" phi  "<<phi<<endl;



	
				if(incVxSmr>=0&&incVySmr>=0&&incVzSmr>=0) {
					vx += G4RandGauss::shoot(0,incVxSmr)*mm;
					vy += G4RandGauss::shoot(0,incVySmr)*mm;
					vz += G4RandGauss::shoot(0,incVzSmr)*mm;
				} else {
					vx += incVxSmr*(2*G4UniformRand()-1.)*mm;
					vy += incVySmr*(2*G4UniformRand()-1.)*mm;
					vz += incVzSmr*(2*G4UniformRand()-1.)*mm;
				}
      }
	 */
//       double dirx = abs(ini_Dir.x());
//       double dirz = abs(ini_Dir.z());
//       double diry = abs(ini_Dir.y());
//       cout<<"PGA after = "<<vx<<" "<<vy<<" "<<vz<<" "<<ini_Dir<<" "<<tmp3v.mag()<<endl;	
//   if(vx<0){
//     ini_Dir.setX(dirx)  ;
//     ini_Dir.setZ(-dirz);
//     // ini_Dir.setY(-diry);
//   }
//   else{
//     ini_Dir.setX(-dirx) ;
// ini_Dir.setZ(-dirz);
// // ini_Dir.setY(-diry);
//   }
      cout<<"PGA after = "<<vx<<" "<<vy<<" "<<vz<<" "<<tmp3v<<" "<<tmp3v.mag()<<endl;	
			//      cout <<"ini_Dir "<<ini_Dir<<endl;	


        particleGun->SetParticleMomentumDirection(tmp3v);
	particleGun->SetParticleMomentum( tmp3v.mag());
	particleGun->SetParticlePosition(G4ThreeVector(vx,vy,vz));

      // particleGun->SetParticleMomentumDirection(G4ThreeVector(-0.0804781,0.00115967,-1.04092));
      // 	  particleGun->SetParticleMomentum( 10662.8);
      // 	  particleGun->SetParticlePosition(G4ThreeVector(328.103,-9798.48,0.0));
	  particleGun->GeneratePrimaryVertex(anEvent);
	  if (count < (int)pAnalysis->ngenmx) {
	    pAnalysis->pidin [count]= particleGun->GetParticleDefinition()->GetPDGEncoding();//_cpid[count];
	    pAnalysis->posxin[count]= particleGun->GetParticlePosition().x();
	    pAnalysis->posyin[count]= particleGun->GetParticlePosition().y();
	    pAnalysis->poszin[count]= particleGun->GetParticlePosition().z();
	    if(particle->GetPDGCharge()==0){
	      pAnalysis->momin[count] = particleGun->GetParticleMomentum()/GeV;
	    } else {
	      pAnalysis->momin[count] = (particleGun->GetParticleMomentum())*(particle->GetPDGCharge())/GeV;
	    }
	    pAnalysis->thein[count] = particleGun->GetParticleMomentumDirection().theta();
	    pAnalysis->phiin[count] = particleGun->GetParticleMomentumDirection().phi();
	  } // if (count < (int)pAnalysis->ngenmx) {
	} // if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
      } // for (int count=0; count<iC_npart; count++) {
    } // if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
    // Corsika Event By Event


	  
	}//if(InputFlag==6)


 else if (InputFlag==7) {
    if (initialiseCor==0) {
      TreeCORSIKA = (TTree*)FileCORSIKA->Get("corsikaTreeAll");
      TreeCORSIKA->SetBranchAddress("iC_nevt",&iC_nevt);
      TreeCORSIKA->SetBranchAddress("iC_eventweight",&iC_eventweight);
      TreeCORSIKA->SetBranchAddress("iC_npart",&iC_npart);
      TreeCORSIKA->SetBranchAddress("iC_cpid",iC_cpid);
      TreeCORSIKA->SetBranchAddress("iC_cvx",iC_cvx);
      TreeCORSIKA->SetBranchAddress("iC_cvy",iC_cvy);
      TreeCORSIKA->SetBranchAddress("iC_cpx",iC_cpx);
      TreeCORSIKA->SetBranchAddress("iC_cpy",iC_cpy);
      TreeCORSIKA->SetBranchAddress("iC_cpz",iC_cpz);
      initialiseCor=1;
    } // if (initialiseCor==0) {
    if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
      TreeCORSIKA->GetEntry(FirstEvt+g_nevt);
      pAnalysis->ievt		= iC_nevt;
      pAnalysis->ievt_wt	= iC_eventweight;
      pAnalysis->ngent = (iC_npart < (int)pAnalysis->ngenmx) ? iC_npart : pAnalysis->ngenmx;
      for (unsigned int count=0; count<pAnalysis->ngent; count++) {
	if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {

	  cout<<"Particle Id: "<<iC_cpid[count]<<endl;
	  particleGun->SetParticlePosition(G4ThreeVector(0,0,0));
	  G4ThreeVector tmp3v(iC_cpx[count]*GeV,iC_cpy[count]*GeV,iC_cpz[count]*GeV);
	  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // if (G4UniformRand()>0.5) partId *= -1;
  // cout<<"partId "<<partId<<endl;
      G4ParticleDefinition* particle = particleTable->FindParticle(partId);
      particleGun->SetParticleDefinition(particle);
        double Pie = acos(-1.0);
	  G4ThreeVector ini_Dir(incDirection);
	  // G4ThreeVector ini_Dir(0,0,1);
	  //  ini_Dir.setTheta((Pie-2.2)*G4UniformRand()+2.2);
	  //  ini_Dir.setPhi(G4UniformRand()*2*Pie);
        cout<<"incDirection "<<ini_Dir<<endl;
     
			vx = incPosition.x()*mm;
			vy = incPosition.y()*mm;
			vz = incPosition.z()*mm;
	        
						cout<<"theta "<<ini_Dir.theta()<<" phi  "<<ini_Dir.phi()<<endl;
				
					
      if (rndmFlag=="on") {
	cout<<"Smear in angle: "<<incThetaSmr<<" "<<incPhiSmr<<endl;

	G4double theta=0;
				if(incThetaSmr>=0    && incDirection.theta() !=0) {
					theta = G4RandGauss::shoot(0,incThetaSmr);
				} else if(incThetaSmr<0 || incDirection.theta()==0) {
					theta = incThetaSmr*(2*G4UniformRand()-1);
					
				}
				cout<<"smeared theta: "<<theta<<endl;
				G4double phi=0;
				if(incPhiSmr>=0    && incDirection.theta() !=0) {
					phi = G4RandGauss::shoot(0,incPhiSmr);
				} else if(incPhiSmr<0 || incDirection.theta()==0) {
					phi = incPhiSmr*(2*G4UniformRand()-1);
				}

				cout<<"smeared phi: "<<phi<<endl;
				ini_Dir.setTheta(ini_Dir.theta()+theta);
				ini_Dir.setPhi(ini_Dir.phi()+phi);

		       
				cout<<"theta "<<theta<<" phi  "<<phi<<endl;



	
				if(incVxSmr>=0&&incVySmr>=0&&incVzSmr>=0) {
					vx += G4RandGauss::shoot(0,incVxSmr)*mm;
					vy += G4RandGauss::shoot(0,incVySmr)*mm;
					vz += G4RandGauss::shoot(0,incVzSmr)*mm;
				} else {
					vx += incVxSmr*(2*G4UniformRand()-1.)*mm;
					vy += incVySmr*(2*G4UniformRand()-1.)*mm;
					vz += incVzSmr*(2*G4UniformRand()-1.)*mm;
				}
      }

      double dirx = abs(ini_Dir.x());
      double dirz = abs(ini_Dir.z());
      double diry = abs(ini_Dir.y());
      // cout<<"PGA after = "<<vx<<" "<<vy<<" "<<vz<<" "<<ini_Dir<<" "<<tmp3v.mag()<<endl;	
 //  if(vx<0){
//     ini_Dir.setX(dirx)  ;
//     ini_Dir.setZ(-dirz);
//     // ini_Dir.setY(-diry);
//   }
//   else{
//     ini_Dir.setX(-dirx) ;
// ini_Dir.setZ(-dirz);
// // ini_Dir.setY(-diry);
//   }
      cout<<"PGA after = "<<vx<<" "<<vy<<" "<<vz<<" "<<ini_Dir<<" "<<tmp3v.mag()<<" "<<ini_Dir.theta()<<" "<<ini_Dir.phi()<<endl;	
			//      cout <<"ini_Dir "<<ini_Dir<<endl;	


        particleGun->SetParticleMomentumDirection(ini_Dir);
	particleGun->SetParticleMomentum( tmp3v.mag());
	particleGun->SetParticlePosition(G4ThreeVector(vx,vy,vz));

      // particleGun->SetParticleMomentumDirection(G4ThreeVector(-0.0456152,-0.304171,-0.951525));
      // particleGun->SetParticleMomentum(2.82895);
      // 	  particleGun->SetParticlePosition(G4ThreeVector(582.008 ,-8133.14 ,-242.995));
	  particleGun->GeneratePrimaryVertex(anEvent);
	  if (count < (int)pAnalysis->ngenmx) {
	    pAnalysis->pidin [count]= particleGun->GetParticleDefinition()->GetPDGEncoding();//_cpid[count];
	    pAnalysis->posxin[count]= particleGun->GetParticlePosition().x();
	    pAnalysis->posyin[count]= particleGun->GetParticlePosition().y();
	    pAnalysis->poszin[count]= particleGun->GetParticlePosition().z();
	    if(particle->GetPDGCharge()==0){
	      pAnalysis->momin[count] = particleGun->GetParticleMomentum()/GeV;
	    } else {
	      pAnalysis->momin[count] = (particleGun->GetParticleMomentum())*(particle->GetPDGCharge())/GeV;
	    }
	    pAnalysis->thein[count] = particleGun->GetParticleMomentumDirection().theta();
	    pAnalysis->phiin[count] = particleGun->GetParticleMomentumDirection().phi();
	  } // if (count < (int)pAnalysis->ngenmx) {
	} // if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
      } // for (int count=0; count<iC_npart; count++) {
    } // if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
    // Corsika Event By Event



 }













  else {
    //Future NuGenerators
  }
  }//InputOutput<3
}

double micalPrimaryGeneratorAction::GetCosmicEnergy(double ELimLow, double ELimUp) {

  double selectedEnergy = 0.001;
  bool energyBool = true;
  double Pmax =0.00364709;//17.7098;// 0.00364709(Alkofer);//168.23;// 0.019372(alkofer);//168.23;
  
  double E1 = ELimLow*0.001;
  double E2 = ELimUp*0.001;

  while(energyBool) {
    double x1 = E1 + (G4UniformRand()*(E2-E1));
    double x2 = G4UniformRand();
    double f1 = energy_func(x1);
    double f2 = x2*Pmax;
    // cout<<"E "<<E1<<" "<<E2<<endl;
    // cout<<"f2 "<<f2<<" "<<x2<<" "<<f1<<" "<<x1<<endl;
    if(f2 < f1) {
      selectedEnergy = x1;
      energyBool = false;
    }
  }
  
  // return 1.5*1000;
  return selectedEnergy*1000.;
  
}

int micalPrimaryGeneratorAction::LinePlaneInt(double* Line, double* Plane, double* Point) { //, double &Dist) {
  double a, b, Dist;
  int ok = 0;
  b = Line[3]*Plane[3]  + Line[4]*Plane[4] + Line[5]*Plane[5];
  ok= (fabs(b) > 1e-10) ? 1 : 0;
  
  if(ok==1){
    a=(Plane[0]-Line[0])*Plane[3] + (Plane[1]-Line[1])*Plane[4] + (Plane[2]-Line[2])*Plane[5];
    Dist = a/b;
    Point[0] = Line[0] + Line[3]*Dist;
    Point[1] = Line[1] + Line[4]*Dist;
    Point[2] = Line[2] + Line[5]*Dist;
  }else{Point[0]=0; Point[1]=0; Point[2]=0;}
  
  return ok;
}

double micalPrimaryGeneratorAction::energy_func(double xx) {
  //double paren[7]={1.44021e-02,4.09983e-01,2.96593e-01,7.46643e-03,-2.20450e+00,1.15441e+00,-4.19923e+00};
  double paren[7]={9.20933e-03,4.57797e-01,1.00738e+00 ,4.82396e-03,-1.54754e+00,1.60996e-01 ,-2.90163e+00}; //Pmax = 0.00364709;
  if (xx<1.5) {
    return paren[0]*(TMath::Gaus(xx, paren[1], paren[2], kTRUE));
  } else if (xx<14.0) {
    return paren[3]*pow(xx, paren[4]);
  } else {
    return paren[5]*pow(xx, paren[6]);
  }
} 

void micalPrimaryGeneratorAction::OpenFileCORSIKA() {
  G4String infile;
  //  infile = CorsikaFileDir;
  //  infile.append(CorsikaFileName);
  //  cout<<"CorsikaFileName "<<CorsikaFileName<<endl;
  //  FileCORSIKA = new TFile(infile,"READ","input file");

  infile = FluxFileDir;
  infile.append(FluxFileName);
  cout<<"FluxFileName "<<FluxFileName<<endl;
  FileFLUX = new TFile(infile,"READ","input file");

}

void micalPrimaryGeneratorAction::CloseFileCORSIKA() {
    
  FileFLUX->Close();
  delete FileFLUX;
  //  FileCORSIKA->Close();
  //  delete FileCORSIKA;
}

void micalPrimaryGeneratorAction::SetInputFlag(G4int p) {
  InputFlag = p;
  // cout<<"Setting Input Flag..."<<InputFlag<<endl;
  
}

void micalPrimaryGeneratorAction::SetIncPosition(G4ThreeVector p) {
  incPosition = p;
  // cout<<"Setting new inc pos "<<incPosition<<endl;
}
