#include "TGeoManager.h"

void gdml(){

  cout<<"Check0"<<endl;
  TGeoManager* gGeoManager = new TGeoManager("myGeo","");
    cout<<"Check0"<<endl;
    gGeoManager->Import("geo_mical_world20rpc.gdml");
    cout<<"Check"<<endl;
  int steps = 36;
  double dz = 0;
  double epsilon = 1.e-6;
  for(int i=0;i<steps;i++){
    gGeoManager->InitTrack(0,-973.2,-185.712-dz,0,0,-1);

    gGeoManager->FindNextBoundary();
  dz+=epsilon+ gGeoManager->GetStep();
  cout<<"Volume: "<< gGeoManager->GetCurrentVolume()->GetName()<<" Stepsize: "<<dz<<" "<<gGeoManager->GetStep()<<endl;
  

  


    }
  delete gGeoManager;

}
