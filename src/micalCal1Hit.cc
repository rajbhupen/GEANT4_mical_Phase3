
#include "micalCal1Hit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include <iomanip>
using namespace std;

G4Allocator<micalcal1Hit> micalcal1HitAllocator;

micalcal1Hit::micalcal1Hit() {
  pdgid=-25;
  edep = 0;
  toff = 1000000;
  
}

micalcal1Hit::~micalcal1Hit()
{;}

micalcal1Hit::micalcal1Hit(const micalcal1Hit &right)
  : G4VHit(),
   pos(G4ThreeVector())
{
  pdgid  = right.pdgid;
  edep = right.edep;
  pos = right.pos;
  mom = right.mom;
  toff = right.toff;
  HitId = right.HitId;
}

const micalcal1Hit& micalcal1Hit::operator=(const micalcal1Hit &right) {
  pdgid = right.pdgid;
  edep = right.edep;
  pos = right.pos;
  mom = right.mom;
  toff = right.toff;
  HitId = right.HitId;

  return *this;
}

G4int micalcal1Hit::operator==(const micalcal1Hit &right) const {
  return (this==&right) ? 1 : 0;
}

void micalcal1Hit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(8.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0.,0.,1.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void micalcal1Hit::Print() {
  G4cout<<"hit "<<HitId<<" "<<pos<<endl;
  
  
}


