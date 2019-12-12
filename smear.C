// decay: pi+ -> mu+
// decay: mu+ -> e+
// decay : kaon- -> e- + pi0
// decay: pi- -> gamma + gamma (hadron inelastic)

// TEST MACRO FOR USE WITH OLDER ROOT6.  DOESN"T WORK WHEN CLING KNOWS ABOUT
// THE VARIOUS CLASSES.
#include <TFile.h>
#include "TSystem.h"
#include <TGeoManager.h>
#include <iostream>
#include <TRandom3.h>
#include "TH1.h"
#include "TH2.h"
#include "TGeoTrd2.h"
#include "TGeoTube.h"
#include <TTree.h>
#include "TDatabasePDG.h"

#include "TG4Event.h"

#include <fstream>
#include <sstream>
#include <utility>
#include <functional>
#include <cassert>
#include <map>
#include <cstdlib>
#include <assert.h> 
using std::cout;
using std::endl;
TTree* gEDepSimTree = NULL;
TGeoManager *geo;
TG4Event* event;
double sigmas=200E-6; // m
double B=0.6; 
double x0=2.8; // m 
TRandom3 *ran;
TH2 *herr_dipAngle_stt;
TH2 *herr_pt_stt;
TH2 *herr_p_stt;
TH1 *herr_E_ecal;
TH1 *herr_theta_ecal;
TH1 *herr_phi_ecal;
TH1 *herr_theta_N;
TH1 *herr_phi_N;
TH1 *herr_p_equa_N;
TH1 *herr_p_beta_N;
TH1 *herr_p_pi0;
TH1 *herr_theta_pi0;
TH1 *herr_phi_pi0;

TFile *outf;
TFile *outTreeF;
TTree * tree;

bool isnumuBar_Htarget;

TDatabasePDG *dbpdg;

TH2* hPi0_mom_recotrue;
TH2* hPi0_ang_recotrue;
TH1* hNeutron_ang_sigma_stt;
TH1* hNeutron_ang_sigma_ecal;
TH2* hNeutron_beta_recotrue_stt;
TH2* hNeutron_beta_recotrue_ecal;
//std::ofstream treefile;
std::map<int,std::pair<int,int> > sttMap;
//std::map<int,std::pair<int,int> > ecalMap;
//std::map<int,std::pair<int,int> > sttMap_prim;
//std::map<int,std::pair<int,int> > ecalMap_prim;
std::map<int, std::vector<int> > sttMap_prim2;
std::map<int, std::vector<int> > ecalMap_prim2;


const int kNPmax = 300;
int         brNPar;
int         brNPrim;
double      brRecoP4[kNPmax][4];
double      brTrueP4[kNPmax][4];
int         brPdg   [kNPmax];
int         brTrackId[kNPmax];
int         brParentId[kNPmax];
int         brTopParentId[kNPmax];
double      brLen[kNPmax];
int         brNXhit[kNPmax];
int         brNYhit[kNPmax];
char        brInfo[kNPmax][10];
int iFillPar;
void smearPar(int trackid, std::string name);
void findEvis_forCell(int starthit, int nhit, std::map<int, double> &cellId_Evis);
void findEvis_forCell(std::vector<int> allhits, std::map<int, double> &cellId_Evis);

void cleanBranch(){
  brNPar=0;
  for(int i=0;i<kNPmax;i++){
    brRecoP4[i][0]=-999;
    brRecoP4[i][1]=-999;
    brRecoP4[i][2]=-999;
    brRecoP4[i][3]=-999;
    brTrueP4[i][0]=-999;
    brTrueP4[i][1]=-999;
    brTrueP4[i][2]=-999;
    brTrueP4[i][3]=-999;
    brPdg[i]=-999;
    brTrackId[i]=-999;
    brParentId[i]=-999;
    brTopParentId[i]=-999;
    brLen[i]=-999;
    brNXhit[i]=-999;
    brNYhit[i]=-999;
    // brInfo[i]=0;
    strcpy(brInfo[i],"");
  }  
}

void showAll(){
  std::cout<<"============================================="<<std::endl;
  for (std::vector<TG4Trajectory>::iterator
         t = event->Trajectories.begin();
       t != event->Trajectories.end(); ++t) {
    std::cout << "   Traj " << t->TrackId;
    std::cout << " " << t->ParentId;
    std::cout << " " << t->Name;
    std::cout << " " << t->Points.size();
    std::cout<< " E:"<<t->GetInitialMomentum().E();
    std::cout<<" beginpro:"<<t->Points.begin()->Process<<" "<<t->Points.begin()->Subprocess<<" endpro:"<<(t->Points.end()-1)->Process<<" "<<(t->Points.end()-1)->Subprocess;
    std::cout << std::endl;
  }
  for (auto d = event->SegmentDetectors.begin();
       d != event->SegmentDetectors.end(); ++d) {
    std::cout << "   det " << d->first;
    std::cout << " " << d->second.size();
    int count = 10;
    std::cout << " up to " << count << " segments";
    std::cout << std::endl;
    //    if(d->first!="Straw") continue;
    int i=0;
    for (std::vector<TG4HitSegment>::iterator
	   h = d->second.begin();
	 h != d->second.end();
	 ++h) {
      std::cout << "      "<<i;
      i++;
      std::cout << " P: " << h->PrimaryId<<" "<<h->Contrib[0];
      std::cout << " E: " << h->EnergyDeposit;
      std::cout << " S: " << h->SecondaryDeposit;
      std::cout << " C: " << h->Contrib.size()<<"->";
      for(unsigned long j=0;j<h->Contrib.size();j++){
	std::cout<<" "<<h->Contrib[j];
      }
      //      std::cout<<" name:"<<h->GetVolName();
      //            std::cout << " L: " << h->TrackLength;
      TLorentzVector mid= (h->Start+h->Stop)*0.5;
      TString name=geo->FindNode(mid.X(),mid.Y(),mid.Z())->GetName();
      std::cout<<" "<<name;
      std::cout<<" start:"<<h->Start.X()<<" "<<h->Start.Y()<<" "<<h->Start.Z()<<" "<<h->Start.T()<<" endT:"<<h->Stop.T();
      if((h+1)!= d->second.end() && (h+1)->Start.T()<h->Start.T()) std::cout<<"   !!!!!!! time reverted";
      std::cout<<std::endl;
    }
  }

  std::cout<<"============================================="<<std::endl;
}

int  findTopParent(int trackid){
  TG4Trajectory trk=event->Trajectories[trackid];
  while(trk.ParentId!=-1){
    trk=event->Trajectories[trk.ParentId];
  }
  return trk.TrackId;
}

void fill1Par2tree(double recoPx, double recoPy, double recoPz, int trackid, double len, int nXhit, int nYhit, const char info[10]){
  brRecoP4[iFillPar][0]=recoPx;
  brRecoP4[iFillPar][1]=recoPy;
  brRecoP4[iFillPar][2]=recoPz;
  int pdg=event->Trajectories[trackid].PDGCode;
  double m=dbpdg->GetParticle(pdg)->Mass()*1000;  // MeV
  brRecoP4[iFillPar][3]=sqrt(recoPx*recoPx+recoPy*recoPy+recoPz*recoPz+m*m);

  brTrueP4[iFillPar][0]=event->Trajectories[trackid].InitialMomentum.X();
  brTrueP4[iFillPar][1]=event->Trajectories[trackid].InitialMomentum.Y();
  brTrueP4[iFillPar][2]=event->Trajectories[trackid].InitialMomentum.Z();
  brTrueP4[iFillPar][3]=event->Trajectories[trackid].InitialMomentum.E();
  brTrackId[iFillPar]= trackid;
  brLen[iFillPar]= len;
  brNXhit[iFillPar]= nXhit;
  brNYhit[iFillPar]= nYhit;
  brParentId[iFillPar]= event->Trajectories[trackid].ParentId;
  brPdg[iFillPar]= pdg;
  brTopParentId[iFillPar]= findTopParent(trackid);
  strcpy(brInfo[iFillPar], info);
  iFillPar++;

}


class Node{
public:
  TG4Trajectory *Traj;
  Node *Parent;
  Node *FirstChild;
  Node *RightSibling;
  Node(TG4Trajectory *traj, Node *parent, Node *firstChild, Node *rightSibling): Traj(traj),Parent(parent),FirstChild(firstChild),RightSibling(rightSibling){}
};
Node *root;
std::vector<int> verlines;

Node *findNodeFast(int trackid){
  if(trackid==-1) return root;
  std::vector<int> parents;
  parents.push_back(trackid);
  TG4Trajectory *traj=&event->Trajectories[trackid];
  while(1){
    if(traj->ParentId==-1) break;
    parents.push_back(traj->ParentId);
    traj=&event->Trajectories[traj->ParentId];
  }
  Node *no=root->FirstChild;
  while(no){
    while(no){
      if(no->Traj->TrackId==parents.back()) { parents.pop_back(); break;}	
      no=no->RightSibling;
    }
    if(parents.empty()) return no;
    no=no->FirstChild;
  }
  return 0;
}
Node *findNode(Node *topNode, int trackid){
  if(trackid==-1) { /*std::cout<<"primary particle"<<std::endl; */return topNode;}
  if(topNode->Traj!=0) 
    if(trackid==topNode->Traj->TrackId) 
      return topNode;
  if(topNode->FirstChild==0) return 0;
  Node *next=topNode->FirstChild;
  while(1){    
    Node *m=findNode(next,trackid);
    if(m!=0) return m;
    if(next->RightSibling==0) break;
    next=next->RightSibling;
  }
  return 0;
}

void insertNode(Node *n){

  Node *parent=findNodeFast(n->Traj->ParentId);
  if(parent==0) std::cout<<"cannot find parent, something must be wrong"<<std::endl;
  n->Parent=parent;
  if(parent->FirstChild==0) {  parent->FirstChild=n; return;}
  Node *next=parent->FirstChild;
  while(1){
    if(next->RightSibling==0) { next->RightSibling=n; break;}
    next=next->RightSibling;    
  }
}

void makeTree(){
  root=new Node(0,0,0,0);
  for (unsigned long i=0;i<event->Trajectories.size();i++){
    insertNode(new Node(&(event->Trajectories[i]),0,0,0));
  }
}
bool findit(int i, std::vector<int> nums){
  for(unsigned long j=0;j<nums.size();j++){
    if(nums[j]==i) return true;
  }
  return false;
}


void drawNode(Node *n, int icol){
  if(icol==0) std::cout<<"_";
  else std::cout<<"__"<<n->Traj->Name<<"("<<n->Traj->TrackId<<")";
  if(n->FirstChild==0) return;
  Node *c=n->FirstChild;
  if(c->RightSibling!=0) { verlines.push_back(icol);}
  while(1){
    if(c->Parent->FirstChild!=c &&  c->RightSibling==0) {verlines.pop_back();}
    drawNode(c, icol+ 2+(c->Traj->Name).size()+ std::to_string(c->Traj->TrackId).size()+2);
    //    if(c->RightSibling==0) {std::cout<<"one pop"<<std::endl;verlines.pop_back();}
    if(c->RightSibling==0) break;
    c=c->RightSibling;
    std::cout<<std::endl;
    for(int i=0;i<=verlines.back();i++){
      if(findit(i,verlines)) std::cout<<"|";
      else std::cout<<" ";        
    }
  }

}

void drawTree(){
  verlines.clear();
  drawNode(root, 0);
  std::cout<<std::endl;
}
void dumpTree(){
  makeTree();
  drawTree();
}

void findEvis_forCell(std::vector<int> hitchains, std::map<int, double> &cellId_Evis){
  int id, cellID, planeID, modID;
  cellId_Evis.clear();
  assert(hitchains.size()%2==0);
  //  if(hitchains.size()%2!=0) std::exit(EXIT_FAILURE);
  for(int i=0;i<hitchains.size()/2;i++){
    for(unsigned int j = hitchains[2*i]; j <= hitchains[2*i+1]; j++){
      const TG4HitSegment& h = event->SegmentDetectors["ECAL"].at(j);
      double x = 0.5*(h.Start.X()+h.Stop.X());
      double y = 0.5*(h.Start.Y()+h.Stop.Y());
      double z = 0.5*(h.Start.Z()+h.Stop.Z());
      TGeoNode* node = geo->FindNode(x,y,z);
      TString slabstr = node->GetName();
      TString modstr=geo->GetMother()->GetName();
      //    std::cout<<"slabstr:"<<slabstr<<std::endl;
      if(slabstr.Contains("volECALActiveSlab") == true)
	{
	TObjArray* obj1 = slabstr.Tokenize("_");  //volECALActiveSlab_125_PV_0
	TObjArray* obj2 = modstr.Tokenize("_");  //ECAL_lv_PV_19

	int slabID;
	modID  = ((TObjString*) obj2->At(3))->GetString().Atoi();
	slabID = ((TObjString*) obj1->At(1))->GetString().Atoi();

	//    std::cout<<"modID:"<<modID<<" slabID:"<<slabID<<std::endl;
	delete obj1;
	delete obj2;
	// planeID==0 -> smallest slab
	// planeID==208 -> biggest slab
	planeID = slabID/40;

	if (planeID > 4) planeID = 4;
	double Pmaster[3]={x,y,z};
	double Plocal[3];
	geo->GetCurrentNavigator()->MasterToLocal(Pmaster,Plocal);
	//    std::cout<<"Plocal[0]:"<<Plocal[0]<<std::endl;  // along circular second longest
	//    std::cout<<"Plocal[1]:"<<Plocal[1]<<std::endl;  // along column longest
	//    std::cout<<"Plocal[2]:"<<Plocal[2]<<std::endl; //along radial, short
	TGeoTrd2* trd = (TGeoTrd2*) node->GetVolume()->GetShape();

	double dx1 = trd->GetDx1();  // shorter one along circumferential
	double dx2 = trd->GetDx2();  // longer one along circumferential
	double dz  = trd->GetDz();   // half thickness , along radial
	//	double dy1 = trd->GetDy1();  // along axial/fiber , same as dy2
	//	d1 = dy1 + Plocal[1];
	//	d2 = dy1 - Plocal[1];
	double dx = (dx2 - dx1) / dz * Plocal[2];
	double dis= Plocal[0]>0? (dx1+dx2)/2. + Plocal[0] - dx/2.: (dx1+dx2)/2. + Plocal[0] + dx/2.;
	double cellw = (dx1+dx2) / 12.;
	cellID = dis / cellw;
      }
    else if(slabstr.Contains("vol_endECALActiveSlab") == true)
      {

	TObjArray* obja = slabstr.Tokenize("_");
	int slabID;
	modID  = x>0?30:40;
	slabID = ((TObjString*) obja->At(1))->GetString().Atoi();

	delete obja;

	planeID = (208 - slabID)/40;

	if (planeID > 4) planeID = 4;

        double Pmaster[3]={x,y,z};
	double Plocal[3];
	//      std::cout<<"x:"<<x<<" y:"<<y<<" z:"<<z<<std::endl;
	geo->GetCurrentNavigator()->MasterToLocal(Pmaster,Plocal);

	TGeoTube* tub = (TGeoTube*) node->GetVolume()->GetShape();
	double rmax = tub->GetRmax();
	//      double dz  = tub->GetDz();
	// Plocal[0] : horizontal distance to tube center
	// Plocal[1] : vertical distance
	//	d1 = rmax * TMath::Sin(TMath::ACos(Plocal[0]/rmax)) - Plocal[1]; // d1 is shorter distance here
	//	d2 = rmax * TMath::Sin(TMath::ACos(Plocal[0]/rmax)) + Plocal[1];

	cellID = int((Plocal[0]/rmax + 1.) * 45);
      }
    else { continue;}
    id = cellID + 100 * planeID + 1000 * modID;
    //    std::cout<<"cellID:"<<cellID<<" planeID:"<<planeID<<" modID:"<<modID<<std::endl;
    cellId_Evis[id]= h.EnergyDeposit;
    
    } // for 
  }

}

void smearFirstHitPosition(int starthit, double sigmaX, double sigmaY, double *pos_smear){
  sigmaX*=10.; //cm -> mm
  sigmaY*=10.; //cm --> mm
  //  std::cout<<"start smearFirstHitPosition starthit"<<starthit<<std::endl;
  const TG4HitSegment& h = event->SegmentDetectors["ECAL"].at(starthit);
  double x = 0.5*(h.Start.X()+h.Stop.X());
  double y = 0.5*(h.Start.Y()+h.Stop.Y());
  double z = 0.5*(h.Start.Z()+h.Stop.Z());

  TGeoNode* node = geo->FindNode(x,y,z);
  double Pmaster[3]={x,y,z};
  double Plocal[3];

  geo->GetCurrentNavigator()->MasterToLocal(Pmaster,Plocal);
  ////// barrel ///////
  ////// barrel ///////

  TString slabstr = node->GetName();
  if(slabstr.Contains("volECALActiveSlab") == true) {
    double Plocal_smear[3]= { Plocal[0]+ran->Gaus(0,sigmaX), Plocal[1]+ran->Gaus(0,sigmaY), Plocal[2]};
    //////////////////////// boundary check 
    TGeoTrd2* trd = (TGeoTrd2*) node->GetVolume()->GetShape();
    double dx1 = trd->GetDx1();  // shorter one along circumferential
    double dx2 = trd->GetDx2();  // longer one along circumferential
    double dz  = trd->GetDz();   // half thickness , along radial
    double dy1 = trd->GetDy1();  // along axial/fiber , same as dy2
    //    std::cout<<"dx1:"<<dx1<<" dx2:"<<dx2<<" dz:"<<dz<<" dy1:"<<dy1<<std::endl;
    double maxX= (dx1+dx2)/2. + Plocal[2]*(dx2-dx1)/dz;
    if(Plocal_smear[0] > maxX) Plocal_smear[0]=maxX;
    else if(Plocal_smear[0] < -maxX) Plocal_smear[0]=-maxX;
    if(Plocal_smear[1]> dy1 ) Plocal_smear[1]= dy1;
    else   if(Plocal_smear[1]< -dy1 ) Plocal_smear[1]=-dy1;
    geo->GetCurrentNavigator()->LocalToMaster(Plocal_smear, pos_smear);
  }
  else if(slabstr.Contains("vol_endECALActiveSlab") == true){
    ////// need to know which simga to use for endcap ???????? ////////
    ////// now temperatly use same as barrel, sigma ??????? /////////
    double Plocal_smear[3]= { Plocal[0]+ran->Gaus(0,sigmaX), Plocal[1]+ran->Gaus(0,sigmaY), Plocal[2]};
    // boundary check 
    TGeoTube* tub = (TGeoTube*) node->GetVolume()->GetShape();
    double rmax = tub->GetRmax();
    double width=rmax*2./90.;
    double low=floor(Plocal[0]/width)*width;
    double high=ceil(Plocal[0]/width)*width;
    if(Plocal_smear[0]>high) Plocal_smear[0]=high;
    else     if(Plocal_smear[0]<low) Plocal_smear[0]=low;
    double ymax= sqrt(rmax*rmax+low*low);
    if(ymax<sqrt(rmax*rmax+high*high)) ymax=sqrt(rmax*rmax+high*high);
    if(Plocal_smear[1]>ymax) Plocal_smear[1]=ymax;
    else if(Plocal_smear[1] <-ymax) Plocal_smear[1]=-ymax;
    geo->GetCurrentNavigator()->LocalToMaster(Plocal_smear, pos_smear);
  }

  return;
  
}

bool isthisParent(int lowId, int highId){
  if(lowId==highId) return true;
  int parentId=event->Trajectories[lowId].ParentId;
  while(parentId!=-1){
    if(parentId==highId) return true;    
    parentId=event->Trajectories[parentId].ParentId;
  }
  return false;
}

void bruteforceFindHit_ecal(int trackid, int primid, std::vector<int> &allhits){
  allhits.clear();
  assert(primid!=-1);
  assert(trackid!=primid);
  std::vector<int> hitchains=ecalMap_prim2[primid];
  for(int i=0;i<hitchains.size()/2;i++){
    for(unsigned int j = hitchains[2*i]; j <= hitchains[2*i+1]; j++){
      if( !isthisParent(event->SegmentDetectors["ECAL"].at(j).Contrib[0], trackid)) { continue;}
      //      std::cout<<"j:"<<j<<" check1pass"<<std::endl;
      if(allhits.size()==0) {allhits.push_back(j);allhits.push_back(j);}
      if(i-allhits.back()==1) allhits.back()=j;
      else
	{allhits.push_back(j);allhits.push_back(j);}	
    }
  }
}

bool smearPar_ecal(int trackid, int primaryId=-1){
  //  std::cout<<"start ecal smear"<<std::endl;
  //  5.7%/sqrt(E)  (GeV)
  int starthit;
  //  int nhit;  
  std::vector<int> allhits;
  std::map<int, double> cellId_Evis;
  if(event->Trajectories[trackid].ParentId==-1 || event->Trajectories[event->Trajectories[trackid].ParentId].Name=="pi0"){
    if(ecalMap_prim2.find(trackid)==ecalMap_prim2.end()) return false;
    findEvis_forCell(ecalMap_prim2[trackid], cellId_Evis);
    starthit=ecalMap_prim2[trackid][0];
  }
  else {
    //    std::cout<<"-------> start bruteforceFindHit"<<std::endl;
    if(primaryId==-1) {
      int parentId=event->Trajectories[trackid].ParentId;
      int grandId=(parentId==-1)?-1:event->Trajectories[parentId].ParentId;
      assert(event->Trajectories[parentId].Name=="gamma");
      if(grandId==-1 || event->Trajectories[grandId].Name=="pi0") primaryId=parentId;
      else { assert(event->Trajectories[grandId].Name=="gamma"); primaryId=grandId;}
    }
    bruteforceFindHit_ecal(trackid, primaryId,  allhits);
    //    std::cout<<"-------> after bruteforceFindHit"<<std::endl;
    if(allhits.size()==0) { /* std::cout<<"ecal smear failed, due to zero ecal hit"<<std::endl;*/ return false;}
    findEvis_forCell(allhits, cellId_Evis);
    starthit=allhits[0]; 
  }
  bool haveDetectableHit=false;
  for(auto cellE: cellId_Evis){
    if(cellE.second > 0.1) {
      haveDetectableHit=true;
      break;
    }
  }
  
  if(!haveDetectableHit) { 
    //    std::cout<<"ecal smear failed due to no detectable ecal hit"<<std::endl; 
    return false;}
  double E=event->Trajectories[trackid].InitialMomentum.E();
  //  double P=event->Trajectories[trackid].InitialMomentum.P();
  double de2e= 0.057/sqrt(E/1000.);
  double E_smear= E*ran->Gaus(1,de2e); // MeV
  int pdg=event->Trajectories[trackid].PDGCode;
  double m=dbpdg->GetParticle(pdg)->Mass() * 1000; // MeV
  double P_smear=sqrt(E_smear*E_smear-m*m); // MeV

  double p0                        =    0.0034793;
  double p1                        =   -0.0151348;
  double p2                        =      1.52382;
  double p3                        =      0.57308;  
  double q0                        =     -1.93866;
  double q1                        =      13.5211;
  // sigmaX fitting function : [0]*x*exp([1]*x+[2])+[3]"
  // sigma Z fitting function:  [0]*log(x)+[1]
  //sigmaX ,sigmaZ, plot extracted from paper <<A FLUKA simulation of the KLOE electromagnetic calorimeter>> 
  // fit on my own 
  //  std::cout<<"smear ecal 4"<<std::endl;
  double sigmaX= p0*E*exp(p1*E+p2)+p3;    // cm  along circumferential
  double sigmaZ= q0*log(E)+q1;  // cm  along fiber/axial 
  double Pos_smear[3];
  smearFirstHitPosition( starthit, sigmaX, sigmaZ, Pos_smear); // sigmaZ is actually sigmaY in local coordinate of Ttrd2
  //  std::cout<<"smear ecal 5"<<std::endl;
  TVector3 firstHitPos_smear(Pos_smear);
  TVector3 vtx=event->Primaries[0].Position.Vect(); 
  TVector3 dir=firstHitPos_smear-vtx;  
  TVector3 P3_smear=P_smear*dir.Unit();
  //  std::cout<<"ecal smear succeeded, fill it"<<std::endl;
  double theta=event->Trajectories[trackid].InitialMomentum.Theta();
  double phi=event->Trajectories[trackid].InitialMomentum.Phi();
  double theta_smear=dir.Theta();
  double phi_smear=dir.Phi();
  herr_E_ecal->Fill((E_smear-E)/E*100);
  herr_theta_ecal->Fill((theta_smear-theta)/theta*100);
  herr_phi_ecal->Fill((phi_smear-phi)/phi*100);
  
  fill1Par2tree(P3_smear.X(), P3_smear.Y(), P3_smear.Z(),  trackid, -999, -999, -999, "ecalsmear ");
}

bool smearChargedPar_stt(int trackid){
  //  std::cout<<"stt smearing start --------------"<<trackid<<std::endl;
  // Ptran RMS from the equation, 
  // the angle between Ptran and Px smear by PDG multiscattering-RMS-equation, in which Px is decided.  
  // angle between Py and Pz is smeared by Roberto-provided equation with multiple-scattering(second) term replaced by PDG one
  //  std::cout<<"start smearChargedPar_stt:"<<trackid<<std::endl;
  int nYhit=0;
  int nXhit=0;
  double Lyz=0;
  double L=0;
  double Lx=0;
  //EDEP FACT: for hitsegment, the primary Id will be its top parent's track id, unless if there's "decay" process happens, the decayed daughters will start to become top parent
  //EDEP FACT: same track's higsegment most time are id-connected,form a long chain,if later broke, only a few hits left,  thoese later hits are wierd hits, themselves are not time-connected, space-connected, better disregard.   for the long chain, most time they are time-forwarding, but the last a few hits may time-reverse, these hits are also wierd hits. 

  // EDEP FACT: for pion0->gamma, gamma will use its own trackid as primaryid, but gamma->e+ e-, e+ e- will use gammas trackid as their primary ID

  if(sttMap.find(trackid)==sttMap.end()) {
    //    std::cout<<"zero stt hit, stt smear fail!"<<std::endl;
    return false;
  }

  TLorentzVector prePos, postPos;
  double dy,dz,dx;
  unsigned int ihit=sttMap[trackid].first;
  int nhit=sttMap[trackid].second;
  //  if(ihit>0) { if(event->SegmentDetectors["Straw"].at(ihit-1).Contrib[0]==trackid) { std::cout<<" wrong"; std::exit(EXIT_FAILURE);}}
  //  if(event->SegmentDetectors["Straw"].size()>(ihit+nhit) ) {if(event->SegmentDetectors["Straw"].at(ihit+nhit).Contrib[0]==trackid) {std::cout<<" wrong"; std::exit(EXIT_FAILURE);}}
  
  TG4HitSegment  h= event->SegmentDetectors["Straw"].at(ihit);
  unsigned int i=(ihit+1);
  if( (event->SegmentDetectors["Straw"].at(ihit+1).Start.T()<h.Start.T() || event->SegmentDetectors["Straw"].at(ihit+1).Stop.T()<h.Stop.T() ) && h.Contrib.size()>1) 
    { h= event->SegmentDetectors["Straw"].at(ihit+1); i=(ihit+2);}
  //  const TG4HitSegment& hseg = ev->SegmentDetectors["Straw"].at(j);
  TLorentzVector mid= (h.Start+h.Stop)*0.5;
  prePos=mid;
  TString name=geo->FindNode(mid.X(),mid.Y(),mid.Z())->GetName();
  if(name.Contains("horizontal")) nYhit++;
  else nXhit++;
  for( ;i<(ihit+nhit);i++){

    h=event->SegmentDetectors["Straw"].at(i);
    postPos= (h.Start+h.Stop)*0.5;
    //    if(h.Contrib[0]!=trackid) { std::cout<<" !!!!!!!!!!!!!!!!!!!!! wrong trackid"<<trackid<<" ihit"<<ihit<<" nhit:"<<nhit<<std::endl; showAll();  std::exit(EXIT_FAILURE);}
    //    if(h.Contrib.size()>1) {  std::cout<<" contribution more than 2 tracks ihit:"<<ihit<<" nhit:"<<nhit<<" i:"<<i<<std::endl;  continue;}
    if(postPos.T()< prePos.T()) { 																					
      //      if((i-ihit)*1.0<=0.5*nhit) std::cout<<"trackid:"<<trackid<<" time reverse, cut! ihit:"<<ihit<<" nhit+ihit:"<<nhit+ihit<<" i:"<<i<<std::endl;
      assert((i-ihit)*1.0>0.5*nhit || nhit<9);
      break;}
    if(h.EnergyDeposit<250E-6) continue;
    //    postPos= (h.Start+h.Stop)*0.5;
    dx= postPos.X()-prePos.X();
    dy= postPos.Y()-prePos.Y();
    dz= postPos.Z()-prePos.Z();
    name=geo->FindNode(postPos.X(),postPos.Y(),postPos.Z())->GetName();
    if(name.Contains("horizontal")) nYhit++;
    else nXhit++;
    Lyz+= sqrt(dy*dy+dz*dz);
    L+= sqrt(dx*dx+dy*dy+dz*dz);    
    prePos=postPos;
  }

  if(nYhit<4)  { 
    //    std::cout<<"stt smear failed "<<nYhit<<std::endl;
    return false;
  } // fill1Par2tree(-999,-999,-999, trackid, L, nXhit, nYhit, "sttFail") ; return false;}
  //  std::cout<<"nYhit:"<<nYhit<<" nXhit:"<<nXhit<<" postT:"<<postPos.T()<<std::endl;
  /*
  if(h==event->SegmentDetectors["Straw"].end()) return false;
  h++;
  std::cout<<" broken points time:";
  for( ;h != event->SegmentDetectors["Straw"].end(); h++){
    if(h->Contrib[0] ==trackid ) std::cout<<" "<<h->Start.T();
  }
  std::cout<<std::endl;
  */
  TVector3 initP=event->Trajectories[trackid].InitialMomentum.Vect();
  double Pt=sqrt(initP.Y()*initP.Y()+initP.Z()*initP.Z());
  double Px=initP.X();
  double P=initP.Mag();
  double dipAng=atan(Px/Pt);
  double thetaYZ=atan(initP.Y()/initP.Z());
  //  std::cout<<"Lyz:"<<Lyz<<" P:"<<P<<std::endl;
  Lx=sqrt(L*L-Lyz*Lyz);
  L/=1000;
  Lyz/=1000.; //mm-> m
  Lx/=1000.; //mm-> m
  Pt/=1000.;  // Mev-->GeV
  Px/=1000.; 
  P/=1000.;
  //  double dPt2Pt=sqrt(pow(sigmas*Pt/0.3/B/Lyz/Lyz*sqrt(720./(nYhit+4)),2)+pow(0.045/B/sqrt(Lyz*x0),2));
  double dPt2Pt=sqrt(pow(sigmas*Pt/0.3/B/L/L*sqrt(720./(nYhit+4)),2)+pow(0.045/B/sqrt(L*x0),2));
  double Pt_smear=Pt*ran->Gaus(1, dPt2Pt);
  double sigma_dipAng=13.6E-3/P*sqrt(L/x0)*(1+0.038*log(L/x0)); // from pdg
  //  double sigma_dipAng2=13.6E-3/P*sqrt(Lx/x0)*(1+0.038*log(Lx/x0)); // from pdg
  //  double sigma_dipAng3=13.6E-3/P*sqrt(Lyz/x0)*(1+0.038*log(Lyz/x0)); // from pdg
  double dipAng_smear=dipAng + ran->Gaus(0,sigma_dipAng);
  //  double dipAng_smear2=dipAng + ran->Gaus(0,sigma_dipAng2);
  //  double dipAng_smear3=dipAng + ran->Gaus(0,sigma_dipAng3);
  double Px_smear=Pt*tan(dipAng_smear);
  double sigma_thetaYZ=13.6E-3/P*sqrt(Lyz/x0)*(1+0.038*log(Lyz/x0)); 
  double thetaYZ_smear=thetaYZ + ran->Gaus(0,sigma_thetaYZ);
  double Py_smear=Pt_smear*sin(thetaYZ_smear);
  double Pz_smear=Pt_smear*cos(thetaYZ_smear);
  double P_smear=sqrt(Pt_smear*Pt_smear+Px_smear*Px_smear);
  //  double Px_smear2=Px + Pt*tan(ran->Gaus(0,sigma_thetaYZ));
  double namecode;
  if(event->Trajectories[trackid].Name=="mu+" || event->Trajectories[trackid].Name=="mu-")
    namecode=0;
  else if(event->Trajectories[trackid].Name=="proton" || event->Trajectories[trackid].Name=="anti_proton")
    namecode=1;
  else if(event->Trajectories[trackid].Name=="pi+" || event->Trajectories[trackid].Name=="pi-")
    namecode=2;
  else if(event->Trajectories[trackid].Name=="e+" || event->Trajectories[trackid].Name=="e-")
    namecode=3;
  else if(event->Trajectories[trackid].Name=="kaon+" || event->Trajectories[trackid].Name=="kaon-")
    namecode=4;
  else 
    namecode=5;
  
  herr_dipAngle_stt->Fill(dipAng_smear-dipAng, namecode);
  herr_pt_stt->Fill((Pt_smear-Pt)/Pt*100,namecode);
  herr_p_stt->Fill((P_smear-P)/P*100,namecode);
  //  std::cout<<"sigma_dipAng:"<<sigma_dipAng<<std::endl;
  //  std::cout<<"dPt2Pt:"<<dPt2Pt<<" sigma_thetaX/thetaX:"<<sigma_thetaX/thetaX<<" sigma_thetaYZ/thetaYZ:"<<sigma_thetaYZ/thetaYZ<<std::endl;
  //  std::cout<<"sigma_thetaX:"<<sigma_thetaX<<" sigma_thetaYZ:"<<sigma_thetaYZ<<std::endl;
  //  std::cout<<"Pt:"<<Pt<<" Px:"<<Px<<std::endl;
  //  std::cout<<"Pt_smear:"<<Pt_smear<<" Px_smear:"<<Px_smear<<" Px_smear2:"<<Px_smear2<<" px ratio:"<<Px_smear/Px<<" 2:"<<Px_smear2/Px<<std::endl;
  //  std::cout<<"---------iFillPar:"<<iFillPar<<" Px_smear:"<<Px_smear<<" Py_smear:"<<Py_smear<<std::endl;
  //  std::cout<<"trackid:"<<trackid<<" parentId:"<<event->Trajectories[trackid].ParentId<<" pdg:"<<pdgMap[event->Trajectories[trackid].Name]<<std::endl;
  //  std::cout<<"stt smear succeed, L"<<L<<" nXhit:"<<nXhit<<" nYhit:"<<nYhit<<std::endl;
  fill1Par2tree(Px_smear*1000., Py_smear*1000., Pz_smear*1000., trackid, L*1000, nXhit, nYhit, "sttsmear  "); // always use MeV to fill
  
  return true;
}

bool smearChargedPar(int trackid){  
  bool sttSmear=false;
  bool ecalSmear=false;
  sttSmear=smearChargedPar_stt(trackid);

  if(!sttSmear) {
    //    std::cout<<"stt smear fail, will start ecal smear"<<std::endl;
    ecalSmear=smearPar_ecal(trackid); 
  }
  return (sttSmear || ecalSmear);  
}

void smearRemnantGamma(int trackid, int primaryId=-1){
  
  if(primaryId==-1) {
    int parentId=event->Trajectories[trackid].ParentId;
    int grandId=(parentId==-1)?-1:event->Trajectories[parentId].ParentId;
    //    int grand2Id=(grandId==-1)?-1:event->Trajectories[grandId].ParentId;
    if(parentId==-1) primaryId=trackid;
    else if(event->Trajectories[parentId].Name=="pi0") primaryId=trackid;
    else {std::cout<<"remnant gamma parent becoming primaryId check:"<<event->Trajectories[parentId].Name<<std::endl; primaryId=parentId;}
  }
  smearPar_ecal(trackid, primaryId); // this remnant gamma very very not possible to have any stt hits, so just do ecal smearing

}

void smearGamma(int trackid){
  // common cases:
  //  gamma daughters: 50 e- 51 e+
  //  gamma daughters: 19 e- 20 gamma
  //  gamma daughters:
  //  gamma daughters: 9 e-
  //  gamma daughters: 9 e- 10 e- 11 e+
  // gamma daughters: 45 e- 46 e- 47 e-
  // gamma daughters: 135 e- 136 e-
  // gamma daughters: 22 e- 23 e- 24 e- 25 e-
  // rare cases:
  //    gamma daughters: 85 neutron 86 neutron
  //    gamma daughters: 44 proton 45 neutron 46 neutron
  //    gamma daughters: 5 neutron 6 neutron 7 pi+ 8 neutron 9 e-
  //    gamma daughters: 65 pi- 66 pi+ 67 proton
  //    gamma daughters: 9 gamma
  //    gamma daughters: 322 proton 323 pi- 324 pi+ 342 gamma 343 gamma

  // only the last rare case, the daughters not connected, all other cases, daughters connected.
  //  std::cout<<"start to smear gamma:"<<trackid<<std::endl;
  int parentId=event->Trajectories[trackid].ParentId;
  int grandId=(parentId==-1)?-1:event->Trajectories[parentId].ParentId;
  int grand2Id=(grandId==-1)?-1:event->Trajectories[grandId].ParentId;
  if(parentId!=-1 && grandId!=-1 && event->Trajectories[parentId].Name=="gamma" && event->Trajectories[grandId].Name=="gamma") { 
    if(grand2Id==-1) {smearRemnantGamma(trackid, grandId); return;}
    else { std::cout<<"need to check this!!!: "<<grand2Id<<" -> "<<grandId<<" -> "<<parentId<<" -> "<< trackid<<std::endl; smearRemnantGamma(trackid, grand2Id); return;}
  }
  if(parentId!=-1 && grandId!=-1 && event->Trajectories[parentId].Name=="gamma"&& event->Trajectories[grandId].Name=="pi0") { smearRemnantGamma(trackid, parentId); return;}

  Node *no=findNodeFast(trackid);
  no=no->FirstChild;
  if(no==0) { smearRemnantGamma(trackid); return;}
  
  while(no){
    //    std::cout<<"start to smear 1 gamma daughter: "<<no->Traj->TrackId<<"  "<<no->Traj->Name<<std::endl;
    if(no->Traj->Name!="neutron") 
      smearPar(no->Traj->TrackId, no->Traj->Name);
    no=no->RightSibling;
  }
}

void smearDaughters(int trackid){
  Node *no=findNodeFast(trackid);
  no=no->FirstChild;
  while(no){
    //  std::cout<<"start to smear 1 gamma daughter: "<<no->Traj->TrackId<<"  "<<no->Traj->Name<<std::endl;
    if(no->Traj->Name!="neutron")
      smearChargedPar_stt(no->Traj->TrackId);  // for these daughters, should only do STT smear, if only ecal hits, even not possible to recognize them
    no=no->RightSibling;
  }
}

void smearSTT_or_Daughters(int trackid){ 
  
  bool sttSucceed=smearChargedPar_stt(trackid);
  if(sttSucceed) return;
  else 
    smearDaughters(trackid);

}
void smearPi0_external(int trackid){
  double P=event->Trajectories[trackid].InitialMomentum.P()/1000.;
  double Phi=event->Trajectories[trackid].InitialMomentum.Phi(); // -pi , pi 
  double Theta=event->Trajectories[trackid].InitialMomentum.Theta(); // 0 - pi
  double lowP=0.;
  double highP=4.; // GeV
  int  nPbin=100;
  double lowAng=0;
  double highAng=180;
  int nAngbin=180;
  int iPbin=TMath::CeilNint(P/(highP-lowP)*nPbin);
  int iThetabin=TMath::CeilNint(Theta/(highAng-lowAng)*nAngbin);
  int iPhibin=TMath::CeilNint(abs(Phi)/(highAng-lowAng)*nAngbin);

  if(iPbin>nPbin) std::cout<<"-------> out of range P:"<<P<<std::endl;
  //  std::cout<<" get random P start"<<"iPbin:"<<iPbin<<std::endl;
  double P_smear= P<4?hPi0_mom_recotrue->ProjectionY("", iPbin,iPbin)->GetRandom():P; // GeV
  if(iPbin>37){ // the stats is too little, have to do in  this way, right now, treat the reco mean same as true mean
    double rms;
    if(iPbin<=39)
      rms=hPi0_mom_recotrue->ProjectionY("",38,39)->GetRMS();
    else if(iPbin<=41)
      rms=hPi0_mom_recotrue->ProjectionY("",40,41)->GetRMS();
    else if(iPbin<=45)
      rms=hPi0_mom_recotrue->ProjectionY("",42,45)->GetRMS();
    else if(iPbin<=50)
      rms=hPi0_mom_recotrue->ProjectionY("",46,50)->GetRMS();
    else if(iPbin<=56)
      rms=hPi0_mom_recotrue->ProjectionY("",51,56)->GetRMS();
    else if(iPbin<=69)
      rms=hPi0_mom_recotrue->ProjectionY("",57,69)->GetRMS();
    else
      rms=hPi0_mom_recotrue->ProjectionY("",70,100)->GetRMS();
    P_smear=P+ran->Gaus(0,rms);
  }

  double Theta_smear= hPi0_ang_recotrue->ProjectionY("", iThetabin,iThetabin)->GetRandom();
  //  std::cout<<"get phi"<<std::endl;
  double Phi_smear= hPi0_ang_recotrue->ProjectionY("", iPhibin,iPhibin)->GetRandom();
  //  std::cout<<" get random end ---<"<<std::endl;
  if(Phi<0) Phi_smear*=-1.;

  double Pz_smear=P_smear*cos(Theta_smear);
  double Px_smear=P_smear*sin(Theta_smear)*cos(Phi_smear);
  double Py_smear=P_smear*sin(Theta_smear)*sin(Phi_smear);
  
  fill1Par2tree(Px_smear*1000, Py_smear*1000, Pz_smear*1000 , trackid, -999, -999, -999, "pi0extern ");
  herr_p_pi0->Fill((P_smear-P)/P*100);
  herr_theta_pi0->Fill((Theta_smear-Theta)/Theta*100);
  herr_phi_pi0->Fill((Phi_smear-Phi)/Phi*100);
  
}

void smearPi0(int trackid){
  // only two mode: gamma + gamma / gamma e+ e-, for each mode, no one will miss
  //for the two modes, all the daughters  are always connected, trackid connected.

  // if both gamma convert in ECAL, use external theta true2reco to smear angle, external E true2reco to smear energy
  // if one gamma convert in STT, use charged-track way to reconstruct this gamma, the other gamma convert in ECAL, use 5.7%/sqrt(E) to smear E, use sigmaZ, sigmaX from KLOE ECAL Paper to smear angle
  
  Node *no=findNodeFast(trackid);
  no=no->FirstChild;
  int nChild=0;
  std::string names[3];
  int trackids[3];
  while(no){
    names[nChild]= no->Traj->Name;
    trackids[nChild]=no->Traj->TrackId;
    no=no->RightSibling;
    nChild++;
  }
  if(nChild==2) {
    assert(names[0]=="gamma" && names[1]=="gamma");
    if(sttMap_prim2.find(trackids[0])==sttMap_prim2.end() && sttMap_prim2.find(trackids[1])==sttMap_prim2.end()){ // <--- here need to use primaryId_organized map
      //      std::cout<<"no any hit for pi0 two gammas, use external way to smear pi0"<<std::endl;
      smearPi0_external(trackid);
    }
    else{
      //      std::cout<<" ==pi0==> smear pi0 first gamma: id:"<<trackids[0]<<std::endl;
      smearGamma(trackids[0]);
      //      std::cout<<" ==pi0==> smear pi0 second gamma: id:"<<trackids[1]<<std::endl;
      smearGamma(trackids[1]);
    }
  }
  else if(nChild==3) {
    //    std::cout<<"  ==pi0==> three children for pi0, smear one by one"<<std::endl;
    smearPar(trackids[0], names[0]);
    smearPar(trackids[1], names[1]);
    smearPar(trackids[2], names[2]);
  }
  else std::cout<<" ==pi0==> more than 3 children from pi0, check!!!!"<<std::endl;
}

bool smearN_byEquation(double &psmear){  // only for antinumu events with 1 neutron produced 
  //  std::cout<<"smear neutron by equation"<<std::endl;
  const double mpr = dbpdg->GetParticle(2212)->Mass()*1000;
  //  const double mpip = 0.13957018;
  //  const double mpi0 = 0.1349766;
  const double mmu = dbpdg->GetParticle(13)->Mass()*1000;
  const double mn = dbpdg->GetParticle(2112)->Mass()*1000;
  

  TLorentzVector p4hadreco(0,0,0,0);  
  if(brPdg[0]!=-13) {std::cout<<"the first smeared track is not mu+, something must be wrong, check!!!!"<<std::endl; return false;}
  TLorentzVector p4mureco(brRecoP4[0][0],brRecoP4[0][1],brRecoP4[0][2],brRecoP4[0][3]);
  for(int i=1;i<iFillPar;i++){
    p4hadreco+=TLorentzVector(brRecoP4[i][0],brRecoP4[i][1],brRecoP4[i][2],brRecoP4[i][3]);
  }
  double en = 0.5*( mmu*mmu + pow(p4hadreco.M(),2) + mpr*mpr - mn*mn - 2*mpr*(p4mureco.E() + p4hadreco.E()) +
  		    2*p4mureco*p4hadreco)/(p4mureco.E() + p4hadreco.E() - p4mureco.Pz() - p4hadreco.Pz() - mpr);
  en = en - p4mureco.E() - p4hadreco.E() + mpr;
  //  std::cout<<"en:"<<en<<" mn:"<<mn<<std::endl;
  //  assert(en>mn);
  if(en<mn) return false;
  psmear=sqrt(en*en-mn*mn);

  return true;
}

bool smearNeutron(int trackid){
  if(event->Trajectories[trackid].ParentId!=-1)  {
    std::cout<<" about to smear a non-primary neutron, stop right now, check why it happens"<<std::endl;
    return false;
  }
  bool isSTTdetectable=false;
  bool isECALdetectable=false;
  if(sttMap_prim2.find(trackid)!=sttMap_prim2.end()){
    const std::vector<int> &vec= sttMap_prim2[trackid];
    for(int i=0;i<vec.size()/2;i++){
      for(int j=vec[2*i];j<=vec[2*i+1];j++){
	if(event->SegmentDetectors["Straw"].at(j).EnergyDeposit>250E-6) {isSTTdetectable=true; break;}
      }
      if(isSTTdetectable) break;
    }
    //    if(!isSTTdetectable) std::cout<<"this neutron has straw hits but none of them are detectable, check!!!"<<std::endl;
  }
  if(!isSTTdetectable) {
    if(ecalMap_prim2.find(trackid)!=ecalMap_prim2.end()){
      std::map<int, double> cellId_Evis;
      findEvis_forCell(ecalMap_prim2[trackid], cellId_Evis);
      for(auto cellE: cellId_Evis){
	if(cellE.second > 0.1) { isECALdetectable=true;	  break;}
      }      
      //      if(!isECALdetectable) std::cout<<"this neutron has ecal hits but none of them are detectable, check!!!"<<std::endl;  
      //      else { std::cout<<"this neutron does have ecal hits"<<std::endl;}
    }
  }

  if( !isSTTdetectable && !isECALdetectable) {
    //    std::cout<<"smear neutron fail due to no detectable stt hit and ecal hit"<<std::endl;
    return false;
  }
  

  double Phi=event->Trajectories[trackid].InitialMomentum.Phi(); // -pi , pi
  double Theta=event->Trajectories[trackid].InitialMomentum.Theta(); // 0 - pi
  double Phi_smear=isSTTdetectable?Phi*ran->Gaus(1,hNeutron_ang_sigma_stt->GetRandom()):Phi*ran->Gaus(1,hNeutron_ang_sigma_ecal->GetRandom());
  double Theta_smear=isSTTdetectable?Theta*ran->Gaus(1,hNeutron_ang_sigma_stt->GetRandom()):Theta*ran->Gaus(1,hNeutron_ang_sigma_ecal->GetRandom());

  double P_smear;
  double P=event->Trajectories[trackid].InitialMomentum.P();
  bool PequationSmearSucceed=false;
  if(isnumuBar_Htarget) PequationSmearSucceed=smearN_byEquation(P_smear);
  if(PequationSmearSucceed)   herr_p_equa_N->Fill((P_smear-P)/P*100);
  if(!PequationSmearSucceed || !isnumuBar_Htarget)  {
    double E=event->Trajectories[trackid].InitialMomentum.E();
    double m=event->Trajectories[trackid].InitialMomentum.Mag();
    double beta=P/E;
    int iTrueBetaBin=TMath::CeilNint(beta/(1.-0.)*100);
    double beta_smear=isSTTdetectable?hNeutron_beta_recotrue_stt->ProjectionY("",iTrueBetaBin,iTrueBetaBin)->GetRandom(): hNeutron_beta_recotrue_ecal->ProjectionY("",iTrueBetaBin,iTrueBetaBin)->GetRandom();
    P_smear=m*beta_smear/sqrt(1-beta_smear*beta_smear);    
    //    std::cout<<"isSTTdetectable: "<<isSTTdetectable<<" ---------------------------------------------->P:"<<P<<" P_smear:"<<P_smear<<" beta:"<<beta<<" beta_smear:"<<beta_smear<<" iTrueBetaBin:"<<iTrueBetaBin<<" sttbin entry:"<<hNeutron_beta_recotrue_stt->ProjectionY("",iTrueBetaBin,iTrueBetaBin)->GetEntries()<<" ecal:"<<hNeutron_beta_recotrue_ecal->ProjectionY("",iTrueBetaBin,iTrueBetaBin)->GetEntries()<<std::endl;
    herr_p_beta_N->Fill((P_smear-P)/P*100);
  }
  double Pz_smear=P_smear*cos(Theta_smear);
  double Px_smear=P_smear*sin(Theta_smear)*cos(Phi_smear);
  double Py_smear=P_smear*sin(Theta_smear)*sin(Phi_smear);
  //  std::cout<<"neutron smear succeed"<<std::endl;
  if(PequationSmearSucceed) fill1Par2tree(Px_smear, Py_smear, Pz_smear , trackid, -999, -999, -999, "NsmearEqua");
  else fill1Par2tree(Px_smear, Py_smear, Pz_smear , trackid, -999, -999, -999, "NsmearBeta");

  herr_theta_N->Fill((Theta_smear-Theta)/Theta*100);
  herr_phi_N->Fill((Phi_smear-Phi)/Phi*100);
  

  return true;
}

void organizeHits(){
  sttMap.clear();
  //  ecalMap.clear();
  int pretrackid;
  int posttrackid;
  int nhit;
  int istart;

  if(event->SegmentDetectors["Straw"].size()>0){    
    pretrackid=event->SegmentDetectors["Straw"].begin()->Contrib[0];
    nhit=1;
    istart=0;
    posttrackid=pretrackid;
    for(unsigned int i=1; i<event->SegmentDetectors["Straw"].size(); i++){
      posttrackid=event->SegmentDetectors["Straw"].at(i).Contrib[0];
      //      std::cout<<"posttrackid:"<<posttrackid<<" nhit:"<<nhit<<std::endl;
      if(posttrackid==pretrackid) { nhit++;continue;}
      if(sttMap.find(pretrackid) ==sttMap.end()) 
	sttMap[pretrackid]= std::make_pair(istart,nhit);
      if(sttMap.find(posttrackid) ==sttMap.end()){
	nhit=1;
	istart=i;
      }    
      pretrackid=posttrackid;
    }
    if(sttMap.find(posttrackid) ==sttMap.end())
      sttMap[posttrackid]= std::make_pair(istart,nhit);
  }
  for(auto it = sttMap.begin(); it != sttMap.end(); ) {
    if(it->second.second <4)
      it = sttMap.erase(it);
    else
      ++it;
  }

  /*
  if(event->SegmentDetectors["ECAL"].size()>0){
    pretrackid=event->SegmentDetectors["ECAL"].begin()->Contrib[0];
    nhit=1;
    istart=0;
    posttrackid=pretrackid;
    for(unsigned long i=1; i<event->SegmentDetectors["ECAL"].size(); i++){
      posttrackid=event->SegmentDetectors["ECAL"].at(i).Contrib[0];
      //      std::cout<<"posttrackid:"<<posttrackid<<" nhit:"<<nhit<<std::endl;
      if(posttrackid==pretrackid) { nhit++;continue;}
      if(ecalMap.find(pretrackid) ==ecalMap.end())
        ecalMap[pretrackid]= std::make_pair(istart,nhit);
      if(ecalMap.find(posttrackid) ==ecalMap.end()){
        nhit=1;
        istart=i;
      }
      pretrackid=posttrackid;
    }
    if(ecalMap.find(posttrackid) ==ecalMap.end())
      ecalMap[posttrackid]= std::make_pair(istart,nhit);
  }
  */
}
/*
void organizeHits_prim(){
  sttMap_prim.clear();
  ecalMap_prim.clear();
  int pretrackid;
  int posttrackid;
  int nhit;
  int istart;
  if(event->SegmentDetectors["Straw"].size()>0){
    pretrackid=event->SegmentDetectors["Straw"].begin()->PrimaryId;
    nhit=1;
    istart=0;
    posttrackid=pretrackid;
    for(unsigned int i=1; i<event->SegmentDetectors["Straw"].size(); i++){
      posttrackid=event->SegmentDetectors["Straw"].at(i).PrimaryId;
      //      std::cout<<"posttrackid:"<<posttrackid<<" nhit:"<<nhit<<std::endl;
      if(posttrackid==pretrackid) { nhit++;continue;}
      if(sttMap_prim.find(pretrackid) ==sttMap_prim.end())
        sttMap_prim[pretrackid]= std::make_pair(istart,nhit);
      if(sttMap_prim.find(posttrackid) ==sttMap_prim.end()){
        nhit=1;
        istart=i;
      }
      pretrackid=posttrackid;
    }
    if(sttMap_prim.find(posttrackid) ==sttMap_prim.end())
      sttMap_prim[posttrackid]= std::make_pair(istart,nhit);
  }

  if(event->SegmentDetectors["ECAL"].size()>0){
    pretrackid=event->SegmentDetectors["ECAL"].begin()->PrimaryId;
    nhit=1;
    istart=0;
    posttrackid=pretrackid;
    for(unsigned long i=1; i<event->SegmentDetectors["ECAL"].size(); i++){
      posttrackid=event->SegmentDetectors["ECAL"].at(i).PrimaryId;
      //      std::cout<<"posttrackid:"<<posttrackid<<" nhit:"<<nhit<<std::endl;
      if(posttrackid==pretrackid) { nhit++;continue;}
      if(ecalMap_prim.find(pretrackid) ==ecalMap_prim.end())
        ecalMap_prim[pretrackid]= std::make_pair(istart,nhit);
      if(ecalMap_prim.find(posttrackid) ==ecalMap_prim.end()){
        nhit=1;
        istart=i;
      }
      pretrackid=posttrackid;
    }
    if(ecalMap_prim.find(posttrackid) ==ecalMap_prim.end())
      ecalMap_prim[posttrackid]= std::make_pair(istart,nhit);
  } 
}
*/
void organizeHits_prim2(){
  sttMap_prim2.clear();
  ecalMap_prim2.clear();
  int trackid;
  for(int i=0; i<event->SegmentDetectors["Straw"].size(); i++){
    trackid=event->SegmentDetectors["Straw"].at(i).PrimaryId;
    if(sttMap_prim2.find(trackid)==sttMap_prim2.end()) { sttMap_prim2[trackid].push_back(i); sttMap_prim2[trackid].push_back(i); continue;}
    if(i-sttMap_prim2[trackid].back()==1)
      sttMap_prim2[trackid].back()=i;
    else
      {sttMap_prim2[trackid].push_back(i); sttMap_prim2[trackid].push_back(i);}    
  } // for
  
  for(int i=0; i<event->SegmentDetectors["ECAL"].size(); i++){
    trackid=event->SegmentDetectors["ECAL"].at(i).PrimaryId;
    if(ecalMap_prim2.find(trackid)==ecalMap_prim2.end()) { ecalMap_prim2[trackid].push_back(i); ecalMap_prim2[trackid].push_back(i);continue;}
    if(i-ecalMap_prim2[trackid].back()==1)
      ecalMap_prim2[trackid].back()=i;
    else
      {ecalMap_prim2[trackid].push_back(i); ecalMap_prim2[trackid].push_back(i);}
  } // for

}
void showHitMap(){
  std::cout<<"-------stt map show:"<<std::endl;
  for(auto pp:sttMap){
    std::cout<<pp.first<<" "<<pp.second.first<<" "<<pp.second.second<<std::endl;
  }
  /*
  std::cout<<"-------ecal map show:"<<std::endl;
  for(auto pp:ecalMap){
    std::cout<<pp.first<<" "<<pp.second.first<<" "<<pp.second.second<<std::endl;
  }
  */
  std::cout<<"-------stt prim2 map show:"<<std::endl;
  for(auto pp:sttMap_prim2){
    std::cout<<pp.first<<" ->";
    assert(pp.second.size()%2==0);
    for(auto ii:pp.second)
      std::cout<<" "<<ii;
    std::cout<<std::endl;
  }
  std::cout<<"-------ecal prim2 map show:"<<std::endl;
  for(auto pp:ecalMap_prim2){
    std::cout<<pp.first<<" ->";
    assert(pp.second.size()%2==0);
    for(auto ii:pp.second)
      std::cout<<" "<<ii;
    std::cout<<std::endl;
  }
  
}
void smearPar(int trackid, std::string name){

  if(name=="mu-" || name=="mu+" || name=="proton" || name=="pi+" || name=="pi-") smearChargedPar(trackid);
  else if(name=="e-" || name=="e+") smearChargedPar(trackid);
  else if(name=="kaon+" || name=="kaon-") smearChargedPar(trackid);
  else if(name=="pi0") smearPi0(trackid);
  else if(name=="neutron") smearNeutron(trackid);
  else if(name=="gamma") smearGamma(trackid);  
  else if(name=="kaon0" || name=="anti_kaon0" || name=="kaon0L") smearDaughters(trackid);
  else if(name=="lambda" || name=="anti_lambda")  smearDaughters(trackid);   // lambda travel very short distance (5cm) then decay to pi- and proton
  else if(name=="sigma+" || name=="sigma-") smearSTT_or_Daughters(trackid); // itself can create medium track then decay 
  else if(name=="anti_proton") smearChargedPar(trackid);
  else if(name=="nu_e"|| name=="anti_nu_e" || name=="anti_nu_mu" || name=="nu_mu" || name=="C12" || name=="Ar40" || name=="deuteron") return;
  else if(name=="anti_neutron") smearNeutron(trackid);
  else std::exit(EXIT_FAILURE);

  //  else std::cout<<"--->unknown par:"<<name<<std::endl;
}

bool checkPrimNeutron(bool &NeutronAtLast, int &Neutron_trackId){

  Neutron_trackId=-1;
  NeutronAtLast=false;
  int nPrim= event->Primaries.begin()->Particles.size();
  int i=0;
  for (;i<nPrim;i++){    
    if(event->Primaries.begin()->Particles[i].Name=="neutron")
      { Neutron_trackId=i; NeutronAtLast=(i==(nPrim-1)); return true;}
  }
  return false;

}

void checkit(){
  for (std::vector<TG4PrimaryVertex>::iterator
	 v = event->Primaries.begin();
       v != event->Primaries.end(); ++v) {
    int ipar=0;
    //    std::cout << "  particles " << v->Particles.size();
    int nn=0;
    int npr=0;
    for (std::vector<TG4PrimaryParticle>::iterator
	   p = v->Particles.begin();
	 p != v->Particles.end(); ++p) {
      //      std::cout << " " << p->Name;
      if(p->Name=="proton") npr++;
      else if(p->Name=="neutron") nn++;
      ipar++;
    }
    //    std::cout << endl;
    if(nn>1) std::cout<<"***************************************nn: "<<nn<<"  npr:"<<npr<<std::endl;
  }

}

void smearEvent(int iEvent){
  cleanBranch();
  //  treefile.open(Form("treemap%d.txt",i));
  gEDepSimTree->GetEntry(iEvent);
  //  std::cout<<" ############################################################## new event ####################################  "<<i<<std::endl;
  //  showAll();
  organizeHits();
  //  organizeHits_prim();
  organizeHits_prim2();
  //  showHitMap();
  makeTree();
  //  dumpTree();
  //  treefile.close();
  iFillPar=0;
  int nPrim= event->Primaries.begin()->Particles.size();
  bool havePrimNeutron=false;
  bool NeutronAtLast=false;
  int neutron_trackid;
  if(isnumuBar_Htarget) havePrimNeutron=checkPrimNeutron(NeutronAtLast, neutron_trackid);
  if(isnumuBar_Htarget && havePrimNeutron && (!NeutronAtLast)) {
    for (int i=0;i<nPrim;i++){      
      if(event->Primaries.begin()->Particles[i].TrackId==neutron_trackid) continue;
      //      std::cout<<" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++smearPar: "<<event->Primaries.begin()->Particles[i].Name<<" id: "<<event->Primaries.begin()->Particles[i].TrackId<<std::endl;
      smearPar(event->Primaries.begin()->Particles[i].TrackId , event->Primaries.begin()->Particles[i].Name);
    }
    //    std::cout<<"smear neutron at the very end since its a numubar hydrogen target event"<<std::endl;
    smearNeutron(neutron_trackid);
  }
  else {
    for (int i=0;i<nPrim;i++){
      //      std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++smearPar: "<<event->Primaries.begin()->Particles[i].Name<<" id: "<<event->Primaries.begin()->Particles[i].TrackId<<std::endl;
      smearPar(event->Primaries.begin()->Particles[i].TrackId , event->Primaries.begin()->Particles[i].Name);
    }
  }
  //  std::cout<<"iFillPar:"<<iFillPar<<std::endl;
  if(iFillPar==0) return;
  brNPar=iFillPar;
  brNPrim=nPrim;
  tree->Fill();

}

std::vector<std::string> makefilelist(std::string st,int Nfilelist=0){
  // if Nfilelist is 0, then the default means input all the lines/files  i have.
  std::vector<std::string> files;
  files.clear();

  std::ifstream filelist(st.c_str());
  //count how many lines in the filelistlist
  int num=std::count(std::istreambuf_iterator<char>(filelist), std::istreambuf_iterator<char>(), '\n');

  if(Nfilelist==0)
    Nfilelist=num;
  else if(Nfilelist>num)
    std::cout<<"we don't have that many filelist, please change the  number."<<std::endl;
  //the following close and reopen is important, cause if you didn't close, right now the pointer is at the end of the filelistlist
  filelist.close();
  filelist.open(st.c_str());

  std::string onefile;
  for(int i=0;i<Nfilelist;i++){
    getline(filelist,onefile);
    files.push_back(onefile);
  }
  filelist.close();
  return files;
}//////////////////////////////////////////////////////////////////


int main(int argc, char *argv[]){
  if(argc!=3) { std::cout<<"two arguments needed"<<std::endl;return 0;}
  //  outf=new TFile("outf.root","recreate");
  TFile *fpi0=new TFile("data/Histograms_1pi0_complete.root");
  TFile *fneutron_ang=new TFile("data/plotres.root");
  TFile *fneutron_beta=new TFile("data/RecoVsTrue_Beta_RealCal_pdg2112_20190315_192325_ThrStt2.500000e-07.root");
  hPi0_mom_recotrue=(TH2*)fpi0->Get("h_mom_recotrue");
  hPi0_ang_recotrue=(TH2*)fpi0->Get("h_arctg_recotrue");
  hNeutron_ang_sigma_stt=(TH1*)fneutron_ang->Get("STT Resolution");
  hNeutron_ang_sigma_ecal=(TH1*)fneutron_ang->Get("Calorimeter Resolution");
  hNeutron_beta_recotrue_stt=(TH2*)fneutron_beta->Get("Beta STT");
  hNeutron_beta_recotrue_ecal=(TH2*)fneutron_beta->Get("Beta calorimeter");

  hPi0_mom_recotrue->Smooth();
  hPi0_mom_recotrue->Smooth();
  hPi0_ang_recotrue->Smooth();
  hPi0_ang_recotrue->Smooth();
  hNeutron_ang_sigma_stt->Smooth();
  hNeutron_ang_sigma_ecal->Smooth();
  hNeutron_beta_recotrue_stt->Smooth();
  hNeutron_beta_recotrue_ecal->Smooth();
  
  outTreeF=new TFile(argv[2],"recreate");
  TBranch * brStdHepPdg=0;
  int  StdHepPdg[100];
  TTree *rootrackerTree;

  tree = new TTree("edep_smeared_tree"," reconstructed on edep sim");
  tree->Branch("NPar", &brNPar, "NPar/I");
  tree->Branch("NPrim", &brNPrim, "NPrim/I");
  tree->Branch("RecoP4",    brRecoP4,        "RecoP4[NPar][4]/D"); 
  tree->Branch("TrueP4",    brTrueP4,        "TrueP4[NPar][4]/D");
  tree->Branch("Pdg",    brPdg,        "Pdg[NPar]/I");
  tree->Branch("TrackId",    brTrackId,        "TrackId[NPar]/I");
  tree->Branch("ParentId",    brParentId,        "ParentId[NPar]/I");
  tree->Branch("TopParentId",    brTopParentId,        "TopParentId[NPar]/I");
  tree->Branch("Len",    brLen,        "Len[NPar]/D");
  tree->Branch("NXhit",    brNXhit,        "NXhit[NPar]/I");
  tree->Branch("NYhit",    brNYhit,        "NYhit[NPar]/I");
  tree->Branch("Info",    brInfo,        "Info[NPar][10]/C");

  herr_dipAngle_stt=new TH2F("herr_dipAngle_stt","",200,-0.05,0.05,10,0,10); // rad
  herr_pt_stt=new TH2F("herr_pt_stt","",100, -30,30,10,0,10); // percent
  herr_p_stt=new TH2F("herr_p_stt","",100, -30,30,10,0,10); //percent
  herr_E_ecal=new TH1F("herr_E_ecal","",100,-30,30); 
  herr_theta_ecal=new TH1F("herr_theta_ecal","",100,-30,30);
  herr_phi_ecal=new TH1F("herr_phi_ecal","",100,-30,30);
  herr_theta_N=new TH1F("herr_theta_N","",100,-30,30);
  herr_phi_N=new TH1F("herr_phi_N","",100,-30,30);
  herr_p_equa_N=new TH1F("herr_p_equa_N","",100,-50,50);
  herr_p_beta_N=new TH1F("herr_p_beta_N","",100,-50,50);
  herr_p_pi0=new TH1F("herr_p_pi0","",100,-50,50);
  herr_theta_pi0=new TH1F("herr_theta_pi0","",100,-50,50);
  herr_phi_pi0=new TH1F("herr_phi_pi0","",100,-50,50);

  ran=new TRandom3(0); // 0 will always give different result when you recreate it
  dbpdg = new TDatabasePDG();
  gSystem->Load("libGeom");  


  gFile=new TFile(argv[1]);
  geo = (TGeoManager*) gFile->Get("EDepSimGeometry");
  gEDepSimTree = (TTree*) gFile->Get("EDepSimEvents");
  gEDepSimTree->SetBranchAddress("Event",&event);    
  rootrackerTree=(TTree*) gFile->Get("DetSimPassThru/gRooTracker");
  brStdHepPdg    = rootrackerTree -> GetBranch ("StdHepPdg");
  brStdHepPdg -> SetAddress (StdHepPdg);
  

  outTreeF->cd();

  int nEntry=gEDepSimTree->GetEntries();
  int rootrackerEntry=rootrackerTree->GetEntries();
  if(nEntry!=rootrackerEntry) {std::cout<<"----->----->not same entries"<<std::endl; return 1;}
  for(int i=0;i<nEntry;i++){
    isnumuBar_Htarget=false;
    if(i%100==0) std::cout<<"ientry:"<<i<<std::endl;
    rootrackerTree->GetEntry(i);
    if (StdHepPdg[0]==-14 && StdHepPdg[1]==2212) 
      { isnumuBar_Htarget=true;} //  std::cout<<"^^^^^^^^^^^^^^^^^^^^^^isnumuBar_Htarget true"<<std::endl;}
    smearEvent(i);
    
  }    
  std::cout<<"close files"<<std::endl;

  tree->Write();
  herr_dipAngle_stt->Write();
  herr_pt_stt->Write();
  herr_p_stt->Write();
  herr_E_ecal->Write();
  herr_theta_ecal->Write();
  herr_phi_ecal->Write();
  herr_theta_N->Write();
  herr_phi_N->Write();
  herr_p_equa_N->Write();
  herr_p_beta_N->Write();
  herr_p_pi0->Write();
  herr_theta_pi0->Write();
  herr_phi_pi0->Write();
  //  outTreeF->Write();
  outTreeF->Close();

}
