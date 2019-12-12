{
  //  TFile *file=new TFile("plotres.root");
  //  TH1 *h1=(TH1*)file->Get("STT Resolution");
  //  h1->Draw();
  //  TFile *file=new TFile("RecoVsTrue_Beta_RealCal_pdg2112_20190315_192325_ThrStt2.500000e-07.root");
  TFile *file=new TFile("Histograms_1pi0_complete.root");
  TH2 *h1=(TH2*)file->Get("h_arctg_recotrue");
  TH2 *h2=(TH2*)file->Get("h_mom_recotrue");
  
  int nxbin=h1->GetNbinsX();
  int nybin=h1->GetNbinsY();
  for(int i=0;i<nxbin;i++){
    std::cout<<"i+1:"<<i+1<<" "<<h1->ProjectionY("",i+1,i+1)->GetEntries()<<std::endl;
  }
  
  /*

  std::cout<<h2->ProjectionY("",38,39)->GetEntries()<<std::endl;
  std::cout<<h2->ProjectionY("",40,41)->GetEntries()<<std::endl;
  std::cout<<h2->ProjectionY("",42,45)->GetEntries()<<std::endl;
  std::cout<<h2->ProjectionY("",46,50)->GetEntries()<<std::endl;
  std::cout<<h2->ProjectionY("",51,56)->GetEntries()<<std::endl;
  std::cout<<h2->ProjectionY("",57,69)->GetEntries()<<std::endl;
  std::cout<<h2->ProjectionY("",70,100)->GetEntries()<<std::endl;

  int l[7]={38,40,42,46,51,57,70};
  int r[7]={39,41,45,50,56,69,100};
  
  TCanvas *can=new TCanvas("can","",1);
  can->Print("plot.pdf(");
  for(int i=0;i<7;i++){
    h2->ProjectionY("",l[i],r[i])->Draw();
    can->Print("plot.pdf");
    std::cout<<"rms:"<<h2->ProjectionY("",l[i],r[i])->GetRMS()<<std::endl;
  }
  can->Print("plot.pdf)");
  */

}
