///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///  Shayne Reichard, 2020 May 26, shayne@physik.uzh.ch
///
///  This script determines the maximum likelihood estimators of the electroweak parameters of the elastic neutrino-electron scattering
///  process. You input a histogram and loop of the range of values.
///
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






#include <TEventList.h>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstring>
#include "TChain.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TCut.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TLine.h"
#include "TGraph.h"
#include "TFormula.h"
#include "TLegend.h"
#include "TAxis.h"
#include <math.h>
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TColor.h"
#include "TRandom3.h"
//#include "nuparentdata.h"
#include "nuFunc.h"
#include "TROOT.h"



const int Nbins=300;
const int hist_LB=0, hist_UB=3000;
const int nBands = 3;
const double tonne=1000.;  //kg
const double year=365.25;    //days

double reduction=100.;


//Load the unsmeared histograms
TFile *fTot = new TFile("./inputs/nudata_total.root"); TH1F *hTot = (TH1F*)fTot->Get("h");
TFile *fMat = new TFile("./inputs/nudata_Materials.root"); TH1F *hMat = (TH1F*)fMat->Get("hEdTotal");
TFile *fRn = new TFile("./inputs/nudata_Rn222.root"); TH1F *hRn = (TH1F*)fRn->Get("hEd_Rn222");
TFile *fKr = new TFile("./inputs/nudata_Kr85.root"); TH1F *hKr = (TH1F*)fKr->Get("hEd_Kr85");
TFile *fDB = new TFile("./inputs/nudata_DB.root"); TH1F *hDB = (TH1F*)fDB->Get("h1_Xe_1");

TFile *fSN = new TFile("./inputs/pp_recoil_spectrum_stepped.root"); TH1F *hSN = (TH1F*)fSN->Get("h_pp");
TFile *fBe7 = new TFile("./inputs/Be7_main_recoil_spectrum_stepped.root"); TH1F *hBe7 = (TH1F*)fBe7->Get("h_Be7");
TFile *fBe7_other = new TFile("./inputs/Be7_other_recoil_spectrum_stepped.root"); TH1F *hBe7_other = (TH1F*)fBe7_other->Get("h_Be7_other");
TFile *fpep = new TFile("./inputs/pep_recoil_spectrum_stepped.root"); TH1F *hpep = (TH1F*)fpep->Get("h_pep");
TFile *fN13 = new TFile("./inputs/N13_recoil_spectrum_stepped.root"); TH1F *hN13 = (TH1F*)fN13->Get("h_N");
TFile *fO15 = new TFile("./inputs/O15_recoil_spectrum_stepped.root"); TH1F *hO15 = (TH1F*)fO15->Get("h_O");

//Load the smeared histograms
TFile *fSmear = new TFile("./inputs/smear_data.root");
TFile *fPatricia = new TFile("./inputs/shayne_materials.root");

//TH1F *SMat = (TH1F*)fSmear->Get("sMat");
TH1F *SMat = (TH1F*)fPatricia->Get("sMat");
TH1F *SRn = (TH1F*)fSmear->Get("sRn");
TH1F *SKr = (TH1F*)fSmear->Get("sKr");
TH1F *SDB = (TH1F*)fSmear->Get("sDB");
TH1F *SSN = (TH1F*)fSmear->Get("sSN");
TH1F *SBe7 = (TH1F*)fSmear->Get("sBe7");
TH1F *SBe7_other = (TH1F*)fSmear->Get("sBe7_other");

//Define fit functions with unsmeared histograms
Double_t fitfMat(Double_t *v, Double_t *par) { Double_t fitval = par[0]*hMat->Interpolate(v[0]);  return fitval; }
Double_t fitfRn(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hRn->Interpolate(v[0]);  return fitval; }
Double_t fitfKr(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hKr->Interpolate(v[0]);  return fitval; }
Double_t fitfDB(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hDB->Interpolate(v[0]);  return fitval; }
Double_t fitfSN(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hSN->Interpolate(v[0]);  return fitval; }

Double_t EnergyRes(double E) {
  double a=0.313;  double b=0.002;
  return a*TMath::Power(E,0.5)+b*E;
}
//Define fit functions with smeared histograms
Double_t fitfMatSmear(Double_t *v, Double_t *par) { Double_t fitval=par[0]*SMat->Interpolate(v[0]);  return fitval; }
Double_t fitfRnSmear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SRn->Interpolate(v[0]);  return fitval; }
Double_t fitfKrSmear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SKr->Interpolate(v[0]);  return fitval; }
Double_t fitfDBSmear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SDB->Interpolate(v[0]);  return fitval; }
Double_t fitfSNSmear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SSN->Interpolate(v[0]);  return fitval; }

Double_t fitfTotSmear(Double_t *v, Double_t *par) {
  Double_t fitval = 
    par[0]*SMat->Interpolate(v[0])
    +par[1]*SRn->Interpolate(v[0])
    +par[2]*SKr->Interpolate(v[0])
    +par[3]*SDB->Interpolate(v[0])
    +par[4]*SSN->Interpolate(v[0]);
  return fitval;
}

std::string double_to_string(double inputnumber) {
  string s;
  stringstream out;
  out << inputnumber;
  s = out.str();
  return s;
}

double ConvolveGauss(TH1F* hist, double E0, int N=Nbins*10) {
  double width=(hist_UB-hist_LB)/N;
  double sum=0;
  for (int i=0; i<N; i++) {
    double E=hist_LB+width*(1.*i+0.5);
    double sigma=EnergyRes(E);
    double f=hist->Interpolate(E);
    double g=1./(sigma*sqrt(2.*TMath::Pi()))*TMath::Exp(-0.5*TMath::Power((E0-E)/sigma,2));
    sum+=f*g*width;
  }
  return sum;
}

TH1F* GetConvolvedHist(TH1F* hist, const char* name="h_tmp", int bins=300, double l=0, double u=3000) {
  double w=(u-l)/bins;
  TH1F *h_smear = new TH1F(name,name,bins,l,u);
  for (int j=1; j<bins+1; j++) {
    double E=w*(j-0.5); //cout<<"E="<<E<<endl;
    h_smear->SetBinContent(j,ConvolveGauss(hist,E));  //cout<<ConvolveGauss(hist,E)<<endl;
  }
  return h_smear;
}



TH1F* HistOfNuRate(double sin2_theta, double elec_surv=0.55, double Q=420., double flux=5.97e10, int bins=300, double l=0, double u=3000, bool mono=0) {
  TH1F* h_tmp = MakeHisto("h_pp_tmp",Q,sin2_theta,elec_surv,flux,l,u,bins,mono);
  TH1F* h_tmp_smear = GetConvolvedHist(h_tmp,"h_pp_tmp_smear",bins,l,u); h_tmp_smear->Scale(1/tonne/year);
  TH1F* h_ret = new TH1F("h_nupars","h_nupars",bins,l,u); 
  for(int k=1; k<bins+1; k++) {
    double fval =
      SMat->GetBinContent(k)
      +SRn->GetBinContent(k)
      +SKr->GetBinContent(k)
      +SDB->GetBinContent(k)
      +h_tmp_smear->GetBinContent(k);
    h_ret->SetBinContent(k,fval); //cout<<"Setting bin "<<k<<": "<<fval<<endl;
    }
  delete gROOT->FindObject("h_pp_tmp");
  delete gROOT->FindObject("h_pp_tmp_smear");
  return h_ret;
}

double GetStats(TH1D* hist, int flag) {
  double prob[nBands] = {0.1587, 0.5, 0.8413};
  double quantiles[nBands];
  hist->GetQuantiles(nBands, quantiles, prob);
  return quantiles[flag];
}


double LogFactorial(int n) {
  double retval = 0;
  for (int i=1; i<n+1; i++) retval += log(1.0*i);
  return retval;
}


double FindLL(vector<int> toy, double mult, double mu_pp=1, double sin2_pp=0.2387, double surv_pp=0.55, double mu_Be=1, double sin2_Be=0.2387, double surv_Be=0.52, int bins=300, double l=0, double u=3000) {
  //give experiment in dru

  double Q_pp=420., flux_pp=5.98e10;
  double Q_Be=862., flux_Be=5.00e9, branch=0.9;
  double Q_Be_other=384.;

  TH1F* h_tmp_pp = MakeHisto("h_pp_tmp",Q_pp,sin2_pp,surv_pp,flux_pp,l,u,bins,0);
  TH1F* h_tmp_pp_smear = GetConvolvedHist(h_tmp_pp,"h_pp_tmp_smear",bins,l,u);                                    
  TH1F* h_tmp_Be = MakeHisto("h_Be_tmp",Q_Be,sin2_Be,surv_Be,branch*flux_Be,l,u,bins,1);
  TH1F* h_tmp_Be_smear = GetConvolvedHist(h_tmp_Be,"h_Be_tmp_smear",bins,l,u);                               
  TH1F* h_tmp_Be_other = MakeHisto("h_Be_other_tmp",Q_Be_other,sin2_pp,surv_pp,(1.-branch)*flux_Be,l,u,bins,1);  //using "pp" values here is NOT a mistake
  TH1F* h_tmp_Be_other_smear = GetConvolvedHist(h_tmp_Be_other,"h_Be_other_tmp_smear",bins,l,u);

  double LLmax = 0;

  for (int j=1; j<Nbins+1; j++) {
//  for (int j=1; j<25+1; j++) {
    int n=toy[j-1]; //cout<<"  n="<<n<<endl;
    double s = mu_pp*h_tmp_pp_smear->GetBinContent(j)+mu_Be*(h_tmp_Be_smear->GetBinContent(j)+h_tmp_Be_other_smear->GetBinContent(j)); 
    //cout<<"s="<<s<<endl; //already in tonne-years
    double b = SMat->GetBinContent(j)+SRn->GetBinContent(j)+SKr->GetBinContent(j)+SDB->GetBinContent(j); //cout<<"b="<<b<<endl;
    double t = (s+b*tonne*year)*mult; //cout<<"  t="<<t<<endl;
    double inc = n*log(t)-t-LogFactorial(n);
    LLmax += inc;
//    cout<<"   The likelihood of this observation is "<<TMath::Exp(inc)<<endl;
//    cout<<j<<", "<<inc<<" ("<<LLmax<<")"<<endl;
  }

  delete gROOT->FindObject("h_pp_tmp");
  delete gROOT->FindObject("h_pp_tmp_smear");
  delete gROOT->FindObject("h_Be_tmp");
  delete gROOT->FindObject("h_Be_tmp_smear");
  delete gROOT->FindObject("h_Be_other_tmp");
  delete gROOT->FindObject("h_Be_other_tmp_smear");

  return LLmax;
}






//------------------------------------------------------------------------------------------------------------------------






void NuLikelihood(int index=0, double mult=30.0) { cout<<"Investigating an exposure of "<<mult<<" tonne-years..."<<endl;

  hRn->Scale(10./reduction);
  hKr->Scale(1./reduction); 
  hMat->Scale(1./reduction);
  hDB->Scale(1./100.);
  SRn->Scale(10./reduction);
  SKr->Scale(1./reduction); 
  SMat->Scale(1./reduction);
  SDB->Scale(1./100.);

  double exposure=mult*tonne*year;
  double FitLB=0, FitUB=3000;

  gRandom = new TRandom3();
  gRandom->SetSeed(0);

  //Define canvas for original histogram
  TCanvas *c = new TCanvas("c", "c",265,50,1009,588);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c->Range(-231.8436,-8.828753,2086.592,-0.5412262);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetLogy();
  c->SetRightMargin(0.0373494);
  c->SetTopMargin(0.05535714);
  c->SetFrameBorderMode(0);

  //hTot->Scale(exposure);  hMat->Scale(exposure);  hRn->Scale(exposure);  hKr->Scale(exposure);  hDB->Scale(exposure);  hSN->Scale(exposure);
  //SMat->Scale(exposure);  SRn->Scale(exposure);  SKr->Scale(exposure);  SDB->Scale(exposure);  SSN->Scale(exposure); //scale for fitting purposes
  hTot->GetXaxis()->SetRange(1,Nbins); hTot->SetMinimum(1.e-7); hTot->SetMaximum(1e0); hTot->SetLineColor(1);
  hTot->GetXaxis()->SetTitle("Energy  (keV)"); hTot->GetYaxis()->SetTitle("events/kg/day/10keV");
  hTot->GetXaxis()->CenterTitle(); hTot->GetYaxis()->CenterTitle();
  hSN->SetLineColor(kGreen); hSN->SetLineWidth(2);  hDB->SetLineColor(kOrange+1);
  hTot->Draw("");  hMat->Draw("same");  hRn->Draw("same");  hKr->Draw("same");  hDB->Draw("same");  hSN->Draw("same");
//  TH1F* hTest = HistOfNuPars(0.2223); hTest->SetLineColor(kCyan); hTest->Draw("same");

  TLegend *leg = new TLegend(0.48,0.15,0.63,0.41,NULL,"brNDC");
  leg->SetFillColor(0);
  leg->AddEntry("hEdTotal","Materials","l");
  leg->AddEntry("hEd_Rn222","Radon","l");
  leg->AddEntry("hEd_Kr85","Krypton","l");
  leg->AddEntry("h1_Xe_1","Double Beta","l");
  leg->AddEntry("h_pp","Solar Neutrino","l");
  leg->Draw("same");
  gPad->Modified();

  //------------------------------------------------------------------------------------

  //Fit functions (unsmeared)
  TF1 *funcMat = new TF1("funcMat",fitfMat,FitLB,FitUB,1);  funcMat->SetNpx(10*Nbins);  funcMat->SetParameter(0,1);  funcMat->SetParName(0,"Materials");
  TF1 *funcRn = new TF1("funcRn",fitfRn,FitLB,FitUB,1);  funcRn->SetNpx(10*Nbins);  funcRn->SetParameter(0,1);  funcRn->SetParName(0,"Radon");
  TF1 *funcKr = new TF1("funcKr",fitfKr,FitLB,FitUB,1);  funcKr->SetNpx(10*Nbins);  funcKr->SetParameter(0,1);  funcKr->SetParName(0,"Krypton");
  TF1 *funcDB = new TF1("funcDB",fitfDB,FitLB,FitUB,1);  funcDB->SetNpx(10*Nbins);  funcDB->SetParameter(0,1);  funcDB->SetParName(0,"Double Beta");
  TF1 *funcSN = new TF1("funcSN",fitfSN,FitLB,FitUB,1);  funcSN->SetNpx(10*Nbins);  funcSN->SetParameter(0,1);  funcSN->SetParName(0,"Solar");

  //Fit functions (smeared)
  TF1 *funcMat_smear = new TF1("funcMat_smear",fitfMatSmear,FitLB,FitUB,1);  funcMat_smear->SetNpx(10*Nbins);  funcMat_smear->SetParameter(0,1);  funcMat_smear->SetParName(0,"sMaterials");
  TF1 *funcRn_smear = new TF1("funcRn_smear",fitfRnSmear,FitLB,FitUB,1);  funcRn_smear->SetNpx(10*Nbins);  funcRn_smear->SetParameter(0,1);  funcRn_smear->SetParName(0,"sRadon");
  TF1 *funcKr_smear = new TF1("funcKr_smear",fitfKrSmear,FitLB,FitUB,1);  funcKr_smear->SetNpx(10*Nbins);  funcKr_smear->SetParameter(0,1);  funcKr_smear->SetParName(0,"sKrypton");
  TF1 *funcDB_smear = new TF1("funcDB_smear",fitfDBSmear,FitLB,FitUB,1);  funcDB_smear->SetNpx(10*Nbins);  funcDB_smear->SetParameter(0,1);  funcDB_smear->SetParName(0,"sDouble Beta");
  TF1 *funcSN_smear = new TF1("funcSN_smear",fitfSNSmear,FitLB,FitUB,1);  funcSN_smear->SetNpx(10*Nbins);  funcSN_smear->SetParameter(0,1);  funcSN_smear->SetParName(0,"sSolar");

  //TH1F *rTot = new TH1F("rTot","",Nbins,hist_LB,hist_UB);  //rTot->Sumw2();

  TH1F *rTot_smear = new TH1F("rTot_smear","",Nbins,hist_LB,hist_UB); rTot_smear->SetLineColor(1);
  rTot_smear->Sumw2();
  double fillval=0;

  //Random sampling for Materials
  double N_Mat=0;
  for (int i=1; i<=Nbins+1; i++) N_Mat+=SMat->GetBinContent(i)*exposure;
  cout<<"N_Mat="<<N_Mat<<endl;
  TH1F *rMat = new TH1F("rMat","",Nbins,hist_LB,hist_UB);
  TH1F *rMat_smear = new TH1F("rMat_smear","",Nbins,hist_LB,hist_UB);
  rMat->Sumw2();  rMat_smear->Sumw2();
  rMat->SetLineColor(kMagenta);  rMat_smear->SetLineColor(kMagenta);
  for (int i=0; i<N_Mat; i++) {rMat->Fill(funcMat->GetRandom()); fillval=funcMat_smear->GetRandom(); rMat_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rMat->Scale(1./exposure); rMat_smear->Scale(1./exposure);

  //Random sampling for Rn
  double N_Rn=0;
  for (int i=1; i<=Nbins+1; i++) N_Rn+=SRn->GetBinContent(i)*exposure;
  cout<<"N_Rn="<<N_Rn<<endl;
  TH1F *rRn = new TH1F("rRn","",Nbins,hist_LB,hist_UB);
  TH1F *rRn_smear = new TH1F("rRn_smear","",Nbins,hist_LB,hist_UB);
  rRn->Sumw2();  rRn_smear->Sumw2();
  rRn->SetLineColor(kRed);  rRn_smear->SetLineColor(kRed);
  for (int i=0; i<N_Rn; i++) {rRn->Fill(funcRn->GetRandom()); fillval=funcRn_smear->GetRandom(); rRn_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rRn->Scale(1./exposure); rRn_smear->Scale(1./exposure);

  //Random sampling for Kr
  double N_Kr=0; 
  for (int i=1; i<=Nbins+1; i++) N_Kr+=SKr->GetBinContent(i)*exposure;
  cout<<"N_Kr="<<N_Kr<<endl;
  TH1F *rKr = new TH1F("rKr","",Nbins,hist_LB,hist_UB);
  TH1F *rKr_smear = new TH1F("rKr_smear","",Nbins,hist_LB,hist_UB);
  rKr->Sumw2();  rKr_smear->Sumw2();
  rKr->SetLineColor(kBlue);  rKr_smear->SetLineColor(kBlue);
  for (int i=0; i<N_Kr; i++) {rKr->Fill(funcKr->GetRandom()); fillval=funcKr_smear->GetRandom(); rKr_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rKr->Scale(1./exposure); rKr_smear->Scale(1./exposure);

  //Random sampling for DB
  double N_DB=0; 
  for (int i=1; i<=Nbins+1; i++) N_DB+=SDB->GetBinContent(i)*exposure;
  cout<<"N_DB="<<N_DB<<endl;
  TH1F *rDB = new TH1F("rDB","",Nbins,hist_LB,hist_UB);
  TH1F *rDB_smear = new TH1F("rDB_smear","",Nbins,hist_LB,hist_UB);
  rDB->Sumw2();  rDB_smear->Sumw2();
  rDB->SetLineColor(kOrange+2); rDB_smear->SetLineColor(kOrange+2);
  for (int i=0; i<N_DB; i++) {rDB->Fill(funcDB->GetRandom()); fillval=funcDB_smear->GetRandom(); rDB_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rDB->Scale(1./exposure); rDB_smear->Scale(1./exposure);

  //Random sampling for SN
  double N_SN=0; 
  for (int i=1; i<=Nbins+1; i++) N_SN+=SSN->GetBinContent(i)*exposure;
  cout<<"N_SN="<<N_SN<<endl;
  TH1F *rSN = new TH1F("rSN","",Nbins,hist_LB,hist_UB);
  TH1F *rSN_smear = new TH1F("rSN_smear","",Nbins,hist_LB,hist_UB);
  rSN->Sumw2();  rSN_smear->Sumw2();
  rSN->SetLineColor(kGreen);  rSN_smear->SetLineColor(kGreen);
  for (int i=0; i<N_SN; i++) {rSN->Fill(funcSN->GetRandom()); fillval=funcSN_smear->GetRandom(); rSN_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rSN->Scale(1./exposure); rSN_smear->Scale(1./exposure);

  rTot_smear->Scale(1./exposure); cout<<endl;

  //Fit unsmeared data
  rMat->Fit(funcMat,"NLQ","",FitLB,FitUB);  //funcMat->Draw("same");
  rRn->Fit(funcRn,"NLQ","",FitLB,FitUB);  //funcRn->Draw("same");
  rKr->Fit(funcKr,"NLQ","",FitLB,FitUB);  //funcKr->Draw("same");
  rDB->Fit(funcDB,"NLQ","",FitLB,FitUB);  //funcDB->Draw("same");
  rSN->Fit(funcSN,"NLQ","",FitLB,FitUB);  //funcSN->Draw("same");

  //Define canvas for sampled histograms
  TCanvas *s = new TCanvas("s", "s",265,50,1009,588);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  s->Range(-231.8436,-8.828753,2086.592,-0.5412262);
  s->SetFillColor(0);
  s->SetBorderMode(0);
  s->SetBorderSize(2);
  s->SetLogy();
  s->SetRightMargin(0.0373494);
  s->SetTopMargin(0.05535714);
  s->SetFrameBorderMode(0);

  rTot_smear->Draw("");  rMat_smear->Draw("same");  rRn_smear->Draw("same");  rKr_smear->Draw("same");  rDB_smear->Draw("same");  rSN_smear->Draw("same");
  //hTest->Draw("same");//SMat->Draw("same");  SRn->Draw("same");  SKr->Draw("same");  SDB->Draw("same");  SSN->Draw("same");

  leg->Draw("same");
  gPad->Modified();

  //Fit smeared data
  rMat_smear->Fit(funcMat_smear,"NLQ","",FitLB,FitUB);  //funcMat_smear->Draw("same");
  rRn_smear->Fit(funcRn_smear,"NLQ","",FitLB,FitUB);  //funcRn_smear->Draw("same");
  rKr_smear->Fit(funcKr_smear,"NLQ","",FitLB,FitUB);  //funcKr_smear->Draw("same");
  rDB_smear->Fit(funcDB_smear,"NLQ","",FitLB,FitUB);  //funcDB_smear->Draw("same");
  rSN_smear->Fit(funcSN_smear,"NLQ","",FitLB,FitUB);  //funcSN_smear->Draw("same");
//  return;

  //Run the campaign----------------------------------------------------------------------------------

  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////                        TOY EXPERIMENTS
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    //const int nToybins = Nbins;
//    TFile* file = new TFile("/home/atp/shayne/SolarNeutrinos_bak/toyMCs/MC_ppBe_40ty_Rn0p1_Kr2ppq_Mat1percent.root","read");
    TFile *file = new TFile(TString::Format("/home/directory/exposure_%ity.root",int(mult)),"read");

    TTree* tree = (TTree*)file->Get("tree");
    int *toyMC = new int[300];
    tree->SetBranchAddress("toyMC",toyMC);
    tree->GetEntry(index);
    

    double sin2_crit=0, surv_p_crit=0, maxLL=0;
    double d_sin2=1, d_surv=1;
    double sin2_prev=0, surv_prev=0;

    double mu_pp=1, mu_Be=1;
	
    cout<<"----------------------------------------------------------------------------"<<endl;
    double true_sin2_pp = 0.2387, true_surv_pp = 0.550;
    double true_sin2_Be = 0.2387, true_surv_Be = 0.520;
    sin2_crit=true_sin2_pp;
    surv_p_crit=true_surv_pp;
    std::vector<int> vtoyMC(toyMC,toyMC+300);
    maxLL=FindLL(vtoyMC, mult, mu_pp, true_sin2_pp, true_surv_pp, mu_Be, true_sin2_Be, true_surv_Be);
    cout<<"Initial Values:  "<<sin2_crit<<", "<<surv_p_crit<<", "<<maxLL<<endl;
    cout<<"----------------------------------------------------------------------------"<<endl;

    double N=100; //, N_B=200;
    double step_sin2=0.0001, step_surv=0.001;    
    double rejection_box_X_A=sin2_crit, rejection_box_X_B=sin2_crit, rejection_box_Y_A=surv_p_crit, rejection_box_Y_B=surv_p_crit;
    
    bool newCP=1, coarse=1, coarse1=0, allowed=1;
    int type=-1, possible_type=-1;
    double possible_sin2_crit=sin2_crit, possible_surv_p_crit=surv_p_crit;
    double LL = maxLL;
    const int n_slopes = 14;
    double dydx = step_surv/step_sin2;
    double slope[n_slopes] = {-dydx,-dydx/2.,-dydx/4.,-dydx/8.,-2.*dydx,-4.*dydx,-8*dydx, dydx,dydx/2.,dydx/4.,dydx/8.,2.*dydx,4.*dydx,8*dydx};
    step_sin2 *= N/1.; step_surv *= N/1.;
    N=2;

    while (newCP || coarse || coarse1) {
    
      newCP=0;
      type=-1, possible_type=-1;

      if (coarse) cout<<"Beginning coarse loop..."<<endl;
      if (coarse1) cout<<"Beginning finer loop..."<<endl;

      for (int m=0; m<N+1; m++) {  /*if (m%10==0)*/ //cout<<" m="<<m<<endl;
        if (m==1) continue;  /// this corresponds to the critical point for which we have already calculated the LL value; no need to waste time on redundant calculations
        cout<<" m="<<m<<endl;
        double x_sin2 = sin2_crit+(m-0.5*N)*step_sin2;
        double x_surv = surv_p_crit+(m-0.5*N)*step_surv;
        if (x_surv<0 || x_surv>1) continue;
        if (x_sin2<0 || x_sin2>1) continue;
        for (int i=0; i<n_slopes; i++) {
//          allowed=1;
          double v_slope = slope[i];
          double y = v_slope*(x_sin2-sin2_crit)+surv_p_crit;
//          if (x_surv<0 || x_surv>1) allowed=0;
          if (x_sin2<rejection_box_X_A || x_sin2>rejection_box_X_B || y<rejection_box_Y_A || y>rejection_box_Y_B || allowed) {
            LL = FindLL(vtoyMC, mult, mu_pp, x_sin2, y, mu_Be, true_sin2_Be, true_surv_Be);  	
            if (LL>maxLL) {
              possible_sin2_crit=x_sin2;
              possible_surv_p_crit=y;  
              maxLL=LL; 
              newCP=1;
              possible_type=i+2;   
            }
          }
        }
        if (x_sin2<rejection_box_X_A || x_sin2>rejection_box_X_B || surv_p_crit<rejection_box_Y_A || surv_p_crit>rejection_box_Y_B) {
          LL = FindLL(vtoyMC, mult, mu_pp, x_sin2, surv_p_crit, mu_Be, true_sin2_Be, true_surv_Be);  	
          if (LL>maxLL) {
            possible_sin2_crit=x_sin2;
            possible_surv_p_crit=surv_p_crit;
            maxLL=LL; 
            newCP=1;  
            possible_type=1; 
          }
        }
        if (x_surv<rejection_box_Y_A || x_surv>rejection_box_Y_B || sin2_crit<rejection_box_X_A || sin2_crit>rejection_box_X_B) {
          LL = FindLL(vtoyMC, mult, mu_pp, sin2_crit, x_surv, mu_Be, true_sin2_Be, true_surv_Be);  	
          if (LL>maxLL) {
            possible_sin2_crit=sin2_crit;
            possible_surv_p_crit=x_surv;
            maxLL=LL; 
            newCP=1; 
            possible_type=0;     
          }
        }

      } ///end loop over m
      rejection_box_X_A = sin2_crit-0.5*N*step_sin2;
      rejection_box_X_B = sin2_crit+0.5*N*step_sin2;
      rejection_box_Y_A = surv_p_crit-0.5*N*step_surv;
      rejection_box_Y_B = surv_p_crit+0.5*N*step_surv;
      if (coarse1 && !newCP) {
        rejection_box_X_A = true_sin2_pp; rejection_box_X_B = rejection_box_X_A;
        rejection_box_Y_A = true_surv_pp; rejection_box_Y_B = rejection_box_Y_A;
        N*=1;  step_sin2/=10.;  step_surv/=10.;
        coarse1=0; newCP=1; cout<<"   Repeating with level 2 finer iterations..."<<endl;
      } else if (coarse && !newCP) { //possible_sin2_crit>sin2_crit-0.5*N*step_sin2 && possible_sin2_crit<sin2_crit+0.5*N*step_sin2 && possible_surv_p_crit>surv_p_crit-0.5*N*step_surv && possible_surv_p_crit<surv_p_crit+0.5*N*step_surv) {
          rejection_box_X_A = true_sin2_pp; rejection_box_X_B = rejection_box_X_A;
          rejection_box_Y_A = true_surv_pp; rejection_box_Y_B = rejection_box_Y_A;
          N*=1;  step_sin2/=10.;  step_surv/=10.;
          coarse=0; coarse1=1; newCP=1; cout<<"   Repeating with finer iterations..."<<endl;
        }
      if(newCP) {
        sin2_prev=sin2_crit;
        surv_prev=surv_p_crit;
        sin2_crit=possible_sin2_crit;
        surv_p_crit=possible_surv_p_crit;
        type=possible_type;
        d_sin2 = sin2_crit-sin2_prev;
        d_surv = surv_p_crit-surv_prev;
        cout<<endl<<"The new setpoint is at [ "<<sin2_crit<<" , "<<surv_p_crit<<" ]  ("<<maxLL<<")    ["<<type<<"]"<<endl<<endl;
        cout<<"Repeating maximization..."<<endl;
        cout<<"  ...with rejection box  ["<<rejection_box_X_A<<","<<rejection_box_X_B<<"]  ["<<rejection_box_Y_A<<","<<rejection_box_Y_B<<"]"<<endl<<endl;

      }
    }

    cout<<endl<<"The maximum LL is at [ "<<sin2_crit<<" , "<<surv_p_crit<<" ]  ("<<maxLL<<")    ["<<type<<"]"<<endl<<endl;

    cout<<endl<<"Finally, we settle on a maximum LL at [ "<<sin2_crit<<" , "<<surv_p_crit<<" ]"<<endl<<endl;

    ofstream myfile(TString::Format("/home/directory/exposure_%ity/filename.txt",int(mult),int(mult),index));
    myfile << sin2_crit << endl;
    myfile << surv_p_crit << endl;
    myfile.close();

  return;
}
