///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///  Shayne Reichard, 2020 May 26, shayne@physik.uzh.ch
///
///  This script simulates the detection of solar neutrinos in liquid xenon through the electron scattering process. You set the background levels
///  and then run the script for a given exposure ("mult" in tonne-years).
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
#include "TROOT.h"
#include "TVirtualFitter.h"



const int Nbins=300;
const int hist_LB=0, hist_UB=3000;
const int nBands = 7;
const double tonne=1000.;  //kg
const double year=365.25;  //days

const double factor=100.;  //unity for XENON1T, 10 for XENONnT, 100 for DARWIN

//const double sigma_Rn=0.1, sigma_Kr=0.2, sigma_DB=0.15;
const double sigma_Rn=0, sigma_Kr=0, sigma_DB=0;

const double multiplier_Mat=1./factor, multiplier_Rn=10./factor*(1+sigma_Rn), multiplier_Kr=1./factor*(1+sigma_Kr), multiplier_DB=1.;// /factor*(1+sigma_DB);//*1e-20;
const double multiplier_SN=1., multiplier_Be7=1., multiplier_Be7_other=multiplier_Be7, multiplier_pep=1., multiplier_N13=1., multiplier_O15=1.;
//const double multiplier_SN=6.03/5.98, multiplier_Be7=4.50/4.93, multiplier_Be7_other=multiplier_Be7, multiplier_pep=1.46/1.44, multiplier_N13=2.04/2.78, multiplier_O15=1.44/2.05;



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
TFile *fNC = new TFile("./inputs/NC_abundance_212.root"); //NC abundance Xe131 21.2%
TFile *fPatricia = new TFile("./inputs/shayne_materials.root");

TH1F *SMat = (TH1F*)fPatricia->Get("sMat");
TH1F *SRn = (TH1F*)fSmear->Get("sRn");
TH1F *SKr = (TH1F*)fSmear->Get("sKr");
TH1F *SDB = (TH1F*)fSmear->Get("sDB");
TH1F *SSN = (TH1F*)fSmear->Get("sSN");
TH1F *SBe7 = (TH1F*)fSmear->Get("sBe7");
TH1F *SBe7_other = (TH1F*)fSmear->Get("sBe7_other");
TH1F *Spep = (TH1F*)fSmear->Get("spep");
TH1F *SN13 = (TH1F*)fSmear->Get("sN13");
TH1F *SO15 = (TH1F*)fSmear->Get("sO15");
TH1F *SNC = (TH1F*)fNC->Get("h_NC");

//Define fit functions with unsmeared histograms
Double_t fitfMat(Double_t *v, Double_t *par) { Double_t fitval = par[0]*hMat->Interpolate(v[0]);  return fitval; }
Double_t fitfRn(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hRn->Interpolate(v[0]);  return fitval; }
Double_t fitfKr(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hKr->Interpolate(v[0]);  return fitval; }
Double_t fitfDB(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hDB->Interpolate(v[0]);  return fitval; }
Double_t fitfSN(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hSN->Interpolate(v[0]);  return fitval; }
Double_t fitfBe7(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*(hBe7->Interpolate(v[0])+hBe7_other->Interpolate(v[0]));  return fitval; }
Double_t fitfpep(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hpep->Interpolate(v[0]);  return fitval; }
Double_t fitfN13(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hN13->Interpolate(v[0]);  return fitval; }
Double_t fitfO15(Double_t *v, Double_t *par) {  Double_t fitval = par[0]*hO15->Interpolate(v[0]);  return fitval; }

//Define fit functions with smeared histograms.
Double_t fitfMatSmear(Double_t *v, Double_t *par) { Double_t fitval=par[0]*SMat->Interpolate(v[0]);  return fitval; }
Double_t fitfRnSmear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SRn->Interpolate(v[0]);  return fitval; }
Double_t fitfKrSmear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SKr->Interpolate(v[0]);  return fitval; }
Double_t fitfDBSmear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SDB->Interpolate(v[0]);  return fitval; }
Double_t fitfSNSmear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SSN->Interpolate(v[0]);  return fitval; }
Double_t fitfBe7Smear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*(SBe7->Interpolate(v[0])+SBe7_other->Interpolate(v[0]));  return fitval; }
Double_t fitfpepSmear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*Spep->Interpolate(v[0]);  return fitval; }
Double_t fitfN13Smear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SN13->Interpolate(v[0]);  return fitval; }
Double_t fitfO15Smear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SO15->Interpolate(v[0]);  return fitval; }
Double_t fitfNCSmear(Double_t *v, Double_t *par) {  Double_t fitval=par[0]*SNC->Interpolate(v[0]);  return fitval; }


Double_t fitfTotSmear(Double_t *v, Double_t *par) {
  Double_t fitval = 
    par[0]*SMat->Interpolate(v[0])
    +par[1]*SRn->Interpolate(v[0])
    +par[2]*SKr->Interpolate(v[0])
    +par[3]*SDB->Interpolate(v[0])
    +par[4]*SSN->Interpolate(v[0])
//    +par[5]*SBe7->Interpolate(v[0])
    +par[5]*(SBe7->Interpolate(v[0])+SBe7_other->Interpolate(v[0]))
    +par[6]*Spep->Interpolate(v[0])
    +par[7]*SN13->Interpolate(v[0])
    +par[8]*SO15->Interpolate(v[0])
    +par[9]*SNC->Interpolate(v[0]);
  return fitval;
}


std::string double_to_string(double inputnumber) {
  string s;
  stringstream out;
  out << inputnumber;
  s = out.str();
  return s;
}


Double_t EnergyRes(double E) {
  double a=0.313;  double b=0.002;
  return a*TMath::Power(E,0.5)+b*E;
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


double GetStats(TH1D* hist, int flag) {
  double prob[nBands] = {2.87e-7, 0.00135, 0.15866, 0.5, 0.84134, 0.99865, 0.999999713};
  double quantiles[nBands];
  hist->GetQuantiles(nBands, quantiles, prob);
  return quantiles[flag];
}



void GetPercentError(TH1D* hist) {
  double median = GetStats(hist,3);
  double lowerr = GetStats(hist,2), low3sig = GetStats(hist,1), low5sig = GetStats(hist,0);
  double upperr = GetStats(hist,4), upp3sig = GetStats(hist,5), upp5sig = GetStats(hist,6);
  double error = 0.5*(upperr-lowerr);
  cout<<"         "<<low5sig<<endl;
  cout<<"       "<<low3sig<<endl;
  cout<<"     "<<lowerr<<endl;
  cout<<"   "<<median<<endl;
  cout<<"     "<<upperr<<endl;
  cout<<"       "<<upp3sig<<endl;
  cout<<"         "<<upp5sig<<endl;

  cout<<"err1 = "<<error<<endl;
  cout<<"err3 = "<<0.5*(upp3sig-low3sig)/3.<<endl;
  //return error;
}


//------------------------------------------------------------------------------------------------------------------------



void GenerateNuData(double mult=1.0) {

  hMat->Scale(multiplier_Mat); hRn->Scale(multiplier_Rn); hKr->Scale(multiplier_Kr); hDB->Scale(multiplier_DB);
  SMat->Scale(multiplier_Mat);
  SRn->Scale(multiplier_Rn); SKr->Scale(multiplier_Kr); SDB->Scale(multiplier_DB);
  hSN->Scale(multiplier_SN); hBe7->Scale(multiplier_Be7); hBe7_other->Scale(multiplier_Be7_other); hpep->Scale(multiplier_pep); hN13->Scale(multiplier_N13); hO15->Scale(multiplier_O15);
  SSN->Scale(multiplier_SN); SBe7->Scale(multiplier_Be7); SBe7_other->Scale(multiplier_Be7_other); Spep->Scale(multiplier_pep); SN13->Scale(multiplier_N13); SO15->Scale(multiplier_O15);

  double exposure=mult*tonne*year; cout<<"exposure="<<exposure<<endl;
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
  hTot->GetXaxis()->SetRange(1,Nbins); hTot->SetMinimum(1.e-9); hTot->SetMaximum(1e0); hTot->SetLineColor(0);
  hTot->GetXaxis()->SetTitle("Energy [keV]"); hTot->GetYaxis()->SetTitle("count/kg/day/10keV");
  hTot->GetXaxis()->CenterTitle(); hTot->GetYaxis()->CenterTitle();

  hSN->SetLineColor(kGreen); hSN->SetLineWidth(2);  hDB->SetLineColor(kOrange+1);
  hBe7->SetLineColor(kCyan); hBe7->SetLineWidth(2);
  hBe7_other->SetLineColor(kCyan); hBe7_other->SetLineWidth(2);
  hpep->SetLineColor(kYellow+1); hpep->SetLineWidth(2);
  hN13->SetLineColor(kViolet); hN13->SetLineWidth(2);
  hO15->SetLineColor(kGray+1); hO15->SetLineWidth(2);

  hTot->Draw("");  hMat->Draw("same");  hRn->Draw("same");  hKr->Draw("same");  hDB->Draw("same");
  hSN->Draw("same"); hBe7->Draw("same"); hBe7_other->Draw("same"); hpep->Draw("same"); hN13->Draw("same"); hO15->Draw("same");

  TLegend *leg = new TLegend(0.75,0.60,0.90,0.95,NULL,"brNDC");
  leg->SetFillColor(0);
  leg->AddEntry("hEdTotal","Materials","l");
  leg->AddEntry("hEd_Rn222","Radon","l");
  leg->AddEntry("hEd_Kr85","Krypton","l");
  leg->AddEntry("h1_Xe_1","Double Beta","l");
  leg->AddEntry("h_pp","pp","l");
  leg->AddEntry("h_Be7","Be7","l");
  leg->AddEntry("h_pep","pep","l");
  leg->AddEntry("h_N","^{13}N","l");
  leg->AddEntry("h_O","^{15}O","l");

  double top_margin=0.02;
  double bottom_margin=0.11;
  double right_margin=0.03;
  double left_margin=0.10;
  double title_size=0.06;
  double label_size=0.06;
  double x_title_offset=0.89;
  double y_title_offset=0.86;

  gPad->SetTopMargin(top_margin);
  gPad->SetBottomMargin(bottom_margin);
  gPad->SetRightMargin(right_margin);
  gPad->SetLeftMargin(left_margin);

  hTot->GetXaxis()->SetTitleSize(title_size);
  hTot->GetYaxis()->SetTitleSize(title_size);
  hTot->GetXaxis()->SetLabelSize(title_size);
  hTot->GetYaxis()->SetLabelSize(title_size);
  hTot->GetXaxis()->SetTitleOffset(x_title_offset);
  hTot->GetYaxis()->SetTitleOffset(y_title_offset);


  leg->Draw("same");
  gPad->Modified();

  //------------------------------------------------------------------------------------

  //Fit functions (unsmeared)
  TF1 *funcMat = new TF1("funcMat",fitfMat,FitLB,FitUB,1);  funcMat->SetNpx(10*Nbins);  funcMat->SetParameter(0,1);  funcMat->SetParName(0,"Materials");
  TF1 *funcRn = new TF1("funcRn",fitfRn,FitLB,FitUB,1);  funcRn->SetNpx(10*Nbins);  funcRn->SetParameter(0,1);  funcRn->SetParName(0,"Radon");
  TF1 *funcKr = new TF1("funcKr",fitfKr,FitLB,FitUB,1);  funcKr->SetNpx(10*Nbins);  funcKr->SetParameter(0,1);  funcKr->SetParName(0,"Krypton");
  TF1 *funcDB = new TF1("funcDB",fitfDB,FitLB,FitUB,1);  funcDB->SetNpx(10*Nbins);  funcDB->SetParameter(0,1);  funcDB->SetParName(0,"Double Beta");
  TF1 *funcSN = new TF1("funcSN",fitfSN,FitLB,FitUB,1);  funcSN->SetNpx(10*Nbins);  funcSN->SetParameter(0,1);  funcSN->SetParName(0,"pp");
  TF1 *funcBe7 = new TF1("funcBe7",fitfBe7,FitLB,FitUB,1);  funcBe7->SetNpx(10*Nbins);  funcBe7->SetParameter(0,1);  funcBe7->SetParName(0,"Be7");
  TF1 *funcpep = new TF1("funcpep",fitfpep,FitLB,FitUB,1);  funcpep->SetNpx(10*Nbins);  funcpep->SetParameter(0,1);  funcpep->SetParName(0,"pep");
  TF1 *funcN13 = new TF1("funcN13",fitfN13,FitLB,FitUB,1);  funcN13->SetNpx(10*Nbins);  funcN13->SetParameter(0,1);  funcN13->SetParName(0,"N13");
  TF1 *funcO15 = new TF1("funcO15",fitfO15,FitLB,FitUB,1);  funcO15->SetNpx(10*Nbins);  funcO15->SetParameter(0,1);  funcO15->SetParName(0,"O15");

  //Fit functions (smeared)
  TF1 *funcMat_smear = new TF1("funcMat_smear",fitfMatSmear,FitLB,FitUB,1);  funcMat_smear->SetNpx(1000*Nbins);  funcMat_smear->SetParameter(0,1);  funcMat_smear->SetParName(0,"sMaterials");
  TF1 *funcRn_smear = new TF1("funcRn_smear",fitfRnSmear,FitLB,FitUB,1);  funcRn_smear->SetNpx(1000*Nbins);  funcRn_smear->SetParameter(0,1);  funcRn_smear->SetParName(0,"sRadon");
  TF1 *funcKr_smear = new TF1("funcKr_smear",fitfKrSmear,FitLB,FitUB,1);  funcKr_smear->SetNpx(1000*Nbins);  funcKr_smear->SetParameter(0,1);  funcKr_smear->SetParName(0,"sKrypton");
  TF1 *funcDB_smear = new TF1("funcDB_smear",fitfDBSmear,FitLB,FitUB,1);  funcDB_smear->SetNpx(1000*Nbins);  funcDB_smear->SetParameter(0,1);  funcDB_smear->SetParName(0,"sDouble Beta");
  TF1 *funcSN_smear = new TF1("funcSN_smear",fitfSNSmear,FitLB,FitUB,1);  funcSN_smear->SetNpx(1000*Nbins);  funcSN_smear->SetParameter(0,1);  funcSN_smear->SetParName(0,"spp");
  TF1 *funcBe7_smear = new TF1("funcBe7_smear",fitfBe7Smear,FitLB,FitUB,1);  funcBe7_smear->SetNpx(1000*Nbins);  funcBe7_smear->SetParameter(0,1);  funcBe7_smear->SetParName(0,"sBe7");
  TF1 *funcpep_smear = new TF1("funcpep_smear",fitfpepSmear,FitLB,FitUB,1);  funcpep_smear->SetNpx(1000*Nbins);  funcpep_smear->SetParameter(0,1);  funcpep_smear->SetParName(0,"spep");
  TF1 *funcN13_smear = new TF1("funcN13_smear",fitfN13Smear,FitLB,FitUB,1);  funcN13_smear->SetNpx(1000*Nbins);  funcN13_smear->SetParameter(0,1);  funcN13_smear->SetParName(0,"sN13");
  TF1 *funcO15_smear = new TF1("funcO15_smear",fitfO15Smear,FitLB,FitUB,1);  funcO15_smear->SetNpx(1000*Nbins);  funcO15_smear->SetParameter(0,1);  funcO15_smear->SetParName(0,"sO15");
  TF1 *funcNC_smear = new TF1("funcNC_smear",fitfNCSmear,FitLB,FitUB,1);  funcNC_smear->SetNpx(1000*Nbins);  funcNC_smear->SetParameter(0,1);  funcNC_smear->SetParName(0,"sNC");


  TH1F *rTot_smear = new TH1F("rTot_smear","",Nbins,hist_LB,hist_UB); rTot_smear->SetLineColor(1);
  rTot_smear->Sumw2();
  double fillval=0;

  //Random sampling for Materials
  double N_Mat=0, N_fnc_Mat=0;
  for (int i=1; i<=Nbins+1; i++) N_Mat+=SMat->GetBinContent(i)*exposure;
  cout<<"N_Mat="<<N_Mat<<endl;
  N_fnc_Mat=funcMat_smear->Integral(0,3000)*exposure/10.;
  cout<<"N_fnc_Mat="<<N_fnc_Mat<<endl;
  N_Mat=N_fnc_Mat;

  TH1F *rMat = new TH1F("rMat","",Nbins,hist_LB,hist_UB);
  TH1F *rMat_smear = new TH1F("rMat_smear","",Nbins,hist_LB,hist_UB);
  rMat->Sumw2();  rMat_smear->Sumw2();
  rMat->SetLineColor(kMagenta);  rMat_smear->SetLineColor(kMagenta); rMat_smear->SetLineWidth(2);
  for (int i=0; i<N_Mat; i++) {rMat->Fill(funcMat->GetRandom()); fillval=funcMat_smear->GetRandom(); rMat_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rMat->Scale(1./exposure); rMat_smear->Scale(1./exposure);

  //Random sampling for Rn
  double N_Rn=0, N_fnc_Rn=0;;
  for (int i=1; i<=Nbins+1; i++) N_Rn+=SRn->GetBinContent(i)*exposure;
  cout<<"N_Rn="<<N_Rn<<endl;
  N_fnc_Rn=funcRn_smear->Integral(0,3000)*exposure/10.;
  cout<<"N_fnc_Rn="<<N_fnc_Rn<<endl;
  N_Rn=N_fnc_Rn;

  TH1F *rRn = new TH1F("rRn","",Nbins,hist_LB,hist_UB);
  TH1F *rRn_smear = new TH1F("rRn_smear","",Nbins,hist_LB,hist_UB);
  rRn->Sumw2();  rRn_smear->Sumw2();
  rRn->SetLineColor(kRed);  rRn_smear->SetLineColor(kRed); rRn_smear->SetLineWidth(2);
  for (int i=0; i<N_Rn; i++) {rRn->Fill(funcRn->GetRandom()); fillval=funcRn_smear->GetRandom(); rRn_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rRn->Scale(1./exposure); rRn_smear->Scale(1./exposure);

  //Random sampling for Kr
  double N_Kr=0, N_fnc_Kr=0; 
  for (int i=1; i<=Nbins+1; i++) N_Kr+=SKr->GetBinContent(i)*exposure;
  cout<<"N_Kr="<<N_Kr<<endl;
  N_fnc_Kr=funcKr_smear->Integral(0,3000)*exposure/10.;
  cout<<"N_fnc_Kr="<<N_fnc_Kr<<endl;
  N_Kr=N_fnc_Kr;

  TH1F *rKr = new TH1F("rKr","",Nbins,hist_LB,hist_UB);
  TH1F *rKr_smear = new TH1F("rKr_smear","",Nbins,hist_LB,hist_UB);
  rKr->Sumw2();  rKr_smear->Sumw2();
  rKr->SetLineColor(kBlue);  rKr_smear->SetLineColor(kBlue); rKr_smear->SetLineWidth(2);
  for (int i=0; i<N_Kr; i++) {rKr->Fill(funcKr->GetRandom()); fillval=funcKr_smear->GetRandom(); rKr_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rKr->Scale(1./exposure); rKr_smear->Scale(1./exposure);

  //Random sampling for DB
  double N_DB=0, N_fnc_DB=0; 
  for (int i=1; i<=Nbins+1; i++) N_DB+=SDB->GetBinContent(i)*exposure;
  cout<<"N_DB="<<N_DB<<endl;
  N_fnc_DB=funcDB_smear->Integral(0,3000)*exposure/10.;
  cout<<"N_fnc_DB="<<N_fnc_DB<<endl;
  N_DB=N_fnc_DB;

  TH1F *rDB = new TH1F("rDB","",Nbins,hist_LB,hist_UB);
  TH1F *rDB_smear = new TH1F("rDB_smear","",Nbins,hist_LB,hist_UB);
  rDB->Sumw2();  rDB_smear->Sumw2();
  rDB->SetLineColor(kOrange+2); rDB_smear->SetLineColor(kOrange+2); rDB_smear->SetLineWidth(2);
  for (int i=0; i<N_DB; i++) {rDB->Fill(funcDB->GetRandom()); fillval=funcDB_smear->GetRandom(); rDB_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rDB->Scale(1./exposure); rDB_smear->Scale(1./exposure);

  //Random sampling for SN
  double N_SN=0, N_fnc_SN=0; 
  for (int i=1; i<=Nbins+1; i++) N_SN+=SSN->GetBinContent(i)*exposure;
  cout<<"N_SN="<<N_SN<<endl;
  N_fnc_SN=funcSN_smear->Integral(0,3000)*exposure/10.;
  cout<<"N_fnc_SN="<<N_fnc_SN<<endl;
  N_SN=N_fnc_SN;

  TH1F *rSN = new TH1F("rSN","",Nbins,hist_LB,hist_UB);
  TH1F *rSN_smear = new TH1F("rSN_smear","",Nbins,hist_LB,hist_UB);
  rSN->Sumw2();  rSN_smear->Sumw2();
  rSN->SetLineColor(kGreen);  rSN_smear->SetLineColor(kGreen); rSN_smear->SetLineWidth(2);
  for (int i=0; i<N_SN; i++) {rSN->Fill(funcSN->GetRandom()); fillval=funcSN_smear->GetRandom(); rSN_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rSN->Scale(1./exposure); rSN_smear->Scale(1./exposure);

  //Random sampling for Be7
  double N_Be7=0, N_fnc_Be7=0;;
  for (int i=1; i<=Nbins+1; i++) N_Be7+=(SBe7->GetBinContent(i)+SBe7_other->GetBinContent(i))*exposure;
  cout<<"N_Be7="<<N_Be7<<endl;
  N_fnc_Be7=funcBe7_smear->Integral(0,3000)*exposure/10.;
  cout<<"N_fnc_Be7="<<N_fnc_Be7<<endl;
  N_Be7=N_fnc_Be7;

  TH1F *rBe7 = new TH1F("rBe7","",Nbins,hist_LB,hist_UB);
  TH1F *rBe7_smear = new TH1F("rBe7_smear","",Nbins,hist_LB,hist_UB);
  rBe7->Sumw2();  rBe7_smear->Sumw2();
  rBe7->SetLineColor(kCyan);  rBe7_smear->SetLineColor(kCyan); rBe7_smear->SetLineWidth(2);
  for (int i=0; i<N_Be7; i++) {rBe7->Fill(funcBe7->GetRandom()); fillval=funcBe7_smear->GetRandom(); rBe7_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rBe7->Scale(1./exposure); rBe7_smear->Scale(1./exposure);

  //Random sampling for pep
  double N_pep=0, N_fnc_pep=0;
  for (int i=1; i<=Nbins+1; i++) N_pep+=Spep->GetBinContent(i)*exposure;
  cout<<"N_pep="<<N_pep<<endl;
  N_fnc_pep=funcpep_smear->Integral(0,3000)*exposure/10.;
  cout<<"N_fnc_pep="<<N_fnc_pep<<endl;
  N_pep=N_fnc_pep;

  TH1F *rpep = new TH1F("rpep","",Nbins,hist_LB,hist_UB);
  TH1F *rpep_smear = new TH1F("rpep_smear","",Nbins,hist_LB,hist_UB);
  rpep->Sumw2();  rpep_smear->Sumw2();
  rpep->SetLineColor(kYellow+1);  rpep_smear->SetLineColor(kYellow+1); rpep_smear->SetLineWidth(2);
  for (int i=0; i<N_pep; i++) {rpep->Fill(funcpep->GetRandom()); fillval=funcpep_smear->GetRandom(); rpep_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rpep->Scale(1./exposure); rpep_smear->Scale(1./exposure);

  //Random sampling for N13
  double N_N13=0, N_fnc_N13=0;
  for (int i=1; i<=Nbins+1; i++) N_N13+=SN13->GetBinContent(i)*exposure;
  cout<<"N_N13="<<N_N13<<endl;
  N_fnc_N13=funcN13_smear->Integral(0,3000)*exposure/10.;
  cout<<"N_fnc_N13="<<N_fnc_N13<<endl;
  N_N13=N_fnc_N13;

  TH1F *rN13 = new TH1F("rN13","",Nbins,hist_LB,hist_UB);
  TH1F *rN13_smear = new TH1F("rN13_smear","",Nbins,hist_LB,hist_UB);
  rN13->Sumw2();  rN13_smear->Sumw2();
  rN13->SetLineColor(kViolet);  rN13_smear->SetLineColor(kViolet); rN13_smear->SetLineWidth(2);
  for (int i=0; i<N_N13; i++) {rN13->Fill(funcN13->GetRandom()); fillval=funcN13_smear->GetRandom(); rN13_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rN13->Scale(1./exposure); rN13_smear->Scale(1./exposure);

  //Random sampling for O15
  double N_O15=0, N_fnc_O15=0;
  for (int i=1; i<=Nbins+1; i++) N_O15+=SO15->GetBinContent(i)*exposure;
  cout<<"N_O15="<<N_O15<<endl;
  N_fnc_O15=funcO15_smear->Integral(0,3000)*exposure/10.;
  cout<<"N_fnc_O15="<<N_fnc_O15<<endl;
  N_O15=N_fnc_O15;

  TH1F *rO15 = new TH1F("rO15","",Nbins,hist_LB,hist_UB);
  TH1F *rO15_smear = new TH1F("rO15_smear","",Nbins,hist_LB,hist_UB);
  rO15->Sumw2();  rO15_smear->Sumw2();
  rO15->SetLineColor(kGray+1);  rO15_smear->SetLineColor(kGray+1); rO15_smear->SetLineWidth(2);
  for (int i=0; i<N_O15; i++) {rO15->Fill(funcO15->GetRandom()); fillval=funcO15_smear->GetRandom(); rO15_smear->Fill(fillval); rTot_smear->Fill(fillval);}
  rO15->Scale(1./exposure); rO15_smear->Scale(1./exposure);

  //Random sampling for NC
  double N_NC=0, N_fnc_NC=0;
  for (int i=1; i<Nbins+1; i++) N_NC+=SNC->GetBinContent(i)*exposure;
  cout<<"N_NC="<<N_NC<<endl;

  N_fnc_NC=funcNC_smear->Integral(0,3000)*exposure/10.;
  cout<<"N_fnc_NC="<<N_fnc_NC<<endl;
  N_NC=N_fnc_NC;

  TH1F *rNC = new TH1F("rNC","",Nbins,hist_LB,hist_UB);
  TH1F *rNC_smear = new TH1F("rNC_smear","",Nbins,hist_LB,hist_UB);
  rNC->Sumw2();  rNC_smear->Sumw2();
  rNC_smear->SetLineColor(kCyan);
  rNC_smear->SetLineWidth(2);
  for (int i=0; i<N_NC; i++) {
    fillval=funcNC_smear->GetRandom(); rNC_smear->Fill(fillval); rTot_smear->Fill(fillval);
  }
  rNC_smear->Scale(1./exposure);



  //Run the campaign----------------------------------------------------------------------------------

  double fitvalueMat=0, fitvalueRn=0, fitvalueKr=0, fitvalueDB=0, fitvalueSN=0, fitvalueBe7=0, fitvaluepep=0, fitvalueN13=0, fitvalueO15=0, fitvalueNC=0; // fitvalueF17=0;
  double LB=-1, UB=3;
  int nbins=(UB-LB)*10000;
  TH1D *hmatA = new TH1D("hmatA",";Materials",nbins,0.98,1.02); hmatA->Sumw2();
  TH1D *hrnA = new TH1D("hrnA",";Radon",nbins,LB,UB); hrnA->Sumw2();
  TH1D *hkrA = new TH1D("hkrA",";Krypton",nbins,LB,UB); hkrA->Sumw2();
  TH1D *hdbA = new TH1D("hdbA",";Double Beta",nbins,0.99,1.01); hdbA->Sumw2();
  TH1D *hnuA = new TH1D("hnuA",";pp",nbins,LB,UB); hnuA->Sumw2();
  TH1D *hBe7A = new TH1D("hBe7A",";Be7",nbins,LB,UB); hBe7A->Sumw2();
  TH1D *hpepA = new TH1D("hpepA",";pep",nbins,LB,UB); hpepA->Sumw2();
  TH1D *hN13A = new TH1D("hN13A",";N13",nbins,LB,UB); hN13A->Sumw2();
  TH1D *hO15A = new TH1D("hO15A",";O15",nbins,LB,UB); hO15A->Sumw2();
  TH1D *hNC = new TH1D("hNC",";nuCap",nbins,LB,UB); hNC->Sumw2();
  TH1D *hsumA = new TH1D("hsumA",";Sum",nbins,0,18); hsumA->Sumw2();


  /////Define the function to be used in fitting during each iteration of the Toy Experiments
  TF1 *funcMC = new TF1("funcMC",fitfTotSmear,FitLB,FitUB,10);  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  funcMC->SetNpx(1000*Nbins);
  funcMC->SetParName(0,"Materials");  funcMC->SetParName(1,"Radon");  funcMC->SetParName(2,"Krypton");  funcMC->SetParName(3,"Double Beta");
  funcMC->SetParName(4,"pp"); funcMC->SetParName(5,"Be7"); funcMC->SetParName(6,"pep"); funcMC->SetParName(7,"N13"); funcMC->SetParName(8,"O15");
  funcMC->SetParName(9,"NC");

  funcMC->SetParameter(0,1.);  funcMC->SetParameter(1,1.);  funcMC->SetParameter(2,1.);  funcMC->SetParameter(3,1);
  funcMC->SetParameter(4,1.);  funcMC->SetParameter(5,1);  funcMC->SetParameter(6,1); funcMC->SetParameter(7,1);  funcMC->SetParameter(8,1);
  funcMC->SetParameter(9,1);

//  return;

  TH1D *dummyMC = new TH1D("dummyMC","",Nbins,hist_LB,hist_UB); dummyMC->Sumw2();
  TH1D *dummyMCmat = new TH1D("dummyMCmat","",Nbins,hist_LB,hist_UB); dummyMCmat->Sumw2();
  TH1D *dummyMCrn = new TH1D("dummyMCrn","",Nbins,hist_LB,hist_UB); dummyMCrn->Sumw2();
  TH1D *dummyMCkr = new TH1D("dummyMCkr","",Nbins,hist_LB,hist_UB); dummyMCkr->Sumw2();
  TH1D *dummyMCdb = new TH1D("dummyMCdb","",Nbins,hist_LB,hist_UB); dummyMCdb->Sumw2();
  TH1D *dummyMCnu = new TH1D("dummyMCnu","",Nbins,hist_LB,hist_UB); dummyMCnu->Sumw2();
  TH1D *dummyMCBe7 = new TH1D("dummyMCBe7","",Nbins,hist_LB,hist_UB); dummyMCBe7->Sumw2();
  TH1D *dummyMCpep = new TH1D("dummyMCpep","",Nbins,hist_LB,hist_UB); dummyMCpep->Sumw2();
  TH1D *dummyMCN13 = new TH1D("dummyMCN13","",Nbins,hist_LB,hist_UB); dummyMCN13->Sumw2();
  TH1D *dummyMCO15 = new TH1D("dummyMCO15","",Nbins,hist_LB,hist_UB); dummyMCO15->Sumw2();
  TH1D *dummyMCNC = new TH1D("dummyMCNC","",Nbins,hist_LB,hist_UB); dummyMCNC->Sumw2();

  bool alreadydone=false;

  double r_N_Mat=0, r_N_Rn=0, r_N_Kr=0, r_N_DB=0, r_N_SN=0, r_N_Be7=0, r_N_pep=0, r_N_N13=0, r_N_O15=0, r_N_NC=0;
  double g_N_SN=0, g_N_Be7=0, g_N_pep=0, g_N_N13=0, g_N_O15=0, g_N_NC=0;
  double nsig = 0;

  int nombre=1e5;
  for (int l=0; l<nombre; l++) { if (l%1000==0) cout<<endl<<"l="<<l<<endl<<endl;
    if (l%1000==0 && l>0) {

      cout<<"MATERIALS"<<endl;
      GetPercentError(hmatA);
      cout<<"RADON"<<endl;
      GetPercentError(hrnA);
      cout<<"KRYPTON"<<endl;
      GetPercentError(hkrA);
      cout<<"DOUBLE BETA"<<endl;
      GetPercentError(hdbA);

      cout<<"pp"<<endl;
      GetPercentError(hnuA);
      cout<<"Be7"<<endl;
      GetPercentError(hBe7A);
      cout<<"N13"<<endl;
      GetPercentError(hN13A);
      cout<<"O15"<<endl;
      GetPercentError(hO15A);
      cout<<"pep"<<endl;
      GetPercentError(hpepA);
      cout<<"NC"<<endl;
      GetPercentError(hNC);
    }
    //reset the histogram
    dummyMC->Reset(); dummyMCmat->Reset(); dummyMCrn->Reset(); dummyMCkr->Reset(); dummyMCdb->Reset(); dummyMCnu->Reset();
    dummyMCBe7->Reset(); dummyMCpep->Reset(); dummyMCN13->Reset(); dummyMCO15->Reset(); dummyMCNC->Reset();


    //account for HZ model uncertainties
//    g_N_SN=gRandom->Gaus(0,0.006*N_SN); g_N_Be7=gRandom->Gaus(0,0.06*N_Be7); g_N_pep=gRandom->Gaus(0,0.01*N_pep); g_N_N13=gRandom->Gaus(0.15*N_N13); g_N_O15=gRandom->Gaus(0,0.17*N_O15);
//    g_N_SN=0; g_N_Be7=0; g_N_pep=0; g_N_N13=0; g_N_O15=0; g_N_NC=0;

    //account for LZ model uncertainties
//    g_N_SN=gRandom->Gaus(0,0.005*N_SN); g_N_Be7=gRandom->Gaus(0,0.06*N_Be7); g_N_pep=gRandom->Gaus(0,0.01*N_pep); g_N_N13=gRandom->Gaus(0.14*N_N13); g_N_O15=gRandom->Gaus(0,0.16*N_O15);


    //randomly sample the total count of each ER source
    r_N_Mat=gRandom->Poisson(N_Mat); r_N_Rn=gRandom->Poisson(N_Rn); r_N_Kr=gRandom->Poisson(N_Kr); r_N_DB=gRandom->Poisson(N_DB);
//    r_N_SN=gRandom->Poisson(N_SN); r_N_Be7=gRandom->Poisson(N_Be7); r_N_pep=gRandom->Poisson(N_pep); r_N_N13=gRandom->Poisson(N_N13); r_N_O15=gRandom->Poisson(N_O15); //r_N_F17=gRandom->Poisson(N_F17);

    /// use to incorporate HZ or LZ model uncertainties
    r_N_SN=gRandom->Poisson(N_SN+g_N_SN); r_N_Be7=gRandom->Poisson(N_Be7+g_N_Be7); r_N_pep=gRandom->Poisson(N_pep+g_N_pep); r_N_N13=gRandom->Poisson(N_N13+g_N_N13); r_N_O15=gRandom->Poisson(N_O15+g_N_O15); 
    /// use if you want to take the arbitrary n-sigma values of the HZ model
/*    r_N_SN=gRandom->Poisson(N_SN*(1.+nsig*0.006));
    r_N_Be7=gRandom->Poisson(N_Be7*(1.+nsig*0.06));
    r_N_pep=gRandom->Poisson(N_pep*(1.+nsig*0.01));
    r_N_N13=gRandom->Poisson(N_N13*(1.+nsig*0.15)); 
    r_N_O15=gRandom->Poisson(N_O15*(1.+nsig*0.17));
    r_N_NC=gRandom->Poisson(N_NC*(1.+nsig*0.0));*/


    //randomly sample
    for (int i=0; i<r_N_Mat; i++) {double v=funcMat_smear->GetRandom(); dummyMC->Fill(v); dummyMCmat->Fill(v);}
    for (int i=0; i<r_N_Rn; i++) { double v=funcRn_smear->GetRandom();  dummyMC->Fill(v); dummyMCrn->Fill(v);}
    for (int i=0; i<r_N_Kr; i++) { double v=funcKr_smear->GetRandom();  dummyMC->Fill(v); dummyMCkr->Fill(v);}
    for (int i=0; i<r_N_DB; i++) { double v=funcDB_smear->GetRandom();  dummyMC->Fill(v); dummyMCdb->Fill(v);}

    for (int i=0; i<r_N_SN; i++) { double v=funcSN_smear->GetRandom();  dummyMC->Fill(v); dummyMCnu->Fill(v);}
    for (int i=0; i<r_N_Be7; i++) { double v=funcBe7_smear->GetRandom();  dummyMC->Fill(v); dummyMCBe7->Fill(v);}
    for (int i=0; i<r_N_pep; i++) { double v=funcpep_smear->GetRandom();  dummyMC->Fill(v); dummyMCpep->Fill(v);}
    for (int i=0; i<r_N_N13; i++) { double v=funcN13_smear->GetRandom();  dummyMC->Fill(v); dummyMCN13->Fill(v);}
    for (int i=0; i<r_N_O15; i++) { double v=funcO15_smear->GetRandom();  dummyMC->Fill(v); dummyMCO15->Fill(v);}
    for (int i=0; i<r_N_NC;  i++) { double v=funcNC_smear->GetRandom();   dummyMC->Fill(v); dummyMCNC->Fill(v);}

    dummyMC->Scale(1./exposure);  dummyMCmat->Scale(1./exposure); dummyMCrn->Scale(1./exposure); 
    dummyMCkr->Scale(1./exposure); dummyMCdb->Scale(1./exposure);
    dummyMCnu->Scale(1./exposure); dummyMCBe7->Scale(1./exposure); dummyMCpep->Scale(1./exposure); dummyMCN13->Scale(1./exposure); dummyMCO15->Scale(1./exposure); //dummyMCF17->Scale(1./exposure);
    dummyMCNC->Scale(1./exposure);

    //Fit
    funcMC->SetParameter(0,1.);  funcMC->SetParameter(1,1.);  funcMC->SetParameter(2,1.);  funcMC->SetParameter(3,1.);
    funcMC->SetParameter(4,1.);  funcMC->SetParameter(5,1.);
    funcMC->SetParameter(6,1.); funcMC->SetParameter(7,1.);  funcMC->SetParameter(8,1.);
    funcMC->SetParameter(9,1.);

//    funcMC->FixParameter(0,0);//    funcMC->FixParameter(1,0);//    funcMC->FixParameter(2,0);
//    funcMC->FixParameter(3,0);
//    funcMC->FixParameter(4,0);//    funcMC->FixParameter(5,0);//    funcMC->FixParameter(6,0);//    funcMC->FixParameter(7,0);    //    funcMC->FixParameter(8,0);
//      funcMC->FixParameter(9,0);

    dummyMC->Fit(funcMC,"WL I N Q","",0,3000);   ///WL I N Q
    fitvalueMat=funcMC->GetParameter(0); fitvalueRn=funcMC->GetParameter(1); fitvalueKr=funcMC->GetParameter(2); fitvalueDB=funcMC->GetParameter(3);
    fitvalueSN=funcMC->GetParameter(4); fitvalueBe7=funcMC->GetParameter(5); fitvaluepep=funcMC->GetParameter(6); fitvalueN13=funcMC->GetParameter(7); fitvalueO15=funcMC->GetParameter(8);
    fitvalueNC=funcMC->GetParameter(9);

//    TCanvas *test = new TCanvas("test", "test",265,50,1009,588); test->SetLogy();
//    dummyMC->Draw();

    hmatA->Fill(fitvalueMat); hrnA->Fill(fitvalueRn); hkrA->Fill(fitvalueKr); hdbA->Fill(fitvalueDB);
    hnuA->Fill(fitvalueSN); hBe7A->Fill(fitvalueBe7); hpepA->Fill(fitvaluepep); hN13A->Fill(fitvalueN13); hO15A->Fill(fitvalueO15); //hF17A->Fill(fitvalueF17);
    hNC->Fill(fitvalueNC);
    hsumA->Fill(fitvalueMat+fitvalueRn+fitvalueKr+fitvalueDB+fitvalueSN+fitvalueBe7+fitvaluepep+fitvalueN13+fitvalueO15);

  }


//---------------------------------------------------------------------------------------------------------------------------------------


  cout<<"----------------------------------------------------------------------------------------------------------------------------"<<endl;

  cout<<"MATERIALS"<<endl;
//  for (int j=0; j<nBands; j++) cout<<"  "<<GetStats(hmatA,j)<<endl; //cout<<endl; 
  GetPercentError(hmatA);
  cout<<"RADON"<<endl;
//  for (int j=0; j<nBands; j++) cout<<"  "<<GetStats(hrnA,j)<<endl; //cout<<endl;
  GetPercentError(hrnA);
  cout<<"KRYPTON"<<endl;
//  for (int j=0; j<nBands; j++) cout<<"  "<<GetStats(hkrA,j)<<endl; //cout<<endl;
  GetPercentError(hkrA);
  cout<<"DOUBLE BETA"<<endl;
//  for (int j=0; j<nBands; j++) cout<<"  "<<GetStats(hdbA,j)<<endl; //cout<<endl;
  GetPercentError(hdbA);
  cout<<"SOLAR NU"<<endl;
//  for (int j=0; j<nBands; j++) cout<<"  "<<GetStats(hnuA,j)<<endl; //cout<<endl;
  GetPercentError(hnuA);
  cout<<"Be7"<<endl;
//  for (int j=0; j<nBands; j++) cout<<"  "<<GetStats(hBe7A,j)<<endl; //cout<<endl;
  GetPercentError(hBe7A);
  cout<<"pep"<<endl;
//  for (int j=0; j<nBands; j++) cout<<"  "<<GetStats(hpepA,j)<<endl; //cout<<endl;
  GetPercentError(hpepA);
  cout<<"N13"<<endl;
//  for (int j=0; j<nBands; j++) cout<<"  "<<GetStats(hN13A,j)<<endl; //cout<<endl;
  GetPercentError(hN13A);
  cout<<"O15"<<endl;
//  for (int j=0; j<nBands; j++) cout<<"  "<<GetStats(hO15A,j)<<endl; //cout<<endl;
  GetPercentError(hO15A);
  cout<<"NC"<<endl;
//  for (int j=0; j<nBands; j++) cout<<"  "<<GetStats(hNC,j)<<endl; //cout<<endl;
  GetPercentError(hNC);

  cout<<"Determined at "<<mult<<" tonne-years with "<<nombre<<" trials."<<endl;


//  TFile *MyFile = new TFile(TString::Format("/home/directory/exposure_%ity.root",int(mult)),"recreate");

//  hnuA->Write(); hBe7A->Write(); hpepA->Write(); hN13A->Write(); hO15A->Write();
//  hmatA->Write(); hrnA->Write(); hkrA->Write(); hdbA->Write();

//  MyFile->Close();

  return;
}
