#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TH2F.h>
#include <TGaxis.h>

#endif

enum {kVsRadius, kVsN};

//____________________________________________________________________________________________
void PlotResult_DBSCAN(TString filenametxt = "input.txt");
Bool_t ReadFile(TString filenametxt, vector<double> &nTotMCvec, vector<double> &nPionMCvec, vector<double> &nHeliumMCvec, vector<double> &nTotRecovec, vector<double> &nTotGoodRecovec, vector<double> &nPionGoodRecovec, vector<double> &nHeliumGoodRecovec, vector<double> &radiusvec, vector<double> &nPtsMinvec);
void SetStyle();

//____________________________________________________________________________________________
void PlotResult_DBSCAN(TString filenametxt) {
  
  SetStyle();
  
  vector<double> nTotMC;
  vector<double> nPionMC;
  vector<double> nHeliumMC;
  vector<double> nTotReco;
  vector<double> nTotGoodReco;
  vector<double> nPionGoodReco;
  vector<double> nHeliumGoodReco;
  vector<double> radius;
  vector<double> nPtsMin;
  
  Bool_t read = ReadFile(filenametxt,nTotMC,nPionMC,nHeliumMC,nTotReco,nTotGoodReco,nPionGoodReco,nHeliumGoodReco,radius,nPtsMin);
  if(!read) {return;}
  
  if(nTotMC.size()<2) {
    cerr << "Only one point in file, exit." << endl;
    return;
  }
  
  Int_t xaxis = kVsRadius;
  TString xaxistit = "#varepsilon (#mum) ";
  if(radius[0] == radius[1]) {
    xaxis = kVsN;
    xaxistit = "nMinPts";
  }
  
  TGraphErrors* gEffTot = new TGraphErrors(0);
  gEffTot->GetXaxis()->SetTitle(xaxistit.Data());
  gEffTot->GetYaxis()->SetTitle("Efficiency");
  gEffTot->SetMarkerStyle(kFullCircle);
  gEffTot->SetMarkerSize(0.5);
  gEffTot->SetMarkerColor(kBlack);
  gEffTot->SetLineColor(kBlack);
  gEffTot->SetLineWidth(2);
  TGraphErrors* gEffTotGood = new TGraphErrors(0); 
  gEffTotGood->GetXaxis()->SetTitle(xaxistit.Data());
  gEffTotGood->GetYaxis()->SetTitle("Efficiency");
  gEffTotGood->SetMarkerStyle(kFullSquare);
  gEffTotGood->SetMarkerSize(0.5);
  gEffTotGood->SetMarkerColor(kRed+1);
  gEffTotGood->SetLineColor(kRed+1);
  gEffTotGood->SetLineWidth(2);
  TGraphErrors* gEffPionGood = new TGraphErrors(0); 
  gEffPionGood->GetXaxis()->SetTitle(xaxistit.Data());
  gEffPionGood->GetYaxis()->SetTitle("Efficiency");
  gEffPionGood->SetMarkerStyle(kFullDiamond);
  gEffPionGood->SetMarkerSize(0.8);
  gEffPionGood->SetMarkerColor(kBlue+1);
  gEffPionGood->SetLineColor(kBlue+1);
  gEffPionGood->SetLineWidth(2);
  TGraphErrors* gEffHeliumGood = new TGraphErrors(0); 
  gEffHeliumGood->GetXaxis()->SetTitle(xaxistit.Data());
  gEffHeliumGood->GetYaxis()->SetTitle("Efficiency");
  gEffHeliumGood->SetMarkerStyle(kFullTriangleUp);
  gEffHeliumGood->SetMarkerSize(0.5);
  gEffHeliumGood->SetMarkerColor(kGreen+2);
  gEffHeliumGood->SetLineColor(kGreen+2);
  gEffHeliumGood->SetLineWidth(2);

  for(int iEntry=0; iEntry<nTotReco.size(); iEntry++) {
    double efftot = nTotReco[iEntry]/nTotMC[iEntry];
    double effgoodtot = nTotGoodReco[iEntry]/nTotMC[iEntry];
    double effpion = nPionGoodReco[iEntry]/nPionMC[iEntry];
    double effhelium = nHeliumGoodReco[iEntry]/nHeliumMC[iEntry];
    
    double efftoterr = std::sqrt(efftot/nTotMC[iEntry]*(1-efftot));
    double effgoodtoterr = std::sqrt(effgoodtot/nTotMC[iEntry]*(1-effgoodtot));
    double effpionerr = std::sqrt(effpion/nPionMC[iEntry]*(1-effpion));
    double effheliumerr = std::sqrt(effhelium/nHeliumMC[iEntry]*(1-effhelium));
    
    if(xaxis == kVsRadius) {
      gEffTot->SetPoint(iEntry,radius[iEntry]*10000,efftot);
      gEffTot->SetPointError(iEntry,0.,efftoterr);
      gEffTotGood->SetPoint(iEntry,radius[iEntry]*10000,effgoodtot);
      gEffTotGood->SetPointError(iEntry,0.,effgoodtoterr);
      gEffPionGood->SetPoint(iEntry,radius[iEntry]*10000,effpion);
      gEffPionGood->SetPointError(iEntry,0.,effpionerr);
      gEffHeliumGood->SetPoint(iEntry,radius[iEntry]*10000,effhelium);
      gEffHeliumGood->SetPointError(iEntry,0.,effheliumerr);
    }
    else {
      gEffTot->SetPoint(iEntry,nPtsMin[iEntry],efftot);
      gEffTot->SetPointError(iEntry,0.,efftoterr);
      gEffTotGood->SetPoint(iEntry,nPtsMin[iEntry],effgoodtot);
      gEffTotGood->SetPointError(iEntry,0.,effgoodtoterr);
      gEffPionGood->SetPoint(iEntry,nPtsMin[iEntry],effpion);
      gEffPionGood->SetPointError(iEntry,0.,effpionerr);
      gEffHeliumGood->SetPoint(iEntry,nPtsMin[iEntry],effhelium);
      gEffHeliumGood->SetPointError(iEntry,0.,effheliumerr);
    }
  }
  
  TH1F* hFrameTotEff = 0x0;
  TH1F* hFramePID = 0x0;
  if(xaxis==kVsRadius) {
    hFrameTotEff = new TH1F("hFrameTotEff",Form("Cluster reconstruction efficiency;%s;Efficiency",xaxistit.Data()),10000,14,radius[radius.size()-1]*10010);
    hFramePID = new TH1F("hFramePID",Form("Good reco clusters / true clusters;%s;Efficiency",xaxistit.Data()),10000,14,radius[radius.size()-1]*10010);
  }
  else {
    hFrameTotEff = new TH1F("hFrameTotEff",Form("Cluster reconstruction efficiency;%s;Efficiency",xaxistit.Data()),10000,0,nPtsMin[nPtsMin.size()-1]+2);
    hFramePID = new TH1F("hFramePID",Form("Good reco clusters / true clusters;%s;Efficiency",xaxistit.Data()),10000,0,nPtsMin[radius.size()-1]+2);
  }
  
  TLegend* legTotEff = new TLegend(0.3,0.75,0.7,0.85);
  legTotEff->SetFillStyle(0);
  legTotEff->SetTextSize(0.035);
  legTotEff->AddEntry(gEffTot,"reco clusters / true clusters","pe");
  legTotEff->AddEntry(gEffTotGood,"good reco clusters / true clusters","pe");
    
  TLegend* legTotPID = new TLegend(0.65,0.70,0.85,0.85);
  legTotPID->SetFillStyle(0);
  legTotPID->SetTextSize(0.035);
  legTotPID->AddEntry(gEffTotGood,"all clusters","pe");
  legTotPID->AddEntry(gEffPionGood,"#pi^{#pm}","pe");
  legTotPID->AddEntry(gEffHeliumGood,"^{3}He","pe");

  TCanvas* cTotEff = new TCanvas("cTotEff","",800,800);
  //cTotEff->SetLogx();
  hFrameTotEff->GetYaxis()->SetRangeUser(0.,1.5);
  hFrameTotEff->Draw();
  gEffTot->Draw("PZ");
  gEffTotGood->Draw("PZ");
  legTotEff->Draw("same");
  
  TCanvas* cEffPID = new TCanvas("cEffPID","",800,800);
  //cEffPID->SetLogx();
  hFramePID->GetYaxis()->SetRangeUser(0.,1.2);
  hFramePID->Draw();
  gEffTotGood->Draw("P");
  gEffPionGood->Draw("P");
  gEffHeliumGood->Draw("P");
  legTotPID->Draw("same");
  
  TString outfilename = filenametxt;
  outfilename.ReplaceAll(".txt","_TotEff.pdf");
  cTotEff->SaveAs(outfilename.Data());
  outfilename.ReplaceAll("_TotEff","_PIDEff");
  cEffPID->SaveAs(outfilename.Data());
}

//_____________________________________________________________________________________________
Bool_t ReadFile(TString filenametxt, vector<double> &nTotMCvec, vector<double> &nPionMCvec, vector<double> &nHeliumMCvec, vector<double> &nTotRecovec, vector<double> &nTotGoodRecovec, vector<double> &nPionGoodRecovec, vector<double> &nHeliumGoodRecovec, vector<double> &radiusvec, vector<double> &nPtsMinvec)
{
  ifstream inSet(filenametxt.Data());
  if(!inSet) {
    cerr<<"File "<<filenametxt.Data() <<" does not exists. "<<endl;
    return kFALSE;
  }

  double nTotMC;
  double nPionMC;
  double nHeliumMC;
  double nTotReco;
  double nTotGoodReco;
  double nPionGoodReco;
  double nHeliumGoodReco;
  double radius;
  double nPtsMin;
  while(inSet>>nTotMC>>nPionMC>>nHeliumMC>>nTotReco>>nTotGoodReco>>nPionGoodReco>>nHeliumGoodReco>>radius>>nPtsMin) {
    nTotMCvec.push_back(nTotMC);
    nPionMCvec.push_back(nPionMC);
    nHeliumMCvec.push_back(nHeliumMC);
    nTotRecovec.push_back(nTotReco);
    nTotGoodRecovec.push_back(nTotGoodReco);
    nPionGoodRecovec.push_back(nPionGoodReco);
    nHeliumGoodRecovec.push_back(nHeliumGoodReco);
    radiusvec.push_back(radius);
    nPtsMinvec.push_back(nPtsMin);
  }

  inSet.close();

  return kTRUE;
}

//__________________________________________________________________________________________________________________
void SetStyle() {
  cout << "Setting drawing style!" << endl;
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetTitleSize(0.045,"xyzt");
  gStyle->SetTitleSize(0.045,"t");
  gStyle->SetLabelSize(0.04,"xyzt");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetTitleFont(42,"xy");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  TGaxis::SetMaxDigits(3);
}

