#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include <TRandom1.h>
#include <vector>
#include <TLorentzVector.h>
#include "THnSparse.h"
#include "TRandom3.h"
#include <cstring>
#include <ctime>
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include "TVector3.h"
#include "TGraph.h"
#include <vector>
#include <map>
#include "TFrame.h"
#include "TH1F.h"
#include "TBenchmark.h"
#include "TSystem.h"
#include "./setSPRACEColors.C"

void LambdaMain(double Rf0, double Rf0stat, double Rf0syst, double d0, double d0stat, double d0syst) {
  
const int LALpoints=1;

// ====================== Other Experiments or Theory ====================== 
// ====================== STAR datapoints ====================== 
double LL_STAR_Rf0[LALpoints]={-1.10};
double LL_STAR_Rf0_stat[LALpoints]={0.37};
double LL_STAR_Rf0_sysl[LALpoints]={0.68};
double LL_STAR_Rf0_sysh[LALpoints]={0.08};

double LL_STAR_d0[LALpoints]={8.52};
double LL_STAR_d0_stat[LALpoints]={2.56};
double LL_STAR_d0_sysl[LALpoints]={0.74};
double LL_STAR_d0_sysh[LALpoints]={2.09};

// ====================== STAR Reanalysis  ====================== 
double LL_STAR_Rf0_reana[LALpoints]={0.9};
double LL_STAR_Rf0_stat_reana[LALpoints]={0.35};

double LL_STAR_d0_reana[LALpoints]={5.25};
double LL_STAR_d0_stat_reana[LALpoints]={1.75};

// ====================== NAGARA events ====================== 
double LL_FG_Rf0[LALpoints]={0.77};
double LL_FG_d0[LALpoints]={6.59};

double LL_HKMYY_Rf0[LALpoints]={0.575};
double LL_HKMYY_d0[LALpoints]={6.45};

//  ====================== Lattice  ====================== 
double LL_LQCD_Rf0[LALpoints]={0.69};
double LL_LQCD_Rf0_stat[LALpoints]={0.1};

double LL_LQCD_d0[LALpoints]={5.16};
double LL_LQCD_d0_stat[LALpoints]={0.82};

//  ====================== CMS PbPb ====================== 
double LL_CMSPbPb_Rf0[LALpoints]={0.74};
double LL_CMSPbPb_Rf0_stat[LALpoints]={0.59};
double LL_CMSPbPb_Rf0_sys[LALpoints]={0.33};

double LL_CMSPbPb_d0[LALpoints]={4.2};
double LL_CMSPbPb_d0_stat[LALpoints]={5.7};
double LL_CMSPbPb_d0_sys[LALpoints]={2.9};

// ====================== Results here ====================== 
//  ====================== CMS pPb ====================== 
//integrated
double LL_CMS_Rf0[LALpoints]={Rf0};
double LL_CMS_Rf0_stat[LALpoints]={Rf0stat};
double LL_CMS_Rf0_sys[LALpoints]={Rf0*Rf0syst};

double LL_CMS_d0[LALpoints]={d0};
double LL_CMS_d0_stat[LALpoints]={d0stat};
double LL_CMS_d0_sys[LALpoints]={d0*d0syst};


auto gr_CMS_d0 = new TGraphErrors(LALpoints,LL_CMS_Rf0,LL_CMS_d0,LL_CMS_Rf0_stat,LL_CMS_d0_stat);
gr_CMS_d0->SetMarkerStyle(20);
gr_CMS_d0->SetMarkerColor(SPcolors[SPred]);
gr_CMS_d0->SetMarkerSize(1.5);
gr_CMS_d0->SetLineColor(SPcolors[SPred]);
auto gr_CMS_d0_sys = new TGraphErrors(LALpoints,LL_CMS_Rf0,LL_CMS_d0,LL_CMS_Rf0_sys,LL_CMS_d0_sys);
gr_CMS_d0_sys->SetMarkerColor(0);
gr_CMS_d0_sys->SetFillColor(0);
gr_CMS_d0_sys->SetFillStyle(1);
gr_CMS_d0_sys->SetLineColor(SPcolors[SPred]);

auto gr_CMS_d0_sysb = new TGraphErrors(LALpoints,LL_CMS_Rf0,LL_CMS_d0,LL_CMS_Rf0_sys,LL_CMS_d0_sys);
gr_CMS_d0_sysb->SetMarkerColor(0);
gr_CMS_d0_sysb->SetFillColor(0);
gr_CMS_d0_sysb->SetFillStyle(1);
gr_CMS_d0_sysb->SetLineColor(1);

auto gr_STAR_d0 = new TGraphErrors(LALpoints,LL_STAR_Rf0,LL_STAR_d0,LL_STAR_Rf0_stat,LL_STAR_d0_stat);
gr_STAR_d0->SetMarkerStyle(30);
gr_STAR_d0->SetMarkerColor(1);
gr_STAR_d0->SetMarkerSize(1.7);

auto gr_STAR_sys = new TGraphAsymmErrors(LALpoints,LL_STAR_Rf0,LL_STAR_d0,LL_STAR_Rf0_sysh,LL_STAR_Rf0_sysl,LL_STAR_d0_sysl,LL_STAR_d0_sysh);
gr_CMS_d0_sysb->SetMarkerColor(0);
gr_CMS_d0_sysb->SetFillColor(0);
gr_CMS_d0_sysb->SetFillStyle(0);
gr_CMS_d0_sysb->SetLineColor(1);

auto gr_STARreana_d0 = new TGraphErrors(LALpoints,LL_STAR_Rf0_reana,LL_STAR_d0_reana,LL_STAR_Rf0_stat_reana,LL_STAR_d0_stat_reana);
gr_STARreana_d0->SetMarkerColor(0);
gr_STARreana_d0->SetFillColor(SPcolors[SPsgrey]);
gr_STARreana_d0->SetFillStyle(1001);
gr_STARreana_d0->SetLineColor(0);

auto gr_FG_d0 = new TGraphErrors(LALpoints,LL_FG_Rf0,LL_FG_d0,0,0);
gr_FG_d0->SetMarkerStyle(55);
gr_FG_d0->SetMarkerColor(kBlue);
gr_FG_d0->SetMarkerSize(1.7);
gr_FG_d0->SetLineColor(kBlue);

auto gr_HKMYY_d0 = new TGraphErrors(LALpoints,LL_HKMYY_Rf0,LL_HKMYY_d0,0,0);
gr_HKMYY_d0->SetMarkerStyle(55);
gr_HKMYY_d0->SetMarkerColor(kBlue);
gr_HKMYY_d0->SetMarkerSize(1.7);
gr_HKMYY_d0->SetLineColor(kBlue);

// LQCD
auto gr_LQCD_d0 = new TGraphErrors(LALpoints,LL_LQCD_Rf0,LL_LQCD_d0,LL_LQCD_Rf0_stat,LL_LQCD_d0_stat);
gr_LQCD_d0->SetMarkerStyle(44);
gr_LQCD_d0->SetMarkerColor(SPcolors[SPpurple]);
gr_LQCD_d0->SetMarkerSize(1.7);
gr_LQCD_d0->SetLineColor(SPcolors[SPpurple]);

TMultiGraph *mg2 = new TMultiGraph();
mg2->Add(gr_STAR_d0,"P");

TCanvas *c2 = new TCanvas("c2","",1000,700);
mg2->Draw("a");
gr_STAR_sys->Draw("5");
gr_STAR_d0->Draw("PSAME");
gr_STARreana_d0->Draw("5");
gr_FG_d0->Draw("PSAME");
gr_HKMYY_d0->Draw("PSAME");
gr_LQCD_d0->Draw("PSAME");
gr_CMS_d0->Draw("PSAME");
gr_CMS_d0_sys->Draw("5");

mg2->GetYaxis()->SetRangeUser(0.0,12.);
mg2->GetXaxis()->SetLimits(-2,2);
mg2->GetYaxis()->SetTitle("#scale[1.6]{d_{0} [fm]}");
mg2->GetXaxis()->SetTitle("#scale[1.6]{f_{0} [fm]}");
mg2->GetXaxis()->SetNdivisions(506);
mg2->GetYaxis()->SetNdivisions(506);
mg2->GetXaxis()->CenterTitle(1);
mg2->GetYaxis()->CenterTitle(1);
mg2->GetXaxis()->SetTickSize(0.03);
mg2->GetYaxis()->SetTickSize(0.03);
mg2->GetXaxis()->SetLabelSize(0.05);
mg2->GetYaxis()->SetLabelSize(0.05);
mg2->GetYaxis()->SetLabelOffset(0.0025);
mg2->GetXaxis()->SetLabelOffset(0.0025);
mg2->GetYaxis()->SetTitleOffset(1.3);
mg2->GetXaxis()->SetTitleOffset(1.3);  

TLegend *legend2 = new TLegend(0.15,0.15,0.5,0.45);
legend2->SetBorderSize(0);
legend2->AddEntry(gr_CMS_d0,"CMS pPb at #sqrt{s_{NN}} = 8.16 TeV","P");
legend2->AddEntry(gr_STAR_d0,"STAR AuAu at #sqrt{s_{NN}} = 200 GeV","P");
legend2->AddEntry(gr_STARreana_d0,"STAR reanalysis","F");
legend2->AddEntry(gr_FG_d0,"NAGARA","P");
legend2->AddEntry(gr_LQCD_d0,"LQCD","P");
legend2->AddEntry(gr_CMS_d0_sysb,"Systematic Uncertainties","F");
legend2->Draw();


std::ostringstream os2; 
TLatex* texCMS2 = new TLatex;
texCMS2->SetNDC();
texCMS2->SetTextFont(42);
texCMS2->SetTextSize(0.05);
texCMS2->SetTextAlign(12);
texCMS2->DrawLatex(.1,.93," #font[61]{CMS} #it{Work in Progress}");
texCMS2->DrawLatex(.67,.8,"#scale[1.2]{#font[61]{#Lambda#Lambda#oplus#bar{#Lambda}#bar{#Lambda}}}");
texCMS2->DrawLatex(0.6,0.931, "#scale[0.95]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");   
c2->Update();
c2->Print("Rf0vsd0.pdf");
c2->Update();

}

