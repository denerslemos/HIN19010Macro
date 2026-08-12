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
#include <vector>
#include <map>
#include "TFrame.h"
#include "TH1D.h"
#include "TBenchmark.h"
#include "TSystem.h"
#include "TMinuit.h"
#include "./setSPRACEColors.C"

#include "TStyle.h"

using namespace std;

void normalize(double qminN, double qmaxN, TH1D* aux_histo_num, TH1D* aux_histo_den){

	Int_t bin_for_normInt_min = aux_histo_num->GetXaxis()->FindBin(qminN);
	Int_t bin_for_normInt_max = aux_histo_num->GetXaxis()->FindBin(qmaxN);
	Double_t int_num_controlRegion = aux_histo_num->Integral(bin_for_normInt_min,bin_for_normInt_max);
	Double_t int_den_controlRegion = aux_histo_den->Integral(bin_for_normInt_min,bin_for_normInt_max);
	aux_histo_den->Scale(int_num_controlRegion/int_den_controlRegion);

}

void sideband(TH1D* aux_histo_ss, TH1D* aux_histo_sb, TH1D* aux_histo_bb, double fss, double fsb, double fbb){

	aux_histo_ss->Scale(1./aux_histo_ss->Integral());
	aux_histo_sb->Scale(fsb/aux_histo_sb->Integral());
	aux_histo_bb->Scale(fbb/aux_histo_bb->Integral());
	aux_histo_ss->Add(aux_histo_sb,-1);
	aux_histo_ss->Add(aux_histo_bb,-1);
	aux_histo_ss->Scale(1/fss);

}

void Standard_N(int syst, int multbin, int ktbin) {

TString V0 = "_V0_"; //just a label

//qinv range
Double_t qmin = 0.0;
Double_t qmax = 5.999999;

//qinv normalization range
Double_t qmin_Norm = qmin;
Double_t qmax_Norm = qmax;

Double_t qminM = 0.0; //max qinv for mass fits (only for usemass = true)
Double_t qmaxM = 2.999999; //max qinv for mass fits (only for usemass = true)

//kT bin
Double_t ktmin, ktmax;
//int ktbin = 0; // == 0 means integrated
if(ktbin==0){ 
	ktmin = 0.0;
	ktmax = 1.999999;
}else if(ktbin==1){ // for K0sK0s only
}

Double_t ntkmin, ntkmax;
if(ktbin != 0) multbin = 0;
//ntrk bin
if(multbin==0){ // == 0 means integrated (LL+ALAL and LAL)
	ntkmin = 0.0;
	ntkmax = 399.9999999;
}else if(multbin==1){ // MB (LAL)
	ntkmin = 0.0;
	ntkmax = 119.9999999;
}else if(multbin==2){ // HM (LAL)
	ntkmin = 185.0;
	ntkmax = 249.9999999;
}else if(multbin==3){ // additional ones for K0sK0s only
}

//root stuff
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

//reading the input file
TString filename; 
filename = Form("./Merged_V0HBT_%i.root",syst);

TFile* f = new TFile(Form("%s",filename.Data()),"READ");

cout << "File Read" << endl;

//reading mass histograms
THnSparseD* h_K0s_Mass=(THnSparseD*)f->Get("HBT_K0s/hMass_K0s_K0s_T");
THnSparseD* h_Lam_Mass=(THnSparseD*)f->Get("HBT_Lam/hMass_Lam_Lam_T");
THnSparseD* h_ALam_Mass=(THnSparseD*)f->Get("HBT_ALam/hMass_ALam_ALam_T");
THnSparseD* h_Lam_Mass_LAL=(THnSparseD*)f->Get("HBT_LAL/hMass_Lam_T");
THnSparseD* h_ALam_Mass_LAL=(THnSparseD*)f->Get("HBT_LAL/hMass_ALam_T");
THnSparseD* h_KL_Mass=(THnSparseD*)f->Get("HBT_K0sLam/hMass_K0sL_T");
THnSparseD* h_KAL_Mass=(THnSparseD*)f->Get("HBT_K0sALam/hMass_K0sAL_T");
THnSparseD* h_LK_Mass=(THnSparseD*)f->Get("HBT_K0sLam/hMass_LamK_T");
THnSparseD* h_ALK_Mass=(THnSparseD*)f->Get("HBT_K0sALam/hMass_ALamK_T");


//Find the bins and set the ranges

Int_t aux_kt_min = h_K0s_Mass->GetAxis(1)->FindBin(ktmin+0.0001); 
Int_t aux_kt_max = h_K0s_Mass->GetAxis(1)->FindBin(ktmax-0.0001); 
Int_t aux_ntk_min = h_K0s_Mass->GetAxis(2)->FindBin(ntkmin+0.0001); 
Int_t aux_ntk_max = h_K0s_Mass->GetAxis(2)->FindBin(ntkmax-0.0001); 

h_K0s_Mass->GetAxis(1)->SetRange(aux_kt_min,aux_kt_max);
h_K0s_Mass->GetAxis(2)->SetRange(aux_ntk_min,aux_ntk_max);
h_Lam_Mass->GetAxis(1)->SetRange(aux_kt_min,aux_kt_max);
h_Lam_Mass->GetAxis(2)->SetRange(aux_ntk_min,aux_ntk_max);
h_ALam_Mass->GetAxis(1)->SetRange(aux_kt_min,aux_kt_max);
h_ALam_Mass->GetAxis(2)->SetRange(aux_ntk_min,aux_ntk_max);
h_Lam_Mass_LAL->GetAxis(1)->SetRange(aux_kt_min,aux_kt_max);
h_Lam_Mass_LAL->GetAxis(2)->SetRange(aux_ntk_min,aux_ntk_max);
h_ALam_Mass_LAL->GetAxis(1)->SetRange(aux_kt_min,aux_kt_max);
h_ALam_Mass_LAL->GetAxis(2)->SetRange(aux_ntk_min,aux_ntk_max);
h_KL_Mass->GetAxis(1)->SetRange(aux_kt_min,aux_kt_max);
h_KL_Mass->GetAxis(2)->SetRange(aux_ntk_min,aux_ntk_max);
h_KAL_Mass->GetAxis(1)->SetRange(aux_kt_min,aux_kt_max);
h_KAL_Mass->GetAxis(2)->SetRange(aux_ntk_min,aux_ntk_max);
h_LK_Mass->GetAxis(1)->SetRange(aux_kt_min,aux_kt_max);
h_LK_Mass->GetAxis(2)->SetRange(aux_ntk_min,aux_ntk_max);
h_ALK_Mass->GetAxis(1)->SetRange(aux_kt_min,aux_kt_max);
h_ALK_Mass->GetAxis(2)->SetRange(aux_ntk_min,aux_ntk_max);

Int_t aux_q_min = h_K0s_Mass->GetAxis(3)->FindBin(qminM); 
Int_t aux_q_max = h_K0s_Mass->GetAxis(3)->FindBin(qmaxM); 
h_K0s_Mass->GetAxis(3)->SetRange(aux_q_min,aux_q_max);
h_Lam_Mass->GetAxis(3)->SetRange(aux_q_min,aux_q_max);
h_ALam_Mass->GetAxis(3)->SetRange(aux_q_min,aux_q_max);
h_Lam_Mass_LAL->GetAxis(3)->SetRange(aux_q_min,aux_q_max);
h_ALam_Mass_LAL->GetAxis(3)->SetRange(aux_q_min,aux_q_max);
h_KL_Mass->GetAxis(3)->SetRange(aux_q_min,aux_q_max);
h_KAL_Mass->GetAxis(3)->SetRange(aux_q_min,aux_q_max);
h_LK_Mass->GetAxis(3)->SetRange(aux_q_min,aux_q_max);
h_ALK_Mass->GetAxis(3)->SetRange(aux_q_min,aux_q_max);

//Make 1D projections

TH1D *h_K0s_Mass_1D=(TH1D*)h_K0s_Mass->Projection(0,"e");
TH1D *h_Lam_Mass_1D=(TH1D*)h_Lam_Mass->Projection(0,"e");
TH1D *h_ALam_Mass_1D=(TH1D*)h_ALam_Mass->Projection(0,"e");
TH1D *h_Lam_Mass_LAL_1D=(TH1D*)h_Lam_Mass_LAL->Projection(0,"e");
TH1D *h_ALam_Mass_LAL_1D=(TH1D*)h_ALam_Mass_LAL->Projection(0,"e");
TH1D *h_KL_Mass_1D=(TH1D*)h_KL_Mass->Projection(0,"e");
TH1D *h_KAL_Mass_1D=(TH1D*)h_KAL_Mass->Projection(0,"e");
TH1D *h_LK_Mass_1D=(TH1D*)h_LK_Mass->Projection(0,"e");
TH1D *h_ALK_Mass_1D=(TH1D*)h_ALK_Mass->Projection(0,"e");


h_K0s_Mass->Reset("ICESM");
h_Lam_Mass->Reset("ICESM");
h_ALam_Mass->Reset("ICESM");
h_Lam_Mass_LAL->Reset("ICESM");
h_ALam_Mass_LAL->Reset("ICESM");
h_KL_Mass->Reset("ICESM");
h_KAL_Mass->Reset("ICESM");
h_LK_Mass->Reset("ICESM");
h_ALK_Mass->Reset("ICESM");

//=================Fits=================

cout << "Start the Mass Fits" << endl;

const std::string gaus_func = "[2]*TMath::Gaus(x, [0], [1])";
const std::string poly_bkg0= "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x";//only for plots
const std::string poly_bkg = "[5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x";//background
const std::string double_gaussian = "[3]*TMath::Gaus(x, [0], [1]) + [4]*TMath::Gaus(x, [0], [2])"; //signal
const std::string massfunc = double_gaussian + " + " + poly_bkg;

//K0s stuff
int iter_ks = 5;
Int_t Nbins_ks = 270;
std::pair<double, double>  full_range_ks(0.450, 0.550);
std::pair<double, double>  sb1_ks(0.450, 0.482);
std::pair<double, double>  sb2_ks(0.513, 0.550);
double mass_ks = 0.497648;
double sigma_ks = 0.01;
double Dmass_ks = 0.005;

//K0s used in K0sK0s
cout << "==================K0s Mass for  K0sK0s==================" << endl;
TH1D* massks = (TH1D*)h_K0s_Mass_1D->Clone();
//massks->Sumw2();
TF1 *func = new TF1("func", massfunc.c_str(), full_range_ks.first, full_range_ks.second);

//sideband fit
int cnt = 0;
for ( int i = 0; i < massks->GetNbinsX(); i++ ) {
   double center = massks->GetBinCenter(i+1);
   if ( (center > sb1_ks.first and center < sb1_ks.second) or (center > sb2_ks.first and center < sb2_ks.second) ) {
   cnt++;
   }
}
TGraphErrors gr(cnt);
cnt = 0;
for ( int i = 0; i < massks->GetNbinsX(); i++ ) {
    double center = massks->GetBinCenter(i+1);
    if ( (center > sb1_ks.first and center < sb1_ks.second) or (center > sb2_ks.first and center < sb2_ks.second) ) {
       gr.GetX()[cnt] = center;
       gr.GetY()[cnt] = massks->GetBinContent(i+1);
       gr.GetEY()[cnt] = massks->GetBinError(i+1);
       cnt++;
    }
}

//sideband fit parameters
func->SetParameter(5, 10.);
func->SetParameter(6, 1.);
func->SetParameter(7, 1.);
func->SetParameter(8, 1.);
func->SetParameter(9, 1.);

//fix the parameters in order to get only bkg
func->FixParameter(0, mass_ks);
func->FixParameter(1, sigma_ks);
func->FixParameter(2, sigma_ks);
func->FixParameter(3, 0.);
func->FixParameter(4, 0.);

//make sideband fit
gr.Fit(func, "q0");
gr.Fit(func, "q0");
gr.Fit(func, "q0");

//Set sideband parameters to full fit
func->FixParameter(5, func->GetParameter(5));
func->FixParameter(6, func->GetParameter(6));
func->FixParameter(7, func->GetParameter(7));
func->FixParameter(8, func->GetParameter(8));
func->FixParameter(9, func->GetParameter(9));

//release the signal parameters
func->ReleaseParameter(0);
func->ReleaseParameter(1);
func->ReleaseParameter(2);
func->ReleaseParameter(3);
func->ReleaseParameter(4);
func->SetParLimits(0, mass_ks-Dmass_ks, mass_ks+Dmass_ks);
func->SetParLimits(1, 0, sigma_ks);
func->SetParLimits(2, 0, sigma_ks);
func->SetParLimits(3, 0, 0.8*massks->GetBinContent(massks->FindBin(mass_ks)));
func->SetParLimits(4, 0, 0.8*massks->GetBinContent(massks->FindBin(mass_ks)));

//make a full fit
massks->Fit(func, "q0", "", full_range_ks.first, full_range_ks.second);
massks->Fit(func, "q0", "", full_range_ks.first, full_range_ks.second);
massks->Fit(func, "q0", "", full_range_ks.first, full_range_ks.second);


//now release all the parameters
func->ReleaseParameter(0);
func->ReleaseParameter(1);
func->ReleaseParameter(2);
func->ReleaseParameter(3);
func->ReleaseParameter(4);
func->ReleaseParameter(5);
func->ReleaseParameter(6);
func->ReleaseParameter(7);
func->ReleaseParameter(8);
func->ReleaseParameter(9);
func->SetParLimits(0, mass_ks-Dmass_ks, mass_ks+Dmass_ks);
func->SetParLimits(1, 0, sigma_ks);
func->SetParLimits(2, 0, sigma_ks);
func->SetParLimits(3, 0, 0.8*massks->GetBinContent(massks->FindBin(mass_ks)));
func->SetParLimits(4, 0, 0.8*massks->GetBinContent(massks->FindBin(mass_ks)));



//make final fits
for ( int i = 0; i < iter_ks; i++ ) {
   massks->Fit(func, "q0", "", full_range_ks.first, full_range_ks.second);
}

//make signal and bkg functions for plots
TF1 * func_signal   = 0;
TF1 * func_bckgnd   = 0;
TF1 * gaus1_signal_ks   = 0;
TF1 * gaus2_signal_ks   = 0;
func_signal  = new TF1("func_signal" , double_gaussian.c_str(), full_range_ks.first, full_range_ks.second);
func_bckgnd = new TF1("func_bckgnd" , poly_bkg0.c_str(), full_range_ks.first, full_range_ks.second);
gaus1_signal_ks = new TF1("gaus1_signal_ks", gaus_func.c_str(), full_range_ks.first, full_range_ks.second);
gaus2_signal_ks = new TF1("gaus2_signal_ks", gaus_func.c_str(), full_range_ks.first, full_range_ks.second);

//set parameters for signal and bkg
func_signal->SetParameter(0, func->GetParameter(0));
func_signal->SetParameter(1, func->GetParameter(1));
func_signal->SetParameter(2, func->GetParameter(2));
func_signal->SetParameter(3, func->GetParameter(3));
func_signal->SetParameter(4, func->GetParameter(4));
func_bckgnd->SetParameter(0, func->GetParameter(5));
func_bckgnd->SetParameter(1, func->GetParameter(6));
func_bckgnd->SetParameter(2, func->GetParameter(7));
func_bckgnd->SetParameter(3, func->GetParameter(8));
func_bckgnd->SetParameter(4, func->GetParameter(9));

//set parameters for gaussians
gaus1_signal_ks->SetParameter(0, func->GetParameter(0));
gaus1_signal_ks->SetParameter(1, func->GetParameter(1));
gaus1_signal_ks->SetParameter(2, func->GetParameter(3));
gaus2_signal_ks->SetParameter(0, func->GetParameter(0));
gaus2_signal_ks->SetParameter(1, func->GetParameter(2));
gaus2_signal_ks->SetParameter(2, func->GetParameter(4));

//make pull and ratio
double vny[400];
double xvn[400];
double pullvn[400];
double ratiovn[400];
double err_pullvn[400];
double err_ratiovn[400];
for(int k=0;k<massks->GetNbinsX();k++)
{
            vny[k] = massks->GetBinContent(k+1);
            xvn[k] = massks->GetBinCenter(k+1);
        	if(massks->GetBinContent(k+1)!=0){
            pullvn[k] = (vny[k] - func->Eval(xvn[k]))/massks->GetBinError(k+1);
            ratiovn[k] = vny[k]/func->Eval(xvn[k]);
            err_pullvn[k] = massks->GetBinError(k+1);
            err_ratiovn[k] = (vny[k]/func->Eval(xvn[k]))*(massks->GetBinError(k+1)/vny[k]);
            }else{
            pullvn[k] = 1.0;
            ratiovn[k] = 1.0;
            err_pullvn[k]=0.0;
            err_ratiovn[k]=0.0;
            }
}

//make the plots
TH1D* ksmass = (TH1D*)massks->Clone();
//ksmass->Sumw2();
TGraphErrors* TGpullx = new TGraphErrors(400,xvn,pullvn,0,0);
TGraphErrors* TGratiox= new TGraphErrors(400,xvn,ratiovn,0,err_ratiovn);
TCanvas* c1 = new TCanvas("c1","c1",600,600);
c1->cd();
TPad* pad1x = new TPad("pad1x","top pad",0.0,0.4,1.0,1.0);
TPad* pad2x = new TPad("pad3x","middle pad",0.0,0.0,1.0,0.25);
TPad* pad3x = new TPad("pad2x","bottom pad",0.0,0.25,1.0,0.4);
pad1x->SetTopMargin(0.1);
pad1x->SetBottomMargin(0.0);
pad1x->SetRightMargin(0.038);
pad1x->SetLeftMargin(0.15);
pad1x->Draw();
pad2x->SetFrameFillStyle(4000);
pad2x->SetTopMargin(0.0);
pad2x->SetBottomMargin(0.37);
pad2x->SetRightMargin(0.038);
pad2x->SetLeftMargin(0.15);
pad2x->Draw();
pad3x->SetFrameFillStyle(4000);
pad3x->SetTopMargin(0.0);
pad3x->SetBottomMargin(0.0);
pad3x->SetRightMargin(0.038);
pad3x->SetLeftMargin(0.15);
pad3x->Draw();
pad1x->cd();
gPad->SetTickx();
gPad->SetTicky();
ksmass->SetTitle("");
ksmass->Draw("E");
TGaxis::SetMaxDigits(3);
ksmass->GetYaxis()->SetRangeUser(0.1,1.05*ksmass->GetMaximum());
ksmass->GetXaxis()->SetRangeUser(full_range_ks.first,full_range_ks.second);
ksmass->SetMarkerStyle(24);
ksmass->SetMarkerSize(1.5);
ksmass->GetYaxis()->SetTitle("#scale[1.7]{Candidates / 0.5 MeV}");
ksmass->GetYaxis()->CenterTitle(1);
ksmass->GetYaxis()->SetNdivisions(810);
ksmass->GetXaxis()->SetTickSize(0.03);
ksmass->GetYaxis()->SetTickSize(0.02);
ksmass->GetXaxis()->SetLabelSize(0.055);
ksmass->GetYaxis()->SetLabelSize(0.058);
ksmass->GetYaxis()->SetLabelOffset(0.01);
ksmass->GetYaxis()->SetTitleOffset(1.8);
func->Draw("LSAME");
func->SetLineColor(SPcolors[SPdred]);
func_bckgnd->Draw("LSAME");
func_bckgnd->SetLineColor(SPcolors[SPdred]);
func_bckgnd->SetLineStyle(2);
double meanks = func->GetParameter(0);
double sigma1ks = func->GetParameter(1);
double sigma2ks = func->GetParameter(2);
double gaus1_yield_ks = gaus1_signal_ks->Integral(meanks - 2*sigma1ks, meanks + 2*sigma1ks)/massks->GetBinWidth(1);
double gaus2_yield_ks = gaus2_signal_ks->Integral(meanks - 2*sigma2ks, meanks + 2*sigma2ks)/massks->GetBinWidth(1);
double gausTot_yield_ks = gaus1_yield_ks + gaus2_yield_ks;
double rms_gaus1_sig_ks = gaus1_yield_ks/gausTot_yield_ks;
double rms_gaus2_sig_ks = gaus2_yield_ks/gausTot_yield_ks;
double rms_true_ks = TMath::Sqrt(rms_gaus1_sig_ks*sigma1ks*sigma1ks + rms_gaus2_sig_ks*sigma2ks*sigma2ks);
std::pair<double, double>  signal_range_ks(meanks - 2*rms_true_ks, meanks + 2*rms_true_ks);
std::pair<double, double>  side_range_ks(meanks - 3*rms_true_ks, meanks + 3*rms_true_ks);
double S_ks = func_signal->Integral(signal_range_ks.first, signal_range_ks.second)/massks->GetBinWidth(1);
double B_ks = func_bckgnd->Integral(signal_range_ks.first, signal_range_ks.second)/massks->GetBinWidth(1);
double Fsig_ks = S_ks/(S_ks+B_ks);
double sig_ks = S_ks/sqrt(S_ks+B_ks);
double chi2ks = func->GetChisquare()/func->GetNDF();
double probks = TMath::Prob(func->GetChisquare(),func->GetNDF());
std::ostringstream os; 
TLatex* texCMS = new TLatex;
texCMS->SetNDC();
texCMS->SetTextFont(42);
texCMS->SetTextSize(0.05);
texCMS->SetTextAlign(12);
texCMS->DrawLatex(0.25,.95," #font[61]{CMS} #it{Work in Progress}");
texCMS->DrawLatex(0.70,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS->DrawLatex(0.26,0.4, "#scale[1.7]{K^{0}_{S}}");
// texCMS->DrawLatex(0.225,0.2, "#scale[1.0]{Used in K^{0}_{S}K^{0}_{S}}"); 
// texCMS->DrawLatex(0.23,0.15, "#scale[1.0]{correlations}"); 
// texCMS->DrawLatex(0.75,0.15, "This Thesis"); 
os << "#scale[1.0]{s/#sqrt{s+b}: " << std::setprecision(4) << sig_ks << "}" << std::setprecision(2);           
texCMS->DrawLatex(0.68,0.60, os.str(  ).c_str(  )); 
os.str(std::string());
os << "#scale[1.0]{s/(s+b): " << std::setprecision(4) << Fsig_ks*100 << " %}" << std::setprecision(2); 
texCMS->DrawLatex(0.68,0.52, os.str(  ).c_str(  ));  
os.str(std::string()); 
os << "#scale[1.0]{Mean: " << std::setprecision(4) << meanks << " GeV}" << std::setprecision(6);   
texCMS->DrawLatex(0.68,0.46, os.str(  ).c_str(  ));
os.str(std::string()); 
os << "#scale[1.0]{#sigma: " << std::setprecision(4) << rms_true_ks << " }" << std::setprecision(2); 
texCMS->DrawLatex(0.68,0.40, os.str(  ).c_str(  )); 
os.str(std::string()); 
os << "#scale[1.0]{Yield: " << std::setprecision(4) << S_ks << " }" << std::setprecision(2); 
texCMS->DrawLatex(0.68,0.33, os.str(  ).c_str(  )); 
os.str(std::string()); 
os << "#scale[1.0]{#chi^{2}/ndf: " << std::setprecision(4) << chi2ks << " }" << std::setprecision(4); 
texCMS->DrawLatex(0.68,0.26, os.str(  ).c_str(  ));  
os.str(std::string()); 
// os << "#scale[1.0]{Prob: " << std::setprecision(4) << probks << "}" << std::setprecision(2); 
// texCMS->DrawLatex(0.68,0.20, os.str(  ).c_str(  )); 
// os.str(std::string()); 
os << "#scale[1.]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS->DrawLatex(0.68,0.75, os.str(  ).c_str(  )); 
os.str(std::string()); 
os << "#scale[1.]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS->DrawLatex(0.68,0.67, os.str(  ).c_str(  ));  
os.str(std::string()); 
os << "#scale[1.]{ " << std::setprecision(4) << qminM << std::setprecision(4) <<" < q_{inv} < "  << qmaxM << " GeV}" << std::setprecision(2); 
texCMS->DrawLatex(0.68,0.83, os.str(  ).c_str(  ));  
os.str(std::string()); 
Float_t yyymax = ksmass->GetMaximum();
TLine *linexx = new TLine(signal_range_ks.first,0.1,signal_range_ks.first,yyymax);
linexx->SetLineColor(SPcolors[SPdblue]);
linexx->Draw();
TLine *line2xx = new TLine(signal_range_ks.second,0.1,signal_range_ks.second,yyymax);
line2xx->SetLineColor(SPcolors[SPdblue]);
line2xx->Draw();
TLegend* leg = new TLegend(0.16,0.6,0.43,0.85);
leg->AddEntry(ksmass,"Data");
leg->AddEntry(func,"Signal+Background","L");
leg->AddEntry(func_bckgnd,"Background","L");
leg->AddEntry(linexx,"Peak Region","L");
leg->SetBorderSize(0);
leg->Draw();
c1->Update();
pad2x->cd();
gPad->SetTickx();
gPad->SetTicky();
double lowx = -49.1;
double highx = 49.1;
TH1F* frame_ks = pad2x->cd()->DrawFrame(full_range_ks.first,lowx,full_range_ks.second,highx);
frame_ks->GetXaxis()->SetTitle("#scale[4.3]{#pi^{+}#pi^{-} invariant mass [GeV]}");
frame_ks->GetYaxis()->SetTitle("#scale[4.0]{Pull}");
frame_ks->GetXaxis()->CenterTitle(1);
frame_ks->GetYaxis()->CenterTitle(1);
frame_ks->GetXaxis()->SetNdivisions(510);
frame_ks->GetYaxis()->SetNdivisions(406);
frame_ks->GetXaxis()->SetTickSize(0.06);
frame_ks->GetYaxis()->SetTickSize(0.02);
frame_ks->GetXaxis()->SetTitleOffset(4.5);
frame_ks->GetYaxis()->SetTitleOffset(1.8);
frame_ks->GetXaxis()->SetLabelSize(0.14);
frame_ks->GetYaxis()->SetLabelSize(0.14);
frame_ks->GetXaxis()->SetLabelOffset(0.01);
frame_ks->GetYaxis()->SetLabelOffset(0.013);
TGpullx->SetMarkerStyle(25);
TGpullx->SetMarkerSize(1.3);
TGpullx->SetMarkerColor(kBlack);
TGpullx->Draw("P");
TLine* liney = new TLine(full_range_ks.first, 0, full_range_ks.second, 0);
liney->SetLineStyle(2);
liney->Draw("same");
pad3x->cd();
gPad->SetTickx();
gPad->SetTicky();
double low2x = 0.1;
double high2x = 1.899999;
frame_ks = pad3x->cd()->DrawFrame(full_range_ks.first,low2x,full_range_ks.second,high2x);
frame_ks->GetYaxis()->SetTitle("#scale[6.3]{#frac{Data}{Fit}}");
frame_ks->GetXaxis()->CenterTitle(1);
frame_ks->GetYaxis()->CenterTitle(1);
frame_ks->GetXaxis()->SetNdivisions(510);
frame_ks->GetYaxis()->SetNdivisions(406);
frame_ks->GetXaxis()->SetTickSize(0.1);
frame_ks->GetYaxis()->SetTickSize(0.01);
frame_ks->GetYaxis()->SetTitleOffset(1.8);
frame_ks->GetXaxis()->SetLabelSize(0.12);
frame_ks->GetYaxis()->SetLabelSize(0.22);
frame_ks->GetYaxis()->SetLabelOffset(0.01);
TGratiox->SetMarkerStyle(25);
TGratiox->SetMarkerSize(1.3);
TGratiox->SetMarkerColor(kBlack);
TGratiox->Draw("P");
TLine* line2y = new TLine(full_range_ks.first, 1, full_range_ks.second, 1);
line2y->SetLineStyle(2);
line2y->Draw("same");
c1->Print("results/Mass"+V0+"KSKS_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");
c1->Update();


//K0s used in K0sLambda
cout << "==================K0s Mass for  K0sL==================" << endl;
TH1D* massksL = (TH1D*)h_KL_Mass_1D->Clone();
//massksL->Sumw2();
TF1 *funcKL = new TF1("funcKL", massfunc.c_str(), full_range_ks.first, full_range_ks.second);

//sideband fit
int cntKL = 0;
for ( int i = 0; i < massksL->GetNbinsX(); i++ ) {
    double center = massksL->GetBinCenter(i+1);
    if ( (center > sb1_ks.first and center < sb1_ks.second) or (center > sb2_ks.first and center < sb2_ks.second) ) {
       cntKL++;
    }
}
TGraphErrors grKL(cntKL);
cntKL = 0;
for ( int i = 0; i < massksL->GetNbinsX(); i++ ) {
    double center = massksL->GetBinCenter(i+1);
    if ( (center > sb1_ks.first and center < sb1_ks.second) or (center > sb2_ks.first and center < sb2_ks.second) ) {
       grKL.GetX()[cntKL] = center;
       grKL.GetY()[cntKL] = massksL->GetBinContent(i+1);
       grKL.GetEY()[cntKL] = massksL->GetBinError(i+1);
       cntKL++;
    }
}

//sideband fit parameters
funcKL->SetParameter(5, 10.);
funcKL->SetParameter(6, 1.);
funcKL->SetParameter(7, 1.);
funcKL->SetParameter(8, 1.);
funcKL->SetParameter(9, 1.);

//fix the parameters in order to get only bkg
funcKL->FixParameter(0, mass_ks);
funcKL->FixParameter(1, sigma_ks);
funcKL->FixParameter(2, sigma_ks);
funcKL->FixParameter(3, 0.);
funcKL->FixParameter(4, 0.);

//make sideband fit
grKL.Fit(funcKL, "q0");
grKL.Fit(funcKL, "q0");
grKL.Fit(funcKL, "q0");

//Set sideband parameters to full fit
funcKL->FixParameter(5, funcKL->GetParameter(5));
funcKL->FixParameter(6, funcKL->GetParameter(6));
funcKL->FixParameter(7, funcKL->GetParameter(7));
funcKL->FixParameter(8, funcKL->GetParameter(8));
funcKL->FixParameter(9, funcKL->GetParameter(9));

//release the signal parameters
funcKL->ReleaseParameter(0);
funcKL->SetParLimits(0, mass_ks-Dmass_ks, mass_ks+Dmass_ks);
funcKL->ReleaseParameter(1);
funcKL->ReleaseParameter(2);
funcKL->ReleaseParameter(3);
funcKL->ReleaseParameter(4);
funcKL->SetParLimits(1, 0, sigma_ks);
funcKL->SetParLimits(2, 0, sigma_ks);
funcKL->SetParLimits(3, 0, 0.8*massksL->GetBinContent(massksL->FindBin(mass_ks)));
funcKL->SetParLimits(4, 0, 0.8*massksL->GetBinContent(massksL->FindBin(mass_ks)));

//make a full fit
massksL->Fit(funcKL, "q0", "", full_range_ks.first, full_range_ks.second);
massksL->Fit(funcKL, "q0", "", full_range_ks.first, full_range_ks.second);
massksL->Fit(funcKL, "q0", "", full_range_ks.first, full_range_ks.second);

//now release all the parameters
funcKL->ReleaseParameter(0);
funcKL->ReleaseParameter(1);
funcKL->ReleaseParameter(2);
funcKL->ReleaseParameter(3);
funcKL->ReleaseParameter(4);
funcKL->ReleaseParameter(5);
funcKL->ReleaseParameter(6);
funcKL->ReleaseParameter(7);
funcKL->ReleaseParameter(8);
funcKL->ReleaseParameter(9);
funcKL->SetParLimits(0, mass_ks-Dmass_ks, mass_ks+Dmass_ks);
funcKL->SetParLimits(1, 0, sigma_ks);
funcKL->SetParLimits(2, 0, sigma_ks);
funcKL->SetParLimits(3, 0, 0.8*massksL->GetBinContent(massksL->FindBin(mass_ks)));
funcKL->SetParLimits(4, 0, 0.8*massksL->GetBinContent(massksL->FindBin(mass_ks)));

//make final fits
for ( int i = 0; i < iter_ks; i++ ) {
        massksL->Fit(funcKL, "q0", "", full_range_ks.first, full_range_ks.second);
}

//make signal and bkg functions for plots
TF1 * func_signal_KL   = 0;
TF1 * func_bckgnd_KL   = 0;
TF1 * gaus1_signal_ks_KL   = 0;
TF1 * gaus2_signal_ks_KL   = 0;
func_signal_KL  = new TF1("func_signal_KL" , double_gaussian.c_str(), full_range_ks.first, full_range_ks.second);
func_bckgnd_KL = new TF1("func_bckgnd_KL" , poly_bkg0.c_str(), full_range_ks.first, full_range_ks.second);
gaus1_signal_ks_KL = new TF1("gaus1_signal_ks_KL", gaus_func.c_str(), full_range_ks.first, full_range_ks.second);
gaus2_signal_ks_KL = new TF1("gaus2_signal_ks_KL", gaus_func.c_str(), full_range_ks.first, full_range_ks.second);

//set parameters for signal and bkg
func_signal_KL->SetParameter(0, funcKL->GetParameter(0));
func_signal_KL->SetParameter(1, funcKL->GetParameter(1));
func_signal_KL->SetParameter(2, funcKL->GetParameter(2));
func_signal_KL->SetParameter(3, funcKL->GetParameter(3));
func_signal_KL->SetParameter(4, funcKL->GetParameter(4));
func_bckgnd_KL->SetParameter(0, funcKL->GetParameter(5));
func_bckgnd_KL->SetParameter(1, funcKL->GetParameter(6));
func_bckgnd_KL->SetParameter(2, funcKL->GetParameter(7));
func_bckgnd_KL->SetParameter(3, funcKL->GetParameter(8));
func_bckgnd_KL->SetParameter(4, funcKL->GetParameter(9));

//set parameters for gaussians
gaus1_signal_ks_KL->SetParameter(0, funcKL->GetParameter(0));
gaus1_signal_ks_KL->SetParameter(1, funcKL->GetParameter(1));
gaus1_signal_ks_KL->SetParameter(2, funcKL->GetParameter(3));
gaus2_signal_ks_KL->SetParameter(0, funcKL->GetParameter(0));
gaus2_signal_ks_KL->SetParameter(1, funcKL->GetParameter(2));
gaus2_signal_ks_KL->SetParameter(2, funcKL->GetParameter(4));

//make pull and ratio
double vny_KL[400];
double xvn_KL[400];
double pullvn_KL[400];
double ratiovn_KL[400];
double err_pullvn_KL[400];
double err_ratiovn_KL[400];
for(int k=0;k<massksL->GetNbinsX();k++)
{
    vny_KL[k] = massksL->GetBinContent(k+1);
    xvn_KL[k] = massksL->GetBinCenter(k+1);
    if(massksL->GetBinContent(k+1)!=0){
    pullvn_KL[k] = (vny_KL[k] - funcKL->Eval(xvn_KL[k]))/massksL->GetBinError(k+1);
    ratiovn_KL[k] = vny_KL[k]/funcKL->Eval(xvn[k]);
    err_pullvn_KL[k] = massksL->GetBinError(k+1);
    err_ratiovn_KL[k] = (vny[k]/funcKL->Eval(xvn[k]))*(massksL->GetBinError(k+1)/vny[k]);
    }else{
    pullvn_KL[k] = 1.0;
    ratiovn_KL[k] = 1.0;
    err_pullvn_KL[k]=0.0;
    err_ratiovn_KL[k]=0.0;
    }
}

//make the plots
TH1D* ksmass_KL = (TH1D*)massksL->Clone();
//ksmass_KL->Sumw2();
TGraphErrors* TGpullx_KL = new TGraphErrors(400,xvn_KL,pullvn_KL,0,0);
TGraphErrors* TGratiox_KL= new TGraphErrors(400,xvn_KL,ratiovn_KL,0,err_ratiovn_KL);
TCanvas* c1_KL = new TCanvas("c1_KL","c1_KL",600,600);
c1_KL->cd();
TPad* pad1x_KL = new TPad("pad1x_KL","top pad",0.0,0.4,1.0,1.0);
TPad* pad2x_KL = new TPad("pad3x_KL","middle pad",0.0,0.0,1.0,0.25);
TPad* pad3x_KL = new TPad("pad2x_KL","bottom pad",0.0,0.25,1.0,0.4);
pad1x_KL->SetTopMargin(0.1);
pad1x_KL->SetBottomMargin(0.0);
pad1x_KL->SetRightMargin(0.038);
pad1x_KL->SetLeftMargin(0.15);
pad1x_KL->Draw();
pad2x_KL->SetFrameFillStyle(4000);
pad2x_KL->SetTopMargin(0.0);
pad2x_KL->SetBottomMargin(0.37);
pad2x_KL->SetRightMargin(0.038);
pad2x_KL->SetLeftMargin(0.15);
pad2x_KL->Draw();
pad3x_KL->SetFrameFillStyle(4000);
pad3x_KL->SetTopMargin(0.0);
pad3x_KL->SetBottomMargin(0.0);
pad3x_KL->SetRightMargin(0.038);
pad3x_KL->SetLeftMargin(0.15);
pad3x_KL->Draw();
pad1x_KL->cd();
gPad->SetTickx();
gPad->SetTicky();
ksmass_KL->SetTitle("");
ksmass_KL->Draw("E");
TGaxis::SetMaxDigits(3);
ksmass_KL->GetYaxis()->SetRangeUser(0.1,1.05*ksmass_KL->GetMaximum());
ksmass_KL->GetXaxis()->SetRangeUser(full_range_ks.first,full_range_ks.second);
ksmass_KL->SetMarkerStyle(24);
ksmass_KL->SetMarkerSize(1.5);
ksmass_KL->GetYaxis()->SetTitle("#scale[1.7]{Candidates / 0.5 MeV}");
ksmass_KL->GetYaxis()->CenterTitle(1);
ksmass_KL->GetYaxis()->SetNdivisions(810);
ksmass_KL->GetXaxis()->SetTickSize(0.03);
ksmass_KL->GetYaxis()->SetTickSize(0.02);
ksmass_KL->GetXaxis()->SetLabelSize(0.055);
ksmass_KL->GetYaxis()->SetLabelSize(0.058);
ksmass_KL->GetYaxis()->SetLabelOffset(0.01);
ksmass_KL->GetYaxis()->SetTitleOffset(1.8);
funcKL->Draw("LSAME");
funcKL->SetLineColor(SPcolors[SPdred]);
func_bckgnd_KL->Draw("LSAME");
func_bckgnd_KL->SetLineColor(SPcolors[SPdred]);
func_bckgnd_KL->SetLineStyle(2);
double meanks_KL = funcKL->GetParameter(0);
double sigma1ks_KL = funcKL->GetParameter(1);
double sigma2ks_KL = funcKL->GetParameter(2);
double gaus1_yield_ks_KL = gaus1_signal_ks_KL->Integral(meanks_KL - 2*sigma1ks_KL, meanks_KL + 2*sigma1ks_KL)/massksL->GetBinWidth(1);
double gaus2_yield_ks_KL = gaus2_signal_ks_KL->Integral(meanks_KL - 2*sigma2ks_KL, meanks_KL + 2*sigma2ks_KL)/massksL->GetBinWidth(1);
double gausTot_yield_ks_KL = gaus1_yield_ks_KL + gaus2_yield_ks_KL;
double rms_gaus1_sig_ks_KL = gaus1_yield_ks_KL/gausTot_yield_ks_KL;
double rms_gaus2_sig_ks_KL = gaus2_yield_ks_KL/gausTot_yield_ks_KL;
double rms_true_ks_KL = TMath::Sqrt(rms_gaus1_sig_ks_KL*sigma1ks_KL*sigma1ks_KL + rms_gaus2_sig_ks_KL*sigma2ks_KL*sigma2ks_KL);
std::pair<double, double>  signal_range_ks_KL(meanks_KL - 2*rms_true_ks_KL, meanks_KL + 2*rms_true_ks_KL);
std::pair<double, double>  side_range_ks_KL(meanks_KL - 3*rms_true_ks_KL, meanks_KL + 3*rms_true_ks_KL);
double S_ks_KL = func_signal_KL->Integral(signal_range_ks.first, signal_range_ks.second)/massksL->GetBinWidth(1);
double B_ks_KL = func_bckgnd_KL->Integral(signal_range_ks.first, signal_range_ks.second)/massksL->GetBinWidth(1);
double Fsig_ks_KL = S_ks_KL/(S_ks_KL+B_ks_KL);
double sig_ks_KL = S_ks_KL/sqrt(S_ks_KL+B_ks_KL);
double chi2ks_KL = funcKL->GetChisquare()/funcKL->GetNDF();
double probks_KL = funcKL->GetProb();//TMath::Prob(funcKL->GetChisquare(),funcKL->GetNDF());
std::ostringstream osKL; 
TLatex* texCMS_KL = new TLatex;
texCMS_KL->SetNDC();
texCMS_KL->SetTextFont(42);
texCMS_KL->SetTextSize(0.05);
texCMS_KL->SetTextAlign(12);
texCMS_KL->DrawLatex(0.25,.95," #font[61]{CMS} #it{Work in Progress}");
texCMS_KL->DrawLatex(0.70,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS_KL->DrawLatex(0.26,0.4, "#scale[1.7]{K^{0}_{S}}");
texCMS_KL->DrawLatex(0.225,0.2, "#scale[1.0]{Used in K^{0}_{S}#Lambda}"); 
texCMS_KL->DrawLatex(0.23,0.15, "#scale[1.0]{correlations}"); 
// texCMS_KL->DrawLatex(0.75,0.15, "This Thesis"); 
osKL << "#scale[1.0]{s/#sqrt{s+b}: " << std::setprecision(4) << sig_ks_KL << "}" << std::setprecision(2);           
texCMS_KL->DrawLatex(0.68,0.60, osKL.str(  ).c_str(  )); 
osKL.str(std::string());
osKL << "#scale[1.0]{s/(s+b): " << std::setprecision(4) << Fsig_ks_KL*100 << " %}" << std::setprecision(2); 
texCMS_KL->DrawLatex(0.68,0.52, osKL.str(  ).c_str(  ));  
osKL.str(std::string()); 
osKL << "#scale[1.0]{Mean: " << std::setprecision(4) << meanks_KL << " GeV}" << std::setprecision(6);   
texCMS_KL->DrawLatex(0.68,0.46, osKL.str(  ).c_str(  ));
osKL.str(std::string()); 
osKL << "#scale[1.0]{#sigma: " << std::setprecision(4) << rms_true_ks_KL << " }" << std::setprecision(2); 
texCMS_KL->DrawLatex(0.68,0.40, osKL.str(  ).c_str(  )); 
osKL.str(std::string()); 
osKL << "#scale[1.0]{Yield: " << std::setprecision(4) << S_ks_KL << " }" << std::setprecision(2); 
texCMS_KL->DrawLatex(0.68,0.33, osKL.str(  ).c_str(  )); 
osKL.str(std::string()); 
osKL << "#scale[1.0]{#chi^{2}/ndf: " << std::setprecision(4) << chi2ks_KL << " }" << std::setprecision(4); 
texCMS_KL->DrawLatex(0.68,0.26, osKL.str(  ).c_str(  ));  
osKL.str(std::string()); 
// osKL << "#scale[1.0]{Prob: " << std::setprecision(4) << probks_KL << "}" << std::setprecision(2); 
// texCMS_KL->DrawLatex(0.68,0.20, osKL.str(  ).c_str(  )); 
// osKL.str(std::string()); 
osKL << "#scale[1.]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_KL->DrawLatex(0.68,0.75, osKL.str(  ).c_str(  )); 
osKL.str(std::string()); 
osKL << "#scale[1.]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_KL->DrawLatex(0.68,0.67, osKL.str(  ).c_str(  ));  
osKL.str(std::string()); 
osKL << "#scale[1.]{ " << std::setprecision(4) << qminM << std::setprecision(4) <<" < q_{inv} < "  << qmaxM << " GeV}" << std::setprecision(2); 
texCMS_KL->DrawLatex(0.68,0.83, osKL.str(  ).c_str(  ));  
osKL.str(std::string()); 
Float_t yyymaxKL = ksmass_KL->GetMaximum();
TLine *linexxKL = new TLine(signal_range_ks_KL.first,0.1,signal_range_ks_KL.first,yyymaxKL);
linexxKL->SetLineColor(SPcolors[SPdblue]);
linexxKL->Draw();
TLine *line2xxKL = new TLine(signal_range_ks_KL.second,0.1,signal_range_ks_KL.second,yyymaxKL);
line2xxKL->SetLineColor(SPcolors[SPdblue]);
line2xxKL->Draw();
TLegend* legKL = new TLegend(0.16,0.6,0.43,0.85);
legKL->AddEntry(ksmass_KL,"Data");
legKL->AddEntry(funcKL,"Signal+Background","L");
legKL->AddEntry(func_bckgnd_KL,"Background","L");
legKL->AddEntry(linexx,"Peak Region","L");
legKL->SetBorderSize(0);
legKL->Draw();
c1_KL->Update();
pad2x_KL->cd();
gPad->SetTickx();
gPad->SetTicky();
double lowxKL = -49.1;
double highxKL = 49.1;
TH1F* frame_ks_KL = pad2x_KL->cd()->DrawFrame(full_range_ks.first,lowxKL,full_range_ks.second,highxKL);
frame_ks_KL->GetXaxis()->SetTitle("#scale[4.3]{#pi^{+}#pi^{-} invariant mass [GeV]}");
frame_ks_KL->GetYaxis()->SetTitle("#scale[4.0]{Pull}");
frame_ks_KL->GetXaxis()->CenterTitle(1);
frame_ks_KL->GetYaxis()->CenterTitle(1);
frame_ks_KL->GetXaxis()->SetNdivisions(510);
frame_ks_KL->GetYaxis()->SetNdivisions(406);
frame_ks_KL->GetXaxis()->SetTickSize(0.06);
frame_ks_KL->GetYaxis()->SetTickSize(0.02);
frame_ks_KL->GetXaxis()->SetTitleOffset(4.5);
frame_ks_KL->GetYaxis()->SetTitleOffset(1.8);
frame_ks_KL->GetXaxis()->SetLabelSize(0.14);
frame_ks_KL->GetYaxis()->SetLabelSize(0.14);
frame_ks_KL->GetXaxis()->SetLabelOffset(0.01);
frame_ks_KL->GetYaxis()->SetLabelOffset(0.013);
TGpullx_KL->SetMarkerStyle(25);
TGpullx_KL->SetMarkerSize(1.3);
TGpullx_KL->SetMarkerColor(kBlack);
TGpullx_KL->Draw("P");
TLine* lineyKL = new TLine(full_range_ks.first, 0, full_range_ks.second, 0);
lineyKL->SetLineStyle(2);
lineyKL->Draw("same");
pad3x_KL->cd();
gPad->SetTickx();
gPad->SetTicky();
double low2xKL = 0.1;
double high2xKL = 1.899999;
frame_ks_KL = pad3x_KL->cd()->DrawFrame(full_range_ks.first,low2xKL,full_range_ks.second,high2xKL);
frame_ks_KL->GetYaxis()->SetTitle("#scale[6.3]{#frac{Data}{Fit}}");
frame_ks_KL->GetXaxis()->CenterTitle(1);
frame_ks_KL->GetYaxis()->CenterTitle(1);
frame_ks_KL->GetXaxis()->SetNdivisions(510);
frame_ks_KL->GetYaxis()->SetNdivisions(406);
frame_ks_KL->GetXaxis()->SetTickSize(0.1);
frame_ks_KL->GetYaxis()->SetTickSize(0.01);
frame_ks_KL->GetYaxis()->SetTitleOffset(1.8);
frame_ks_KL->GetXaxis()->SetLabelSize(0.12);
frame_ks_KL->GetYaxis()->SetLabelSize(0.22);
frame_ks_KL->GetYaxis()->SetLabelOffset(0.01);
TGratiox_KL->SetMarkerStyle(25);
TGratiox_KL->SetMarkerSize(1.3);
TGratiox_KL->SetMarkerColor(kBlack);
TGratiox_KL->Draw("P");
TLine* line2yKL = new TLine(full_range_ks.first, 1, full_range_ks.second, 1);
line2yKL->SetLineStyle(2);
line2yKL->Draw("same");
c1_KL->Print("results/Mass"+V0+"KSL_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");
c1_KL->Update();

//K0s used in K0sALam

cout << "==================K0s Mass for  K0sAL==================" << endl;
TH1D* massksAL = (TH1D*)h_KAL_Mass_1D->Clone();
//massksAL->Sumw2();
TF1 *funcKAL = new TF1("funcKAL", massfunc.c_str(), full_range_ks.first, full_range_ks.second);

//sideband fit
int cntKAL = 0;
for ( int i = 0; i < massksAL->GetNbinsX(); i++ ) {
	double center = massksAL->GetBinCenter(i+1);
	if ( (center > sb1_ks.first and center < sb1_ks.second) or (center > sb2_ks.first and center < sb2_ks.second) ) {
	cntKAL++;
	}
}
TGraphErrors grKAL(cntKAL);
cntKAL = 0;
for ( int i = 0; i < massksAL->GetNbinsX(); i++ ) {
	double center = massksAL->GetBinCenter(i+1);
	if ( (center > sb1_ks.first and center < sb1_ks.second) or (center > sb2_ks.first and center < sb2_ks.second) ) {
		grKAL.GetX()[cntKAL] = center;
		grKAL.GetY()[cntKAL] = massksAL->GetBinContent(i+1);
		grKAL.GetEY()[cntKAL] = massksAL->GetBinError(i+1);
		cntKAL++;
	}
}

//sideband fit parameters
funcKAL->SetParameter(5, 10.);
funcKAL->SetParameter(6, 1.);
funcKAL->SetParameter(7, 1.);
funcKAL->SetParameter(8, 1.);
funcKAL->SetParameter(9, 1.);

//fix the parameters in order to get only bkg
funcKAL->FixParameter(0, mass_ks);
funcKAL->FixParameter(1, sigma_ks);
funcKAL->FixParameter(2, sigma_ks);
funcKAL->FixParameter(3, 0.);
funcKAL->FixParameter(4, 0.);

//make sideband fit
grKAL.Fit(funcKAL, "q0");
grKAL.Fit(funcKAL, "q0");
grKAL.Fit(funcKAL, "q0");

//Set sideband parameters to full fit
funcKAL->FixParameter(5, funcKAL->GetParameter(5));
funcKAL->FixParameter(6, funcKAL->GetParameter(6));
funcKAL->FixParameter(7, funcKAL->GetParameter(7));
funcKAL->FixParameter(8, funcKAL->GetParameter(8));
funcKAL->FixParameter(9, funcKAL->GetParameter(9));

//release the signal parameters
funcKAL->ReleaseParameter(0);
funcKAL->SetParLimits(0, mass_ks-Dmass_ks, mass_ks+Dmass_ks);
funcKAL->ReleaseParameter(1);
funcKAL->ReleaseParameter(2);
funcKAL->ReleaseParameter(3);
funcKAL->ReleaseParameter(4);
funcKAL->SetParLimits(1, 0, sigma_ks);
funcKAL->SetParLimits(2, 0, sigma_ks);
funcKAL->SetParLimits(3, 0, 0.8*massksAL->GetBinContent(massksAL->FindBin(mass_ks)));
funcKAL->SetParLimits(4, 0, 0.8*massksAL->GetBinContent(massksAL->FindBin(mass_ks)));

//make a full fit
massksAL->Fit(funcKAL, "q0", "", full_range_ks.first, full_range_ks.second);
massksAL->Fit(funcKAL, "q0", "", full_range_ks.first, full_range_ks.second);
massksAL->Fit(funcKAL, "q0", "", full_range_ks.first, full_range_ks.second);

//now release all the parameters
funcKAL->ReleaseParameter(0);
funcKAL->ReleaseParameter(1);
funcKAL->ReleaseParameter(2);
funcKAL->ReleaseParameter(3);
funcKAL->ReleaseParameter(4);
funcKAL->ReleaseParameter(5);
funcKAL->ReleaseParameter(6);
funcKAL->ReleaseParameter(7);
funcKAL->ReleaseParameter(8);
funcKAL->ReleaseParameter(9);
funcKAL->SetParLimits(0, mass_ks-Dmass_ks, mass_ks+Dmass_ks);
funcKAL->SetParLimits(1, 0, sigma_ks);
funcKAL->SetParLimits(2, 0, sigma_ks);
funcKAL->SetParLimits(3, 0, 0.8*massksAL->GetBinContent(massksAL->FindBin(mass_ks)));
funcKAL->SetParLimits(4, 0, 0.8*massksAL->GetBinContent(massksAL->FindBin(mass_ks)));

//make final fits
for ( int i = 0; i < iter_ks; i++ ) {
    massksAL->Fit(funcKAL, "q0", "", full_range_ks.first, full_range_ks.second);
}

//make signal and bkg functions for plots
TF1 * func_signal_KAL   = 0;
TF1 * func_bckgnd_KAL   = 0;
TF1 * gaus1_signal_ks_KAL   = 0;
TF1 * gaus2_signal_ks_KAL   = 0;
func_signal_KAL  = new TF1("func_signal_KAL" , double_gaussian.c_str(), full_range_ks.first, full_range_ks.second);
func_bckgnd_KAL = new TF1("func_bckgnd_KAL" , poly_bkg0.c_str(), full_range_ks.first, full_range_ks.second);
gaus1_signal_ks_KAL = new TF1("gaus1_signal_ks_KAL", gaus_func.c_str(), full_range_ks.first, full_range_ks.second);
gaus2_signal_ks_KAL = new TF1("gaus2_signal_ks_KAL", gaus_func.c_str(), full_range_ks.first, full_range_ks.second);

//set parameters for signal and bkg
func_signal_KAL->SetParameter(0, funcKAL->GetParameter(0));
func_signal_KAL->SetParameter(1, funcKAL->GetParameter(1));
func_signal_KAL->SetParameter(2, funcKAL->GetParameter(2));
func_signal_KAL->SetParameter(3, funcKAL->GetParameter(3));
func_signal_KAL->SetParameter(4, funcKAL->GetParameter(4));
func_bckgnd_KAL->SetParameter(0, funcKAL->GetParameter(5));
func_bckgnd_KAL->SetParameter(1, funcKAL->GetParameter(6));
func_bckgnd_KAL->SetParameter(2, funcKAL->GetParameter(7));
func_bckgnd_KAL->SetParameter(3, funcKAL->GetParameter(8));
func_bckgnd_KAL->SetParameter(4, funcKAL->GetParameter(9));

//set parameters for gaussians
gaus1_signal_ks_KAL->SetParameter(0, funcKAL->GetParameter(0));
gaus1_signal_ks_KAL->SetParameter(1, funcKAL->GetParameter(1));
gaus1_signal_ks_KAL->SetParameter(2, funcKAL->GetParameter(3));
gaus2_signal_ks_KAL->SetParameter(0, funcKAL->GetParameter(0));
gaus2_signal_ks_KAL->SetParameter(1, funcKAL->GetParameter(2));
gaus2_signal_ks_KAL->SetParameter(2, funcKAL->GetParameter(4));

//make pull and ratio
double vny_KAL[400];
double xvn_KAL[400];
double pullvn_KAL[400];
double ratiovn_KAL[400];
double err_pullvn_KAL[400];
double err_ratiovn_KAL[400];
for(int k=0;k<massksAL->GetNbinsX();k++)
        {
            vny_KAL[k] = massksAL->GetBinContent(k+1);
            xvn_KAL[k] = massksAL->GetBinCenter(k+1);
        	if(massksAL->GetBinContent(k+1)!=0){
            pullvn_KAL[k] = (vny_KAL[k] - funcKAL->Eval(xvn_KAL[k]))/massksAL->GetBinError(k+1);
            ratiovn_KAL[k] = vny_KAL[k]/funcKAL->Eval(xvn[k]);
            err_pullvn_KAL[k] = massksAL->GetBinError(k+1);
            err_ratiovn_KAL[k] = (vny[k]/funcKAL->Eval(xvn[k]))*(massksAL->GetBinError(k+1)/vny[k]);
            }else{
            pullvn_KAL[k] = 1.0;
            ratiovn_KAL[k] = 1.0;
            err_pullvn_KAL[k]=0.0;
            err_ratiovn_KAL[k]=0.0;
            }
}

//make the plots
TH1D* ksmass_KAL = (TH1D*)massksAL->Clone();
//ksmass_KAL->Sumw2();
TGraphErrors* TGpullx_KAL = new TGraphErrors(400,xvn_KAL,pullvn_KAL,0,0);
TGraphErrors* TGratiox_KAL= new TGraphErrors(400,xvn_KAL,ratiovn_KAL,0,err_ratiovn_KAL);
TCanvas* c1_KAL = new TCanvas("c1_KAL","c1_KAL",600,600);
c1_KAL->cd();
TPad* pad1x_KAL = new TPad("pad1x_KAL","top pad",0.0,0.4,1.0,1.0);
TPad* pad2x_KAL = new TPad("pad3x_KAL","middle pad",0.0,0.0,1.0,0.25);
TPad* pad3x_KAL = new TPad("pad2x_KAL","bottom pad",0.0,0.25,1.0,0.4);
pad1x_KAL->SetTopMargin(0.1);
pad1x_KAL->SetBottomMargin(0.0);
pad1x_KAL->SetRightMargin(0.038);
pad1x_KAL->SetLeftMargin(0.15);
pad1x_KAL->Draw();
pad2x_KAL->SetFrameFillStyle(4000);
pad2x_KAL->SetTopMargin(0.0);
pad2x_KAL->SetBottomMargin(0.37);
pad2x_KAL->SetRightMargin(0.038);
pad2x_KAL->SetLeftMargin(0.15);
pad2x_KAL->Draw();
pad3x_KAL->SetFrameFillStyle(4000);
pad3x_KAL->SetTopMargin(0.0);
pad3x_KAL->SetBottomMargin(0.0);
pad3x_KAL->SetRightMargin(0.038);
pad3x_KAL->SetLeftMargin(0.15);
pad3x_KAL->Draw();
pad1x_KAL->cd();
gPad->SetTickx();
gPad->SetTicky();
ksmass_KAL->SetTitle("");
ksmass_KAL->Draw("E");
TGaxis::SetMaxDigits(3);
ksmass_KAL->GetYaxis()->SetRangeUser(0.1,1.05*ksmass_KAL->GetMaximum());
ksmass_KAL->GetXaxis()->SetRangeUser(full_range_ks.first,full_range_ks.second);
ksmass_KAL->SetMarkerStyle(24);
ksmass_KAL->SetMarkerSize(1.5);
ksmass_KAL->GetYaxis()->SetTitle("#scale[1.7]{Candidates / 0.5 MeV}");
ksmass_KAL->GetYaxis()->CenterTitle(1);
ksmass_KAL->GetYaxis()->SetNdivisions(810);
ksmass_KAL->GetXaxis()->SetTickSize(0.03);
ksmass_KAL->GetYaxis()->SetTickSize(0.02);
ksmass_KAL->GetXaxis()->SetLabelSize(0.055);
ksmass_KAL->GetYaxis()->SetLabelSize(0.058);
ksmass_KAL->GetYaxis()->SetLabelOffset(0.01);
ksmass_KAL->GetYaxis()->SetTitleOffset(1.8);
funcKAL->Draw("LSAME");
funcKAL->SetLineColor(SPcolors[SPdred]);
func_bckgnd_KAL->Draw("LSAME");
func_bckgnd_KAL->SetLineColor(SPcolors[SPdred]);
func_bckgnd_KAL->SetLineStyle(2);
double meanks_KAL = funcKAL->GetParameter(0);
double sigma1ks_KAL = funcKAL->GetParameter(1);
double sigma2ks_KAL = funcKAL->GetParameter(2);
double gaus1_yield_ks_KAL = gaus1_signal_ks_KAL->Integral(meanks_KAL - 2*sigma1ks_KAL, meanks_KAL + 2*sigma1ks_KAL)/massksAL->GetBinWidth(1);
double gaus2_yield_ks_KAL = gaus2_signal_ks_KAL->Integral(meanks_KAL - 2*sigma2ks_KAL, meanks_KAL + 2*sigma2ks_KAL)/massksAL->GetBinWidth(1);
double gausTot_yield_ks_KAL = gaus1_yield_ks_KAL + gaus2_yield_ks_KAL;
double rms_gaus1_sig_ks_KAL = gaus1_yield_ks_KAL/gausTot_yield_ks_KAL;
double rms_gaus2_sig_ks_KAL = gaus2_yield_ks_KAL/gausTot_yield_ks_KAL;
double rms_true_ks_KAL = TMath::Sqrt(rms_gaus1_sig_ks_KAL*sigma1ks_KAL*sigma1ks_KAL + rms_gaus2_sig_ks_KAL*sigma2ks_KAL*sigma2ks_KAL);
std::pair<double, double>  signal_range_ks_KAL(meanks_KAL - 2*rms_true_ks_KAL, meanks_KAL + 2*rms_true_ks_KAL);
std::pair<double, double>  side_range_ks_KAL(meanks_KAL - 3*rms_true_ks_KAL, meanks_KAL + 3*rms_true_ks_KAL);
double S_ks_KAL = func_signal_KAL->Integral(signal_range_ks.first, signal_range_ks.second)/massksAL->GetBinWidth(1);
double B_ks_KAL = func_bckgnd_KAL->Integral(signal_range_ks.first, signal_range_ks.second)/massksAL->GetBinWidth(1);
double Fsig_ks_KAL = S_ks_KAL/(S_ks_KAL+B_ks_KAL);
double sig_ks_KAL = S_ks_KAL/sqrt(S_ks_KAL+B_ks_KAL);
double chi2ks_KAL = funcKAL->GetChisquare()/funcKAL->GetNDF();
double probks_KAL = TMath::Prob(funcKAL->GetChisquare(),funcKAL->GetNDF());
std::ostringstream osKAL; 
TLatex* texCMS_KAL = new TLatex;
texCMS_KAL->SetNDC();
texCMS_KAL->SetTextFont(42);
texCMS_KAL->SetTextSize(0.05);
texCMS_KAL->SetTextAlign(12);
texCMS_KAL->DrawLatex(0.25,.95," #font[61]{CMS} #it{Work in Progress}");
texCMS_KAL->DrawLatex(0.70,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS_KAL->DrawLatex(0.26,0.4, "#scale[1.7]{K^{0}_{S}}");
texCMS_KAL->DrawLatex(0.225,0.2, "#scale[1.0]{Used in K^{0}_{S}#bar{#Lambda}}"); 
texCMS_KAL->DrawLatex(0.23,0.15, "#scale[1.0]{correlations}"); 
// texCMS_KAL->DrawLatex(0.75,0.15, "This Thesis"); 
osKAL << "#scale[1.0]{s/#sqrt{s+b}: " << std::setprecision(4) << sig_ks_KAL << "}" << std::setprecision(2);           
texCMS_KAL->DrawLatex(0.68,0.60, osKAL.str(  ).c_str(  )); 
osKAL.str(std::string());
osKAL << "#scale[1.0]{s/(s+b): " << std::setprecision(4) << Fsig_ks_KAL*100 << " %}" << std::setprecision(2); 
texCMS_KAL->DrawLatex(0.68,0.52, osKAL.str(  ).c_str(  ));  
osKAL.str(std::string()); 
osKAL << "#scale[1.0]{Mean: " << std::setprecision(4) << meanks_KAL << " GeV}" << std::setprecision(6);   
texCMS_KAL->DrawLatex(0.68,0.46, osKAL.str(  ).c_str(  ));
osKAL.str(std::string()); 
osKAL << "#scale[1.0]{#sigma: " << std::setprecision(4) << rms_true_ks_KAL << " }" << std::setprecision(2); 
texCMS_KAL->DrawLatex(0.68,0.40, osKAL.str(  ).c_str(  )); 
osKAL.str(std::string()); 
osKAL << "#scale[1.0]{Yield: " << std::setprecision(4) << S_ks_KAL << " }" << std::setprecision(2); 
texCMS_KAL->DrawLatex(0.68,0.33, osKAL.str(  ).c_str(  )); 
osKAL.str(std::string()); 
osKAL << "#scale[1.0]{#chi^{2}/ndf: " << std::setprecision(4) << chi2ks_KAL << " }" << std::setprecision(4); 
texCMS_KAL->DrawLatex(0.68,0.26, osKAL.str(  ).c_str(  ));  
osKAL.str(std::string()); 
// osKAL << "#scale[1.0]{Prob: " << std::setprecision(4) << probks_KAL << "}" << std::setprecision(2); 
// texCMS_KAL->DrawLatex(0.68,0.20, osKAL.str(  ).c_str(  )); 
// osKAL.str(std::string()); 
osKAL << "#scale[1.]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_KAL->DrawLatex(0.68,0.75, osKAL.str(  ).c_str(  )); 
osKAL.str(std::string()); 
osKAL << "#scale[1.]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_KAL->DrawLatex(0.68,0.67, osKAL.str(  ).c_str(  ));  
osKAL.str(std::string()); 
osKAL << "#scale[1.]{ " << std::setprecision(4) << qminM << std::setprecision(4) <<" < q_{inv} < "  << qmaxM << " GeV}" << std::setprecision(2); 
texCMS_KAL->DrawLatex(0.68,0.83, osKAL.str(  ).c_str(  ));  
osKAL.str(std::string()); 
Float_t yyymaxKAL = ksmass_KAL->GetMaximum();
TLine *linexxKAL = new TLine(signal_range_ks_KAL.first,0.1,signal_range_ks_KAL.first,yyymaxKAL);
linexxKAL->SetLineColor(SPcolors[SPdblue]);
linexxKAL->Draw();
TLine *line2xxKAL = new TLine(signal_range_ks_KAL.second,0.1,signal_range_ks_KAL.second,yyymaxKAL);
line2xxKAL->SetLineColor(SPcolors[SPdblue]);
line2xxKAL->Draw();
TLegend* legKAL = new TLegend(0.16,0.6,0.43,0.85);
legKAL->AddEntry(ksmass_KAL,"Data");
legKAL->AddEntry(funcKAL,"Signal+Background","L");
legKAL->AddEntry(func_bckgnd_KAL,"Background","L");
legKAL->AddEntry(linexx,"Peak Region","L");
legKAL->SetBorderSize(0);
legKAL->Draw();
c1_KAL->Update();
pad2x_KAL->cd();
gPad->SetTickx();
gPad->SetTicky();
double lowxKAL = -49.1;
double highxKAL = 49.1;
TH1F* frame_ks_KAL = pad2x_KAL->cd()->DrawFrame(full_range_ks.first,lowxKAL,full_range_ks.second,highxKAL);
frame_ks_KAL->GetXaxis()->SetTitle("#scale[4.3]{#pi^{+}#pi^{-} invariant mass [GeV]}");
frame_ks_KAL->GetYaxis()->SetTitle("#scale[4.0]{Pull}");
frame_ks_KAL->GetXaxis()->CenterTitle(1);
frame_ks_KAL->GetYaxis()->CenterTitle(1);
frame_ks_KAL->GetXaxis()->SetNdivisions(510);
frame_ks_KAL->GetYaxis()->SetNdivisions(406);
frame_ks_KAL->GetXaxis()->SetTickSize(0.06);
frame_ks_KAL->GetYaxis()->SetTickSize(0.02);
frame_ks_KAL->GetXaxis()->SetTitleOffset(4.5);
frame_ks_KAL->GetYaxis()->SetTitleOffset(1.8);
frame_ks_KAL->GetXaxis()->SetLabelSize(0.14);
frame_ks_KAL->GetYaxis()->SetLabelSize(0.14);
frame_ks_KAL->GetXaxis()->SetLabelOffset(0.01);
frame_ks_KAL->GetYaxis()->SetLabelOffset(0.013);
TGpullx_KAL->SetMarkerStyle(25);
TGpullx_KAL->SetMarkerSize(1.3);
TGpullx_KAL->SetMarkerColor(kBlack);
TGpullx_KAL->Draw("P");
TLine* lineyKAL = new TLine(full_range_ks.first, 0, full_range_ks.second, 0);
lineyKAL->SetLineStyle(2);
lineyKAL->Draw("same");
pad3x_KAL->cd();
gPad->SetTickx();
gPad->SetTicky();
double low2xKAL = 0.1;
double high2xKAL = 1.899999;
frame_ks_KAL = pad3x_KAL->cd()->DrawFrame(full_range_ks.first,low2xKAL,full_range_ks.second,high2xKAL);
frame_ks_KAL->GetYaxis()->SetTitle("#scale[6.3]{#frac{Data}{Fit}}");
frame_ks_KAL->GetXaxis()->CenterTitle(1);
frame_ks_KAL->GetYaxis()->CenterTitle(1);
frame_ks_KAL->GetXaxis()->SetNdivisions(510);
frame_ks_KAL->GetYaxis()->SetNdivisions(406);
frame_ks_KAL->GetXaxis()->SetTickSize(0.1);
frame_ks_KAL->GetYaxis()->SetTickSize(0.01);
frame_ks_KAL->GetYaxis()->SetTitleOffset(1.8);
frame_ks_KAL->GetXaxis()->SetLabelSize(0.12);
frame_ks_KAL->GetYaxis()->SetLabelSize(0.22);
frame_ks_KAL->GetYaxis()->SetLabelOffset(0.01);
TGratiox_KAL->SetMarkerStyle(25);
TGratiox_KAL->SetMarkerSize(1.3);
TGratiox_KAL->SetMarkerColor(kBlack);
TGratiox_KAL->Draw("P");
TLine* line2yKAL = new TLine(full_range_ks.first, 1, full_range_ks.second, 1);
line2yKAL->SetLineStyle(2);
line2yKAL->Draw("same");
c1_KAL->Print("results/Mass"+V0+"KSAL_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");
c1_KAL->Update();

//Lambda stuff
int iter_la = 5;
Int_t Nbins_la = 270;
std::pair<double, double>  full_range_la(1.095, 1.14);
std::pair<double, double>  sb1_la(1.095, 1.106);
std::pair<double, double>  sb2_la(1.125, 1.14);
double mass_la = 1.115683;
double sigma_la = 0.005;
double Dmass_la = 0.003;

//Lam used in LamLam
cout << "==================Lam Mass for LL==================" << endl;
TH1D* masslam = (TH1D*)h_Lam_Mass_1D->Clone();
//masslam->Sumw2();
TF1 *func2 = new TF1("func2", massfunc.c_str(), full_range_la.first, full_range_la.second);

//sideband fit
int cnt2 = 0;
for ( int i = 0; i < masslam->GetNbinsX(); i++ ) {
    double center = masslam->GetBinCenter(i+1);
    if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
       cnt2++;
    }
}
TGraphErrors gr2(cnt2);
cnt2 = 0;
for ( int i = 0; i < masslam->GetNbinsX(); i++ ) {
   double center = masslam->GetBinCenter(i+1);
   if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
      gr2.GetX()[cnt2] = center;
      gr2.GetY()[cnt2] = masslam->GetBinContent(i+1);
      gr2.GetEY()[cnt2] = masslam->GetBinError(i+1);
      cnt2++;
    }
}

//sideband fit parameters
func2->SetParameter(5, 10.);
func2->SetParameter(6, 1.);
func2->SetParameter(7, 1.);
func2->SetParameter(8, 1.);
func2->SetParameter(9, 1.);

//fix the parameters in order to get only bkg
func2->FixParameter(0, mass_la);
func2->FixParameter(1, sigma_la);
func2->FixParameter(2, sigma_la);
func2->FixParameter(3, 0.);
func2->FixParameter(4, 0.);

//make sideband fit
gr2.Fit(func2, "q0");
gr2.Fit(func2, "q0");
gr2.Fit(func2, "q0");

//Set sideband parameters to full fit
func2->FixParameter(5, func2->GetParameter(5));
func2->FixParameter(6, func2->GetParameter(6));
func2->FixParameter(7, func2->GetParameter(7));
func2->FixParameter(8, func2->GetParameter(8));
func2->FixParameter(9, func2->GetParameter(9));

//release the signal parameters
func2->ReleaseParameter(0);
func2->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2->ReleaseParameter(1);
func2->ReleaseParameter(2);
func2->ReleaseParameter(3);
func2->ReleaseParameter(4);
func2->SetParLimits(1, 0, sigma_la);
func2->SetParLimits(2, 0, sigma_la);
func2->SetParLimits(3, 0, 0.8*masslam->GetBinContent(masslam->FindBin(mass_la)));
func2->SetParLimits(4, 0, 0.8*masslam->GetBinContent(masslam->FindBin(mass_la)));

//make a full fit
masslam->Fit(func2, "q0", "", full_range_la.first, full_range_la.second);
masslam->Fit(func2, "q0", "", full_range_la.first, full_range_la.second);
masslam->Fit(func2, "q0", "", full_range_la.first, full_range_la.second);

//now release all the parameters
func2->ReleaseParameter(0);
func2->ReleaseParameter(1);
func2->ReleaseParameter(2);
func2->ReleaseParameter(3);
func2->ReleaseParameter(4);
func2->ReleaseParameter(5);
func2->ReleaseParameter(6);
func2->ReleaseParameter(7);
func2->ReleaseParameter(8);
func2->ReleaseParameter(9);
func2->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2->SetParLimits(1, 0, sigma_la);
func2->SetParLimits(2, 0, sigma_la);
func2->SetParLimits(3, 0, 0.8*masslam->GetBinContent(masslam->FindBin(mass_la)));
func2->SetParLimits(4, 0, 0.8*masslam->GetBinContent(masslam->FindBin(mass_la)));

//make final fits
for ( int i = 0; i < iter_la-1; i++ ) {
        masslam->Fit(func2, "q0", "", full_range_la.first, full_range_la.second);
}

//make signal and bkg functions for plots
TF1 * func2_signal   = 0;
TF1 * func2_bckgnd   = 0;
TF1 * gaus1_signal_la   = 0;
TF1 * gaus2_signal_la   = 0;
func2_signal  = new TF1("func2_signal" , double_gaussian.c_str(), full_range_la.first, full_range_la.second);
func2_bckgnd = new TF1("func2_bckgnd" , poly_bkg0.c_str(), full_range_la.first, full_range_la.second);
gaus1_signal_la = new TF1("gaus1_signal_la", gaus_func.c_str(), full_range_la.first, full_range_la.second);
gaus2_signal_la = new TF1("gaus2_signal_la", gaus_func.c_str(), full_range_la.first, full_range_la.second);

//set parameters for signal and bkg
func2_signal->SetParameter(0, func2->GetParameter(0));
func2_signal->SetParameter(1, func2->GetParameter(1));
func2_signal->SetParameter(2, func2->GetParameter(2));
func2_signal->SetParameter(3, func2->GetParameter(3));
func2_signal->SetParameter(4, func2->GetParameter(4));
func2_bckgnd->SetParameter(0, func2->GetParameter(5));
func2_bckgnd->SetParameter(1, func2->GetParameter(6));
func2_bckgnd->SetParameter(2, func2->GetParameter(7));
func2_bckgnd->SetParameter(3, func2->GetParameter(8));
func2_bckgnd->SetParameter(4, func2->GetParameter(9));

//set parameters for gaussians
gaus1_signal_la->SetParameter(0, func2->GetParameter(0));
gaus1_signal_la->SetParameter(1, func2->GetParameter(1));
gaus1_signal_la->SetParameter(2, func2->GetParameter(3));
gaus2_signal_la->SetParameter(0, func2->GetParameter(0));
gaus2_signal_la->SetParameter(1, func2->GetParameter(2));
gaus2_signal_la->SetParameter(2, func2->GetParameter(4));

//make pull and ratio
double vny2[400];
double xvn2[400];
double pullvn2[400];
double ratiovn2[400];
double err_pullvn2[400];
double err_ratiovn2[400];
for(int k=0;k<masslam->GetNbinsX();k++)
        {
            vny2[k] = masslam->GetBinContent(k+1);
            xvn2[k] = masslam->GetBinCenter(k+1);
        	if(masslam->GetBinContent(k+1)!=0){
            pullvn2[k] = (vny2[k] - func2->Eval(xvn2[k]))/masslam->GetBinError(k+1);
            ratiovn2[k] = vny2[k]/func2->Eval(xvn2[k]);
            err_pullvn2[k] = masslam->GetBinError(k+1);
            err_ratiovn2[k] = (vny2[k]/func2->Eval(xvn2[k]))*(masslam->GetBinError(k+1)/vny2[k]);
            }else{
            pullvn2[k] = 0.0;
            ratiovn2[k] = 0.0;
            err_pullvn2[k]=0.0;
            err_ratiovn2[k]=0.0;
            }
}

//make the plots
TH1D* lamass = (TH1D*)masslam->Clone();
//lamass->Sumw2();
TGraphErrors* TGpull = new TGraphErrors(400,xvn2,pullvn2,0,0);
TGraphErrors* TGratio= new TGraphErrors(400,xvn2,ratiovn2,0,err_ratiovn2);
TCanvas* c2 = new TCanvas("c2","c2",600,600);
c2->cd();
TPad* pad1 = new TPad("pad1","top pad",0.0,0.4,1.0,1.0);
TPad* pad2 = new TPad("pad3","middle pad",0.0,0.0,1.0,0.25);
TPad* pad3 = new TPad("pad2","bottom pad",0.0,0.25,1.0,0.4);
pad1->SetTopMargin(0.1);
pad1->SetBottomMargin(0.0);
pad1->SetRightMargin(0.038);
pad1->SetLeftMargin(0.15);
pad1->Draw();
pad2->SetFrameFillStyle(4000);
pad2->SetTopMargin(0.0);
pad2->SetBottomMargin(0.37);
pad2->SetRightMargin(0.038);
pad2->SetLeftMargin(0.15);
pad2->Draw();
pad3->SetFrameFillStyle(4000);
pad3->SetTopMargin(0.0);
pad3->SetBottomMargin(0.0);
pad3->SetRightMargin(0.038);
pad3->SetLeftMargin(0.15);
pad3->Draw();
pad1->cd();
gPad->SetTickx();
gPad->SetTicky();
lamass->SetTitle("");
lamass->Draw("E");
TGaxis::SetMaxDigits(3);
lamass->GetXaxis()->SetRangeUser(full_range_la.first,full_range_la.second);
lamass->GetYaxis()->SetRangeUser(0.1,1.05*lamass->GetMaximum());
lamass->SetMarkerStyle(24);
lamass->SetMarkerSize(1.5);
lamass->GetYaxis()->SetTitle("#scale[1.7]{Candidates / 0.5 MeV}");
lamass->GetYaxis()->CenterTitle(1);
lamass->GetYaxis()->CenterTitle(1);
lamass->GetYaxis()->SetNdivisions(810);
lamass->GetXaxis()->SetTickSize(0.03);
lamass->GetYaxis()->SetTickSize(0.02);
lamass->GetXaxis()->SetLabelSize(0.055);
lamass->GetYaxis()->SetLabelSize(0.058);
lamass->GetYaxis()->SetLabelOffset(0.01);
lamass->GetYaxis()->SetTitleOffset(1.8);
func2->Draw("LSAME");
func2->SetLineColor(SPcolors[SPdred]);
func2_bckgnd->Draw("LSAME");
func2_bckgnd->SetLineColor(SPcolors[SPdred]);
func2_bckgnd->SetLineStyle(2);
double meanla = func2->GetParameter(0);
double sigma1la = func2->GetParameter(1);
double sigma2la = func2->GetParameter(2);
double gaus1_yield_la = gaus1_signal_la->Integral(meanla - 2*sigma1la, meanla + 2*sigma1la)/masslam->GetBinWidth(1);
double gaus2_yield_la = gaus2_signal_la->Integral(meanla - 2*sigma2la, meanla + 2*sigma2la)/masslam->GetBinWidth(1);
double gausTot_yield_la = gaus1_yield_la + gaus2_yield_la;
double rms_gaus1_sig_la = gaus1_yield_la/gausTot_yield_la;
double rms_gaus2_sig_la = gaus2_yield_la/gausTot_yield_la;
double rms_true_la = TMath::Sqrt(rms_gaus1_sig_la*sigma1la*sigma1la + rms_gaus2_sig_la*sigma2la*sigma2la);
std::pair<double, double>  signal_range_la(meanla - 2*rms_true_la, meanla + 2*rms_true_la);
double S_la = func2_signal->Integral(signal_range_la.first, signal_range_la.second)/masslam->GetBinWidth(1);
double B_la = func2_bckgnd->Integral(signal_range_la.first, signal_range_la.second)/masslam->GetBinWidth(1);
double Fsig_la = S_la/(S_la+B_la);
double sig_la = S_la/sqrt(S_la+B_la);
double chi2la = func2->GetChisquare()/func2->GetNDF();
double probla = TMath::Prob(func2->GetChisquare(),func2->GetNDF());
std::ostringstream osx; 
TLatex* texCMSx = new TLatex;
texCMSx->SetNDC();
texCMSx->SetTextFont(42);
texCMSx->SetTextSize(0.05);
texCMSx->SetTextAlign(12);
texCMSx->DrawLatex(0.25,.95," #font[61]{CMS} #it{Work in Progress}");
texCMSx->DrawLatex(0.70,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMSx->DrawLatex(0.26,0.4, "#scale[1.7]{#Lambda}");
// texCMSx->DrawLatex(0.225,0.2, "#scale[1.0]{Used in #Lambda#Lambda}"); 
// texCMSx->DrawLatex(0.225,0.15, "#scale[1.0]{correlations}"); 
// texCMSx->DrawLatex(0.75,0.15, "This Thesis"); 
osx << "#scale[1.0]{s/#sqrt{s+b}: " << std::setprecision(4) << sig_la << "}" << std::setprecision(2);           
texCMSx->DrawLatex(0.68,0.60, osx.str(  ).c_str(  )); 
osx.str(std::string());
osx << "#scale[1.0]{s/(s+b): " << std::setprecision(4) << Fsig_la*100 << " %}" << std::setprecision(2); 
texCMSx->DrawLatex(0.68,0.52, osx.str(  ).c_str(  ));  
osx.str(std::string()); 
osx << "#scale[1.0]{Mean: " << std::setprecision(4) << meanla << " GeV}" << std::setprecision(6);   
texCMSx->DrawLatex(0.68,0.46, osx.str(  ).c_str(  ));
osx.str(std::string()); 
osx << "#scale[1.0]{#sigma: " << std::setprecision(4) << rms_true_la << " }" << std::setprecision(2); 
texCMSx->DrawLatex(0.68,0.40, osx.str(  ).c_str(  )); 
osx.str(std::string()); 
osx << "#scale[1.0]{Yield: " << std::setprecision(4) << S_la << " }" << std::setprecision(2); 
texCMSx->DrawLatex(0.68,0.33, osx.str(  ).c_str(  )); 
osx.str(std::string()); 
osx << "#scale[1.0]{#chi^{2}/ndf: " << std::setprecision(4) << chi2la << " }" << std::setprecision(4); 
texCMSx->DrawLatex(0.68,0.26, osx.str(  ).c_str(  ));  
osx.str(std::string()); 
// osx << "#scale[1.0]{Prob: " << std::setprecision(4) << probla << "}" << std::setprecision(2); 
// texCMSx->DrawLatex(0.68,0.20, osx.str(  ).c_str(  )); 
// osx.str(std::string()); 
osx << "#scale[1.]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMSx->DrawLatex(0.68,0.75, osx.str(  ).c_str(  )); 
osx.str(std::string()); 
osx << "#scale[1.]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMSx->DrawLatex(0.68,0.67, osx.str(  ).c_str(  ));  
osx.str(std::string()); 
osx << "#scale[1.]{ " << std::setprecision(4) << qminM << std::setprecision(4) <<" < q_{inv} < "  << qmaxM << " GeV}" << std::setprecision(2); 
texCMSx->DrawLatex(0.68,0.83, osx.str(  ).c_str(  ));  
osx.str(std::string()); 
Float_t yymax = lamass->GetMaximum();
TLine *linex = new TLine(signal_range_la.first,0.,signal_range_la.first,yymax);
linex->SetLineColor(SPcolors[SPdblue]);
linex->Draw();
TLine *line2x = new TLine(signal_range_la.second,0.,signal_range_la.second,yymax);
line2x->SetLineColor(SPcolors[SPdblue]);
line2x->Draw();
TLegend* legLL = new TLegend(0.16,0.6,0.43,0.85);
legLL->AddEntry(lamass,"Data");
legLL->AddEntry(func2,"Signal+Background","L");
legLL->AddEntry(func2_bckgnd,"Background","L");
legLL->AddEntry(linex,"Peak Region","L");
legLL->SetBorderSize(0);
legLL->Draw();
c2->Update();
pad2->cd();
gPad->SetTickx();
gPad->SetTicky();
double low = -49.1;;
double high = 49.1;;
TH1F* frame_la = pad2->cd()->DrawFrame(full_range_la.first,low,full_range_la.second,high);
frame_la->GetXaxis()->SetTitle("#scale[4.3]{p #pi^{-} invariant mass [GeV]}");
frame_la->GetYaxis()->SetTitle("#scale[4.0]{Pull}");
frame_la->GetXaxis()->CenterTitle(1);
frame_la->GetYaxis()->CenterTitle(1);
frame_la->GetXaxis()->SetNdivisions(510);
frame_la->GetYaxis()->SetNdivisions(406);
frame_la->GetXaxis()->SetTickSize(0.06);
frame_la->GetYaxis()->SetTickSize(0.02);
frame_la->GetXaxis()->SetTitleOffset(4.5);
frame_la->GetYaxis()->SetTitleOffset(1.8);
frame_la->GetXaxis()->SetLabelSize(0.14);
frame_la->GetYaxis()->SetLabelSize(0.14);
frame_la->GetXaxis()->SetLabelOffset(0.01);
frame_la->GetYaxis()->SetLabelOffset(0.013);
TGpull->SetMarkerStyle(25);
TGpull->SetMarkerSize(1.3);
TGpull->SetMarkerColor(kBlack);
TGpull->Draw("P");
TLine* line = new TLine(full_range_la.first, 0, full_range_la.second, 0);
line->SetLineStyle(2);
line->Draw("same");
pad3->cd();
gPad->SetTickx();
gPad->SetTicky();
double low2 = 0.1;
double high2 = 1.899999;
frame_la = pad3->cd()->DrawFrame(full_range_la.first,low2,full_range_la.second,high2);
frame_la->GetYaxis()->SetTitle("#scale[6.3]{#frac{Data}{Fit}}");
frame_la->GetXaxis()->CenterTitle(1);
frame_la->GetYaxis()->CenterTitle(1);
frame_la->GetXaxis()->SetNdivisions(510);
frame_la->GetYaxis()->SetNdivisions(406);
frame_la->GetXaxis()->SetTickSize(0.1);
frame_la->GetYaxis()->SetTickSize(0.01);
frame_la->GetYaxis()->SetTitleOffset(1.8);
frame_la->GetXaxis()->SetLabelSize(0.12);
frame_la->GetYaxis()->SetLabelSize(0.22);
frame_la->GetYaxis()->SetLabelOffset(0.01);
TGratio->SetMarkerStyle(25);
TGratio->SetMarkerSize(1.3);
TGratio->SetMarkerColor(kBlack);
TGratio->Draw("P");
TLine* line2 = new TLine(full_range_la.first, 1, full_range_la.second, 1);
line2->SetLineStyle(2);
line2->Draw("same");
c2->Update();
c2->Print("results/Mass"+V0+"LL_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");


//ALam used in ALamALam
cout << "==================ALam Mass for ALAL==================" << endl;
TH1D* masslam_LAL = (TH1D*)h_ALam_Mass_1D->Clone();
//masslam_LAL->Sumw2();
TF1 *func2_LAL = new TF1("func2_LAL", massfunc.c_str(), full_range_la.first, full_range_la.second);

//sideband fit
int cnt2_LAL = 0;
for ( int i = 0; i < masslam_LAL->GetNbinsX(); i++ ) {
   double center = masslam_LAL->GetBinCenter(i+1);
   if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
      cnt2_LAL++;
   }
}
TGraphErrors gr2_LAL(cnt2_LAL);
cnt2_LAL = 0;
for ( int i = 0; i < masslam_LAL->GetNbinsX(); i++ ) {
	double center = masslam_LAL->GetBinCenter(i+1);
    if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
       gr2_LAL.GetX()[cnt2_LAL] = center;
       gr2_LAL.GetY()[cnt2_LAL] = masslam_LAL->GetBinContent(i+1);
       gr2_LAL.GetEY()[cnt2_LAL] = masslam_LAL->GetBinError(i+1);
        cnt2_LAL++;
    }
}

//sideband fit parameters
func2_LAL->SetParameter(5, 10.);
func2_LAL->SetParameter(6, 1.);
func2_LAL->SetParameter(7, 1.);
func2_LAL->SetParameter(8, 1.);
func2_LAL->SetParameter(9, 1.);

//fix the parameters in order to get only bkg
func2_LAL->FixParameter(0, mass_la);
func2_LAL->FixParameter(1, sigma_la);
func2_LAL->FixParameter(2, sigma_la);
func2_LAL->FixParameter(3, 0.);
func2_LAL->FixParameter(4, 0.);

//make sideband fit
gr2_LAL.Fit(func2_LAL, "q0");
gr2_LAL.Fit(func2_LAL, "q0");
gr2_LAL.Fit(func2_LAL, "q0");

//Set sideband parameters to full fit
func2_LAL->FixParameter(5, func2_LAL->GetParameter(5));
func2_LAL->FixParameter(6, func2_LAL->GetParameter(6));
func2_LAL->FixParameter(7, func2_LAL->GetParameter(7));
func2_LAL->FixParameter(8, func2_LAL->GetParameter(8));
func2_LAL->FixParameter(9, func2_LAL->GetParameter(9));

//release the signal parameters
func2_LAL->ReleaseParameter(0);
func2_LAL->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2_LAL->ReleaseParameter(1);
func2_LAL->ReleaseParameter(2);
func2_LAL->ReleaseParameter(3);
func2_LAL->ReleaseParameter(4);
func2_LAL->SetParLimits(1, 0, sigma_la);
func2_LAL->SetParLimits(2, 0, sigma_la);
func2_LAL->SetParLimits(3, 0, 0.8*masslam_LAL->GetBinContent(masslam_LAL->FindBin(mass_la)));
func2_LAL->SetParLimits(4, 0, 0.8*masslam_LAL->GetBinContent(masslam_LAL->FindBin(mass_la)));

//make a full fit
masslam_LAL->Fit(func2_LAL, "q0", "", full_range_la.first, full_range_la.second);
masslam_LAL->Fit(func2_LAL, "q0", "", full_range_la.first, full_range_la.second);
masslam_LAL->Fit(func2_LAL, "q0", "", full_range_la.first, full_range_la.second);

//now release all the parameters
func2_LAL->ReleaseParameter(0);
func2_LAL->ReleaseParameter(1);
func2_LAL->ReleaseParameter(2);
func2_LAL->ReleaseParameter(3);
func2_LAL->ReleaseParameter(4);
func2_LAL->ReleaseParameter(5);
func2_LAL->ReleaseParameter(6);
func2_LAL->ReleaseParameter(7);
func2_LAL->ReleaseParameter(8);
func2_LAL->ReleaseParameter(9);
func2_LAL->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2_LAL->SetParLimits(1, 0, sigma_la);
func2_LAL->SetParLimits(2, 0, sigma_la);
func2_LAL->SetParLimits(3, 0, 0.8*masslam_LAL->GetBinContent(masslam_LAL->FindBin(mass_la)));
func2_LAL->SetParLimits(4, 0, 0.8*masslam_LAL->GetBinContent(masslam_LAL->FindBin(mass_la)));

//make final fits
for ( int i = 0; i < iter_la-1; i++ ) {
    masslam_LAL->Fit(func2_LAL, "q0", "", full_range_la.first, full_range_la.second);

}

//make signal and bkg functions for plots
TF1 * func2_signal_LAL   = 0;
TF1 * func2_bckgnd_LAL   = 0;
TF1 * gaus1_signal_la_LAL   = 0;
TF1 * gaus2_signal_la_LAL   = 0;
func2_signal_LAL  = new TF1("func2_signal_LAL" , double_gaussian.c_str(), full_range_la.first, full_range_la.second);
func2_bckgnd_LAL = new TF1("func2_bckgnd_LAL" , poly_bkg0.c_str(), full_range_la.first, full_range_la.second);
gaus1_signal_la_LAL = new TF1("gaus1_signal_la_LAL", gaus_func.c_str(), full_range_la.first, full_range_la.second);
gaus2_signal_la_LAL = new TF1("gaus2_signal_la_LAL", gaus_func.c_str(), full_range_la.first, full_range_la.second);

//set parameters for signal and bkg
func2_signal_LAL->SetParameter(0, func2_LAL->GetParameter(0));
func2_signal_LAL->SetParameter(1, func2_LAL->GetParameter(1));
func2_signal_LAL->SetParameter(2, func2_LAL->GetParameter(2));
func2_signal_LAL->SetParameter(3, func2_LAL->GetParameter(3));
func2_signal_LAL->SetParameter(4, func2_LAL->GetParameter(4));
func2_bckgnd_LAL->SetParameter(0, func2_LAL->GetParameter(5));
func2_bckgnd_LAL->SetParameter(1, func2_LAL->GetParameter(6));
func2_bckgnd_LAL->SetParameter(2, func2_LAL->GetParameter(7));
func2_bckgnd_LAL->SetParameter(3, func2_LAL->GetParameter(8));
func2_bckgnd_LAL->SetParameter(4, func2_LAL->GetParameter(9));

//set parameters for gaussians
gaus1_signal_la_LAL->SetParameter(0, func2_LAL->GetParameter(0));
gaus1_signal_la_LAL->SetParameter(1, func2_LAL->GetParameter(1));
gaus1_signal_la_LAL->SetParameter(2, func2_LAL->GetParameter(3));
gaus2_signal_la_LAL->SetParameter(0, func2_LAL->GetParameter(0));
gaus2_signal_la_LAL->SetParameter(1, func2_LAL->GetParameter(2));
gaus2_signal_la_LAL->SetParameter(2, func2_LAL->GetParameter(4));

//make pull and ratio
double vny2_LAL[400];
double xvn2_LAL[400];
double pullvn2_LAL[400];
double ratiovn2_LAL[400];
double err_pullvn2_LAL[400];
double err_ratiovn2_LAL[400];
for(int k=0;k<masslam_LAL->GetNbinsX();k++)
            {
            vny2_LAL[k] = masslam_LAL->GetBinContent(k+1);
            xvn2_LAL[k] = masslam_LAL->GetBinCenter(k+1);
        	if(masslam_LAL->GetBinContent(k+1)!=0){
            pullvn2_LAL[k] = (vny2_LAL[k] - func2_LAL->Eval(xvn2_LAL[k]))/masslam_LAL->GetBinError(k+1);
            ratiovn2_LAL[k] = vny2_LAL[k]/func2_LAL->Eval(xvn2_LAL[k]);
            err_pullvn2_LAL[k] = masslam_LAL->GetBinError(k+1);
            err_ratiovn2_LAL[k] = (vny2_LAL[k]/func2_LAL->Eval(xvn2_LAL[k]))*(masslam_LAL->GetBinError(k+1)/vny2_LAL[k]);
            }else{
            pullvn2_LAL[k] = 0.0;
            ratiovn2_LAL[k] = 0.0;
            err_pullvn2_LAL[k]=0.0;
            err_ratiovn2_LAL[k]=0.0;
            }
}

//make the plots
TH1D* lamass_LAL = (TH1D*)masslam_LAL->Clone();
//lamass_LAL->Sumw2();
TGraphErrors* TGpull_LAL = new TGraphErrors(400,xvn2_LAL,pullvn2_LAL,0,0);
TGraphErrors* TGratio_LAL= new TGraphErrors(400,xvn2_LAL,ratiovn2_LAL,0,err_ratiovn2_LAL);
TCanvas* c2_LAL = new TCanvas("c2_LAL","c2_LAL",600,600);
c2_LAL->cd();
TPad* pad1_LAL = new TPad("pad1_LAL","top pad",0.0,0.4,1.0,1.0);
TPad* pad2_LAL = new TPad("pad3_LAL","middle pad",0.0,0.0,1.0,0.25);
TPad* pad3_LAL = new TPad("pad2_LAL","bottom pad",0.0,0.25,1.0,0.4);
pad1_LAL->SetTopMargin(0.1);
pad1_LAL->SetBottomMargin(0.0);
pad1_LAL->SetRightMargin(0.038);
pad1_LAL->SetLeftMargin(0.15);
pad1_LAL->Draw();
pad2_LAL->SetFrameFillStyle(4000);
pad2_LAL->SetTopMargin(0.0);
pad2_LAL->SetBottomMargin(0.37);
pad2_LAL->SetRightMargin(0.038);
pad2_LAL->SetLeftMargin(0.15);
pad2_LAL->Draw();
pad3_LAL->SetFrameFillStyle(4000);
pad3_LAL->SetTopMargin(0.0);
pad3_LAL->SetBottomMargin(0.0);
pad3_LAL->SetRightMargin(0.038);
pad3_LAL->SetLeftMargin(0.15);
pad3_LAL->Draw();
pad1_LAL->cd();
gPad->SetTickx();
gPad->SetTicky();
lamass_LAL->SetTitle("");
lamass_LAL->Draw("E");
TGaxis::SetMaxDigits(3);
lamass_LAL->GetXaxis()->SetRangeUser(full_range_la.first,full_range_la.second);
lamass_LAL->GetYaxis()->SetRangeUser(0.1,1.05*lamass_LAL->GetMaximum());
lamass_LAL->SetMarkerStyle(24);
lamass_LAL->SetMarkerSize(1.5);
lamass_LAL->GetYaxis()->SetTitle("#scale[1.7]{Candidates / 0.5 MeV}");
lamass_LAL->GetYaxis()->CenterTitle(1);
lamass_LAL->GetYaxis()->CenterTitle(1);
lamass_LAL->GetYaxis()->SetNdivisions(810);
lamass_LAL->GetXaxis()->SetTickSize(0.03);
lamass_LAL->GetYaxis()->SetTickSize(0.02);
lamass_LAL->GetXaxis()->SetLabelSize(0.055);
lamass_LAL->GetYaxis()->SetLabelSize(0.058);
lamass_LAL->GetYaxis()->SetLabelOffset(0.01);
lamass_LAL->GetYaxis()->SetTitleOffset(1.8);
func2_LAL->Draw("LSAME");
func2_LAL->SetLineColor(SPcolors[SPdred]);
func2_bckgnd_LAL->Draw("LSAME");
func2_bckgnd_LAL->SetLineColor(SPcolors[SPdred]);
func2_bckgnd_LAL->SetLineStyle(2);
double meanla_LAL  = func2_LAL->GetParameter(0);
double sigma1la_LAL = func2_LAL->GetParameter(1);
double sigma2la_LAL = func2_LAL->GetParameter(2);
double gaus1_yield_la_LAL = gaus1_signal_la_LAL->Integral(meanla_LAL - 2*sigma1la_LAL, meanla_LAL + 2*sigma1la)/masslam_LAL->GetBinWidth(1);
double gaus2_yield_la_LAL = gaus2_signal_la_LAL->Integral(meanla_LAL - 2*sigma2la_LAL, meanla_LAL + 2*sigma2la)/masslam_LAL->GetBinWidth(1);
double gausTot_yield_la_LAL = gaus1_yield_la_LAL + gaus2_yield_la_LAL;
double rms_gaus1_sig_la_LAL = gaus1_yield_la_LAL/gausTot_yield_la_LAL;
double rms_gaus2_sig_la_LAL = gaus2_yield_la_LAL/gausTot_yield_la_LAL;
double rms_true_la_LAL = TMath::Sqrt(rms_gaus1_sig_la_LAL*sigma1la_LAL*sigma1la_LAL + rms_gaus2_sig_la_LAL*sigma2la_LAL*sigma2la_LAL);
std::pair<double, double>  signal_range_la_LAL(meanla_LAL - 2*rms_true_la_LAL, meanla_LAL + 2*rms_true_la_LAL);
double S_la_LAL = func2_signal_LAL->Integral(signal_range_la_LAL.first, signal_range_la_LAL.second)/masslam_LAL->GetBinWidth(1);
double B_la_LAL = func2_bckgnd_LAL->Integral(signal_range_la_LAL.first, signal_range_la_LAL.second)/masslam_LAL->GetBinWidth(1);
double Fsig_la_LAL = S_la_LAL/(S_la_LAL+B_la_LAL);
double sig_la_LAL = S_la_LAL/sqrt(S_la_LAL+B_la_LAL);
double chi2la_LAL = func2_LAL->GetChisquare()/func2_LAL->GetNDF();
double probla_LAL = TMath::Prob(func2_LAL->GetChisquare(),func2_LAL->GetNDF());
std::ostringstream osx_LAL; 
TLatex* texCMSx_LAL = new TLatex;
texCMSx_LAL->SetNDC();
texCMSx_LAL->SetTextFont(42);
texCMSx_LAL->SetTextSize(0.05);
texCMSx_LAL->SetTextAlign(12);
texCMSx_LAL->DrawLatex(0.25,.95," #font[61]{CMS} #it{Work in Progress}");
texCMSx_LAL->DrawLatex(0.70,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMSx_LAL->DrawLatex(0.26,0.4, "#scale[1.7]{#bar{#Lambda}}");
// texCMSx_LAL->DrawLatex(0.225,0.2, "#scale[1.0]{Used in #bar{#Lambda}#bar{#Lambda}}"); 
// texCMSx_LAL->DrawLatex(0.225,0.15, "#scale[1.0]{correlations}"); 
// texCMSx_LAL->DrawLatex(0.75,0.15, "This Thesis"); 
osx_LAL << "#scale[1.0]{s/#sqrt{s+b}: " << std::setprecision(4) << sig_la_LAL << "}" << std::setprecision(2);           
texCMSx_LAL->DrawLatex(0.68,0.60, osx_LAL.str(  ).c_str(  )); 
osx_LAL.str(std::string());
osx_LAL << "#scale[1.0]{s/(s+b): " << std::setprecision(4) << Fsig_la_LAL*100 << " %}" << std::setprecision(2); 
texCMSx_LAL->DrawLatex(0.68,0.52, osx_LAL.str(  ).c_str(  ));  
osx_LAL.str(std::string()); 
osx_LAL << "#scale[1.0]{Mean: " << std::setprecision(4) << meanla_LAL << " GeV}" << std::setprecision(6);   
texCMSx_LAL->DrawLatex(0.68,0.46, osx_LAL.str(  ).c_str(  ));
osx_LAL.str(std::string()); 
osx_LAL << "#scale[1.0]{#sigma: " << std::setprecision(4) << rms_true_la_LAL << " }" << std::setprecision(2); 
texCMSx_LAL->DrawLatex(0.68,0.40, osx_LAL.str(  ).c_str(  )); 
osx_LAL.str(std::string()); 
osx_LAL << "#scale[1.0]{Yield: " << std::setprecision(4) << S_la_LAL << " }" << std::setprecision(2); 
texCMSx_LAL->DrawLatex(0.68,0.33, osx_LAL.str(  ).c_str(  )); 
osx_LAL.str(std::string()); 
osx_LAL << "#scale[1.0]{#chi^{2}/ndf: " << std::setprecision(4) << chi2la_LAL << " }" << std::setprecision(4); 
texCMSx_LAL->DrawLatex(0.68,0.26, osx_LAL.str(  ).c_str(  ));  
osx_LAL.str(std::string()); 
// osx_LAL << "#scale[1.0]{Prob: " << std::setprecision(4) << probla_LAL << "}" << std::setprecision(2); 
// texCMSx_LAL->DrawLatex(0.68,0.20, osx_LAL.str(  ).c_str(  )); 
// osx_LAL.str(std::string()); 
osx_LAL << "#scale[1.]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMSx_LAL->DrawLatex(0.68,0.75, osx_LAL.str(  ).c_str(  )); 
osx_LAL.str(std::string()); 
osx_LAL << "#scale[1.]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMSx_LAL->DrawLatex(0.68,0.67, osx_LAL.str(  ).c_str(  ));  
osx_LAL.str(std::string()); 
osx_LAL << "#scale[1.]{ " << std::setprecision(4) << qminM << std::setprecision(4) <<" < q_{inv} < "  << qmaxM << " GeV}" << std::setprecision(2); 
texCMSx_LAL->DrawLatex(0.68,0.83, osx_LAL.str(  ).c_str(  ));  
osx_LAL.str(std::string()); 
Float_t yymax_LAL = lamass_LAL->GetMaximum();
TLine *linex_LAL = new TLine(signal_range_la_LAL.first,0.,signal_range_la_LAL.first,yymax_LAL);
linex_LAL->SetLineColor(SPcolors[SPdblue]);
linex_LAL->Draw();
TLine *line2x_LAL = new TLine(signal_range_la_LAL.second,0.,signal_range_la_LAL.second,yymax_LAL);
line2x_LAL->SetLineColor(SPcolors[SPdblue]);
line2x_LAL->Draw();
TLegend* legLL_LAL = new TLegend(0.16,0.6,0.43,0.85);
legLL_LAL->AddEntry(lamass_LAL,"Data");
legLL_LAL->AddEntry(func2_LAL,"Signal+Background","L");
legLL_LAL->AddEntry(func2_bckgnd_LAL,"Background","L");
legLL_LAL->AddEntry(linex,"Peak Region","L");
legLL_LAL->SetBorderSize(0);
legLL_LAL->Draw();
c2_LAL->Update();
pad2_LAL->cd();
gPad->SetTickx();
gPad->SetTicky();
double low_LAL = -49.1;;
double high_LAL = 49.1;;
TH1F* frame_la_LAL = pad2_LAL->cd()->DrawFrame(full_range_la.first,low_LAL,full_range_la.second,high_LAL);
frame_la_LAL->GetXaxis()->SetTitle("#scale[4.3]{#bar{p} #pi^{+} invariant mass [GeV]}");
frame_la_LAL->GetYaxis()->SetTitle("#scale[4.0]{Pull}");
frame_la_LAL->GetXaxis()->CenterTitle(1);
frame_la_LAL->GetYaxis()->CenterTitle(1);
frame_la_LAL->GetXaxis()->SetNdivisions(510);
frame_la_LAL->GetYaxis()->SetNdivisions(406);
frame_la_LAL->GetXaxis()->SetTickSize(0.06);
frame_la_LAL->GetYaxis()->SetTickSize(0.02);
frame_la_LAL->GetXaxis()->SetTitleOffset(4.5);
frame_la_LAL->GetYaxis()->SetTitleOffset(1.8);
frame_la_LAL->GetXaxis()->SetLabelSize(0.14);
frame_la_LAL->GetYaxis()->SetLabelSize(0.14);
frame_la_LAL->GetXaxis()->SetLabelOffset(0.01);
frame_la_LAL->GetYaxis()->SetLabelOffset(0.013);
TGpull_LAL->SetMarkerStyle(25);
TGpull_LAL->SetMarkerSize(1.3);
TGpull_LAL->SetMarkerColor(kBlack);
TGpull_LAL->Draw("P");
TLine* line_LAL = new TLine(full_range_la.first, 0, full_range_la.second, 0);
line_LAL->SetLineStyle(2);
line_LAL->Draw("same");
pad3_LAL->cd();
gPad->SetTickx();
gPad->SetTicky();
double low2_LAL = 0.1;
double high2_LAL = 1.8999999;
frame_la_LAL = pad3_LAL->cd()->DrawFrame(full_range_la.first,low2_LAL,full_range_la.second,high2_LAL);
frame_la_LAL->GetYaxis()->SetTitle("#scale[6.3]{#frac{Data}{Fit}}");
frame_la_LAL->GetXaxis()->CenterTitle(1);
frame_la_LAL->GetYaxis()->CenterTitle(1);
frame_la_LAL->GetXaxis()->SetNdivisions(510);
frame_la_LAL->GetYaxis()->SetNdivisions(406);
frame_la_LAL->GetXaxis()->SetTickSize(0.1);
frame_la_LAL->GetYaxis()->SetTickSize(0.01);
frame_la_LAL->GetYaxis()->SetTitleOffset(1.8);
frame_la_LAL->GetXaxis()->SetLabelSize(0.12);
frame_la_LAL->GetYaxis()->SetLabelSize(0.22);
frame_la_LAL->GetYaxis()->SetLabelOffset(0.01);
TGratio_LAL->SetMarkerStyle(25);
TGratio_LAL->SetMarkerSize(1.3);
TGratio_LAL->SetMarkerColor(kBlack);
TGratio_LAL->Draw("P");
TLine* line2_LAL = new TLine(full_range_la.first, 1, full_range_la.second, 1);
line2_LAL->SetLineStyle(2);
line2_LAL->Draw("same");
c2_LAL->Update();
c2_LAL->Print("results/Mass"+V0+"ALAL_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");


//Lam used in LamALam
cout << "==================Lam Mass for LAL==================" << endl;
TH1D* masslam_LLAL = (TH1D*)h_Lam_Mass_LAL_1D->Clone();
//masslam_LLAL->Sumw2();
TF1 *func2_LLAL = new TF1("func2_LLAL", massfunc.c_str(), full_range_la.first, full_range_la.second);
//sideband fit
int cnt2_LLAL = 0;
for ( int i = 0; i < masslam_LLAL->GetNbinsX(); i++ ) {
    double center = masslam_LLAL->GetBinCenter(i+1);
    if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
       cnt2_LLAL++;
    }
}
TGraphErrors gr2_LLAL(cnt2_LLAL);
cnt2_LLAL = 0;
for ( int i = 0; i < masslam_LLAL->GetNbinsX(); i++ ) {
    double center = masslam_LLAL->GetBinCenter(i+1);
    if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
       gr2_LLAL.GetX()[cnt2_LLAL] = center;
       gr2_LLAL.GetY()[cnt2_LLAL] = masslam_LLAL->GetBinContent(i+1);
       gr2_LLAL.GetEY()[cnt2_LLAL] = masslam_LLAL->GetBinError(i+1);
       cnt2_LLAL++;
    }
}

//sideband fit parameters
func2_LLAL->SetParameter(5, 10.);
func2_LLAL->SetParameter(6, 1.);
func2_LLAL->SetParameter(7, 1.);
func2_LLAL->SetParameter(8, 1.);
func2_LLAL->SetParameter(9, 1.);

//fix the parameters in order to get only bkg
func2_LLAL->FixParameter(0, mass_la);
func2_LLAL->FixParameter(1, sigma_la);
func2_LLAL->FixParameter(2, sigma_la);
func2_LLAL->FixParameter(3, 0.);
func2_LLAL->FixParameter(4, 0.);

//make sideband fit
gr2_LLAL.Fit(func2_LLAL, "q0");
gr2_LLAL.Fit(func2_LLAL, "q0");
gr2_LLAL.Fit(func2_LLAL, "q0");

//Set sideband parameters to full fit
func2_LLAL->FixParameter(5, func2_LLAL->GetParameter(5));
func2_LLAL->FixParameter(6, func2_LLAL->GetParameter(6));
func2_LLAL->FixParameter(7, func2_LLAL->GetParameter(7));
func2_LLAL->FixParameter(8, func2_LLAL->GetParameter(8));
func2_LLAL->FixParameter(9, func2_LLAL->GetParameter(9));

//release the signal parameters
func2_LLAL->ReleaseParameter(0);
func2_LLAL->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2_LLAL->ReleaseParameter(1);
func2_LLAL->ReleaseParameter(2);
func2_LLAL->ReleaseParameter(3);
func2_LLAL->ReleaseParameter(4);
func2_LLAL->SetParLimits(1, 0, sigma_la);
func2_LLAL->SetParLimits(2, 0, sigma_la);
func2_LLAL->SetParLimits(3, 0, 0.8*masslam_LLAL->GetBinContent(masslam_LLAL->FindBin(mass_la)));
func2_LLAL->SetParLimits(4, 0, 0.8*masslam_LLAL->GetBinContent(masslam_LLAL->FindBin(mass_la)));

//make a full fit
masslam_LLAL->Fit(func2_LLAL, "q0", "", full_range_la.first, full_range_la.second);
masslam_LLAL->Fit(func2_LLAL, "q0", "", full_range_la.first, full_range_la.second);
masslam_LLAL->Fit(func2_LLAL, "q0", "", full_range_la.first, full_range_la.second);

//now release all the parameters
func2_LLAL->ReleaseParameter(0);
func2_LLAL->ReleaseParameter(1);
func2_LLAL->ReleaseParameter(2);
func2_LLAL->ReleaseParameter(3);
func2_LLAL->ReleaseParameter(4);
func2_LLAL->ReleaseParameter(5);
func2_LLAL->ReleaseParameter(6);
func2_LLAL->ReleaseParameter(7);
func2_LLAL->ReleaseParameter(8);
func2_LLAL->ReleaseParameter(9);
func2_LLAL->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2_LLAL->SetParLimits(1, 0, sigma_la);
func2_LLAL->SetParLimits(2, 0, sigma_la);
func2_LLAL->SetParLimits(3, 0, 0.8*masslam_LLAL->GetBinContent(masslam_LLAL->FindBin(mass_la)));
func2_LLAL->SetParLimits(4, 0, 0.8*masslam_LLAL->GetBinContent(masslam_LLAL->FindBin(mass_la)));

//make final fits
for ( int i = 0; i < iter_la-1; i++ ) {
    masslam_LLAL->Fit(func2_LLAL, "q0", "", full_range_la.first, full_range_la.second);

}

//make signal and bkg functions for plots
TF1 * func2_signal_LLAL   = 0;
TF1 * func2_bckgnd_LLAL   = 0;
TF1 * gaus1_signal_la_LLAL   = 0;
TF1 * gaus2_signal_la_LLAL   = 0;
func2_signal_LLAL  = new TF1("func2_signal_LLAL" , double_gaussian.c_str(), full_range_la.first, full_range_la.second);
func2_bckgnd_LLAL = new TF1("func2_bckgnd_LLAL" , poly_bkg0.c_str(), full_range_la.first, full_range_la.second);
gaus1_signal_la_LLAL = new TF1("gaus1_signal_la_LLAL", gaus_func.c_str(), full_range_la.first, full_range_la.second);
gaus2_signal_la_LLAL = new TF1("gaus2_signal_la_LLAL", gaus_func.c_str(), full_range_la.first, full_range_la.second);

//set parameters for signal and bkg
func2_signal_LLAL->SetParameter(0, func2_LLAL->GetParameter(0));
func2_signal_LLAL->SetParameter(1, func2_LLAL->GetParameter(1));
func2_signal_LLAL->SetParameter(2, func2_LLAL->GetParameter(2));
func2_signal_LLAL->SetParameter(3, func2_LLAL->GetParameter(3));
func2_signal_LLAL->SetParameter(4, func2_LLAL->GetParameter(4));
func2_bckgnd_LLAL->SetParameter(0, func2_LLAL->GetParameter(5));
func2_bckgnd_LLAL->SetParameter(1, func2_LLAL->GetParameter(6));
func2_bckgnd_LLAL->SetParameter(2, func2_LLAL->GetParameter(7));
func2_bckgnd_LLAL->SetParameter(3, func2_LLAL->GetParameter(8));
func2_bckgnd_LLAL->SetParameter(4, func2_LLAL->GetParameter(9));

//set parameters for gaussians
gaus1_signal_la_LLAL->SetParameter(0, func2_LLAL->GetParameter(0));
gaus1_signal_la_LLAL->SetParameter(1, func2_LLAL->GetParameter(1));
gaus1_signal_la_LLAL->SetParameter(2, func2_LLAL->GetParameter(3));
gaus2_signal_la_LLAL->SetParameter(0, func2_LLAL->GetParameter(0));
gaus2_signal_la_LLAL->SetParameter(1, func2_LLAL->GetParameter(2));
gaus2_signal_la_LLAL->SetParameter(2, func2_LLAL->GetParameter(4));

//make pull and ratio
double vny2_LLAL[400];
double xvn2_LLAL[400];
double pullvn2_LLAL[400];
double ratiovn2_LLAL[400];
double err_pullvn2_LLAL[400];
double err_ratiovn2_LLAL[400];
for(int k=0;k<masslam_LLAL->GetNbinsX();k++)
       		{
            vny2_LLAL[k] = masslam_LLAL->GetBinContent(k+1);
            xvn2_LLAL[k] = masslam_LLAL->GetBinCenter(k+1);
        	if(masslam_LLAL->GetBinContent(k+1)!=0){
            pullvn2_LLAL[k] = (vny2_LLAL[k] - func2_LLAL->Eval(xvn2_LLAL[k]))/masslam_LLAL->GetBinError(k+1);
            ratiovn2_LLAL[k] = vny2_LLAL[k]/func2_LLAL->Eval(xvn2_LLAL[k]);
            err_pullvn2_LLAL[k] = masslam_LLAL->GetBinError(k+1);
            err_ratiovn2_LLAL[k] = (vny2_LLAL[k]/func2_LLAL->Eval(xvn2_LLAL[k]))*(masslam_LLAL->GetBinError(k+1)/vny2_LLAL[k]);
            }else{
            pullvn2_LLAL[k] = 0.0;
            ratiovn2_LLAL[k] = 0.0;
            err_pullvn2_LLAL[k]=0.0;
            err_ratiovn2_LLAL[k]=0.0;
            }
        }

//make the plots
TH1D* lamass_LLAL = (TH1D*)masslam_LLAL->Clone();
//lamass_LLAL->Sumw2();
TGraphErrors* TGpull_LLAL = new TGraphErrors(400,xvn2_LLAL,pullvn2_LLAL,0,0);
TGraphErrors* TGratio_LLAL= new TGraphErrors(400,xvn2_LLAL,ratiovn2_LLAL,0,err_ratiovn2_LLAL);
TCanvas* c2_LLAL = new TCanvas("c2_LLAL","c2_LLAL",600,600);
c2_LLAL->cd();
TPad* pad1_LLAL = new TPad("pad1_LLAL","top pad",0.0,0.4,1.0,1.0);
TPad* pad2_LLAL = new TPad("pad3_LLAL","middle pad",0.0,0.0,1.0,0.25);
TPad* pad3_LLAL = new TPad("pad2_LLAL","bottom pad",0.0,0.25,1.0,0.4);
pad1_LLAL->SetTopMargin(0.1);
pad1_LLAL->SetBottomMargin(0.0);
pad1_LLAL->SetRightMargin(0.038);
pad1_LLAL->SetLeftMargin(0.15);
pad1_LLAL->Draw();
pad2_LLAL->SetFrameFillStyle(4000);
pad2_LLAL->SetTopMargin(0.0);
pad2_LLAL->SetBottomMargin(0.37);
pad2_LLAL->SetRightMargin(0.038);
pad2_LLAL->SetLeftMargin(0.15);
pad2_LLAL->Draw();
pad3_LLAL->SetFrameFillStyle(4000);
pad3_LLAL->SetTopMargin(0.0);
pad3_LLAL->SetBottomMargin(0.0);
pad3_LLAL->SetRightMargin(0.038);
pad3_LLAL->SetLeftMargin(0.15);
pad3_LLAL->Draw();
pad1_LLAL->cd();
gPad->SetTickx();
gPad->SetTicky();
lamass_LLAL->SetTitle("");
lamass_LLAL->Draw("E");
TGaxis::SetMaxDigits(3);
lamass_LLAL->GetXaxis()->SetRangeUser(full_range_la.first,full_range_la.second);
lamass_LLAL->GetYaxis()->SetRangeUser(0.1,1.05*lamass_LLAL->GetMaximum());
lamass_LLAL->SetMarkerStyle(24);
lamass_LLAL->SetMarkerSize(1.5);
lamass_LLAL->GetYaxis()->SetTitle("#scale[1.7]{Candidates / 0.5 MeV}");
lamass_LLAL->GetYaxis()->CenterTitle(1);
lamass_LLAL->GetYaxis()->CenterTitle(1);
lamass_LLAL->GetYaxis()->SetNdivisions(810);
lamass_LLAL->GetXaxis()->SetTickSize(0.03);
lamass_LLAL->GetYaxis()->SetTickSize(0.02);
lamass_LLAL->GetXaxis()->SetLabelSize(0.055);
lamass_LLAL->GetYaxis()->SetLabelSize(0.058);
lamass_LLAL->GetYaxis()->SetLabelOffset(0.01);
lamass_LLAL->GetYaxis()->SetTitleOffset(1.8);
func2_LLAL->Draw("LSAME");
func2_LLAL->SetLineColor(SPcolors[SPdred]);
func2_bckgnd_LLAL->Draw("LSAME");
func2_bckgnd_LLAL->SetLineColor(SPcolors[SPdred]);
func2_bckgnd_LLAL->SetLineStyle(2);
double meanla_LLAL  = func2_LLAL->GetParameter(0);
double sigma1la_LLAL = func2_LLAL->GetParameter(1);
double sigma2la_LLAL = func2_LLAL->GetParameter(2);
double gaus1_yield_la_LLAL = gaus1_signal_la_LLAL->Integral(meanla_LLAL - 2*sigma1la_LLAL, meanla_LLAL + 2*sigma1la)/masslam_LLAL->GetBinWidth(1);
double gaus2_yield_la_LLAL = gaus2_signal_la_LLAL->Integral(meanla_LLAL - 2*sigma2la_LLAL, meanla_LLAL + 2*sigma2la)/masslam_LLAL->GetBinWidth(1);
double gausTot_yield_la_LLAL = gaus1_yield_la_LLAL + gaus2_yield_la_LLAL;
double rms_gaus1_sig_la_LLAL = gaus1_yield_la_LLAL/gausTot_yield_la_LLAL;
double rms_gaus2_sig_la_LLAL = gaus2_yield_la_LLAL/gausTot_yield_la_LLAL;
double rms_true_la_LLAL = TMath::Sqrt(rms_gaus1_sig_la_LLAL*sigma1la_LLAL*sigma1la_LLAL + rms_gaus2_sig_la_LLAL*sigma2la_LLAL*sigma2la_LLAL);
std::pair<double, double>  signal_range_la_LLAL(meanla_LLAL - 2*rms_true_la_LLAL, meanla_LLAL + 2*rms_true_la_LLAL);
double S_la_LLAL = func2_signal_LLAL->Integral(signal_range_la_LLAL.first, signal_range_la_LLAL.second)/masslam_LLAL->GetBinWidth(1);
double B_la_LLAL = func2_bckgnd_LLAL->Integral(signal_range_la_LLAL.first, signal_range_la_LLAL.second)/masslam_LLAL->GetBinWidth(1);
double Fsig_la_LLAL = S_la_LLAL/(S_la_LLAL+B_la_LLAL);
double sig_la_LLAL = S_la_LLAL/sqrt(S_la_LLAL+B_la_LLAL);
double chi2la_LLAL = func2_LLAL->GetChisquare()/func2_LLAL->GetNDF();
double probla_LLAL = TMath::Prob(func2_LLAL->GetChisquare(),func2_LLAL->GetNDF());
std::ostringstream osx_LLAL; 
TLatex* texCMSx_LLAL = new TLatex;
texCMSx_LLAL->SetNDC();
texCMSx_LLAL->SetTextFont(42);
texCMSx_LLAL->SetTextSize(0.05);
texCMSx_LLAL->SetTextAlign(12);
texCMSx_LLAL->DrawLatex(0.25,.95," #font[61]{CMS} #it{Work in Progress}");
texCMSx_LLAL->DrawLatex(0.70,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMSx_LLAL->DrawLatex(0.26,0.4, "#scale[1.7]{#Lambda}");
texCMSx_LLAL->DrawLatex(0.225,0.2, "#scale[1.0]{Used in #Lambda#bar{#Lambda}}"); 
texCMSx_LLAL->DrawLatex(0.225,0.15, "#scale[1.0]{correlations}"); 
// texCMSx_LLAL->DrawLatex(0.75,0.15, "This Thesis"); 
osx_LLAL << "#scale[1.0]{s/#sqrt{s+b}: " << std::setprecision(4) << sig_la_LLAL << "}" << std::setprecision(2);           
texCMSx_LLAL->DrawLatex(0.68,0.60, osx_LLAL.str(  ).c_str(  )); 
osx_LLAL.str(std::string());
osx_LLAL << "#scale[1.0]{s/(s+b): " << std::setprecision(4) << Fsig_la_LLAL*100 << " %}" << std::setprecision(2); 
texCMSx_LLAL->DrawLatex(0.68,0.52, osx_LLAL.str(  ).c_str(  ));  
osx_LLAL.str(std::string()); 
osx_LLAL << "#scale[1.0]{Mean: " << std::setprecision(4) << meanla_LLAL << " GeV}" << std::setprecision(6);   
texCMSx_LLAL->DrawLatex(0.68,0.46, osx_LLAL.str(  ).c_str(  ));
osx_LLAL.str(std::string()); 
osx_LLAL << "#scale[1.0]{#sigma: " << std::setprecision(4) << rms_true_la_LLAL << " }" << std::setprecision(2); 
texCMSx_LLAL->DrawLatex(0.68,0.40, osx_LLAL.str(  ).c_str(  )); 
osx_LLAL.str(std::string()); 
osx_LLAL << "#scale[1.0]{Yield: " << std::setprecision(4) << S_la_LLAL << " }" << std::setprecision(2); 
texCMSx_LLAL->DrawLatex(0.68,0.33, osx_LLAL.str(  ).c_str(  )); 
osx_LLAL.str(std::string()); 
osx_LLAL << "#scale[1.0]{#chi^{2}/ndf: " << std::setprecision(4) << chi2la_LLAL << " }" << std::setprecision(4); 
texCMSx_LLAL->DrawLatex(0.68,0.26, osx_LLAL.str(  ).c_str(  ));  
osx_LLAL.str(std::string()); 
// osx_LLAL << "#scale[1.0]{Prob: " << std::setprecision(4) << probla_LLAL << "}" << std::setprecision(2); 
// texCMSx_LLAL->DrawLatex(0.68,0.20, osx_LLAL.str(  ).c_str(  )); 
// osx_LLAL.str(std::string()); 
osx_LLAL << "#scale[1.]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMSx_LLAL->DrawLatex(0.68,0.75, osx_LLAL.str(  ).c_str(  )); 
osx_LLAL.str(std::string()); 
osx_LLAL << "#scale[1.]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMSx_LLAL->DrawLatex(0.68,0.67, osx_LLAL.str(  ).c_str(  ));  
osx_LLAL.str(std::string()); 
osx_LLAL << "#scale[1.]{ " << std::setprecision(4) << qminM << std::setprecision(4) <<" < q_{inv} < "  << qmaxM << " GeV}" << std::setprecision(2); 
texCMSx_LLAL->DrawLatex(0.68,0.83, osx_LLAL.str(  ).c_str(  ));  
osx_LLAL.str(std::string()); 
Float_t yymax_LLAL = lamass_LLAL->GetMaximum();
TLine *linex_LLAL = new TLine(signal_range_la_LLAL.first,0.,signal_range_la_LLAL.first,yymax_LLAL);
linex_LLAL->SetLineColor(SPcolors[SPdblue]);
linex_LLAL->Draw();
TLine *line2x_LLAL = new TLine(signal_range_la_LLAL.second,0.,signal_range_la_LLAL.second,yymax_LLAL);
line2x_LLAL->SetLineColor(SPcolors[SPdblue]);
line2x_LLAL->Draw();
TLegend* legLL_LLAL = new TLegend(0.16,0.6,0.43,0.85);
legLL_LLAL->AddEntry(lamass_LLAL,"Data");
legLL_LLAL->AddEntry(func2_LLAL,"Signal+Background","L");
legLL_LLAL->AddEntry(func2_bckgnd_LLAL,"Background","L");
legLL_LLAL->AddEntry(linex,"Peak Region","L");
legLL_LLAL->SetBorderSize(0);
legLL_LLAL->Draw();
c2_LLAL->Update();
pad2_LLAL->cd();
gPad->SetTickx();
gPad->SetTicky();
double low_LLAL = -49.1;;
double high_LLAL = 49.1;;
TH1F* frame_la_LLAL = pad2_LLAL->cd()->DrawFrame(full_range_la.first,low_LLAL,full_range_la.second,high_LLAL);
frame_la_LLAL->GetXaxis()->SetTitle("#scale[4.3]{p #pi^{-} invariant mass [GeV]}");
frame_la_LLAL->GetYaxis()->SetTitle("#scale[4.0]{Pull}");
frame_la_LLAL->GetXaxis()->CenterTitle(1);
frame_la_LLAL->GetYaxis()->CenterTitle(1);
frame_la_LLAL->GetXaxis()->SetNdivisions(510);
frame_la_LLAL->GetYaxis()->SetNdivisions(406);
frame_la_LLAL->GetXaxis()->SetTickSize(0.06);
frame_la_LLAL->GetYaxis()->SetTickSize(0.02);
frame_la_LLAL->GetXaxis()->SetTitleOffset(4.5);
frame_la_LLAL->GetYaxis()->SetTitleOffset(1.8);
frame_la_LLAL->GetXaxis()->SetLabelSize(0.14);
frame_la_LLAL->GetYaxis()->SetLabelSize(0.14);
frame_la_LLAL->GetXaxis()->SetLabelOffset(0.01);
frame_la_LLAL->GetYaxis()->SetLabelOffset(0.013);
TGpull_LLAL->SetMarkerStyle(25);
TGpull_LLAL->SetMarkerSize(1.3);
TGpull_LLAL->SetMarkerColor(kBlack);
TGpull_LLAL->Draw("P");
TLine* line_LLAL = new TLine(full_range_la.first, 0, full_range_la.second, 0);
line_LLAL->SetLineStyle(2);
line_LLAL->Draw("same");
pad3_LLAL->cd();
gPad->SetTickx();
gPad->SetTicky();
double low2_LLAL = 0.1;
double high2_LLAL = 1.8999999;
frame_la_LLAL = pad3_LLAL->cd()->DrawFrame(full_range_la.first,low2_LLAL,full_range_la.second,high2_LLAL);
frame_la_LLAL->GetYaxis()->SetTitle("#scale[6.3]{#frac{Data}{Fit}}");
frame_la_LLAL->GetXaxis()->CenterTitle(1);
frame_la_LLAL->GetYaxis()->CenterTitle(1);
frame_la_LLAL->GetXaxis()->SetNdivisions(510);
frame_la_LLAL->GetYaxis()->SetNdivisions(406);
frame_la_LLAL->GetXaxis()->SetTickSize(0.1);
frame_la_LLAL->GetYaxis()->SetTickSize(0.01);
frame_la_LLAL->GetYaxis()->SetTitleOffset(1.8);
frame_la_LLAL->GetXaxis()->SetLabelSize(0.12);
frame_la_LLAL->GetYaxis()->SetLabelSize(0.22);
frame_la_LLAL->GetYaxis()->SetLabelOffset(0.01);
TGratio_LLAL->SetMarkerStyle(25);
TGratio_LLAL->SetMarkerSize(1.3);
TGratio_LLAL->SetMarkerColor(kBlack);
TGratio_LLAL->Draw("P");
TLine* line2_LLAL = new TLine(full_range_la.first, 1, full_range_la.second, 1);
line2_LLAL->SetLineStyle(2);
line2_LLAL->Draw("same");
c2_LLAL->Update();
c2_LLAL->Print("results/Mass"+V0+"LLAL_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");


//ALam used in ALamALam

cout << "==================ALam Mass for LAL==================" << endl;
TH1D* masslam_ALAL = (TH1D*)h_ALam_Mass_LAL_1D->Clone();
//masslam_ALAL->Sumw2();
TF1 *func2_ALAL = new TF1("func2_ALAL", massfunc.c_str(), full_range_la.first, full_range_la.second);

//sideband fit
int cnt2_ALAL = 0;
for ( int i = 0; i < masslam_ALAL->GetNbinsX(); i++ ) {
    double center = masslam_ALAL->GetBinCenter(i+1);
    if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
       cnt2_ALAL++;
    }
}
TGraphErrors gr2_ALAL(cnt2_ALAL);
cnt2_ALAL = 0;
for ( int i = 0; i < masslam_ALAL->GetNbinsX(); i++ ) {
    double center = masslam_ALAL->GetBinCenter(i+1);
    if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
       gr2_ALAL.GetX()[cnt2_ALAL] = center;
       gr2_ALAL.GetY()[cnt2_ALAL] = masslam_ALAL->GetBinContent(i+1);
       gr2_ALAL.GetEY()[cnt2_ALAL] = masslam_ALAL->GetBinError(i+1);
       cnt2_ALAL++;
	}
}

//sideband fit parameters
func2_ALAL->SetParameter(5, 10.);
func2_ALAL->SetParameter(6, 1.);
func2_ALAL->SetParameter(7, 1.);
func2_ALAL->SetParameter(8, 1.);
func2_ALAL->SetParameter(9, 1.);

//fix the parameters in order to get only bkg
func2_ALAL->FixParameter(0, mass_la);
func2_ALAL->FixParameter(1, sigma_la);
func2_ALAL->FixParameter(2, sigma_la);
func2_ALAL->FixParameter(3, 0.);
func2_ALAL->FixParameter(4, 0.);

//make sideband fit
gr2_ALAL.Fit(func2_ALAL, "q0");
gr2_ALAL.Fit(func2_ALAL, "q0");
gr2_ALAL.Fit(func2_ALAL, "q0");

//Set sideband parameters to full fit
func2_ALAL->FixParameter(5, func2_ALAL->GetParameter(5));
func2_ALAL->FixParameter(6, func2_ALAL->GetParameter(6));
func2_ALAL->FixParameter(7, func2_ALAL->GetParameter(7));
func2_ALAL->FixParameter(8, func2_ALAL->GetParameter(8));
func2_ALAL->FixParameter(9, func2_ALAL->GetParameter(9));

//release the signal parameters
func2_ALAL->ReleaseParameter(0);
func2_ALAL->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2_ALAL->ReleaseParameter(1);
func2_ALAL->ReleaseParameter(2);
func2_ALAL->ReleaseParameter(3);
func2_ALAL->ReleaseParameter(4);
func2_ALAL->SetParLimits(1, 0, sigma_la);
func2_ALAL->SetParLimits(2, 0, sigma_la);
func2_ALAL->SetParLimits(3, 0, 0.8*masslam_ALAL->GetBinContent(masslam_ALAL->FindBin(mass_la)));
func2_ALAL->SetParLimits(4, 0, 0.8*masslam_ALAL->GetBinContent(masslam_ALAL->FindBin(mass_la)));

//make a full fit
masslam_ALAL->Fit(func2_ALAL, "q0", "", full_range_la.first, full_range_la.second);
masslam_ALAL->Fit(func2_ALAL, "q0", "", full_range_la.first, full_range_la.second);
masslam_ALAL->Fit(func2_ALAL, "q0", "", full_range_la.first, full_range_la.second);

//now release all the parameters
func2_ALAL->ReleaseParameter(0);
func2_ALAL->ReleaseParameter(1);
func2_ALAL->ReleaseParameter(2);
func2_ALAL->ReleaseParameter(3);
func2_ALAL->ReleaseParameter(4);
func2_ALAL->ReleaseParameter(5);
func2_ALAL->ReleaseParameter(6);
func2_ALAL->ReleaseParameter(7);
func2_ALAL->ReleaseParameter(8);
func2_ALAL->ReleaseParameter(9);
func2_ALAL->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2_ALAL->SetParLimits(1, 0, sigma_la);
func2_ALAL->SetParLimits(2, 0, sigma_la);
func2_ALAL->SetParLimits(3, 0, 0.8*masslam_ALAL->GetBinContent(masslam_ALAL->FindBin(mass_la)));
func2_ALAL->SetParLimits(4, 0, 0.8*masslam_ALAL->GetBinContent(masslam_ALAL->FindBin(mass_la)));

//make final fits
for ( int i = 0; i < iter_la-1; i++ ) {
    masslam_ALAL->Fit(func2_ALAL, "q0", "", full_range_la.first, full_range_la.second);
}

//make signal and bkg functions for plots
TF1 * func2_signal_ALAL   = 0;
TF1 * func2_bckgnd_ALAL   = 0;
TF1 * gaus1_signal_la_ALAL   = 0;
TF1 * gaus2_signal_la_ALAL   = 0;
func2_signal_ALAL  = new TF1("func2_signal_ALAL" , double_gaussian.c_str(), full_range_la.first, full_range_la.second);
func2_bckgnd_ALAL = new TF1("func2_bckgnd_ALAL" , poly_bkg0.c_str(), full_range_la.first, full_range_la.second);
gaus1_signal_la_ALAL = new TF1("gaus1_signal_la_ALAL", gaus_func.c_str(), full_range_la.first, full_range_la.second);
gaus2_signal_la_ALAL = new TF1("gaus2_signal_la_ALAL", gaus_func.c_str(), full_range_la.first, full_range_la.second);

//set parameters for signal and bkg
func2_signal_ALAL->SetParameter(0, func2_ALAL->GetParameter(0));
func2_signal_ALAL->SetParameter(1, func2_ALAL->GetParameter(1));
func2_signal_ALAL->SetParameter(2, func2_ALAL->GetParameter(2));
func2_signal_ALAL->SetParameter(3, func2_ALAL->GetParameter(3));
func2_signal_ALAL->SetParameter(4, func2_ALAL->GetParameter(4));
func2_bckgnd_ALAL->SetParameter(0, func2_ALAL->GetParameter(5));
func2_bckgnd_ALAL->SetParameter(1, func2_ALAL->GetParameter(6));
func2_bckgnd_ALAL->SetParameter(2, func2_ALAL->GetParameter(7));
func2_bckgnd_ALAL->SetParameter(3, func2_ALAL->GetParameter(8));
func2_bckgnd_ALAL->SetParameter(4, func2_ALAL->GetParameter(9));

//set parameters for gaussians
gaus1_signal_la_ALAL->SetParameter(0, func2_ALAL->GetParameter(0));
gaus1_signal_la_ALAL->SetParameter(1, func2_ALAL->GetParameter(1));
gaus1_signal_la_ALAL->SetParameter(2, func2_ALAL->GetParameter(3));
gaus2_signal_la_ALAL->SetParameter(0, func2_ALAL->GetParameter(0));
gaus2_signal_la_ALAL->SetParameter(1, func2_ALAL->GetParameter(2));
gaus2_signal_la_ALAL->SetParameter(2, func2_ALAL->GetParameter(4));

//make pull and ratio
double vny2_ALAL[400];
double xvn2_ALAL[400];
double pullvn2_ALAL[400];
double ratiovn2_ALAL[400];
double err_pullvn2_ALAL[400];
double err_ratiovn2_ALAL[400];
for(int k=0;k<masslam_ALAL->GetNbinsX();k++)
       	    {
            vny2_ALAL[k] = masslam_ALAL->GetBinContent(k+1);
            xvn2_ALAL[k] = masslam_ALAL->GetBinCenter(k+1);
        	if(masslam_ALAL->GetBinContent(k+1)!=0){
            pullvn2_ALAL[k] = (vny2_ALAL[k] - func2_ALAL->Eval(xvn2_ALAL[k]))/masslam_ALAL->GetBinError(k+1);
            ratiovn2_ALAL[k] = vny2_ALAL[k]/func2_ALAL->Eval(xvn2_ALAL[k]);
            err_pullvn2_ALAL[k] = masslam_ALAL->GetBinError(k+1);
            err_ratiovn2_ALAL[k] = (vny2_ALAL[k]/func2_ALAL->Eval(xvn2_ALAL[k]))*(masslam_ALAL->GetBinError(k+1)/vny2_ALAL[k]);
            }else{
            pullvn2_ALAL[k] = 0.0;
            ratiovn2_ALAL[k] = 0.0;
            err_pullvn2_ALAL[k]=0.0;
            err_ratiovn2_ALAL[k]=0.0;
            }
}

//make the plots
TH1D* lamass_ALAL = (TH1D*)masslam_ALAL->Clone();
//lamass_ALAL->Sumw2();
TGraphErrors* TGpull_ALAL = new TGraphErrors(400,xvn2_ALAL,pullvn2_ALAL,0,0);
TGraphErrors* TGratio_ALAL= new TGraphErrors(400,xvn2_ALAL,ratiovn2_ALAL,0,err_ratiovn2_ALAL);
TCanvas* c2_ALAL = new TCanvas("c2_ALAL","c2_ALAL",600,600);
c2_ALAL->cd();
TPad* pad1_ALAL = new TPad("pad1_ALAL","top pad",0.0,0.4,1.0,1.0);
TPad* pad2_ALAL = new TPad("pad3_ALAL","middle pad",0.0,0.0,1.0,0.25);
TPad* pad3_ALAL = new TPad("pad2_ALAL","bottom pad",0.0,0.25,1.0,0.4);
pad1_ALAL->SetTopMargin(0.1);
pad1_ALAL->SetBottomMargin(0.0);
pad1_ALAL->SetRightMargin(0.038);
pad1_ALAL->SetLeftMargin(0.15);
pad1_ALAL->Draw();
pad2_ALAL->SetFrameFillStyle(4000);
pad2_ALAL->SetTopMargin(0.0);
pad2_ALAL->SetBottomMargin(0.37);
pad2_ALAL->SetRightMargin(0.038);
pad2_ALAL->SetLeftMargin(0.15);
pad2_ALAL->Draw();
pad3_ALAL->SetFrameFillStyle(4000);
pad3_ALAL->SetTopMargin(0.0);
pad3_ALAL->SetBottomMargin(0.0);
pad3_ALAL->SetRightMargin(0.038);
pad3_ALAL->SetLeftMargin(0.15);
pad3_ALAL->Draw();
pad1_ALAL->cd();
gPad->SetTickx();
gPad->SetTicky();
lamass_ALAL->SetTitle("");
lamass_ALAL->Draw("E");
TGaxis::SetMaxDigits(3);
lamass_ALAL->GetXaxis()->SetRangeUser(full_range_la.first,full_range_la.second);
lamass_ALAL->GetYaxis()->SetRangeUser(0.1,1.05*lamass_ALAL->GetMaximum());
lamass_ALAL->SetMarkerStyle(24);
lamass_ALAL->SetMarkerSize(1.5);
lamass_ALAL->GetYaxis()->SetTitle("#scale[1.7]{Candidates / 0.5 MeV}");
lamass_ALAL->GetYaxis()->CenterTitle(1);
lamass_ALAL->GetYaxis()->CenterTitle(1);
lamass_ALAL->GetYaxis()->SetNdivisions(810);
lamass_ALAL->GetXaxis()->SetTickSize(0.03);
lamass_ALAL->GetYaxis()->SetTickSize(0.02);
lamass_ALAL->GetXaxis()->SetLabelSize(0.055);
lamass_ALAL->GetYaxis()->SetLabelSize(0.058);
lamass_ALAL->GetYaxis()->SetLabelOffset(0.01);
lamass_ALAL->GetYaxis()->SetTitleOffset(1.8);
func2_ALAL->Draw("LSAME");
func2_ALAL->SetLineColor(SPcolors[SPdred]);
func2_bckgnd_ALAL->Draw("LSAME");
func2_bckgnd_ALAL->SetLineColor(SPcolors[SPdred]);
func2_bckgnd_ALAL->SetLineStyle(2);
double meanla_ALAL  = func2_ALAL->GetParameter(0);
double sigma1la_ALAL = func2_ALAL->GetParameter(1);
double sigma2la_ALAL = func2_ALAL->GetParameter(2);
double gaus1_yield_la_ALAL = gaus1_signal_la_ALAL->Integral(meanla_ALAL - 2*sigma1la_ALAL, meanla_ALAL + 2*sigma1la)/masslam_ALAL->GetBinWidth(1);
double gaus2_yield_la_ALAL = gaus2_signal_la_ALAL->Integral(meanla_ALAL - 2*sigma2la_ALAL, meanla_ALAL + 2*sigma2la)/masslam_ALAL->GetBinWidth(1);
double gausTot_yield_la_ALAL = gaus1_yield_la_ALAL + gaus2_yield_la_ALAL;
double rms_gaus1_sig_la_ALAL = gaus1_yield_la_ALAL/gausTot_yield_la_ALAL;
double rms_gaus2_sig_la_ALAL = gaus2_yield_la_ALAL/gausTot_yield_la_ALAL;
double rms_true_la_ALAL = TMath::Sqrt(rms_gaus1_sig_la_ALAL*sigma1la_ALAL*sigma1la_ALAL + rms_gaus2_sig_la_ALAL*sigma2la_ALAL*sigma2la_ALAL);
std::pair<double, double>  signal_range_la_ALAL(meanla_ALAL - 2*rms_true_la_ALAL, meanla_ALAL + 2*rms_true_la_ALAL);
double S_la_ALAL = func2_signal_ALAL->Integral(signal_range_la_ALAL.first, signal_range_la_ALAL.second)/masslam_ALAL->GetBinWidth(1);
double B_la_ALAL = func2_bckgnd_ALAL->Integral(signal_range_la_ALAL.first, signal_range_la_ALAL.second)/masslam_ALAL->GetBinWidth(1);
double Fsig_la_ALAL = S_la_ALAL/(S_la_ALAL+B_la_ALAL);
double sig_la_ALAL = S_la_ALAL/sqrt(S_la_ALAL+B_la_ALAL);
double chi2la_ALAL = func2_ALAL->GetChisquare()/func2_ALAL->GetNDF();
double probla_ALAL = TMath::Prob(func2_ALAL->GetChisquare(),func2_ALAL->GetNDF());
std::ostringstream osx_ALAL; 
TLatex* texCMSx_ALAL = new TLatex;
texCMSx_ALAL->SetNDC();
texCMSx_ALAL->SetTextFont(42);
texCMSx_ALAL->SetTextSize(0.05);
texCMSx_ALAL->SetTextAlign(12);
texCMSx_ALAL->DrawLatex(0.25,.95," #font[61]{CMS} #it{Work in Progress}");
texCMSx_ALAL->DrawLatex(0.70,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMSx_ALAL->DrawLatex(0.26,0.4, "#scale[1.7]{#bar{#Lambda}}");
texCMSx_ALAL->DrawLatex(0.225,0.2, "#scale[1.0]{Used in #Lambda#bar{#Lambda}}"); 
texCMSx_ALAL->DrawLatex(0.225,0.15, "#scale[1.0]{correlations}"); 
// texCMSx_ALAL->DrawLatex(0.75,0.15, "This Thesis"); 
osx_ALAL << "#scale[1.0]{s/#sqrt{s+b}: " << std::setprecision(4) << sig_la_ALAL << "}" << std::setprecision(2);           
texCMSx_ALAL->DrawLatex(0.68,0.60, osx_ALAL.str(  ).c_str(  )); 
osx_ALAL.str(std::string());
osx_ALAL << "#scale[1.0]{s/(s+b): " << std::setprecision(4) << Fsig_la_ALAL*100 << " %}" << std::setprecision(2); 
texCMSx_ALAL->DrawLatex(0.68,0.52, osx_ALAL.str(  ).c_str(  ));  
osx_ALAL.str(std::string()); 
osx_ALAL << "#scale[1.0]{Mean: " << std::setprecision(4) << meanla_ALAL << " GeV}" << std::setprecision(6);   
texCMSx_ALAL->DrawLatex(0.68,0.46, osx_ALAL.str(  ).c_str(  ));
osx_ALAL.str(std::string()); 
osx_ALAL << "#scale[1.0]{#sigma: " << std::setprecision(4) << rms_true_la_ALAL << " }" << std::setprecision(2); 
texCMSx_ALAL->DrawLatex(0.68,0.40, osx_ALAL.str(  ).c_str(  )); 
osx_ALAL.str(std::string()); 
osx_ALAL << "#scale[1.0]{Yield: " << std::setprecision(4) << S_la_ALAL << " }" << std::setprecision(2); 
texCMSx_ALAL->DrawLatex(0.68,0.33, osx_ALAL.str(  ).c_str(  )); 
osx_ALAL.str(std::string()); 
osx_ALAL << "#scale[1.0]{#chi^{2}/ndf: " << std::setprecision(4) << chi2la_ALAL << " }" << std::setprecision(4); 
texCMSx_ALAL->DrawLatex(0.68,0.26, osx_ALAL.str(  ).c_str(  ));  
osx_ALAL.str(std::string()); 
// osx_ALAL << "#scale[1.0]{Prob: " << std::setprecision(4) << probla_ALAL << "}" << std::setprecision(2); 
// texCMSx_ALAL->DrawLatex(0.68,0.20, osx_ALAL.str(  ).c_str(  )); 
// osx_ALAL.str(std::string()); 
osx_ALAL << "#scale[1.]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMSx_ALAL->DrawLatex(0.68,0.75, osx_ALAL.str(  ).c_str(  )); 
osx_ALAL.str(std::string()); 
osx_ALAL << "#scale[1.]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMSx_ALAL->DrawLatex(0.68,0.67, osx_ALAL.str(  ).c_str(  ));  
osx_ALAL.str(std::string()); 
osx_ALAL << "#scale[1.]{ " << std::setprecision(4) << qminM << std::setprecision(4) <<" < q_{inv} < "  << qmaxM << " GeV}" << std::setprecision(2); 
texCMSx_ALAL->DrawLatex(0.68,0.83, osx_ALAL.str(  ).c_str(  ));  
osx_ALAL.str(std::string()); 
Float_t yymax_ALAL = lamass_ALAL->GetMaximum();
TLine *linex_ALAL = new TLine(signal_range_la_ALAL.first,0.,signal_range_la_ALAL.first,yymax_ALAL);
linex_ALAL->SetLineColor(SPcolors[SPdblue]);
linex_ALAL->Draw();
TLine *line2x_ALAL = new TLine(signal_range_la_ALAL.second,0.,signal_range_la_ALAL.second,yymax_ALAL);
line2x_ALAL->SetLineColor(SPcolors[SPdblue]);
line2x_ALAL->Draw();
TLegend* legLL_ALAL = new TLegend(0.16,0.6,0.43,0.85);
legLL_ALAL->AddEntry(lamass_ALAL,"Data");
legLL_ALAL->AddEntry(func2_ALAL,"Signal+Background","L");
legLL_ALAL->AddEntry(func2_bckgnd_ALAL,"Background","L");
legLL_ALAL->AddEntry(linex,"Peak Region","L");
legLL_ALAL->SetBorderSize(0);
legLL_ALAL->Draw();
c2_ALAL->Update();
pad2_ALAL->cd();
gPad->SetTickx();
gPad->SetTicky();
double low_ALAL = -49.1;;
double high_ALAL = 49.1;;
TH1F* frame_la_ALAL = pad2_ALAL->cd()->DrawFrame(full_range_la.first,low_ALAL,full_range_la.second,high_ALAL);
frame_la_ALAL->GetXaxis()->SetTitle("#scale[4.3]{#bar{p} #pi^{+} invariant mass [GeV]}");
frame_la_ALAL->GetYaxis()->SetTitle("#scale[4.0]{Pull}");
frame_la_ALAL->GetXaxis()->CenterTitle(1);
frame_la_ALAL->GetYaxis()->CenterTitle(1);
frame_la_ALAL->GetXaxis()->SetNdivisions(510);
frame_la_ALAL->GetYaxis()->SetNdivisions(406);
frame_la_ALAL->GetXaxis()->SetTickSize(0.06);
frame_la_ALAL->GetYaxis()->SetTickSize(0.02);
frame_la_ALAL->GetXaxis()->SetTitleOffset(4.5);
frame_la_ALAL->GetYaxis()->SetTitleOffset(1.8);
frame_la_ALAL->GetXaxis()->SetLabelSize(0.14);
frame_la_ALAL->GetYaxis()->SetLabelSize(0.14);
frame_la_ALAL->GetXaxis()->SetLabelOffset(0.01);
frame_la_ALAL->GetYaxis()->SetLabelOffset(0.013);
TGpull_ALAL->SetMarkerStyle(25);
TGpull_ALAL->SetMarkerSize(1.3);
TGpull_ALAL->SetMarkerColor(kBlack);
TGpull_ALAL->Draw("P");
TLine* line_ALAL = new TLine(full_range_la.first, 0, full_range_la.second, 0);
line_ALAL->SetLineStyle(2);
line_ALAL->Draw("same");
pad3_ALAL->cd();
gPad->SetTickx();
gPad->SetTicky();
double low2_ALAL = 0.1;
double high2_ALAL = 1.8999999;
frame_la_ALAL = pad3_ALAL->cd()->DrawFrame(full_range_la.first,low2_ALAL,full_range_la.second,high2_ALAL);
frame_la_ALAL->GetYaxis()->SetTitle("#scale[6.3]{#frac{Data}{Fit}}");
frame_la_ALAL->GetXaxis()->CenterTitle(1);
frame_la_ALAL->GetYaxis()->CenterTitle(1);
frame_la_ALAL->GetXaxis()->SetNdivisions(510);
frame_la_ALAL->GetYaxis()->SetNdivisions(406);
frame_la_ALAL->GetXaxis()->SetTickSize(0.1);
frame_la_ALAL->GetYaxis()->SetTickSize(0.01);
frame_la_ALAL->GetYaxis()->SetTitleOffset(1.8);
frame_la_ALAL->GetXaxis()->SetLabelSize(0.12);
frame_la_ALAL->GetYaxis()->SetLabelSize(0.22);
frame_la_ALAL->GetYaxis()->SetLabelOffset(0.01);
TGratio_ALAL->SetMarkerStyle(25);
TGratio_ALAL->SetMarkerSize(1.3);
TGratio_ALAL->SetMarkerColor(kBlack);
TGratio_ALAL->Draw("P");
TLine* line2_ALAL = new TLine(full_range_la.first, 1, full_range_la.second, 1);
line2_ALAL->SetLineStyle(2);
line2_ALAL->Draw("same");
c2_ALAL->Update();
c2_ALAL->Print("results/Mass"+V0+"ALLAL_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");


//Lam used in K0sLam
cout << "==================Lam Mass for K0sL==================" << endl;
TH1D* masslam_LKs = (TH1D*)h_LK_Mass_1D->Clone();
//masslam_LKs->Sumw2();
TF1 *func2_LKs = new TF1("func2_LKs", massfunc.c_str(), full_range_la.first, full_range_la.second);
//sideband fit
int cnt2_LKs = 0;
for ( int i = 0; i < masslam_LKs->GetNbinsX(); i++ ) {
    double center = masslam_LKs->GetBinCenter(i+1);
    if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
       cnt2_LKs++;
    }
}
TGraphErrors gr2_LKs(cnt2_LKs);
cnt2_LKs = 0;
for ( int i = 0; i < masslam_LKs->GetNbinsX(); i++ ) {
    double center = masslam_LKs->GetBinCenter(i+1);
    if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
       gr2_LKs.GetX()[cnt2_LKs] = center;
       gr2_LKs.GetY()[cnt2_LKs] = masslam_LKs->GetBinContent(i+1);
       gr2_LKs.GetEY()[cnt2_LKs] = masslam_LKs->GetBinError(i+1);
       cnt2_LKs++;
    }
}

//sideband fit parameters
func2_LKs->SetParameter(5, 10.);
func2_LKs->SetParameter(6, 1.);
func2_LKs->SetParameter(7, 1.);
func2_LKs->SetParameter(8, 1.);
func2_LKs->SetParameter(9, 1.);

//fix the parameters in order to get only bkg
func2_LKs->FixParameter(0, mass_la);
func2_LKs->FixParameter(1, sigma_la);
func2_LKs->FixParameter(2, sigma_la);
func2_LKs->FixParameter(3, 0.);
func2_LKs->FixParameter(4, 0.);

//make sideband fit
gr2_LKs.Fit(func2_LKs, "q0");
gr2_LKs.Fit(func2_LKs, "q0");
gr2_LKs.Fit(func2_LKs, "q0");

//Set sideband parameters to full fit
func2_LKs->FixParameter(5, func2_LKs->GetParameter(5));
func2_LKs->FixParameter(6, func2_LKs->GetParameter(6));
func2_LKs->FixParameter(7, func2_LKs->GetParameter(7));
func2_LKs->FixParameter(8, func2_LKs->GetParameter(8));
func2_LKs->FixParameter(9, func2_LKs->GetParameter(9));

//release the signal parameters
func2_LKs->ReleaseParameter(0);
func2_LKs->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2_LKs->ReleaseParameter(1);
func2_LKs->ReleaseParameter(2);
func2_LKs->ReleaseParameter(3);
func2_LKs->ReleaseParameter(4);
func2_LKs->SetParLimits(1, 0, sigma_la);
func2_LKs->SetParLimits(2, 0, sigma_la);
func2_LKs->SetParLimits(3, 0, 0.8*masslam_LKs->GetBinContent(masslam_LKs->FindBin(mass_la)));
func2_LKs->SetParLimits(4, 0, 0.8*masslam_LKs->GetBinContent(masslam_LKs->FindBin(mass_la)));

//make a full fit
masslam_LKs->Fit(func2_LKs, "q0", "", full_range_la.first, full_range_la.second);
masslam_LKs->Fit(func2_LKs, "q0", "", full_range_la.first, full_range_la.second);
masslam_LKs->Fit(func2_LKs, "q0", "", full_range_la.first, full_range_la.second);

//now release all the parameters
func2_LKs->ReleaseParameter(0);
func2_LKs->ReleaseParameter(1);
func2_LKs->ReleaseParameter(2);
func2_LKs->ReleaseParameter(3);
func2_LKs->ReleaseParameter(4);
func2_LKs->ReleaseParameter(5);
func2_LKs->ReleaseParameter(6);
func2_LKs->ReleaseParameter(7);
func2_LKs->ReleaseParameter(8);
func2_LKs->ReleaseParameter(9);
func2_LKs->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2_LKs->SetParLimits(1, 0, sigma_la);
func2_LKs->SetParLimits(2, 0, sigma_la);
func2_LKs->SetParLimits(3, 0, 0.8*masslam_LKs->GetBinContent(masslam_LKs->FindBin(mass_la)));
func2_LKs->SetParLimits(4, 0, 0.8*masslam_LKs->GetBinContent(masslam_LKs->FindBin(mass_la)));

//make final fits
for ( int i = 0; i < iter_la-1; i++ ) {
    masslam_LKs->Fit(func2_LKs, "q0", "", full_range_la.first, full_range_la.second);
}

//make signal and bkg functions for plots
TF1 * func2_signal_LKs   = 0;
TF1 * func2_bckgnd_LKs   = 0;
TF1 * gaus1_signal_la_LKs   = 0;
TF1 * gaus2_signal_la_LKs   = 0;
func2_signal_LKs  = new TF1("func2_signal_LKs" , double_gaussian.c_str(), full_range_la.first, full_range_la.second);
func2_bckgnd_LKs = new TF1("func2_bckgnd_LKs" , poly_bkg0.c_str(), full_range_la.first, full_range_la.second);
gaus1_signal_la_LKs = new TF1("gaus1_signal_la_LKs", gaus_func.c_str(), full_range_la.first, full_range_la.second);
gaus2_signal_la_LKs = new TF1("gaus2_signal_la_LKs", gaus_func.c_str(), full_range_la.first, full_range_la.second);

//set parameters for signal and bkg
func2_signal_LKs->SetParameter(0, func2_LKs->GetParameter(0));
func2_signal_LKs->SetParameter(1, func2_LKs->GetParameter(1));
func2_signal_LKs->SetParameter(2, func2_LKs->GetParameter(2));
func2_signal_LKs->SetParameter(3, func2_LKs->GetParameter(3));
func2_signal_LKs->SetParameter(4, func2_LKs->GetParameter(4));
func2_bckgnd_LKs->SetParameter(0, func2_LKs->GetParameter(5));
func2_bckgnd_LKs->SetParameter(1, func2_LKs->GetParameter(6));
func2_bckgnd_LKs->SetParameter(2, func2_LKs->GetParameter(7));
func2_bckgnd_LKs->SetParameter(3, func2_LKs->GetParameter(8));
func2_bckgnd_LKs->SetParameter(4, func2_LKs->GetParameter(9));

//set parameters for gaussians
gaus1_signal_la_LKs->SetParameter(0, func2_LKs->GetParameter(0));
gaus1_signal_la_LKs->SetParameter(1, func2_LKs->GetParameter(1));
gaus1_signal_la_LKs->SetParameter(2, func2_LKs->GetParameter(3));
gaus2_signal_la_LKs->SetParameter(0, func2_LKs->GetParameter(0));
gaus2_signal_la_LKs->SetParameter(1, func2_LKs->GetParameter(2));
gaus2_signal_la_LKs->SetParameter(2, func2_LKs->GetParameter(4));

//make pull and ratio
double vny2_LKs[400];
double xvn2_LKs[400];
double pullvn2_LKs[400];
double ratiovn2_LKs[400];
double err_pullvn2_LKs[400];
double err_ratiovn2_LKs[400];
for(int k=0;k<masslam_LKs->GetNbinsX();k++)
        	{
            vny2_LKs[k] = masslam_LKs->GetBinContent(k+1);
            xvn2_LKs[k] = masslam_LKs->GetBinCenter(k+1);
        	if(masslam_LKs->GetBinContent(k+1)!=0){
            pullvn2_LKs[k] = (vny2_LKs[k] - func2_LKs->Eval(xvn2_LKs[k]))/masslam_LKs->GetBinError(k+1);
            ratiovn2_LKs[k] = vny2_LKs[k]/func2_LKs->Eval(xvn2_LKs[k]);
            err_pullvn2_LKs[k] = masslam_LKs->GetBinError(k+1);
            err_ratiovn2_LKs[k] = (vny2_LKs[k]/func2_LKs->Eval(xvn2_LKs[k]))*(masslam_LKs->GetBinError(k+1)/vny2_LKs[k]);
            }else{
            pullvn2_LKs[k] = 0.0;
            ratiovn2_LKs[k] = 0.0;
            err_pullvn2_LKs[k]=0.0;
            err_ratiovn2_LKs[k]=0.0;
            }
}

//make the plots
TH1D* lamass_LKs = (TH1D*)masslam_LKs->Clone();
//lamass_LKs->Sumw2();
TGraphErrors* TGpull_LKs = new TGraphErrors(400,xvn2_LKs,pullvn2_LKs,0,0);
TGraphErrors* TGratio_LKs= new TGraphErrors(400,xvn2_LKs,ratiovn2_LKs,0,err_ratiovn2_LKs);
TCanvas* c2_LKs = new TCanvas("c2_LKs","c2_LKs",600,600);
c2_LKs->cd();
TPad* pad1_LKs = new TPad("pad1_LKs","top pad",0.0,0.4,1.0,1.0);
TPad* pad2_LKs = new TPad("pad3_LKs","middle pad",0.0,0.0,1.0,0.25);
TPad* pad3_LKs = new TPad("pad2_LKs","bottom pad",0.0,0.25,1.0,0.4);
pad1_LKs->SetTopMargin(0.1);
pad1_LKs->SetBottomMargin(0.0);
pad1_LKs->SetRightMargin(0.038);
pad1_LKs->SetLeftMargin(0.15);
pad1_LKs->Draw();
pad2_LKs->SetFrameFillStyle(4000);
pad2_LKs->SetTopMargin(0.0);
pad2_LKs->SetBottomMargin(0.37);
pad2_LKs->SetRightMargin(0.038);
pad2_LKs->SetLeftMargin(0.15);
pad2_LKs->Draw();
pad3_LKs->SetFrameFillStyle(4000);
pad3_LKs->SetTopMargin(0.0);
pad3_LKs->SetBottomMargin(0.0);
pad3_LKs->SetRightMargin(0.038);
pad3_LKs->SetLeftMargin(0.15);
pad3_LKs->Draw();
pad1_LKs->cd();
gPad->SetTickx();
gPad->SetTicky();
lamass_LKs->SetTitle("");
lamass_LKs->Draw("E");
TGaxis::SetMaxDigits(3);
lamass_LKs->GetXaxis()->SetRangeUser(full_range_la.first,full_range_la.second);
lamass_LKs->GetYaxis()->SetRangeUser(0.1,1.05*lamass_LKs->GetMaximum());
lamass_LKs->SetMarkerStyle(24);
lamass_LKs->SetMarkerSize(1.5);
lamass_LKs->GetYaxis()->SetTitle("#scale[1.7]{Candidates / 0.5 MeV}");
lamass_LKs->GetYaxis()->CenterTitle(1);
lamass_LKs->GetYaxis()->CenterTitle(1);
lamass_LKs->GetYaxis()->SetNdivisions(810);
lamass_LKs->GetXaxis()->SetTickSize(0.03);
lamass_LKs->GetYaxis()->SetTickSize(0.02);
lamass_LKs->GetXaxis()->SetLabelSize(0.055);
lamass_LKs->GetYaxis()->SetLabelSize(0.058);
lamass_LKs->GetYaxis()->SetLabelOffset(0.01);
lamass_LKs->GetYaxis()->SetTitleOffset(1.8);
func2_LKs->Draw("LSAME");
func2_LKs->SetLineColor(SPcolors[SPdred]);
func2_bckgnd_LKs->Draw("LSAME");
func2_bckgnd_LKs->SetLineColor(SPcolors[SPdred]);
func2_bckgnd_LKs->SetLineStyle(2);
double meanla_LKs  = func2_LKs->GetParameter(0);
double sigma1la_LKs = func2_LKs->GetParameter(1);
double sigma2la_LKs = func2_LKs->GetParameter(2);
double gaus1_yield_la_LKs = gaus1_signal_la_LKs->Integral(meanla_LKs - 2*sigma1la_LKs, meanla_LKs + 2*sigma1la)/masslam_LKs->GetBinWidth(1);
double gaus2_yield_la_LKs = gaus2_signal_la_LKs->Integral(meanla_LKs - 2*sigma2la_LKs, meanla_LKs + 2*sigma2la)/masslam_LKs->GetBinWidth(1);
double gausTot_yield_la_LKs = gaus1_yield_la_LKs + gaus2_yield_la_LKs;
double rms_gaus1_sig_la_LKs = gaus1_yield_la_LKs/gausTot_yield_la_LKs;
double rms_gaus2_sig_la_LKs = gaus2_yield_la_LKs/gausTot_yield_la_LKs;
double rms_true_la_LKs = TMath::Sqrt(rms_gaus1_sig_la_LKs*sigma1la_LKs*sigma1la_LKs + rms_gaus2_sig_la_LKs*sigma2la_LKs*sigma2la_LKs);
std::pair<double, double>  signal_range_la_LKs(meanla_LKs - 2*rms_true_la_LKs, meanla_LKs + 2*rms_true_la_LKs);
double S_la_LKs = func2_signal_LKs->Integral(signal_range_la_LKs.first, signal_range_la_LKs.second)/masslam_LKs->GetBinWidth(1);
double B_la_LKs = func2_bckgnd_LKs->Integral(signal_range_la_LKs.first, signal_range_la_LKs.second)/masslam_LKs->GetBinWidth(1);
double Fsig_la_LKs = S_la_LKs/(S_la_LKs+B_la_LKs);
double sig_la_LKs = S_la_LKs/sqrt(S_la_LKs+B_la_LKs);
double chi2la_LKs = func2_LKs->GetChisquare()/func2_LKs->GetNDF();
double probla_LKs = TMath::Prob(func2_LKs->GetChisquare(),func2_LKs->GetNDF());
std::ostringstream osx_LKs; 
TLatex* texCMSx_LKs = new TLatex;
texCMSx_LKs->SetNDC();
texCMSx_LKs->SetTextFont(42);
texCMSx_LKs->SetTextSize(0.05);
texCMSx_LKs->SetTextAlign(12);
texCMSx_LKs->DrawLatex(0.25,.95," #font[61]{CMS} #it{Work in Progress}");
texCMSx_LKs->DrawLatex(0.70,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMSx_LKs->DrawLatex(0.26,0.4, "#scale[1.7]{#Lambda}");
texCMSx_LKs->DrawLatex(0.225,0.2, "#scale[1.0]{Used in K^{0}_{S}#Lambda}"); 
texCMSx_LKs->DrawLatex(0.225,0.15, "#scale[1.0]{correlations}"); 
// texCMSx_LKs->DrawLatex(0.75,0.15, "This Thesis"); 
osx_LKs << "#scale[1.0]{s/#sqrt{s+b}: " << std::setprecision(4) << sig_la_LKs << "}" << std::setprecision(2);           
texCMSx_LKs->DrawLatex(0.68,0.60, osx_LKs.str(  ).c_str(  )); 
osx_LKs.str(std::string());
osx_LKs << "#scale[1.0]{s/(s+b): " << std::setprecision(4) << Fsig_la_LKs*100 << " %}" << std::setprecision(2); 
texCMSx_LKs->DrawLatex(0.68,0.52, osx_LKs.str(  ).c_str(  ));  
osx_LKs.str(std::string()); 
osx_LKs << "#scale[1.0]{Mean: " << std::setprecision(4) << meanla_LKs << " GeV}" << std::setprecision(6);   
texCMSx_LKs->DrawLatex(0.68,0.46, osx_LKs.str(  ).c_str(  ));
osx_LKs.str(std::string()); 
osx_LKs << "#scale[1.0]{#sigma: " << std::setprecision(4) << rms_true_la_LKs << " }" << std::setprecision(2); 
texCMSx_LKs->DrawLatex(0.68,0.40, osx_LKs.str(  ).c_str(  )); 
osx_LKs.str(std::string()); 
osx_LKs << "#scale[1.0]{Yield: " << std::setprecision(4) << S_la_LKs << " }" << std::setprecision(2); 
texCMSx_LKs->DrawLatex(0.68,0.33, osx_LKs.str(  ).c_str(  )); 
osx_LKs.str(std::string()); 
osx_LKs << "#scale[1.0]{#chi^{2}/ndf: " << std::setprecision(4) << chi2la_LKs << " }" << std::setprecision(4); 
texCMSx_LKs->DrawLatex(0.68,0.26, osx_LKs.str(  ).c_str(  ));  
osx_LKs.str(std::string()); 
// osx_LKs << "#scale[1.0]{Prob: " << std::setprecision(4) << probla_LKs << "}" << std::setprecision(2); 
// texCMSx_LKs->DrawLatex(0.68,0.20, osx_LKs.str(  ).c_str(  )); 
// osx_LKs.str(std::string()); 
osx_LKs << "#scale[1.]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMSx_LKs->DrawLatex(0.68,0.75, osx_LKs.str(  ).c_str(  )); 
osx_LKs.str(std::string()); 
osx_LKs << "#scale[1.]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMSx_LKs->DrawLatex(0.68,0.67, osx_LKs.str(  ).c_str(  ));  
osx_LKs.str(std::string()); 
osx_LKs << "#scale[1.]{ " << std::setprecision(4) << qminM << std::setprecision(4) <<" < q_{inv} < "  << qmaxM << " GeV}" << std::setprecision(2); 
texCMSx_LKs->DrawLatex(0.68,0.83, osx_LKs.str(  ).c_str(  ));  
osx_LKs.str(std::string()); 
Float_t yymax_LKs = lamass_LKs->GetMaximum();
TLine *linex_LKs = new TLine(signal_range_la_LKs.first,0.,signal_range_la_LKs.first,yymax_LKs);
linex_LKs->SetLineColor(SPcolors[SPdblue]);
linex_LKs->Draw();
TLine *line2x_LKs = new TLine(signal_range_la_LKs.second,0.,signal_range_la_LKs.second,yymax_LKs);
line2x_LKs->SetLineColor(SPcolors[SPdblue]);
line2x_LKs->Draw();
TLegend* legLL_LKs = new TLegend(0.16,0.6,0.43,0.85);
legLL_LKs->AddEntry(lamass_LKs,"Data");
legLL_LKs->AddEntry(func2_LKs,"Signal+Background","L");
legLL_LKs->AddEntry(func2_bckgnd_LKs,"Background","L");
legLL_LKs->AddEntry(linex,"Peak Region","L");
legLL_LKs->SetBorderSize(0);
legLL_LKs->Draw();
c2_LKs->Update();
pad2_LKs->cd();
gPad->SetTickx();
gPad->SetTicky();
double low_LKs = -49.1;;
double high_LKs = 49.1;;
TH1F* frame_la_LKs = pad2_LKs->cd()->DrawFrame(full_range_la.first,low_LKs,full_range_la.second,high_LKs);
frame_la_LKs->GetXaxis()->SetTitle("#scale[4.3]{p #pi^{-} invariant mass [GeV]}");
frame_la_LKs->GetYaxis()->SetTitle("#scale[4.0]{Pull}");
frame_la_LKs->GetXaxis()->CenterTitle(1);
frame_la_LKs->GetYaxis()->CenterTitle(1);
frame_la_LKs->GetXaxis()->SetNdivisions(510);
frame_la_LKs->GetYaxis()->SetNdivisions(406);
frame_la_LKs->GetXaxis()->SetTickSize(0.06);
frame_la_LKs->GetYaxis()->SetTickSize(0.02);
frame_la_LKs->GetXaxis()->SetTitleOffset(4.5);
frame_la_LKs->GetYaxis()->SetTitleOffset(1.8);
frame_la_LKs->GetXaxis()->SetLabelSize(0.14);
frame_la_LKs->GetYaxis()->SetLabelSize(0.14);
frame_la_LKs->GetXaxis()->SetLabelOffset(0.01);
frame_la_LKs->GetYaxis()->SetLabelOffset(0.013);
TGpull_LKs->SetMarkerStyle(25);
TGpull_LKs->SetMarkerSize(1.3);
TGpull_LKs->SetMarkerColor(kBlack);
TGpull_LKs->Draw("P");
TLine* line_LKs = new TLine(full_range_la.first, 0, full_range_la.second, 0);
line_LKs->SetLineStyle(2);
line_LKs->Draw("same");
pad3_LKs->cd();
gPad->SetTickx();
gPad->SetTicky();
double low2_LKs = 0.1;
double high2_LKs = 1.89999999;
frame_la_LKs = pad3_LKs->cd()->DrawFrame(full_range_la.first,low2_LKs,full_range_la.second,high2_LKs);
frame_la_LKs->GetYaxis()->SetTitle("#scale[6.3]{#frac{Data}{Fit}}");
frame_la_LKs->GetXaxis()->CenterTitle(1);
frame_la_LKs->GetYaxis()->CenterTitle(1);
frame_la_LKs->GetXaxis()->SetNdivisions(510);
frame_la_LKs->GetYaxis()->SetNdivisions(406);
frame_la_LKs->GetXaxis()->SetTickSize(0.1);
frame_la_LKs->GetYaxis()->SetTickSize(0.01);
frame_la_LKs->GetYaxis()->SetTitleOffset(1.8);
frame_la_LKs->GetXaxis()->SetLabelSize(0.12);
frame_la_LKs->GetYaxis()->SetLabelSize(0.22);
frame_la_LKs->GetYaxis()->SetLabelOffset(0.01);
TGratio_LKs->SetMarkerStyle(25);
TGratio_LKs->SetMarkerSize(1.3);
TGratio_LKs->SetMarkerColor(kBlack);
TGratio_LKs->Draw("P");
TLine* line2_LKs = new TLine(full_range_la.first, 1, full_range_la.second, 1);
line2_LKs->SetLineStyle(2);
line2_LKs->Draw("same");
c2_LKs->Update();
c2_LKs->Print("results/Mass"+V0+"LamK_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");

//ALam used in K0sALam
cout << "==================ALam Mass for K0sAL==================" << endl;
TH1D* masslam_ALKs = (TH1D*)h_ALK_Mass_1D->Clone();
//masslam_ALKs->Sumw2();
TF1 *func2_ALKs = new TF1("func2_ALKs", massfunc.c_str(), full_range_la.first, full_range_la.second);

//sideband fit
int cnt2_ALKs = 0;
for ( int i = 0; i < masslam_ALKs->GetNbinsX(); i++ ) {
    double center = masslam_ALKs->GetBinCenter(i+1);
    if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
        cnt2_ALKs++;
    }
}
TGraphErrors gr2_ALKs(cnt2_ALKs);
cnt2_ALKs = 0;
for ( int i = 0; i < masslam_ALKs->GetNbinsX(); i++ ) {
    double center = masslam_ALKs->GetBinCenter(i+1);
    if ( (center > sb1_la.first and center < sb1_la.second) or (center > sb2_la.first and center < sb2_la.second) ) {
       gr2_ALKs.GetX()[cnt2_ALKs] = center;
       gr2_ALKs.GetY()[cnt2_ALKs] = masslam_ALKs->GetBinContent(i+1);
       gr2_ALKs.GetEY()[cnt2_ALKs] = masslam_ALKs->GetBinError(i+1);
       cnt2_ALKs++;
     }
}

//sideband fit parameters
func2_ALKs->SetParameter(5, 10.);
func2_ALKs->SetParameter(6, 1.);
func2_ALKs->SetParameter(7, 1.);
func2_ALKs->SetParameter(8, 1.);
func2_ALKs->SetParameter(9, 1.);

//fix the parameters in order to get only bkg
func2_ALKs->FixParameter(0, mass_la);
func2_ALKs->FixParameter(1, sigma_la);
func2_ALKs->FixParameter(2, sigma_la);
func2_ALKs->FixParameter(3, 0.);
func2_ALKs->FixParameter(4, 0.);

//make sideband fit
gr2_ALKs.Fit(func2_ALKs, "q0");
gr2_ALKs.Fit(func2_ALKs, "q0");
gr2_ALKs.Fit(func2_ALKs, "q0");

//Set sideband parameters to full fit
func2_ALKs->FixParameter(5, func2_ALKs->GetParameter(5));
func2_ALKs->FixParameter(6, func2_ALKs->GetParameter(6));
func2_ALKs->FixParameter(7, func2_ALKs->GetParameter(7));
func2_ALKs->FixParameter(8, func2_ALKs->GetParameter(8));
func2_ALKs->FixParameter(9, func2_ALKs->GetParameter(9));

//release the signal parameters
func2_ALKs->ReleaseParameter(0);
func2_ALKs->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2_ALKs->ReleaseParameter(1);
func2_ALKs->ReleaseParameter(2);
func2_ALKs->ReleaseParameter(3);
func2_ALKs->ReleaseParameter(4);
func2_ALKs->SetParLimits(1, 0, sigma_la);
func2_ALKs->SetParLimits(2, 0, sigma_la);
func2_ALKs->SetParLimits(3, 0, 0.8*masslam_ALKs->GetBinContent(masslam_ALKs->FindBin(mass_la)));
func2_ALKs->SetParLimits(4, 0, 0.8*masslam_ALKs->GetBinContent(masslam_ALKs->FindBin(mass_la)));

//make a full fit
masslam_ALKs->Fit(func2_ALKs, "q0", "", full_range_la.first, full_range_la.second);
masslam_ALKs->Fit(func2_ALKs, "q0", "", full_range_la.first, full_range_la.second);
masslam_ALKs->Fit(func2_ALKs, "q0", "", full_range_la.first, full_range_la.second);

//now release all the parameters
func2_ALKs->ReleaseParameter(0);
func2_ALKs->ReleaseParameter(1);
func2_ALKs->ReleaseParameter(2);
func2_ALKs->ReleaseParameter(3);
func2_ALKs->ReleaseParameter(4);
func2_ALKs->ReleaseParameter(5);
func2_ALKs->ReleaseParameter(6);
func2_ALKs->ReleaseParameter(7);
func2_ALKs->ReleaseParameter(8);
func2_ALKs->ReleaseParameter(9);
func2_ALKs->SetParLimits(0, mass_la-Dmass_la, mass_la+Dmass_la);
func2_ALKs->SetParLimits(1, 0, sigma_la);
func2_ALKs->SetParLimits(2, 0, sigma_la);
func2_ALKs->SetParLimits(3, 0, 0.8*masslam_ALKs->GetBinContent(masslam_ALKs->FindBin(mass_la)));
func2_ALKs->SetParLimits(4, 0, 0.8*masslam_ALKs->GetBinContent(masslam_ALKs->FindBin(mass_la)));

//make final fits
for ( int i = 0; i < iter_la-1; i++ ) {
    masslam_ALKs->Fit(func2_ALKs, "q0", "", full_range_la.first, full_range_la.second);
}

//make signal and bkg functions for plots
TF1 * func2_signal_ALKs   = 0;
TF1 * func2_bckgnd_ALKs   = 0;
TF1 * gaus1_signal_la_ALKs   = 0;
TF1 * gaus2_signal_la_ALKs   = 0;
func2_signal_ALKs  = new TF1("func2_signal_ALKs" , double_gaussian.c_str(), full_range_la.first, full_range_la.second);
func2_bckgnd_ALKs = new TF1("func2_bckgnd_ALKs" , poly_bkg0.c_str(), full_range_la.first, full_range_la.second);
gaus1_signal_la_ALKs = new TF1("gaus1_signal_la_ALKs", gaus_func.c_str(), full_range_la.first, full_range_la.second);
gaus2_signal_la_ALKs = new TF1("gaus2_signal_la_ALKs", gaus_func.c_str(), full_range_la.first, full_range_la.second);

//set parameters for signal and bkg
func2_signal_ALKs->SetParameter(0, func2_ALKs->GetParameter(0));
func2_signal_ALKs->SetParameter(1, func2_ALKs->GetParameter(1));
func2_signal_ALKs->SetParameter(2, func2_ALKs->GetParameter(2));
func2_signal_ALKs->SetParameter(3, func2_ALKs->GetParameter(3));
func2_signal_ALKs->SetParameter(4, func2_ALKs->GetParameter(4));
func2_bckgnd_ALKs->SetParameter(0, func2_ALKs->GetParameter(5));
func2_bckgnd_ALKs->SetParameter(1, func2_ALKs->GetParameter(6));
func2_bckgnd_ALKs->SetParameter(2, func2_ALKs->GetParameter(7));
func2_bckgnd_ALKs->SetParameter(3, func2_ALKs->GetParameter(8));
func2_bckgnd_ALKs->SetParameter(4, func2_ALKs->GetParameter(9));

//set parameters for gaussians
gaus1_signal_la_ALKs->SetParameter(0, func2_ALKs->GetParameter(0));
gaus1_signal_la_ALKs->SetParameter(1, func2_ALKs->GetParameter(1));
gaus1_signal_la_ALKs->SetParameter(2, func2_ALKs->GetParameter(3));
gaus2_signal_la_ALKs->SetParameter(0, func2_ALKs->GetParameter(0));
gaus2_signal_la_ALKs->SetParameter(1, func2_ALKs->GetParameter(2));
gaus2_signal_la_ALKs->SetParameter(2, func2_ALKs->GetParameter(4));

//make pull and ratio
double vny2_ALKs[400];
double xvn2_ALKs[400];
double pullvn2_ALKs[400];
double ratiovn2_ALKs[400];
double err_pullvn2_ALKs[400];
double err_ratiovn2_ALKs[400];
for(int k=0;k<masslam_ALKs->GetNbinsX();k++)
        	{
            vny2_ALKs[k] = masslam_ALKs->GetBinContent(k+1);
            xvn2_ALKs[k] = masslam_ALKs->GetBinCenter(k+1);
        	if(masslam_ALKs->GetBinContent(k+1)!=0){
            pullvn2_ALKs[k] = (vny2_ALKs[k] - func2_ALKs->Eval(xvn2_ALKs[k]))/masslam_ALKs->GetBinError(k+1);
            ratiovn2_ALKs[k] = vny2_ALKs[k]/func2_ALKs->Eval(xvn2_ALKs[k]);
            err_pullvn2_ALKs[k] = masslam_ALKs->GetBinError(k+1);
            err_ratiovn2_ALKs[k] = (vny2_ALKs[k]/func2_ALKs->Eval(xvn2_ALKs[k]))*(masslam_ALKs->GetBinError(k+1)/vny2_ALKs[k]);
            }else{
            pullvn2_ALKs[k] = 0.0;
            ratiovn2_ALKs[k] = 0.0;
            err_pullvn2_ALKs[k]=0.0;
            err_ratiovn2_ALKs[k]=0.0;
            }
}

//make the plots
TH1D* lamass_ALKs = (TH1D*)masslam_ALKs->Clone();
//lamass_ALKs->Sumw2();
TGraphErrors* TGpull_ALKs = new TGraphErrors(400,xvn2_ALKs,pullvn2_ALKs,0,0);
TGraphErrors* TGratio_ALKs= new TGraphErrors(400,xvn2_ALKs,ratiovn2_ALKs,0,err_ratiovn2_ALKs);
TCanvas* c2_ALKs = new TCanvas("c2_ALKs","c2_ALKs",600,600);
c2_ALKs->cd();
TPad* pad1_ALKs = new TPad("pad1_ALKs","top pad",0.0,0.4,1.0,1.0);
TPad* pad2_ALKs = new TPad("pad3_ALKs","middle pad",0.0,0.0,1.0,0.25);
TPad* pad3_ALKs = new TPad("pad2_ALKs","bottom pad",0.0,0.25,1.0,0.4);
pad1_ALKs->SetTopMargin(0.1);
pad1_ALKs->SetBottomMargin(0.0);
pad1_ALKs->SetRightMargin(0.038);
pad1_ALKs->SetLeftMargin(0.15);
pad1_ALKs->Draw();
pad2_ALKs->SetFrameFillStyle(4000);
pad2_ALKs->SetTopMargin(0.0);
pad2_ALKs->SetBottomMargin(0.37);
pad2_ALKs->SetRightMargin(0.038);
pad2_ALKs->SetLeftMargin(0.15);
pad2_ALKs->Draw();
pad3_ALKs->SetFrameFillStyle(4000);
pad3_ALKs->SetTopMargin(0.0);
pad3_ALKs->SetBottomMargin(0.0);
pad3_ALKs->SetRightMargin(0.038);
pad3_ALKs->SetLeftMargin(0.15);
pad3_ALKs->Draw();
pad1_ALKs->cd();
gPad->SetTickx();
gPad->SetTicky();
lamass_ALKs->SetTitle("");
lamass_ALKs->Draw("E");
TGaxis::SetMaxDigits(3);
lamass_ALKs->GetXaxis()->SetRangeUser(full_range_la.first,full_range_la.second);
lamass_ALKs->GetYaxis()->SetRangeUser(0.1,1.05*lamass_ALKs->GetMaximum());
lamass_ALKs->SetMarkerStyle(24);
lamass_ALKs->SetMarkerSize(1.5);
lamass_ALKs->GetYaxis()->SetTitle("#scale[1.7]{Candidates / 0.5 MeV}");
lamass_ALKs->GetYaxis()->CenterTitle(1);
lamass_ALKs->GetYaxis()->CenterTitle(1);
lamass_ALKs->GetYaxis()->SetNdivisions(810);
lamass_ALKs->GetXaxis()->SetTickSize(0.03);
lamass_ALKs->GetYaxis()->SetTickSize(0.02);
lamass_ALKs->GetXaxis()->SetLabelSize(0.055);
lamass_ALKs->GetYaxis()->SetLabelSize(0.058);
lamass_ALKs->GetYaxis()->SetLabelOffset(0.01);
lamass_ALKs->GetYaxis()->SetTitleOffset(1.8);
func2_ALKs->Draw("LSAME");
func2_ALKs->SetLineColor(SPcolors[SPdred]);
func2_bckgnd_ALKs->Draw("LSAME");
func2_bckgnd_ALKs->SetLineColor(SPcolors[SPdred]);
func2_bckgnd_ALKs->SetLineStyle(2);
double meanla_ALKs  = func2_ALKs->GetParameter(0);
double sigma1la_ALKs = func2_ALKs->GetParameter(1);
double sigma2la_ALKs = func2_ALKs->GetParameter(2);
double gaus1_yield_la_ALKs = gaus1_signal_la_ALKs->Integral(meanla_ALKs - 2*sigma1la_ALKs, meanla_ALKs + 2*sigma1la)/masslam_ALKs->GetBinWidth(1);
double gaus2_yield_la_ALKs = gaus2_signal_la_ALKs->Integral(meanla_ALKs - 2*sigma2la_ALKs, meanla_ALKs + 2*sigma2la)/masslam_ALKs->GetBinWidth(1);
double gausTot_yield_la_ALKs = gaus1_yield_la_ALKs + gaus2_yield_la_ALKs;
double rms_gaus1_sig_la_ALKs = gaus1_yield_la_ALKs/gausTot_yield_la_ALKs;
double rms_gaus2_sig_la_ALKs = gaus2_yield_la_ALKs/gausTot_yield_la_ALKs;
double rms_true_la_ALKs = TMath::Sqrt(rms_gaus1_sig_la_ALKs*sigma1la_ALKs*sigma1la_ALKs + rms_gaus2_sig_la_ALKs*sigma2la_ALKs*sigma2la_ALKs);
std::pair<double, double>  signal_range_la_ALKs(meanla_ALKs - 2*rms_true_la_ALKs, meanla_ALKs + 2*rms_true_la_ALKs);
double S_la_ALKs = func2_signal_ALKs->Integral(signal_range_la_ALKs.first, signal_range_la_ALKs.second)/masslam_ALKs->GetBinWidth(1);
double B_la_ALKs = func2_bckgnd_ALKs->Integral(signal_range_la_ALKs.first, signal_range_la_ALKs.second)/masslam_ALKs->GetBinWidth(1);
double Fsig_la_ALKs = S_la_ALKs/(S_la_ALKs+B_la_ALKs);
double sig_la_ALKs = S_la_ALKs/sqrt(S_la_ALKs+B_la_ALKs);
double chi2la_ALKs = func2_ALKs->GetChisquare()/func2_ALKs->GetNDF();
double probla_ALKs = TMath::Prob(func2_ALKs->GetChisquare(),func2_ALKs->GetNDF());
std::ostringstream osx_ALKs; 
TLatex* texCMSx_ALKs = new TLatex;
texCMSx_ALKs->SetNDC();
texCMSx_ALKs->SetTextFont(42);
texCMSx_ALKs->SetTextSize(0.05);
texCMSx_ALKs->SetTextAlign(12);
texCMSx_ALKs->DrawLatex(0.25,.95," #font[61]{CMS} #it{Work in Progress}");
texCMSx_ALKs->DrawLatex(0.70,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMSx_ALKs->DrawLatex(0.26,0.4, "#scale[1.7]{#bar{#Lambda}}");
texCMSx_ALKs->DrawLatex(0.225,0.2, "#scale[1.0]{Used in K^{0}_{S}#bar{#Lambda}}"); 
texCMSx_ALKs->DrawLatex(0.225,0.15, "#scale[1.0]{correlations}"); 
// texCMSx_ALKs->DrawLatex(0.75,0.15, "This Thesis"); 
osx_ALKs << "#scale[1.0]{s/#sqrt{s+b}: " << std::setprecision(4) << sig_la_ALKs << "}" << std::setprecision(2);           
texCMSx_ALKs->DrawLatex(0.68,0.60, osx_ALKs.str(  ).c_str(  )); 
osx_ALKs.str(std::string());
osx_ALKs << "#scale[1.0]{s/(s+b): " << std::setprecision(4) << Fsig_la_ALKs*100 << " %}" << std::setprecision(2); 
texCMSx_ALKs->DrawLatex(0.68,0.52, osx_ALKs.str(  ).c_str(  ));  
osx_ALKs.str(std::string()); 
osx_ALKs << "#scale[1.0]{Mean: " << std::setprecision(4) << meanla_ALKs << " GeV}" << std::setprecision(6);   
texCMSx_ALKs->DrawLatex(0.68,0.46, osx_ALKs.str(  ).c_str(  ));
osx_ALKs.str(std::string()); 
osx_ALKs << "#scale[1.0]{#sigma: " << std::setprecision(4) << rms_true_la_ALKs << " }" << std::setprecision(2); 
texCMSx_ALKs->DrawLatex(0.68,0.40, osx_ALKs.str(  ).c_str(  )); 
osx_ALKs.str(std::string()); 
osx_ALKs << "#scale[1.0]{Yield: " << std::setprecision(4) << S_la_ALKs << " }" << std::setprecision(2); 
texCMSx_ALKs->DrawLatex(0.68,0.33, osx_ALKs.str(  ).c_str(  )); 
osx_ALKs.str(std::string()); 
osx_ALKs << "#scale[1.0]{#chi^{2}/ndf: " << std::setprecision(4) << chi2la_ALKs << " }" << std::setprecision(4); 
texCMSx_ALKs->DrawLatex(0.68,0.26, osx_ALKs.str(  ).c_str(  ));  
osx_ALKs.str(std::string()); 
// osx_ALKs << "#scale[1.0]{Prob: " << std::setprecision(4) << probla_ALKs << "}" << std::setprecision(2); 
// texCMSx_ALKs->DrawLatex(0.68,0.20, osx_ALKs.str(  ).c_str(  )); 
// osx_ALKs.str(std::string()); 
osx_ALKs << "#scale[1.]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMSx_ALKs->DrawLatex(0.68,0.75, osx_ALKs.str(  ).c_str(  )); 
osx_ALKs.str(std::string()); 
osx_ALKs << "#scale[1.]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMSx_ALKs->DrawLatex(0.68,0.67, osx_ALKs.str(  ).c_str(  ));  
osx_ALKs.str(std::string()); 
osx_ALKs << "#scale[1.]{ " << std::setprecision(4) << qminM << std::setprecision(4) <<" < q_{inv} < "  << qmaxM << " GeV}" << std::setprecision(2); 
texCMSx_ALKs->DrawLatex(0.68,0.83, osx_ALKs.str(  ).c_str(  ));  
osx_ALKs.str(std::string()); 
Float_t yymax_ALKs = lamass_ALKs->GetMaximum();
TLine *linex_ALKs = new TLine(signal_range_la_ALKs.first,0.,signal_range_la_ALKs.first,yymax_ALKs);
linex_ALKs->SetLineColor(SPcolors[SPdblue]);
linex_ALKs->Draw();
TLine *line2x_ALKs = new TLine(signal_range_la_ALKs.second,0.,signal_range_la_ALKs.second,yymax_ALKs);
line2x_ALKs->SetLineColor(SPcolors[SPdblue]);
line2x_ALKs->Draw();
TLegend* legLL_ALKs = new TLegend(0.16,0.6,0.43,0.85);
legLL_ALKs->AddEntry(lamass_ALKs,"Data");
legLL_ALKs->AddEntry(func2_ALKs,"Signal+Background","L");
legLL_ALKs->AddEntry(func2_bckgnd_ALKs,"Background","L");
legLL_ALKs->AddEntry(linex,"Peak Region","L");
legLL_ALKs->SetBorderSize(0);
legLL_ALKs->Draw();
c2_ALKs->Update();
pad2_ALKs->cd();
gPad->SetTickx();
gPad->SetTicky();
double low_ALKs = -49.1;
double high_ALKs = 49.1;
TH1F* frame_la_ALKs = pad2_ALKs->cd()->DrawFrame(full_range_la.first,low_ALKs,full_range_la.second,high_ALKs);
frame_la_ALKs->GetXaxis()->SetTitle("#scale[4.3]{#bar{p} #pi^{+} invariant mass [GeV]}");
frame_la_ALKs->GetYaxis()->SetTitle("#scale[4.0]{Pull}");
frame_la_ALKs->GetXaxis()->CenterTitle(1);
frame_la_ALKs->GetYaxis()->CenterTitle(1);
frame_la_ALKs->GetXaxis()->SetNdivisions(510);
frame_la_ALKs->GetYaxis()->SetNdivisions(406);
frame_la_ALKs->GetXaxis()->SetTickSize(0.06);
frame_la_ALKs->GetYaxis()->SetTickSize(0.02);
frame_la_ALKs->GetXaxis()->SetTitleOffset(4.5);
frame_la_ALKs->GetYaxis()->SetTitleOffset(1.8);
frame_la_ALKs->GetXaxis()->SetLabelSize(0.14);
frame_la_ALKs->GetYaxis()->SetLabelSize(0.14);
frame_la_ALKs->GetXaxis()->SetLabelOffset(0.01);
frame_la_ALKs->GetYaxis()->SetLabelOffset(0.013);
TGpull_ALKs->SetMarkerStyle(25);
TGpull_ALKs->SetMarkerSize(1.3);
TGpull_ALKs->SetMarkerColor(kBlack);
TGpull_ALKs->Draw("P");
TLine* line_ALKs = new TLine(full_range_la.first, 0, full_range_la.second, 0);
line_ALKs->SetLineStyle(2);
line_ALKs->Draw("same");
pad3_ALKs->cd();
gPad->SetTickx();
gPad->SetTicky();
double low2_ALKs = 0.1;
double high2_ALKs = 1.89999999;
frame_la_ALKs = pad3_ALKs->cd()->DrawFrame(full_range_la.first,low2_ALKs,full_range_la.second,high2_ALKs);
frame_la_ALKs->GetYaxis()->SetTitle("#scale[6.3]{#frac{Data}{Fit}}");
frame_la_ALKs->GetXaxis()->CenterTitle(1);
frame_la_ALKs->GetYaxis()->CenterTitle(1);
frame_la_ALKs->GetXaxis()->SetNdivisions(510);
frame_la_ALKs->GetYaxis()->SetNdivisions(406);
frame_la_ALKs->GetXaxis()->SetTickSize(0.1);
frame_la_ALKs->GetYaxis()->SetTickSize(0.01);
frame_la_ALKs->GetYaxis()->SetTitleOffset(1.8);
frame_la_ALKs->GetXaxis()->SetLabelSize(0.12);
frame_la_ALKs->GetYaxis()->SetLabelSize(0.22);
frame_la_ALKs->GetYaxis()->SetLabelOffset(0.01);
TGratio_ALKs->SetMarkerStyle(25);
TGratio_ALKs->SetMarkerSize(1.3);
TGratio_ALKs->SetMarkerColor(kBlack);
TGratio_ALKs->Draw("P");
TLine* line2_ALKs = new TLine(full_range_la.first, 1, full_range_la.second, 1);
line2_ALKs->SetLineStyle(2);
line2_ALKs->Draw("same");
c2_ALKs->Update();
c2_ALKs->Print("results/Mass"+V0+"ALamK_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");


//end of mass fits 

//Time to extract the parameters for the sideband method

//K0sK0s

double tott =(double) (S_ks+B_ks);
double sigg =(double) S_ks;
double bkgg =(double) B_ks;

double ftot =(double) (0.5*(pow(tott,2) - tott));//TMath::Factorial(tott);//(2*(TMath::Factorial(tott-2)));    
double fss = (double) ((0.5*(pow(sigg,2) - sigg))/ftot);//TMath::Factorial(sigg);//(2*(TMath::Factorial(sigg-2)));
double fbb =(double) ((0.5*(pow(bkgg,2) - bkgg))/ftot);//TMath::Factorial(bkgg);//(2*(TMath::Factorial(bkgg-2)));
double fsb =(double) (ftot/ftot - fss - fbb);
double FF = fss + fbb + fsb;


//LambdaLambda

double tott_LL =(double) (S_la+B_la);
double sigg_LL =(double) S_la;
double bkgg_LL =(double) B_la;

double ftot_LL =(double) (0.5*(pow(tott_LL,2) - tott_LL));//TMath::Factorial(tott);//(2*(TMath::Factorial(tott-2)));    
double fss_LL = (double) ((0.5*(pow(sigg_LL,2) - sigg_LL))/ftot_LL);//TMath::Factorial(sigg);//(2*(TMath::Factorial(sigg-2)));
double fbb_LL =(double) ((0.5*(pow(bkgg_LL,2) - bkgg_LL))/ftot_LL);//TMath::Factorial(bkgg);//(2*(TMath::Factorial(bkgg-2)));
double fsb_LL =(double) (ftot_LL/ftot_LL - fss_LL - fbb_LL);
double FF_LL = fss_LL + fbb_LL + fsb_LL;

//AntiLambdaAntiLambda

double tott_ALAL =(double) (S_la_LAL+B_la_LAL);
double sigg_ALAL =(double) S_la_LAL;
double bkgg_ALAL =(double) B_la_LAL;

double ftot_ALAL =(double) (0.5*(pow(tott_ALAL,2) - tott_ALAL));//TMath::Factorial(tott);//(2*(TMath::Factorial(tott-2)));    
double fss_ALAL = (double) ((0.5*(pow(sigg_ALAL,2) - sigg_ALAL))/ftot_ALAL);//TMath::Factorial(sigg);//(2*(TMath::Factorial(sigg-2)));
double fbb_ALAL =(double) ((0.5*(pow(bkgg_ALAL,2) - bkgg_ALAL))/ftot_ALAL);//TMath::Factorial(bkgg);//(2*(TMath::Factorial(bkgg-2)));
double fsb_ALAL =(double) (ftot_ALAL/ftot_ALAL - fss_ALAL - fbb_ALAL);
double FF_ALAL = fss_ALAL + fbb_ALAL + fsb_ALAL;


//LambdaLambda+AntiLambdaAntiLambda

double tott_LLALAL =(double) (S_la+S_la_LAL+B_la_LAL+B_la);
double sigg_LLALAL =(double) S_la_LAL+S_la;
double bkgg_LLALAL =(double) B_la_LAL+B_la;

double ftot_LLALAL =(double) (0.5*(pow(tott_LLALAL,2) - tott_LLALAL));//TMath::Factorial(tott);//(2*(TMath::Factorial(tott-2)));    
double fss_LLALAL = (double) ((0.5*(pow(sigg_LLALAL,2) - sigg_LLALAL))/ftot_LLALAL);//TMath::Factorial(sigg);//(2*(TMath::Factorial(sigg-2)));
double fbb_LLALAL =(double) ((0.5*(pow(bkgg_LLALAL,2) - bkgg_LLALAL))/ftot_LLALAL);//TMath::Factorial(bkgg);//(2*(TMath::Factorial(bkgg-2)));
double fsb_LLALAL =(double) (ftot_LLALAL/ftot_LLALAL - fss_LLALAL - fbb_LLALAL);
double FF_LLALAL = fss_LLALAL + fbb_LLALAL + fsb_LLALAL;

//LambdaAntiLambda

double tott_LAL =(double) (S_la_LLAL+S_la_ALAL+B_la_LLAL+B_la_ALAL);
double sigg_LAL =(double) S_la_LLAL+S_la_ALAL;
double bkgg_LAL =(double) B_la_LLAL+B_la_ALAL;

double ftot_LAL =(double) (0.5*(pow(tott_LAL,2) - tott_LAL));//TMath::Factorial(tott);//(2*(TMath::Factorial(tott-2)));    
double fss_LAL = (double) ((0.5*(pow(sigg_LAL,2) - sigg_LAL))/ftot_LAL);//TMath::Factorial(sigg);//(2*(TMath::Factorial(sigg-2)));
double fbb_LAL =(double) ((0.5*(pow(bkgg_LAL,2) - bkgg_LAL))/ftot_LAL);//TMath::Factorial(bkgg);//(2*(TMath::Factorial(bkgg-2)));
double fsb_LAL =(double) (ftot_LAL/ftot_LAL - fss_LAL - fbb_LAL);
double FF_LAL = fss_LAL + fbb_LAL + fsb_LAL;


//K0sLambda

double tott_KL =(double) (S_ks_KL+S_la_LKs+B_ks_KL+B_la_LKs);
double sigg_KL =(double) S_ks_KL+S_la_LKs;
double bkgg_KL =(double) B_ks_KL+B_la_LKs;

double ftot_KL =(double) (0.5*(pow(tott_KL,2) - tott_KL));//TMath::Factorial(tott);//(2*(TMath::Factorial(tott-2)));    
double fss_KL = (double) ((0.5*(pow(sigg_KL,2) - sigg_KL))/ftot_KL);//TMath::Factorial(sigg);//(2*(TMath::Factorial(sigg-2)));
double fbb_KL =(double) ((0.5*(pow(bkgg_KL,2) - bkgg_KL))/ftot_KL);//TMath::Factorial(bkgg);//(2*(TMath::Factorial(bkgg-2)));
double fsb_KL =(double) (ftot_KL/ftot_KL - fss_KL - fbb_KL);
double FF_KL = fss_KL + fbb_KL + fsb_KL;

//K0sALambda

double tott_KAL =(double) (S_ks_KAL+S_la_ALKs+B_ks_KAL+B_la_ALKs);
double sigg_KAL =(double) S_ks_KAL+S_la_ALKs;
double bkgg_KAL =(double) B_ks_KAL+B_la_ALKs;

double ftot_KAL =(double) (0.5*(pow(tott_KAL,2) - tott_KAL));//TMath::Factorial(tott);//(2*(TMath::Factorial(tott-2)));    
double fss_KAL = (double) ((0.5*(pow(sigg_KAL,2) - sigg_KAL))/ftot_KAL);//TMath::Factorial(sigg);//(2*(TMath::Factorial(sigg-2)));
double fbb_KAL =(double) ((0.5*(pow(bkgg_KAL,2) - bkgg_KAL))/ftot_KAL);//TMath::Factorial(bkgg);//(2*(TMath::Factorial(bkgg-2)));
double fsb_KAL =(double) (ftot_KAL/ftot_KAL - fss_KAL - fbb_KAL);
double FF_KAL = fss_KAL + fbb_KAL + fsb_KAL;

//K0sLambda + K0sALambda

double tott_KLAL =(double) (S_ks_KL+S_la_LKs+S_ks_KAL+S_la_ALKs+B_ks_KL+B_la_LKs+B_ks_KAL+B_la_ALKs);
double sigg_KLAL =(double) S_ks_KL+S_la_LKs+S_ks_KAL+S_la_ALKs;
double bkgg_KLAL =(double) B_ks_KL+B_la_LKs+B_ks_KAL+B_la_ALKs;

double ftot_KLAL =(double) (0.5*(pow(tott_KLAL,2) - tott_KLAL));//TMath::Factorial(tott);//(2*(TMath::Factorial(tott-2)));    
double fss_KLAL = (double) ((0.5*(pow(sigg_KLAL,2) - sigg_KLAL))/ftot_KLAL);//TMath::Factorial(sigg);//(2*(TMath::Factorial(sigg-2)));
double fbb_KLAL =(double) ((0.5*(pow(bkgg_KLAL,2) - bkgg_KLAL))/ftot_KLAL);//TMath::Factorial(bkgg);//(2*(TMath::Factorial(bkgg-2)));
double fsb_KLAL =(double) (ftot_KLAL/ftot_KLAL - fss_KLAL - fbb_KLAL);
double FF_KLAL = fss_KLAL + fbb_KLAL + fsb_KLAL;


cout << "=============Fraction Extraction=============" << endl;

cout << "=============K0sK0s=============" << endl;
cout << "fss = " << fss << endl;
cout << "fbb = " << fbb << endl; 
cout << "fsb = " << fsb << endl;
cout << "fss + fbb + fsb = " << FF << endl;
cout << endl;

cout << "=============LL=============" << endl;
cout << "fss = " << fss_LL << endl;
cout << "fbb = " << fbb_LL << endl; 
cout << "fsb = " << fsb_LL << endl;
cout << "fss + fbb + fsb = " << FF_LL << endl;
cout << endl;

cout << "=============ALAL=============" << endl;
cout << "fss = " << fss_ALAL << endl;
cout << "fbb = " << fbb_ALAL << endl; 
cout << "fsb = " << fsb_ALAL << endl;
cout << "fss + fbb + fsb = " << FF_ALAL << endl;
cout << endl;

cout << "=============LL+ALAL=============" << endl;
cout << "fss = " << fss_LLALAL << endl;
cout << "fbb = " << fbb_LLALAL << endl; 
cout << "fsb = " << fsb_LLALAL << endl;
cout << "fss + fbb + fsb = " << FF_LLALAL << endl;
cout << endl;

cout << "=============LAL=============" << endl;
cout << "fss = " << fss_LAL << endl;
cout << "fbb = " << fbb_LAL << endl; 
cout << "fsb = " << fsb_LAL << endl;
cout << "fss + fbb + fsb = " << FF_LAL << endl;
cout << endl;

cout << "=============KL=============" << endl;
cout << "fss = " << fss_KL << endl;
cout << "fbb = " << fbb_KL << endl; 
cout << "fsb = " << fsb_KL << endl;
cout << "fss + fbb + fsb = " << FF_KL << endl;
cout << endl;

cout << "=============KAL=============" << endl;
cout << "fss = " << fss_KAL << endl;
cout << "fbb = " << fbb_KAL << endl; 
cout << "fsb = " << fsb_KAL << endl;
cout << "fss + fbb + fsb = " << FF_KAL << endl;
cout << endl;

cout << "=============KLAL=============" << endl;
cout << "fss = " << fss_KLAL << endl;
cout << "fbb = " << fbb_KLAL << endl; 
cout << "fsb = " << fsb_KLAL << endl;
cout << "fss + fbb + fsb = " << FF_KLAL << endl;
cout << "==============================" << endl;
cout << endl;

cout << endl;
    
h_K0s_Mass->Reset("ICESM");
h_Lam_Mass->Reset("ICESM");
h_ALam_Mass->Reset("ICESM");
h_Lam_Mass_LAL->Reset("ICESM");
h_ALam_Mass_LAL->Reset("ICESM");
h_KL_Mass->Reset("ICESM");
h_KAL_Mass->Reset("ICESM");
h_LK_Mass->Reset("ICESM");
h_ALK_Mass->Reset("ICESM");


//until here, everything is ok    
    
//K0sK0s
THnSparseD* h_K0s_qinv_ss=(THnSparseD*)f->Get("HBT_K0s/hpeak_K0s_K0s_T");
THnSparseD* h_K0s_qinv_mix_ss=(THnSparseD*)f->Get("HBT_K0s/hpeak_K0s_K0s_T_mix");
THnSparseD* h_K0s_qinv_sb=(THnSparseD*)f->Get("HBT_K0s/hpeakside_K0s_K0s_T");
THnSparseD* h_K0s_qinv_mix_sb=(THnSparseD*)f->Get("HBT_K0s/hpeakside_K0s_K0s_T_mix");
THnSparseD* h_K0s_qinv_bb=(THnSparseD*)f->Get("HBT_K0s/hside_K0s_K0s_T");
THnSparseD* h_K0s_qinv_mix_bb=(THnSparseD*)f->Get("HBT_K0s/hside_K0s_K0s_T_mix");

//ALamALam
THnSparseD* h_ALam_qinv_ss=(THnSparseD*)f->Get("HBT_ALam/hpeak_ALam_ALam_T");
THnSparseD* h_ALam_qinv_mix_ss=(THnSparseD*)f->Get("HBT_ALam/hpeak_ALam_ALam_T_mix");
THnSparseD* h_ALam_qinv_sb=(THnSparseD*)f->Get("HBT_ALam/hpeakside_ALam_ALam_T");
THnSparseD* h_ALam_qinv_mix_sb=(THnSparseD*)f->Get("HBT_ALam/hpeakside_ALam_ALam_T_mix");
THnSparseD* h_ALam_qinv_bb=(THnSparseD*)f->Get("HBT_ALam/hside_ALam_ALam_T");
THnSparseD* h_ALam_qinv_mix_bb=(THnSparseD*)f->Get("HBT_ALam/hside_ALam_ALam_T_mix");

//LamLam
THnSparseD* h_Lam_qinv_ss=(THnSparseD*)f->Get("HBT_Lam/hpeak_Lam_Lam_T");
THnSparseD* h_Lam_qinv_mix_ss=(THnSparseD*)f->Get("HBT_Lam/hpeak_Lam_Lam_T_mix");
THnSparseD* h_Lam_qinv_sb=(THnSparseD*)f->Get("HBT_Lam/hpeakside_Lam_Lam_T");
THnSparseD* h_Lam_qinv_mix_sb=(THnSparseD*)f->Get("HBT_Lam/hpeakside_Lam_Lam_T_mix");
THnSparseD* h_Lam_qinv_bb=(THnSparseD*)f->Get("HBT_Lam/hside_Lam_Lam_T");
THnSparseD* h_Lam_qinv_mix_bb=(THnSparseD*)f->Get("HBT_Lam/hside_Lam_Lam_T_mix");

//LamALam
THnSparseD* h_LAL_qinv_ss=(THnSparseD*)f->Get("HBT_LAL/hpeak_Lam_ALam_T");
THnSparseD* h_LAL_qinv_mix_ss=(THnSparseD*)f->Get("HBT_LAL/hpeak_Lam_ALam_T_mix");
THnSparseD* h_LAL_qinv_sb=(THnSparseD*)f->Get("HBT_LAL/hpeakside_Lam_ALam_T");
THnSparseD* h_LAL_qinv_mix_sb=(THnSparseD*)f->Get("HBT_LAL/hpeakside_Lam_ALam_T_mix");
THnSparseD* h_LAL_qinv_bb=(THnSparseD*)f->Get("HBT_LAL/hside_Lam_ALam_T");
THnSparseD* h_LAL_qinv_mix_bb=(THnSparseD*)f->Get("HBT_LAL/hside_Lam_ALam_T_mix");


//K0sLam
THnSparseD* h_KL_qinv_ss=(THnSparseD*)f->Get("HBT_K0sLam/hpeak_K0s_Lam_T");
THnSparseD* h_KL_qinv_mix_ss=(THnSparseD*)f->Get("HBT_K0sLam/hpeak_K0s_Lam_T_mix");
THnSparseD* h_KL_qinv_sb=(THnSparseD*)f->Get("HBT_K0sLam/hpeakside_K0s_Lam_T");
THnSparseD* h_KL_qinv_mix_sb=(THnSparseD*)f->Get("HBT_K0sLam/hpeakside_K0s_Lam_T_mix");
THnSparseD* h_KL_qinv_bb=(THnSparseD*)f->Get("HBT_K0sLam/hside_K0s_Lam_T");
THnSparseD* h_KL_qinv_mix_bb=(THnSparseD*)f->Get("HBT_K0sLam/hside_K0s_Lam_T_mix");

//K0sALam
THnSparseD* h_KAL_qinv_ss=(THnSparseD*)f->Get("HBT_K0sALam/hpeak_K0s_ALam_T");
THnSparseD* h_KAL_qinv_mix_ss=(THnSparseD*)f->Get("HBT_K0sALam/hpeak_K0s_ALam_T_mix");
THnSparseD* h_KAL_qinv_sb=(THnSparseD*)f->Get("HBT_K0sALam/hpeakside_K0s_ALam_T");
THnSparseD* h_KAL_qinv_mix_sb=(THnSparseD*)f->Get("HBT_K0sALam/hpeakside_K0s_ALam_T_mix");
THnSparseD* h_KAL_qinv_bb=(THnSparseD*)f->Get("HBT_K0sALam/hside_K0s_ALam_T");
THnSparseD* h_KAL_qinv_mix_bb=(THnSparseD*)f->Get("HBT_K0sALam/hside_K0s_ALam_T_mix");




//setup ranges
Int_t aux_kt_minx = h_K0s_qinv_ss->GetAxis(1)->FindBin(ktmin+0.0001); 
Int_t aux_kt_maxx = h_K0s_qinv_ss->GetAxis(1)->FindBin(ktmax-0.0001); 

Int_t aux_ntk_minx = h_K0s_qinv_ss->GetAxis(2)->FindBin(ntkmin+0.0001); 
Int_t aux_ntk_maxx = h_K0s_qinv_ss->GetAxis(2)->FindBin(ntkmax-0.0001); 

Int_t aux_q_minx = h_K0s_qinv_ss->GetAxis(0)->FindBin(qmin+0.0001); 
Int_t aux_q_maxx = h_K0s_qinv_ss->GetAxis(0)->FindBin(qmax-0.0001); 


h_K0s_qinv_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_K0s_qinv_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_K0s_qinv_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_K0s_qinv_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_K0s_qinv_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_K0s_qinv_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_K0s_qinv_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_K0s_qinv_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_K0s_qinv_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_K0s_qinv_mix_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_K0s_qinv_mix_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_K0s_qinv_mix_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_K0s_qinv_mix_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_K0s_qinv_mix_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_K0s_qinv_mix_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_K0s_qinv_mix_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_K0s_qinv_mix_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_K0s_qinv_mix_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);


h_Lam_qinv_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_Lam_qinv_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_Lam_qinv_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_Lam_qinv_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_Lam_qinv_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_Lam_qinv_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_Lam_qinv_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_Lam_qinv_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_Lam_qinv_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_Lam_qinv_mix_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_Lam_qinv_mix_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_Lam_qinv_mix_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_Lam_qinv_mix_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_Lam_qinv_mix_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_Lam_qinv_mix_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_Lam_qinv_mix_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_Lam_qinv_mix_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_Lam_qinv_mix_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);


h_ALam_qinv_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_ALam_qinv_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_ALam_qinv_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_ALam_qinv_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_ALam_qinv_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_ALam_qinv_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_ALam_qinv_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_ALam_qinv_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_ALam_qinv_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_ALam_qinv_mix_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_ALam_qinv_mix_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_ALam_qinv_mix_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_ALam_qinv_mix_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_ALam_qinv_mix_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_ALam_qinv_mix_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_ALam_qinv_mix_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_ALam_qinv_mix_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_ALam_qinv_mix_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);


h_LAL_qinv_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_LAL_qinv_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_LAL_qinv_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_LAL_qinv_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_LAL_qinv_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_LAL_qinv_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_LAL_qinv_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_LAL_qinv_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_LAL_qinv_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_LAL_qinv_mix_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_LAL_qinv_mix_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_LAL_qinv_mix_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_LAL_qinv_mix_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_LAL_qinv_mix_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_LAL_qinv_mix_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_LAL_qinv_mix_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_LAL_qinv_mix_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_LAL_qinv_mix_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);

h_KL_qinv_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KL_qinv_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KL_qinv_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_KL_qinv_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KL_qinv_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KL_qinv_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_KL_qinv_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KL_qinv_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KL_qinv_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_KL_qinv_mix_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KL_qinv_mix_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KL_qinv_mix_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_KL_qinv_mix_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KL_qinv_mix_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KL_qinv_mix_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_KL_qinv_mix_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KL_qinv_mix_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KL_qinv_mix_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);

h_KAL_qinv_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KAL_qinv_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KAL_qinv_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_KAL_qinv_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KAL_qinv_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KAL_qinv_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_KAL_qinv_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KAL_qinv_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KAL_qinv_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_KAL_qinv_mix_ss->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KAL_qinv_mix_ss->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KAL_qinv_mix_ss->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_KAL_qinv_mix_sb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KAL_qinv_mix_sb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KAL_qinv_mix_sb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);
h_KAL_qinv_mix_bb->GetAxis(0)->SetRange(aux_q_minx,aux_q_maxx);
h_KAL_qinv_mix_bb->GetAxis(1)->SetRange(aux_kt_minx,aux_kt_maxx);
h_KAL_qinv_mix_bb->GetAxis(2)->SetRange(aux_ntk_minx,aux_ntk_maxx);

//make projections

TH1D *h_K0s_qinv_ss_1D=(TH1D*)h_K0s_qinv_ss->Projection(0,"e");
TH1D *h_K0s_qinv_sb_1D=(TH1D*)h_K0s_qinv_sb->Projection(0,"e");
TH1D *h_K0s_qinv_bb_1D=(TH1D*)h_K0s_qinv_bb->Projection(0,"e");
TH1D *h_K0s_qinv_mix_ss_1D=(TH1D*)h_K0s_qinv_mix_ss->Projection(0,"e");
TH1D *h_K0s_qinv_mix_sb_1D=(TH1D*)h_K0s_qinv_mix_sb->Projection(0,"e");
TH1D *h_K0s_qinv_mix_bb_1D=(TH1D*)h_K0s_qinv_mix_bb->Projection(0,"e");

TH1D *h_ALam_qinv_ss_1D=(TH1D*)h_ALam_qinv_ss->Projection(0,"e");
TH1D *h_ALam_qinv_sb_1D=(TH1D*)h_ALam_qinv_sb->Projection(0,"e");
TH1D *h_ALam_qinv_bb_1D=(TH1D*)h_ALam_qinv_bb->Projection(0,"e");
TH1D *h_ALam_qinv_mix_ss_1D=(TH1D*)h_ALam_qinv_mix_ss->Projection(0,"e");
TH1D *h_ALam_qinv_mix_sb_1D=(TH1D*)h_ALam_qinv_mix_sb->Projection(0,"e");
TH1D *h_ALam_qinv_mix_bb_1D=(TH1D*)h_ALam_qinv_mix_bb->Projection(0,"e");

TH1D *h_Lam_qinv_ss_1D=(TH1D*)h_Lam_qinv_ss->Projection(0,"e");
TH1D *h_Lam_qinv_sb_1D=(TH1D*)h_Lam_qinv_sb->Projection(0,"e");
TH1D *h_Lam_qinv_bb_1D=(TH1D*)h_Lam_qinv_bb->Projection(0,"e");
TH1D *h_Lam_qinv_mix_ss_1D=(TH1D*)h_Lam_qinv_mix_ss->Projection(0,"e");
TH1D *h_Lam_qinv_mix_sb_1D=(TH1D*)h_Lam_qinv_mix_sb->Projection(0,"e");
TH1D *h_Lam_qinv_mix_bb_1D=(TH1D*)h_Lam_qinv_mix_bb->Projection(0,"e");

TH1D *h_LAL_qinv_ss_1D=(TH1D*)h_LAL_qinv_ss->Projection(0,"e");
TH1D *h_LAL_qinv_sb_1D=(TH1D*)h_LAL_qinv_sb->Projection(0,"e");
TH1D *h_LAL_qinv_bb_1D=(TH1D*)h_LAL_qinv_bb->Projection(0,"e");
TH1D *h_LAL_qinv_mix_ss_1D=(TH1D*)h_LAL_qinv_mix_ss->Projection(0,"e");
TH1D *h_LAL_qinv_mix_sb_1D=(TH1D*)h_LAL_qinv_mix_sb->Projection(0,"e");
TH1D *h_LAL_qinv_mix_bb_1D=(TH1D*)h_LAL_qinv_mix_bb->Projection(0,"e");

TH1D *h_KL_qinv_ss_1D=(TH1D*)h_KL_qinv_ss->Projection(0,"e");
TH1D *h_KL_qinv_sb_1D=(TH1D*)h_KL_qinv_sb->Projection(0,"e");
TH1D *h_KL_qinv_bb_1D=(TH1D*)h_KL_qinv_bb->Projection(0,"e");
TH1D *h_KL_qinv_mix_ss_1D=(TH1D*)h_KL_qinv_mix_ss->Projection(0,"e");
TH1D *h_KL_qinv_mix_sb_1D=(TH1D*)h_KL_qinv_mix_sb->Projection(0,"e");
TH1D *h_KL_qinv_mix_bb_1D=(TH1D*)h_KL_qinv_mix_bb->Projection(0,"e");

TH1D *h_KAL_qinv_ss_1D=(TH1D*)h_KAL_qinv_ss->Projection(0,"e");
TH1D *h_KAL_qinv_sb_1D=(TH1D*)h_KAL_qinv_sb->Projection(0,"e");
TH1D *h_KAL_qinv_bb_1D=(TH1D*)h_KAL_qinv_bb->Projection(0,"e");
TH1D *h_KAL_qinv_mix_ss_1D=(TH1D*)h_KAL_qinv_mix_ss->Projection(0,"e");
TH1D *h_KAL_qinv_mix_sb_1D=(TH1D*)h_KAL_qinv_mix_sb->Projection(0,"e");
TH1D *h_KAL_qinv_mix_bb_1D=(TH1D*)h_KAL_qinv_mix_bb->Projection(0,"e");




TH1D *h_K0s_qinv_ss_1D_kt=(TH1D*)h_K0s_qinv_ss->Projection(1,"e");
TH1D *h_ALam_qinv_ss_1D_kt=(TH1D*)h_ALam_qinv_ss->Projection(1,"e");
TH1D *h_Lam_qinv_ss_1D_kt=(TH1D*)h_Lam_qinv_ss->Projection(1,"e");
TH1D *h_LAL_qinv_ss_1D_kt=(TH1D*)h_LAL_qinv_ss->Projection(1,"e");
TH1D *h_KL_qinv_ss_1D_kt=(TH1D*)h_KL_qinv_ss->Projection(1,"e");
TH1D *h_KAL_qinv_ss_1D_kt=(TH1D*)h_KAL_qinv_ss->Projection(1,"e");
h_Lam_qinv_ss_1D_kt->Add(h_ALam_qinv_ss_1D_kt);
h_KL_qinv_ss_1D_kt->Add(h_KAL_qinv_ss_1D_kt);

TH1D *h_K0s_qinv_ss_1D_ntk=(TH1D*)h_K0s_qinv_ss->Projection(2,"e");
TH1D *h_ALam_qinv_ss_1D_ntk=(TH1D*)h_ALam_qinv_ss->Projection(2,"e");
TH1D *h_Lam_qinv_ss_1D_ntk=(TH1D*)h_Lam_qinv_ss->Projection(2,"e");
TH1D *h_LAL_qinv_ss_1D_ntk=(TH1D*)h_LAL_qinv_ss->Projection(2,"e");
TH1D *h_KL_qinv_ss_1D_ntk=(TH1D*)h_KL_qinv_ss->Projection(2,"e");
TH1D *h_KAL_qinv_ss_1D_ntk=(TH1D*)h_KAL_qinv_ss->Projection(2,"e");
h_Lam_qinv_ss_1D_ntk->Add(h_ALam_qinv_ss_1D_ntk);
h_KL_qinv_ss_1D_ntk->Add(h_KAL_qinv_ss_1D_ntk);

double ntk_K0sK0s = h_K0s_qinv_ss_1D_ntk->GetMean();
double ntk_K0sK0s_err = h_K0s_qinv_ss_1D_ntk->GetMeanError();
double kt_K0sK0s = h_K0s_qinv_ss_1D_kt->GetMean();
double kt_K0sK0s_err = h_K0s_qinv_ss_1D_kt->GetMeanError();

double ntk_K0sL = h_KL_qinv_ss_1D_ntk->GetMean();
double ntk_K0sL_err = h_KL_qinv_ss_1D_ntk->GetMeanError();
double kt_K0sL = h_KL_qinv_ss_1D_kt->GetMean();
double kt_K0sL_err = h_KL_qinv_ss_1D_kt->GetMeanError();

double ntk_LL = h_Lam_qinv_ss_1D_ntk->GetMean();
double ntk_LL_err = h_Lam_qinv_ss_1D_ntk->GetMeanError();
double kt_LL = h_Lam_qinv_ss_1D_kt->GetMean();
double kt_LL_err = h_Lam_qinv_ss_1D_kt->GetMeanError();

double ntk_LAL = h_LAL_qinv_ss_1D_ntk->GetMean();
double ntk_LAL_err = h_LAL_qinv_ss_1D_ntk->GetMeanError();
double kt_LAL = h_LAL_qinv_ss_1D_kt->GetMean();
double kt_LAL_err = h_LAL_qinv_ss_1D_kt->GetMeanError();

double a=1.02, b=0.32, bbeta = 0.49, gamma = 0.18;

double Rpar_K0sK0s = pow(pow(a,2) + pow(b*pow(ntk_K0sK0s,bbeta),2),0.5)*pow((0.2/kt_K0sK0s),gamma);//following Sandra/Ferenc paper
double Rpar_K0sL = pow(pow(a,2) + pow(b*pow(ntk_K0sL,bbeta),2),0.5)*pow((0.2/kt_K0sL),gamma);//following Sandra/Ferenc paper
double Rpar_LL = pow(pow(a,2) + pow(b*pow(ntk_LL,bbeta),2),0.5)*pow((0.2/kt_LL),gamma);//following Sandra/Ferenc paper
double Rpar_LAL = pow(pow(a,2) + pow(b*pow(ntk_LAL,bbeta),2),0.5)*pow((0.2/kt_LAL),gamma);//following Sandra/Ferenc paper

cout << "Ntrk LAL: " << ntk_LAL << " +/- " << ntk_LAL_err << endl;
cout << "kT LAL: " << kt_LAL<< " +/- " << kt_LAL_err << endl;
cout << "Rpar LAL: " << Rpar_LAL << endl;
cout << "Rpar K0sK0s: " << Rpar_K0sK0s << endl;
cout << "Rpar K0sL: " << Rpar_K0sL << endl;
cout << "Rpar LL: " << Rpar_LL << endl;

//rebin -> 40 MeV

h_K0s_qinv_ss_1D->Rebin();
h_K0s_qinv_sb_1D->Rebin();
h_K0s_qinv_bb_1D->Rebin();
h_K0s_qinv_mix_ss_1D->Rebin();
h_K0s_qinv_mix_sb_1D->Rebin();
h_K0s_qinv_mix_bb_1D->Rebin();

h_ALam_qinv_ss_1D->Rebin();
h_ALam_qinv_sb_1D->Rebin();
h_ALam_qinv_bb_1D->Rebin();
h_ALam_qinv_mix_ss_1D->Rebin();
h_ALam_qinv_mix_sb_1D->Rebin();
h_ALam_qinv_mix_bb_1D->Rebin();

h_Lam_qinv_ss_1D->Rebin();
h_Lam_qinv_sb_1D->Rebin();
h_Lam_qinv_bb_1D->Rebin();
h_Lam_qinv_mix_ss_1D->Rebin();
h_Lam_qinv_mix_sb_1D->Rebin();
h_Lam_qinv_mix_bb_1D->Rebin();

h_LAL_qinv_ss_1D->Rebin();
h_LAL_qinv_sb_1D->Rebin();
h_LAL_qinv_bb_1D->Rebin();
h_LAL_qinv_mix_ss_1D->Rebin();
h_LAL_qinv_mix_sb_1D->Rebin();
h_LAL_qinv_mix_bb_1D->Rebin();

h_KL_qinv_ss_1D->Rebin();
h_KL_qinv_sb_1D->Rebin();
h_KL_qinv_bb_1D->Rebin();
h_KL_qinv_mix_ss_1D->Rebin();
h_KL_qinv_mix_sb_1D->Rebin();
h_KL_qinv_mix_bb_1D->Rebin();

h_KAL_qinv_ss_1D->Rebin();
h_KAL_qinv_sb_1D->Rebin();
h_KAL_qinv_bb_1D->Rebin();
h_KAL_qinv_mix_ss_1D->Rebin();
h_KAL_qinv_mix_sb_1D->Rebin();
h_KAL_qinv_mix_bb_1D->Rebin();


//make CF histograms after the

TH1D* K0s_hist_data = (TH1D*) h_K0s_qinv_ss_1D->Clone();
TH1D* K0s_hist_mix_data = (TH1D*) h_K0s_qinv_mix_ss_1D->Clone();
normalize(qmin_Norm,qmax_Norm,K0s_hist_data,K0s_hist_mix_data);
K0s_hist_data->Divide(K0s_hist_mix_data);

TH1D* LAL_hist_data = (TH1D*) h_LAL_qinv_ss_1D->Clone();
TH1D* LAL_hist_mix_data = (TH1D*) h_LAL_qinv_mix_ss_1D->Clone();
normalize(qmin_Norm,qmax_Norm,LAL_hist_data,LAL_hist_mix_data);
LAL_hist_data->Divide(LAL_hist_mix_data);

TH1D* K0sLam_hist_data = (TH1D*) h_KL_qinv_ss_1D->Clone();
TH1D* K0sLam_hist_mix_data = (TH1D*) h_KL_qinv_mix_ss_1D->Clone();
TH1D* K0sALam_hist_data = (TH1D*) h_KAL_qinv_ss_1D->Clone();
TH1D* K0sALam_hist_mix_data = (TH1D*) h_KAL_qinv_mix_ss_1D->Clone();
K0sLam_hist_data->Add(K0sALam_hist_data);
K0sLam_hist_mix_data->Add(K0sALam_hist_mix_data);
normalize(qmin_Norm,qmax_Norm,K0sLam_hist_data,K0sLam_hist_mix_data);
K0sLam_hist_data->Divide(K0sLam_hist_mix_data);

TH1D* LL_hist_data = (TH1D*) h_Lam_qinv_ss_1D->Clone();
TH1D* LL_hist_mix_data = (TH1D*) h_Lam_qinv_mix_ss_1D->Clone();
TH1D* ALAL_hist_data = (TH1D*) h_ALam_qinv_ss_1D->Clone();
TH1D* ALAL_hist_mix_data = (TH1D*) h_ALam_qinv_mix_ss_1D->Clone();
LL_hist_data->Add(ALAL_hist_data);
LL_hist_mix_data->Add(ALAL_hist_mix_data);
normalize(qmin_Norm,qmax_Norm,LL_hist_data,LL_hist_mix_data);
LL_hist_data->Divide(LL_hist_mix_data);

//K0sK0s
TH1D* K0s_hist = (TH1D*) h_K0s_qinv_ss_1D->Clone();
TH1D* K0s_ss = (TH1D*) h_K0s_qinv_ss_1D->Clone();
TH1D* K0s_sb = (TH1D*) h_K0s_qinv_sb_1D->Clone();
TH1D* K0s_bb = (TH1D*) h_K0s_qinv_bb_1D->Clone();

TH1D* K0s_hist_mix = (TH1D*) h_K0s_qinv_mix_ss_1D->Clone();
TH1D* K0s_ss_mix = (TH1D*) h_K0s_qinv_mix_ss_1D->Clone();
TH1D* K0s_sb_mix = (TH1D*) h_K0s_qinv_mix_sb_1D->Clone();
TH1D* K0s_bb_mix = (TH1D*) h_K0s_qinv_mix_bb_1D->Clone();



K0s_hist->Scale(1.0/K0s_hist->Integral());
K0s_ss->Scale(fss/K0s_ss->Integral());
K0s_sb->Scale(fsb/K0s_sb->Integral());
K0s_bb->Scale(fbb/K0s_bb->Integral());

K0s_hist_mix->Scale(1.0/K0s_hist_mix->Integral());
K0s_ss_mix->Scale(fss/K0s_ss_mix->Integral());
K0s_sb_mix->Scale(fsb/K0s_sb_mix->Integral());
K0s_bb_mix->Scale(fbb/K0s_bb_mix->Integral());

//K0sK0s -> sideband plot

TCanvas *c_KsKs_sidemet = new TCanvas("c_KsKs_sidemet","",800,600);
gPad->SetTickx();
gPad->SetTicky();
std::ostringstream oKsKssidemet; 
oKsKssidemet << 1000*K0s_hist->GetBinWidth(1) << std::setprecision(2);
K0s_hist->Draw("e");
K0s_hist->GetXaxis()->SetRangeUser(0.0,2.9999999);
K0s_hist->GetYaxis()->SetRangeUser(50e-10,0.5e-1);
K0s_hist->SetMarkerStyle(24);
K0s_hist->SetMarkerSize(1.5);
K0s_hist->SetMarkerColor(kBlack);
K0s_hist->SetLineColor(kBlack);
K0s_ss->Draw("histsame");
K0s_ss->SetLineColor(SPcolors[SPdred]);
K0s_ss->SetLineWidth(4);
K0s_sb->Draw("histsame");
K0s_sb->SetLineColor(SPcolors[SPgreen]);
K0s_sb->SetLineStyle(2);
K0s_sb->SetLineWidth(4);
K0s_bb->Draw("histsame");
K0s_bb->SetLineColor(SPcolors[SPblue]);
K0s_bb->SetLineStyle(9);
K0s_bb->SetLineWidth(4);
K0s_hist->GetYaxis()->SetTitle(Form("#scale[1.2]{Normalized Distributions/ %s MeV}",oKsKssidemet.str(  ).c_str(  )));
K0s_hist->GetYaxis()->CenterTitle(1);
K0s_hist->GetYaxis()->SetNdivisions(810);
K0s_hist->GetXaxis()->SetTickSize(0.03);
K0s_hist->GetYaxis()->SetTickSize(0.03);
K0s_hist->GetXaxis()->SetLabelSize(0.042);
K0s_hist->GetYaxis()->SetLabelSize(0.042);
K0s_hist->GetYaxis()->SetLabelOffset(0.01);
K0s_hist->GetYaxis()->SetTitleOffset(1.4);
K0s_hist->GetXaxis()->SetTitleOffset(1.3);
K0s_hist->GetXaxis()->SetTitle("#scale[1.25]{q_{inv} [GeV]}");
K0s_hist->GetXaxis()->CenterTitle(1);

TLegend* legKs_sig = new TLegend(0.55,0.15,0.86,0.4);
legKs_sig->SetBorderSize(0);
legKs_sig->SetTextSize(0.04);
legKs_sig->AddEntry(K0s_hist,"D(q_{inv}^{meas})");
legKs_sig->AddEntry(K0s_ss,"f_{ss}D(q_{inv}^{ss})","L");
legKs_sig->AddEntry(K0s_sb,"f_{sb}D(q_{inv}^{sb})","L");
legKs_sig->AddEntry(K0s_bb,"f_{bb}D(q_{inv}^{bb})","L");
legKs_sig->Draw();

std::ostringstream os_KsKssidemet; 
TLatex* texCMS_yKsKssidemet = new TLatex;
texCMS_yKsKssidemet->SetNDC();
texCMS_yKsKssidemet->SetTextFont(42);
texCMS_yKsKssidemet->SetTextSize(0.05);
texCMS_yKsKssidemet->SetTextAlign(12);
texCMS_yKsKssidemet->DrawLatex(.12,.95," #font[61]{CMS} #it{Work in Progress}");
texCMS_yKsKssidemet->DrawLatex(0.60,0.95, "#scale[0.9]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS_yKsKssidemet->DrawLatex(0.15,0.83, "#scale[1.2]{K^{0}_{S}K^{0}_{S}}");  
//texCMS_yKsKssidemet->DrawLatex(0.15,0.83, "#scale[1.0]{Signal}");  
os_KsKssidemet << "#scale[0.9]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_yKsKssidemet->DrawLatex(0.25,0.30, os_KsKssidemet.str(  ).c_str(  )); 
os_KsKssidemet.str(std::string()); 
os_KsKssidemet << "#scale[0.9]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_yKsKssidemet->DrawLatex(0.25,0.24, os_KsKssidemet.str(  ).c_str(  ));  
os_KsKssidemet.str(std::string()); 
gPad->SetLogy();
c_KsKs_sidemet->Update();
c_KsKs_sidemet->Print("results/sidemet_log_sig_KsKs_"+V0+"_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");


TH1D* h_sig_K0sK0s_befpur_ss = (TH1D*)h_K0s_qinv_ss_1D->Clone();
TH1D* h_sig_K0sK0s_befpur_sb = (TH1D*)h_K0s_qinv_sb_1D->Clone();
TH1D* h_sig_K0sK0s_befpur_bb = (TH1D*)h_K0s_qinv_bb_1D->Clone();
TH1D* h_mix_K0sK0s_befpur_ss = (TH1D*)h_K0s_qinv_mix_ss_1D->Clone();
TH1D* h_mix_K0sK0s_befpur_sb = (TH1D*)h_K0s_qinv_mix_sb_1D->Clone();
TH1D* h_mix_K0sK0s_befpur_bb = (TH1D*)h_K0s_qinv_mix_bb_1D->Clone();

sideband(h_K0s_qinv_ss_1D, h_K0s_qinv_sb_1D, h_K0s_qinv_bb_1D, fss, fsb, fbb);
sideband(h_K0s_qinv_mix_ss_1D, h_K0s_qinv_mix_sb_1D, h_K0s_qinv_mix_bb_1D, fss, fsb, fbb);

TH1D* h_sig_K0sK0s = (TH1D*)h_K0s_qinv_ss_1D->Clone();
TH1D* h_mix_K0sK0s = (TH1D*)h_K0s_qinv_mix_ss_1D->Clone();

normalize(qmin_Norm,qmax_Norm,h_K0s_qinv_ss_1D,h_K0s_qinv_mix_ss_1D);
TGraphAsymmErrors h_K0s_qinv_ss_graph;
h_K0s_qinv_ss_graph.Divide(h_K0s_qinv_ss_1D, h_K0s_qinv_mix_ss_1D, "pois");
h_K0s_qinv_ss_1D->Divide(h_K0s_qinv_mix_ss_1D);


//LamLam+ALamALam
TH1D* h_Lam_plus_ALam_qinv_ss_1D = (TH1D*)h_Lam_qinv_ss_1D->Clone();
TH1D* h_Lam_plus_ALam_qinv_sb_1D = (TH1D*)h_Lam_qinv_sb_1D->Clone();
TH1D* h_Lam_plus_ALam_qinv_bb_1D = (TH1D*)h_Lam_qinv_bb_1D->Clone();
TH1D* h_Lam_plus_ALam_qinv_ss_1D_mix = (TH1D*)h_Lam_qinv_mix_ss_1D->Clone();
TH1D* h_Lam_plus_ALam_qinv_sb_1D_mix = (TH1D*)h_Lam_qinv_mix_sb_1D->Clone();
TH1D* h_Lam_plus_ALam_qinv_bb_1D_mix = (TH1D*)h_Lam_qinv_mix_bb_1D->Clone();
h_Lam_plus_ALam_qinv_ss_1D->Add(h_ALam_qinv_ss_1D);
h_Lam_plus_ALam_qinv_sb_1D->Add(h_ALam_qinv_sb_1D);
h_Lam_plus_ALam_qinv_bb_1D->Add(h_ALam_qinv_bb_1D);
h_Lam_plus_ALam_qinv_ss_1D_mix->Add(h_ALam_qinv_mix_ss_1D);
h_Lam_plus_ALam_qinv_sb_1D_mix->Add(h_ALam_qinv_mix_sb_1D);
h_Lam_plus_ALam_qinv_bb_1D_mix->Add(h_ALam_qinv_mix_bb_1D);

//LLLL -> sideband plot
TH1D* LL_hist = (TH1D*) h_Lam_plus_ALam_qinv_ss_1D->Clone();
TH1D* LL_ss = (TH1D*) h_Lam_plus_ALam_qinv_ss_1D->Clone();
TH1D* LL_sb = (TH1D*) h_Lam_plus_ALam_qinv_sb_1D->Clone();
TH1D* LL_bb = (TH1D*) h_Lam_plus_ALam_qinv_bb_1D->Clone();

TH1D* LL_hist_mix = (TH1D*) h_Lam_plus_ALam_qinv_ss_1D_mix->Clone();
TH1D* LL_ss_mix = (TH1D*) h_Lam_plus_ALam_qinv_ss_1D_mix->Clone();
TH1D* LL_sb_mix = (TH1D*) h_Lam_plus_ALam_qinv_sb_1D_mix->Clone();
TH1D* LL_bb_mix = (TH1D*) h_Lam_plus_ALam_qinv_bb_1D_mix->Clone();

LL_hist->Scale(1.0/LL_hist->Integral());
LL_ss->Scale(fss_LLALAL/LL_ss->Integral());
LL_sb->Scale(fsb_LLALAL/LL_sb->Integral());
LL_bb->Scale(fbb_LLALAL/LL_bb->Integral());

LL_hist_mix->Scale(1.0/LL_hist_mix->Integral());
LL_ss_mix->Scale(fss_LLALAL/LL_ss_mix->Integral());
LL_sb_mix->Scale(fsb_LLALAL/LL_sb_mix->Integral());
LL_bb_mix->Scale(fbb_LLALAL/LL_bb_mix->Integral());

//LLLL

TCanvas *c_LL_sidemet = new TCanvas("c_LL_sidemet","",800,600);
gPad->SetTickx();
gPad->SetTicky();
std::ostringstream oLLsidemet; 
oLLsidemet << 1000*LL_hist->GetBinWidth(1) << std::setprecision(2);
LL_hist->Draw("e");
LL_hist->GetXaxis()->SetRangeUser(0.0,2.9999999);
LL_hist->GetYaxis()->SetRangeUser(50e-10,0.5e-1);
LL_hist->SetMarkerStyle(24);
LL_hist->SetMarkerSize(1.5);
LL_hist->SetMarkerColor(kBlack);
LL_hist->SetLineColor(kBlack);
LL_ss->Draw("histsame");
LL_ss->SetLineColor(SPcolors[SPdred]);
LL_ss->SetLineWidth(4);
LL_sb->Draw("histsame");
LL_sb->SetLineColor(SPcolors[SPgreen]);
LL_sb->SetLineStyle(2);
LL_sb->SetLineWidth(4);
LL_bb->Draw("histsame");
LL_bb->SetLineColor(SPcolors[SPblue]);
LL_bb->SetLineStyle(9);
LL_bb->SetLineWidth(4);
LL_hist->GetYaxis()->SetTitle(Form("#scale[1.2]{Normalized Distributions/ %s MeV}",oLLsidemet.str(  ).c_str(  )));
LL_hist->GetYaxis()->CenterTitle(1);
LL_hist->GetYaxis()->SetNdivisions(810);
LL_hist->GetXaxis()->SetTickSize(0.03);
LL_hist->GetYaxis()->SetTickSize(0.03);
LL_hist->GetXaxis()->SetLabelSize(0.042);
LL_hist->GetYaxis()->SetLabelSize(0.042);
LL_hist->GetYaxis()->SetLabelOffset(0.01);
LL_hist->GetYaxis()->SetTitleOffset(1.4);
LL_hist->GetXaxis()->SetTitleOffset(1.3);
LL_hist->GetXaxis()->SetTitle("#scale[1.25]{q_{inv} [GeV]}");
LL_hist->GetXaxis()->CenterTitle(1);

TLegend* legLL_sig = new TLegend(0.55,0.15,0.86,0.4);
legLL_sig->SetBorderSize(0);
legLL_sig->SetTextSize(0.04);
legLL_sig->AddEntry(LL_hist,"D(q_{inv}^{meas})");
legLL_sig->AddEntry(LL_ss,"f_{ss}D(q_{inv}^{ss})","L");
legLL_sig->AddEntry(LL_sb,"f_{sb}D(q_{inv}^{sb})","L");
legLL_sig->AddEntry(LL_bb,"f_{bb}D(q_{inv}^{bb})","L");
legLL_sig->Draw();

std::ostringstream os_LLsidemet; 
TLatex* texCMS_yLLsidemet = new TLatex;
texCMS_yLLsidemet->SetNDC();
texCMS_yLLsidemet->SetTextFont(42);
texCMS_yLLsidemet->SetTextSize(0.05);
texCMS_yLLsidemet->SetTextAlign(12);
texCMS_yLLsidemet->DrawLatex(.12,.95," #font[61]{CMS} #it{Work in Progress}");
texCMS_yLLsidemet->DrawLatex(0.60,0.95, "#scale[0.9]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS_yLLsidemet->DrawLatex(0.15,0.83, "#scale[1.2]{#Lambda#Lambda #oplus #bar{#Lambda}#bar{#Lambda}}");  
//texCMS_yLLsidemet->DrawLatex(0.15,0.83, "#scale[1.0]{Signal}");  
os_LLsidemet << "#scale[0.9]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_yLLsidemet->DrawLatex(0.25,0.30, os_LLsidemet.str(  ).c_str(  )); 
os_LLsidemet.str(std::string()); 
os_LLsidemet << "#scale[0.9]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_yLLsidemet->DrawLatex(0.25,0.24, os_LLsidemet.str(  ).c_str(  ));  
os_LLsidemet.str(std::string()); 
gPad->SetLogy();
c_LL_sidemet->Update();
c_LL_sidemet->Print("results/sidemet_log_sig_LL_"+V0+"_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");

TH1D* h_sig_LLALAL_befpur_ss = (TH1D*)h_Lam_plus_ALam_qinv_ss_1D->Clone();
TH1D* h_sig_LLALAL_befpur_sb = (TH1D*)h_Lam_plus_ALam_qinv_sb_1D->Clone();
TH1D* h_sig_LLALAL_befpur_bb = (TH1D*)h_Lam_plus_ALam_qinv_bb_1D->Clone();
TH1D* h_mix_LLALAL_befpur_ss = (TH1D*)h_Lam_plus_ALam_qinv_ss_1D_mix->Clone();
TH1D* h_mix_LLALAL_befpur_sb = (TH1D*)h_Lam_plus_ALam_qinv_sb_1D_mix->Clone();
TH1D* h_mix_LLALAL_befpur_bb = (TH1D*)h_Lam_plus_ALam_qinv_bb_1D_mix->Clone();

sideband(h_Lam_plus_ALam_qinv_ss_1D, h_Lam_plus_ALam_qinv_sb_1D, h_Lam_plus_ALam_qinv_bb_1D, fss_LLALAL, fsb_LLALAL, fbb_LLALAL);
sideband(h_Lam_plus_ALam_qinv_ss_1D_mix, h_Lam_plus_ALam_qinv_sb_1D_mix, h_Lam_plus_ALam_qinv_bb_1D_mix, fss_LLALAL, fsb_LLALAL, fbb_LLALAL);

TH1D* h_sig_LLALAL = (TH1D*)h_Lam_plus_ALam_qinv_ss_1D->Clone();
TH1D* h_mix_LLALAL = (TH1D*)h_Lam_plus_ALam_qinv_ss_1D_mix->Clone();

normalize(qmin_Norm,qmax_Norm,h_Lam_plus_ALam_qinv_ss_1D,h_Lam_plus_ALam_qinv_ss_1D_mix);
TGraphAsymmErrors h_Lam_qinv_ss_graph;
h_Lam_qinv_ss_graph.Divide(h_Lam_plus_ALam_qinv_ss_1D, h_Lam_plus_ALam_qinv_ss_1D_mix, "pois");
h_Lam_plus_ALam_qinv_ss_1D->Divide(h_Lam_plus_ALam_qinv_ss_1D_mix);

//LamLam
sideband(h_Lam_qinv_ss_1D, h_Lam_qinv_sb_1D, h_Lam_qinv_bb_1D, fss_LL, fsb_LL, fbb_LL);
sideband(h_Lam_qinv_mix_ss_1D, h_Lam_qinv_mix_sb_1D, h_Lam_qinv_mix_bb_1D, fss_LL, fsb_LL, fbb_LL);
normalize(qmin_Norm,qmax_Norm,h_Lam_qinv_ss_1D,h_Lam_qinv_mix_ss_1D);
h_Lam_qinv_ss_1D->Divide(h_Lam_qinv_mix_ss_1D);

//ALamALam
sideband(h_ALam_qinv_ss_1D, h_ALam_qinv_sb_1D, h_ALam_qinv_bb_1D, fss_ALAL, fsb_ALAL, fbb_ALAL);
sideband(h_ALam_qinv_mix_ss_1D, h_ALam_qinv_mix_sb_1D, h_ALam_qinv_mix_bb_1D, fss_ALAL, fsb_ALAL, fbb_ALAL);
normalize(qmin_Norm,qmax_Norm,h_ALam_qinv_ss_1D,h_ALam_qinv_mix_ss_1D);
h_ALam_qinv_ss_1D->Divide(h_ALam_qinv_mix_ss_1D);


//=================================================================
//Lam ALam cross checks
//=================================================================
TCanvas *c_x = new TCanvas("c_x","",600,600);
TPad* pad1_x = new TPad("pad1_x","top pad",0.0,0.3,1.0,1.0);
TPad* pad2_x = new TPad("pad2_x","bottom pad",0.0,0.0,1.0,0.3);
pad1_x->SetTopMargin(0.1);
pad1_x->SetBottomMargin(0.0);
pad1_x->SetRightMargin(0.038);
pad1_x->SetLeftMargin(0.15);
pad1_x->Draw();
pad2_x->SetFrameFillStyle(4000);
pad2_x->SetTopMargin(0.0);
pad2_x->SetBottomMargin(0.37);
pad2_x->SetRightMargin(0.038);
pad2_x->SetLeftMargin(0.15);
pad2_x->Draw();
pad1_x->cd();
gPad->SetTickx();
gPad->SetTicky();
std::ostringstream oLAL; 
oLAL << 1000*h_Lam_qinv_ss_1D->GetBinWidth(1) << std::setprecision(2);
h_Lam_qinv_ss_1D->Draw("e");
h_Lam_qinv_ss_1D->GetXaxis()->SetRangeUser(0.0,2.9999999);
h_Lam_qinv_ss_1D->SetMarkerStyle(24);
h_Lam_qinv_ss_1D->SetMarkerSize(1.5);
h_Lam_qinv_ss_1D->SetMarkerColor(SPcolors[SPdred]);
h_Lam_qinv_ss_1D->SetLineColor(SPcolors[SPdred]);
h_Lam_qinv_ss_1D->GetYaxis()->SetTitle(Form("#scale[1.7]{Single Ratio / %s MeV}",oLAL.str(  ).c_str(  )));
h_Lam_qinv_ss_1D->GetYaxis()->CenterTitle(1);
h_Lam_qinv_ss_1D->GetYaxis()->SetNdivisions(810);
h_Lam_qinv_ss_1D->GetXaxis()->SetTickSize(0.03);
h_Lam_qinv_ss_1D->GetYaxis()->SetTickSize(0.02);
h_Lam_qinv_ss_1D->GetXaxis()->SetLabelSize(0.055);
h_Lam_qinv_ss_1D->GetYaxis()->SetLabelSize(0.055);
h_Lam_qinv_ss_1D->GetYaxis()->SetLabelOffset(0.01);
h_Lam_qinv_ss_1D->GetYaxis()->SetTitleOffset(1.8);
h_ALam_qinv_ss_1D->Draw("same");
h_ALam_qinv_ss_1D->SetMarkerStyle(25);
h_ALam_qinv_ss_1D->SetMarkerColor(1);
h_ALam_qinv_ss_1D->SetLineColor(1);
h_ALam_qinv_ss_1D->SetMarkerSize(1.5);
TLegend* leg0 = new TLegend(0.70,0.30,0.90,0.45);
leg0->SetBorderSize(0);
leg0->AddEntry(h_Lam_qinv_ss_1D,"#Lambda#Lambda");
leg0->AddEntry(h_ALam_qinv_ss_1D,"#bar{#Lambda}#bar{#Lambda}");
leg0->Draw();
std::ostringstream os_x; 
TLatex* texCMS_x = new TLatex;
texCMS_x->SetNDC();
texCMS_x->SetTextFont(42);
texCMS_x->SetTextSize(0.05);
texCMS_x->SetTextAlign(12);
texCMS_x->DrawLatex(.20,.95," #font[61]{CMS} #it{Work in Progress}");
texCMS_x->DrawLatex(0.60,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
os_x << "#scale[1.0]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_x->DrawLatex(0.68,0.60, os_x.str(  ).c_str(  )); 
os_x.str(std::string()); 
os_x << "#scale[1.0]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_x->DrawLatex(0.68,0.52, os_x.str(  ).c_str(  ));  
os_x.str(std::string()); 
TLine* line_xx = new TLine(0.0, 1, 2.99999, 1);
line_xx->SetLineStyle(2);
line_xx->Draw("same");
c_x->Update();
pad2_x->cd();
gPad->SetTickx();
gPad->SetTicky();
TH1D* ratio = (TH1D*) h_Lam_qinv_ss_1D->Clone();
ratio->Divide(h_Lam_qinv_ss_1D,h_ALam_qinv_ss_1D,1,1);
double low_x = 0.49999;
double high_x = 1.499999999;
TH1F* frame_ks_x = pad2_x->cd()->DrawFrame(0.00000,low_x,2.99999,high_x);
frame_ks_x->GetXaxis()->SetTitle("#scale[4.3]{q_{inv} [GeV]}");
frame_ks_x->GetYaxis()->SetTitle("#scale[4.0]{#Lambda#Lambda/#bar{#Lambda}#bar{#Lambda}}");
frame_ks_x->GetXaxis()->CenterTitle(1);
frame_ks_x->GetYaxis()->CenterTitle(1);
frame_ks_x->GetXaxis()->SetNdivisions(510);
frame_ks_x->GetYaxis()->SetNdivisions(406);
frame_ks_x->GetXaxis()->SetTickSize(0.06);
frame_ks_x->GetYaxis()->SetTickSize(0.02);
frame_ks_x->GetXaxis()->SetTitleOffset(4.5);
frame_ks_x->GetYaxis()->SetTitleOffset(1.8);
frame_ks_x->GetXaxis()->SetLabelSize(0.14);
frame_ks_x->GetYaxis()->SetLabelSize(0.14);
frame_ks_x->GetXaxis()->SetLabelOffset(0.01);
frame_ks_x->GetYaxis()->SetLabelOffset(0.013);
ratio->Draw("esame");
ratio->SetMarkerStyle(20);
ratio->SetMarkerSize(1.5);
ratio->SetMarkerColor(kBlack);
ratio->SetLineColor(kBlack);
TLine* line_x = new TLine(0.0, 1, 2.99999, 1);
line_x->SetLineStyle(2);
line_x->Draw("same");
c_x->Update();
c_x->Print("results/xcheck_LAL_"+V0+"_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");
//=================================================================


//LALLAL
TH1D* LAL_hist = (TH1D*) h_LAL_qinv_ss_1D->Clone();
TH1D* LAL_ss = (TH1D*) h_LAL_qinv_ss_1D->Clone();
TH1D* LAL_sb = (TH1D*) h_LAL_qinv_sb_1D->Clone();
TH1D* LAL_bb = (TH1D*) h_LAL_qinv_bb_1D->Clone();

TH1D* LAL_hist_mix = (TH1D*) h_LAL_qinv_mix_ss_1D->Clone();
TH1D* LAL_ss_mix = (TH1D*) h_LAL_qinv_mix_ss_1D->Clone();
TH1D* LAL_sb_mix = (TH1D*) h_LAL_qinv_mix_sb_1D->Clone();
TH1D* LAL_bb_mix = (TH1D*) h_LAL_qinv_mix_bb_1D->Clone();

LAL_hist->Scale(1.0/LAL_hist->Integral());
LAL_ss->Scale(fss_LAL/LAL_ss->Integral());
LAL_sb->Scale(fsb_LAL/LAL_sb->Integral());
LAL_bb->Scale(fbb_LAL/LAL_bb->Integral());

LAL_hist_mix->Scale(1.0/LAL_hist_mix->Integral());
LAL_ss_mix->Scale(fss_LAL/LAL_ss_mix->Integral());
LAL_sb_mix->Scale(fsb_LAL/LAL_sb_mix->Integral());
LAL_bb_mix->Scale(fbb_LAL/LAL_bb_mix->Integral());

//LALLAL

TCanvas *c_LAL_sidemet = new TCanvas("c_LAL_sidemet","",800,600);
gPad->SetTickx();
gPad->SetTicky();
std::ostringstream oLALsidemet; 
oLALsidemet << 1000*LAL_hist->GetBinWidth(1) << std::setprecision(2);
LAL_hist->Draw("e");
LAL_hist->GetXaxis()->SetRangeUser(0.0,2.9999999);
LAL_hist->GetYaxis()->SetRangeUser(50e-10,0.5e-1);
LAL_hist->SetMarkerStyle(24);
LAL_hist->SetMarkerSize(1.5);
LAL_hist->SetMarkerColor(kBlack);
LAL_hist->SetLineColor(kBlack);
LAL_ss->Draw("histsame");
LAL_ss->SetLineColor(SPcolors[SPdred]);
LAL_ss->SetLineWidth(4);
LAL_sb->Draw("histsame");
LAL_sb->SetLineColor(SPcolors[SPgreen]);
LAL_sb->SetLineStyle(2);
LAL_sb->SetLineWidth(4);
LAL_bb->Draw("histsame");
LAL_bb->SetLineColor(SPcolors[SPblue]);
LAL_bb->SetLineStyle(9);
LAL_bb->SetLineWidth(4);
LAL_hist->GetYaxis()->SetTitle(Form("#scale[1.2]{Normalized Distributions/ %s MeV}",oLALsidemet.str(  ).c_str(  )));
LAL_hist->GetYaxis()->CenterTitle(1);
LAL_hist->GetYaxis()->SetNdivisions(810);
LAL_hist->GetXaxis()->SetTickSize(0.03);
LAL_hist->GetYaxis()->SetTickSize(0.03);
LAL_hist->GetXaxis()->SetLabelSize(0.042);
LAL_hist->GetYaxis()->SetLabelSize(0.042);
LAL_hist->GetYaxis()->SetLabelOffset(0.01);
LAL_hist->GetYaxis()->SetTitleOffset(1.4);
LAL_hist->GetXaxis()->SetTitleOffset(1.3);
LAL_hist->GetXaxis()->SetTitle("#scale[1.25]{q_{inv} [GeV]}");
LAL_hist->GetXaxis()->CenterTitle(1);

TLegend* legLAL_sig = new TLegend(0.55,0.15,0.86,0.4);
legLAL_sig->SetBorderSize(0);
legLAL_sig->SetTextSize(0.04);
legLAL_sig->AddEntry(LAL_hist,"D(q_{inv}^{meas})");
legLAL_sig->AddEntry(LAL_ss,"f_{ss}D(q_{inv}^{ss})","L");
legLAL_sig->AddEntry(LAL_sb,"f_{sb}D(q_{inv}^{sb})","L");
legLAL_sig->AddEntry(LAL_bb,"f_{bb}D(q_{inv}^{bb})","L");
legLAL_sig->Draw();

std::ostringstream os_LALsidemet; 
TLatex* texCMS_yLALsidemet = new TLatex;
texCMS_yLALsidemet->SetNDC();
texCMS_yLALsidemet->SetTextFont(42);
texCMS_yLALsidemet->SetTextSize(0.05);
texCMS_yLALsidemet->SetTextAlign(12);
texCMS_yLALsidemet->DrawLatex(.12,.95," #font[61]{CMS} #it{Work in Progress}");
texCMS_yLALsidemet->DrawLatex(0.60,0.95, "#scale[0.9]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS_yLALsidemet->DrawLatex(0.15,0.83, "#scale[1.2]{#Lambda#bar{#Lambda}}");  
//texCMS_yLALsidemet->DrawLatex(0.15,0.83, "#scale[1.0]{Signal}");  
os_LALsidemet << "#scale[0.9]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_yLALsidemet->DrawLatex(0.25,0.30, os_LALsidemet.str(  ).c_str(  )); 
os_LALsidemet.str(std::string()); 
os_LALsidemet << "#scale[0.9]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_yLALsidemet->DrawLatex(0.25,0.24, os_LALsidemet.str(  ).c_str(  ));  
os_LALsidemet.str(std::string()); 
gPad->SetLogy();
c_LAL_sidemet->Update();
c_LAL_sidemet->Print("results/sidemet_log_sig_LAL_"+V0+"_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");

//LAL
double LAL_sig_int = h_LAL_qinv_ss_1D->Integral();
double LAL_bkg_int = h_LAL_qinv_mix_ss_1D->Integral();

TH1D* h_sig_LAL_befpur_ss = (TH1D*)h_LAL_qinv_ss_1D->Clone();
TH1D* h_sig_LAL_befpur_sb = (TH1D*)h_LAL_qinv_sb_1D->Clone();
TH1D* h_sig_LAL_befpur_bb = (TH1D*)h_LAL_qinv_bb_1D->Clone();
TH1D* h_mix_LAL_befpur_ss = (TH1D*)h_LAL_qinv_mix_ss_1D->Clone();
TH1D* h_mix_LAL_befpur_sb = (TH1D*)h_LAL_qinv_mix_sb_1D->Clone();
TH1D* h_mix_LAL_befpur_bb = (TH1D*)h_LAL_qinv_mix_bb_1D->Clone();

sideband(h_LAL_qinv_ss_1D, h_LAL_qinv_sb_1D, h_LAL_qinv_bb_1D, fss_LAL, fsb_LAL, fbb_LAL);
sideband(h_LAL_qinv_mix_ss_1D, h_LAL_qinv_mix_sb_1D, h_LAL_qinv_mix_bb_1D, fss_LAL, fsb_LAL, fbb_LAL);

TH1D* h_sig_LAL = (TH1D*)h_LAL_qinv_ss_1D->Clone();
TH1D* h_mix_LAL = (TH1D*)h_LAL_qinv_mix_ss_1D->Clone();

normalize(qmin_Norm,qmax_Norm,h_LAL_qinv_ss_1D,h_LAL_qinv_mix_ss_1D);
TH1D* A_LAL = (TH1D*)h_LAL_qinv_ss_1D->Clone();
TH1D* B_LAL = (TH1D*)h_LAL_qinv_mix_ss_1D->Clone();
A_LAL->Scale(LAL_sig_int);
B_LAL->Scale(LAL_bkg_int);
TGraphAsymmErrors h_LAL_qinv_ss_graph;
h_LAL_qinv_ss_graph.Divide(h_LAL_qinv_ss_1D, h_LAL_qinv_mix_ss_1D, "pois");
h_LAL_qinv_ss_1D->Divide(h_LAL_qinv_mix_ss_1D);



//KLam+KALam

TH1D* h_K_plus_LAL_qinv_ss_1D = (TH1D*)h_KL_qinv_ss_1D->Clone();
TH1D* h_K_plus_LAL_qinv_sb_1D = (TH1D*)h_KL_qinv_sb_1D->Clone();
TH1D* h_K_plus_LAL_qinv_bb_1D = (TH1D*)h_KL_qinv_bb_1D->Clone();
TH1D* h_K_plus_LAL_qinv_ss_1D_mix = (TH1D*)h_KL_qinv_mix_ss_1D->Clone();
TH1D* h_K_plus_LAL_qinv_sb_1D_mix = (TH1D*)h_KL_qinv_mix_sb_1D->Clone();
TH1D* h_K_plus_LAL_qinv_bb_1D_mix = (TH1D*)h_KL_qinv_mix_bb_1D->Clone();
h_K_plus_LAL_qinv_ss_1D->Add(h_KAL_qinv_ss_1D);
h_K_plus_LAL_qinv_sb_1D->Add(h_KAL_qinv_sb_1D);
h_K_plus_LAL_qinv_bb_1D->Add(h_KAL_qinv_bb_1D);
h_K_plus_LAL_qinv_ss_1D_mix->Add(h_KAL_qinv_mix_ss_1D);
h_K_plus_LAL_qinv_sb_1D_mix->Add(h_KAL_qinv_mix_sb_1D);
h_K_plus_LAL_qinv_bb_1D_mix->Add(h_KAL_qinv_mix_bb_1D);


//KLKL
TH1D* KL_hist = (TH1D*) h_K_plus_LAL_qinv_ss_1D->Clone();
TH1D* KL_ss = (TH1D*) h_K_plus_LAL_qinv_ss_1D->Clone();
TH1D* KL_sb = (TH1D*) h_K_plus_LAL_qinv_sb_1D->Clone();
TH1D* KL_bb = (TH1D*) h_K_plus_LAL_qinv_bb_1D->Clone();

TH1D* KL_hist_mix = (TH1D*) h_K_plus_LAL_qinv_ss_1D_mix->Clone();
TH1D* KL_ss_mix = (TH1D*) h_K_plus_LAL_qinv_ss_1D_mix->Clone();
TH1D* KL_sb_mix = (TH1D*) h_K_plus_LAL_qinv_sb_1D_mix->Clone();
TH1D* KL_bb_mix = (TH1D*) h_K_plus_LAL_qinv_bb_1D_mix->Clone();

KL_hist->Scale(1.0/KL_hist->Integral());
KL_ss->Scale(fss_KLAL/KL_ss->Integral());
KL_sb->Scale(fsb_KLAL/KL_sb->Integral());
KL_bb->Scale(fbb_KLAL/KL_bb->Integral());

KL_hist_mix->Scale(1.0/KL_hist_mix->Integral());
KL_ss_mix->Scale(fss_KLAL/KL_ss_mix->Integral());
KL_sb_mix->Scale(fsb_KLAL/KL_sb_mix->Integral());
KL_bb_mix->Scale(fbb_KLAL/KL_bb_mix->Integral());

//KLKL

TCanvas *c_KL_sidemet = new TCanvas("c_KL_sidemet","",800,600);
gPad->SetTickx();
gPad->SetTicky();
std::ostringstream oKLsidemet; 
oKLsidemet << 1000*KL_hist->GetBinWidth(1) << std::setprecision(2);
KL_hist->Draw("e");
KL_hist->GetXaxis()->SetRangeUser(0.0,2.9999999);
KL_hist->GetYaxis()->SetRangeUser(50e-10,0.5e-1);
KL_hist->SetMarkerStyle(24);
KL_hist->SetMarkerSize(1.5);
KL_hist->SetMarkerColor(kBlack);
KL_hist->SetLineColor(kBlack);
KL_ss->Draw("histsame");
KL_ss->SetLineColor(SPcolors[SPdred]);
KL_ss->SetLineWidth(4);
KL_sb->Draw("histsame");
KL_sb->SetLineColor(SPcolors[SPgreen]);
KL_sb->SetLineStyle(2);
KL_sb->SetLineWidth(4);
KL_bb->Draw("histsame");
KL_bb->SetLineColor(SPcolors[SPblue]);
KL_bb->SetLineStyle(9);
KL_bb->SetLineWidth(4);
KL_hist->GetYaxis()->SetTitle(Form("#scale[1.2]{Normalized Distributions/ %s MeV}",oKLsidemet.str(  ).c_str(  )));
KL_hist->GetYaxis()->CenterTitle(1);
KL_hist->GetYaxis()->SetNdivisions(810);
KL_hist->GetXaxis()->SetTickSize(0.03);
KL_hist->GetYaxis()->SetTickSize(0.03);
KL_hist->GetXaxis()->SetLabelSize(0.042);
KL_hist->GetYaxis()->SetLabelSize(0.042);
KL_hist->GetYaxis()->SetLabelOffset(0.01);
KL_hist->GetYaxis()->SetTitleOffset(1.4);
KL_hist->GetXaxis()->SetTitleOffset(1.3);
KL_hist->GetXaxis()->SetTitle("#scale[1.25]{q_{inv} [GeV]}");
KL_hist->GetXaxis()->CenterTitle(1);

TLegend* legKL_sig = new TLegend(0.55,0.15,0.86,0.4);
legKL_sig->SetBorderSize(0);
legKL_sig->SetTextSize(0.04);
legKL_sig->AddEntry(KL_hist,"D(q_{inv}^{meas})");
legKL_sig->AddEntry(KL_ss,"f_{ss}D(q_{inv}^{ss})","L");
legKL_sig->AddEntry(KL_sb,"f_{sb}D(q_{inv}^{sb})","L");
legKL_sig->AddEntry(KL_bb,"f_{bb}D(q_{inv}^{bb})","L");
legKL_sig->Draw();

std::ostringstream os_KLsidemet; 
TLatex* texCMS_yKLsidemet = new TLatex;
texCMS_yKLsidemet->SetNDC();
texCMS_yKLsidemet->SetTextFont(42);
texCMS_yKLsidemet->SetTextSize(0.05);
texCMS_yKLsidemet->SetTextAlign(12);
texCMS_yKLsidemet->DrawLatex(.12,.95," #font[61]{CMS} #it{Work in Progress}");
texCMS_yKLsidemet->DrawLatex(0.60,0.95, "#scale[0.9]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS_yKLsidemet->DrawLatex(0.15,0.83, "#scale[1.2]{K^{0}_{S}#Lambda #oplus K^{0}_{S}#bar{#Lambda}}");  
//texCMS_yKLsidemet->DrawLatex(0.15,0.83, "#scale[1.0]{Signal}");  
os_KLsidemet << "#scale[0.9]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_yKLsidemet->DrawLatex(0.25,0.30, os_KLsidemet.str(  ).c_str(  )); 
os_KLsidemet.str(std::string()); 
os_KLsidemet << "#scale[0.9]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_yKLsidemet->DrawLatex(0.25,0.24, os_KLsidemet.str(  ).c_str(  ));  
os_KLsidemet.str(std::string()); 
gPad->SetLogy();
c_KL_sidemet->Update();
c_KL_sidemet->Print("results/sidemet_log_sig_KL_"+V0+"_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");

TH1D* h_sig_K0sLAL_befpur_ss = (TH1D*)h_K_plus_LAL_qinv_ss_1D->Clone();
TH1D* h_sig_K0sLAL_befpur_sb = (TH1D*)h_K_plus_LAL_qinv_sb_1D->Clone();
TH1D* h_sig_K0sLAL_befpur_bb = (TH1D*)h_K_plus_LAL_qinv_bb_1D->Clone();
TH1D* h_mix_K0sLAL_befpur_ss = (TH1D*)h_K_plus_LAL_qinv_ss_1D_mix->Clone();
TH1D* h_mix_K0sLAL_befpur_sb = (TH1D*)h_K_plus_LAL_qinv_sb_1D_mix->Clone();
TH1D* h_mix_K0sLAL_befpur_bb = (TH1D*)h_K_plus_LAL_qinv_bb_1D_mix->Clone();

sideband(h_K_plus_LAL_qinv_ss_1D, h_K_plus_LAL_qinv_sb_1D, h_K_plus_LAL_qinv_bb_1D, fss_KLAL, fsb_KLAL, fbb_KLAL);
sideband(h_K_plus_LAL_qinv_ss_1D_mix, h_K_plus_LAL_qinv_sb_1D_mix, h_K_plus_LAL_qinv_bb_1D_mix, fss_KLAL, fsb_KLAL, fbb_KLAL);

TH1D* h_sig_K0sLAL = (TH1D*)h_K_plus_LAL_qinv_ss_1D->Clone();
TH1D* h_mix_K0sLAL = (TH1D*)h_K_plus_LAL_qinv_ss_1D_mix->Clone();

normalize(qmin_Norm,qmax_Norm,h_K_plus_LAL_qinv_ss_1D,h_K_plus_LAL_qinv_ss_1D_mix);
TGraphAsymmErrors h_K_plus_LAL_qinv_ss_graph;
h_K_plus_LAL_qinv_ss_graph.Divide(h_K_plus_LAL_qinv_ss_1D, h_K_plus_LAL_qinv_ss_1D_mix, "pois");
h_K_plus_LAL_qinv_ss_1D->Divide(h_K_plus_LAL_qinv_ss_1D_mix);

//KLam

sideband(h_KL_qinv_ss_1D, h_KL_qinv_sb_1D, h_KL_qinv_bb_1D, fss_KL, fsb_KL, fbb_KL);
sideband(h_KL_qinv_mix_ss_1D, h_KL_qinv_mix_sb_1D, h_KL_qinv_mix_bb_1D, fss_KL, fsb_KL, fbb_KL);
normalize(qmin_Norm,qmax_Norm,h_KL_qinv_ss_1D,h_KL_qinv_mix_ss_1D);
h_KL_qinv_ss_1D->Divide(h_KL_qinv_mix_ss_1D);

//KALam

sideband(h_KAL_qinv_ss_1D, h_KAL_qinv_sb_1D, h_KAL_qinv_bb_1D, fss_KAL, fsb_KAL, fbb_KAL);
sideband(h_KAL_qinv_mix_ss_1D, h_KAL_qinv_mix_sb_1D, h_KAL_qinv_mix_bb_1D, fss_KAL, fsb_KAL, fbb_KAL);
normalize(qmin_Norm,qmax_Norm,h_KAL_qinv_ss_1D,h_KAL_qinv_mix_ss_1D);
h_KAL_qinv_ss_1D->Divide(h_KAL_qinv_mix_ss_1D);

//=================================================================
//K0s Lam cross checks
//=================================================================
TCanvas *c_y = new TCanvas("c_y","",600,600);
TPad* pad1_y = new TPad("pad1_y","top pad",0.0,0.3,1.0,1.0);
TPad* pad2_y = new TPad("pad2_y","bottom pad",0.0,0.0,1.0,0.3);
pad1_y->SetTopMargin(0.1);
pad1_y->SetBottomMargin(0.0);
pad1_y->SetRightMargin(0.038);
pad1_y->SetLeftMargin(0.15);
pad1_y->Draw();
pad2_y->SetFrameFillStyle(4000);
pad2_y->SetTopMargin(0.0);
pad2_y->SetBottomMargin(0.37);
pad2_y->SetRightMargin(0.038);
pad2_y->SetLeftMargin(0.15);
pad2_y->Draw();
pad1_y->cd();
gPad->SetTickx();
gPad->SetTicky();
std::ostringstream oKL; 
oKL << 1000*h_KL_qinv_ss_1D->GetBinWidth(1) << std::setprecision(2);
h_KL_qinv_ss_1D->Draw("e");
h_KL_qinv_ss_1D->GetXaxis()->SetRangeUser(0.0,2.9999999);
h_KL_qinv_ss_1D->SetMarkerStyle(24);
h_KL_qinv_ss_1D->SetMarkerSize(1.5);
h_KL_qinv_ss_1D->SetMarkerColor(SPcolors[SPdred]);
h_KL_qinv_ss_1D->SetLineColor(SPcolors[SPdred]);
h_KL_qinv_ss_1D->GetYaxis()->SetTitle(Form("#scale[1.7]{Single Ratio / %s MeV}",oLAL.str(  ).c_str(  )));
h_KL_qinv_ss_1D->GetYaxis()->CenterTitle(1);
h_KL_qinv_ss_1D->GetYaxis()->SetNdivisions(810);
h_KL_qinv_ss_1D->GetXaxis()->SetTickSize(0.03);
h_KL_qinv_ss_1D->GetYaxis()->SetTickSize(0.02);
h_KL_qinv_ss_1D->GetXaxis()->SetLabelSize(0.055);
h_KL_qinv_ss_1D->GetYaxis()->SetLabelSize(0.055);
h_KL_qinv_ss_1D->GetYaxis()->SetLabelOffset(0.01);
h_KL_qinv_ss_1D->GetYaxis()->SetTitleOffset(1.8);
h_KAL_qinv_ss_1D->Draw("same");
h_KAL_qinv_ss_1D->SetMarkerStyle(25);
h_KAL_qinv_ss_1D->SetMarkerColor(1);
h_KAL_qinv_ss_1D->SetLineColor(1);
h_KAL_qinv_ss_1D->SetMarkerSize(1.5);
TLegend* leg00 = new TLegend(0.70,0.32,0.90,0.45);
leg00->SetBorderSize(0);
leg00->AddEntry(h_KL_qinv_ss_1D,"K^{0}_{S}#Lambda");
leg00->AddEntry(h_KAL_qinv_ss_1D,"K^{0}_{S}#bar{#Lambda}");
leg00->Draw();
std::ostringstream os_y; 
TLatex* texCMS_y = new TLatex;
texCMS_y->SetNDC();
texCMS_y->SetTextFont(42);
texCMS_y->SetTextSize(0.05);
texCMS_y->SetTextAlign(12);
texCMS_y->DrawLatex(.20,.95," #font[61]{CMS} #it{Work in Progress}");
texCMS_y->DrawLatex(0.60,0.95, "#scale[1.0]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
os_y << "#scale[1.0]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_y->DrawLatex(0.68,0.60, os_y.str(  ).c_str(  )); 
os_y.str(std::string()); 
os_y << "#scale[1.0]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_y->DrawLatex(0.68,0.52, os_y.str(  ).c_str(  ));  
os_y.str(std::string()); 
TLine* line_yx = new TLine(0.0, 1, 2.99999, 1);
line_yx->SetLineStyle(2);
line_yx->Draw("same");
c_y->Update();
pad2_y->cd();
gPad->SetTickx();
gPad->SetTicky();
TH1D* ratio_K = (TH1D*) h_KL_qinv_ss_1D->Clone();
ratio_K->Divide(h_KL_qinv_ss_1D,h_KAL_qinv_ss_1D,1,1);
double low_y = 0.49999;
double high_y = 1.499999999;
TH1F* frame_ks_y = pad2_y->cd()->DrawFrame(0.00000,low_y,2.99999,high_y);
frame_ks_y->GetXaxis()->SetTitle("#scale[4.3]{q_{inv} [GeV]}");
frame_ks_y->GetYaxis()->SetTitle("#scale[4.0]{K^{0}_{S}#Lambda/K^{0}_{S}#bar{#Lambda}}");
frame_ks_y->GetXaxis()->CenterTitle(1);
frame_ks_y->GetYaxis()->CenterTitle(1);
frame_ks_y->GetXaxis()->SetNdivisions(510);
frame_ks_y->GetYaxis()->SetNdivisions(406);
frame_ks_y->GetXaxis()->SetTickSize(0.06);
frame_ks_y->GetYaxis()->SetTickSize(0.02);
frame_ks_y->GetXaxis()->SetTitleOffset(4.5);
frame_ks_y->GetYaxis()->SetTitleOffset(1.8);
frame_ks_y->GetXaxis()->SetLabelSize(0.14);
frame_ks_y->GetYaxis()->SetLabelSize(0.14);
frame_ks_y->GetXaxis()->SetLabelOffset(0.01);
frame_ks_y->GetYaxis()->SetLabelOffset(0.013);
ratio_K->Draw("esame");
ratio_K->SetMarkerStyle(20);
ratio_K->SetMarkerSize(1.5);
ratio_K->SetMarkerColor(kBlack);
ratio_K->SetLineColor(kBlack);
TLine* line_y = new TLine(0.0, 1, 2.99999, 1);
line_y->SetLineStyle(2);
line_y->Draw("same");
c_y->Update();
c_y->Print("results/xcheck_KLKAL_"+V0+"_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");
//=================================================================


//=================================================================
//plot single ratios
//=================================================================

//K0sLam+K0sALam

TCanvas *c_KL_SR = new TCanvas("c_KL_SR","",800,600);
gPad->SetTickx();
gPad->SetTicky();
std::ostringstream oKLSR;
oKLSR << 1000*h_K_plus_LAL_qinv_ss_1D->GetBinWidth(1) << std::setprecision(2);
h_K_plus_LAL_qinv_ss_1D->Draw("e");
h_K_plus_LAL_qinv_ss_1D->GetXaxis()->SetRangeUser(0.0,2.9999999);
h_K_plus_LAL_qinv_ss_1D->GetYaxis()->SetRangeUser(0.901,1.14999999);
h_K_plus_LAL_qinv_ss_1D->SetMarkerStyle(26);
h_K_plus_LAL_qinv_ss_1D->SetMarkerSize(1.6);
h_K_plus_LAL_qinv_ss_1D->SetMarkerColor(1);
h_K_plus_LAL_qinv_ss_1D->SetMarkerColorAlpha(1,1.2);
h_K_plus_LAL_qinv_ss_1D->SetLineColor(1);
h_K_plus_LAL_qinv_ss_1D->GetYaxis()->SetTitle(Form("#scale[1.3]{Single Ratio / %s MeV}",oLAL.str(  ).c_str(  )));
h_K_plus_LAL_qinv_ss_1D->GetYaxis()->CenterTitle(1);
h_K_plus_LAL_qinv_ss_1D->GetYaxis()->SetNdivisions(810);
h_K_plus_LAL_qinv_ss_1D->GetXaxis()->SetTickSize(0.03);
h_K_plus_LAL_qinv_ss_1D->GetYaxis()->SetTickSize(0.03);
h_K_plus_LAL_qinv_ss_1D->GetXaxis()->SetLabelSize(0.042);
h_K_plus_LAL_qinv_ss_1D->GetYaxis()->SetLabelSize(0.042);
h_K_plus_LAL_qinv_ss_1D->GetYaxis()->SetLabelOffset(0.01);
h_K_plus_LAL_qinv_ss_1D->GetYaxis()->SetTitleOffset(1.3);
h_K_plus_LAL_qinv_ss_1D->GetXaxis()->SetTitleOffset(1.3);
h_K_plus_LAL_qinv_ss_1D->GetXaxis()->SetTitle("#scale[1.3]{q_{inv} [GeV]}");
h_K_plus_LAL_qinv_ss_1D->GetXaxis()->CenterTitle(1);
std::ostringstream os_KLSR; 
TLatex* texCMS_yKLSR = new TLatex;
texCMS_yKLSR->SetNDC();
texCMS_yKLSR->SetTextFont(42);
texCMS_yKLSR->SetTextSize(0.05);
texCMS_yKLSR->SetTextAlign(12);
texCMS_yKLSR->DrawLatex(0.13,0.93," #font[61]{CMS} #it{Work in Progress}");
texCMS_yKLSR->DrawLatex(0.57,0.93, "#scale[0.9]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS_yKLSR->DrawLatex(0.50,0.75, "#scale[1.2]{K^{0}_{S}#Lambda#oplusK^{0}_{S}#bar{#Lambda}}");  
os_KLSR << "#scale[1.0]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_yKLSR->DrawLatex(0.48,0.60, os_KLSR.str(  ).c_str(  )); 
os_KLSR.str(std::string()); 
os_KLSR << "#scale[1.0]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_yKLSR->DrawLatex(0.48,0.54, os_KLSR.str(  ).c_str(  ));  
os_KLSR.str(std::string()); 
TLine* line_yxKLSR = new TLine(0.0, 1, 2.99999, 1);
line_yxKLSR->SetLineStyle(2);
line_yxKLSR->Draw("same");
texCMS_yKLSR->DrawLatex(0.22,0.78, "#scale[0.9]{#Xi(1820)}");  
TArrow *ar_0 = new TArrow(0.57,1.088,0.57,1.105,0.01,"<");
ar_0->SetAngle(40);
ar_0->SetLineWidth(2);
ar_0->Draw();
c_KL_SR->Update();
c_KL_SR->Print("results/SR_KLKAL_"+V0+"_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");

//LamLam+ALamALam

TCanvas *c_LL_SR = new TCanvas("c_LL_SR","",800,600);
gPad->SetTickx();
gPad->SetTicky();
std::ostringstream oLLSR; 
oLLSR << 1000*h_Lam_plus_ALam_qinv_ss_1D->GetBinWidth(1) << std::setprecision(2);
h_Lam_plus_ALam_qinv_ss_1D->Draw("e");
h_Lam_plus_ALam_qinv_ss_1D->GetXaxis()->SetRangeUser(0.0,2.9999999);
h_Lam_plus_ALam_qinv_ss_1D->GetYaxis()->SetRangeUser(0.0,1.14999999);
h_Lam_plus_ALam_qinv_ss_1D->SetMarkerStyle(25);
h_Lam_plus_ALam_qinv_ss_1D->SetMarkerSize(1.6);
h_Lam_plus_ALam_qinv_ss_1D->SetMarkerColor(1);
h_Lam_plus_ALam_qinv_ss_1D->SetLineColor(1);
h_Lam_plus_ALam_qinv_ss_1D->SetMarkerColorAlpha(1,1.0);
h_Lam_plus_ALam_qinv_ss_1D->GetYaxis()->SetTitle(Form("#scale[1.3]{Single Ratio / %s MeV}",oLAL.str(  ).c_str(  )));
h_Lam_plus_ALam_qinv_ss_1D->GetYaxis()->CenterTitle(1);
h_Lam_plus_ALam_qinv_ss_1D->GetYaxis()->SetNdivisions(810);
h_Lam_plus_ALam_qinv_ss_1D->GetXaxis()->SetTickSize(0.03);
h_Lam_plus_ALam_qinv_ss_1D->GetYaxis()->SetTickSize(0.03);
h_Lam_plus_ALam_qinv_ss_1D->GetXaxis()->SetLabelSize(0.042);
h_Lam_plus_ALam_qinv_ss_1D->GetYaxis()->SetLabelSize(0.042);
h_Lam_plus_ALam_qinv_ss_1D->GetYaxis()->SetLabelOffset(0.01);
h_Lam_plus_ALam_qinv_ss_1D->GetYaxis()->SetTitleOffset(1.3);
h_Lam_plus_ALam_qinv_ss_1D->GetXaxis()->SetTitleOffset(1.3);
h_Lam_plus_ALam_qinv_ss_1D->GetXaxis()->SetTitle("#scale[1.3]{q_{inv} [GeV]}");
h_Lam_plus_ALam_qinv_ss_1D->GetXaxis()->CenterTitle(1);
std::ostringstream os_LLSR; 
TLatex* texCMS_yLLSR = new TLatex;
texCMS_yLLSR->SetNDC();
texCMS_yLLSR->SetTextFont(42);
texCMS_yLLSR->SetTextSize(0.05);
texCMS_yLLSR->SetTextAlign(12);
texCMS_yLLSR->DrawLatex(0.13,0.93," #font[61]{CMS} #it{Work in Progress}");
texCMS_yLLSR->DrawLatex(0.57,0.93, "#scale[0.9]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS_yLLSR->DrawLatex(0.50,0.65, "#scale[1.2]{#Lambda#Lambda#oplus#bar{#Lambda}#bar{#Lambda}}");  
os_LLSR << "#scale[1.0]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_yLLSR->DrawLatex(0.48,0.50, os_LLSR.str(  ).c_str(  )); 
os_LLSR.str(std::string()); 
os_LLSR << "#scale[1.0]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_yLLSR->DrawLatex(0.48,0.44, os_LLSR.str(  ).c_str(  ));  
os_LLSR.str(std::string()); 
TLine* line_yxLLSR = new TLine(0.0, 1, 2.99999, 1);
line_yxLLSR->SetLineStyle(2);
line_yxLLSR->Draw("same");
c_LL_SR->Update();
c_LL_SR->Print("results/SR_LL_"+V0+"_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");

//LamALam

TCanvas *c_LAL_SR = new TCanvas("c_LAL_SR","",800,600);
gPad->SetTickx();
gPad->SetTicky();
std::ostringstream oLALSR; 
oLALSR << 1000*h_LAL_qinv_ss_1D->GetBinWidth(1) << std::setprecision(2);
h_LAL_qinv_ss_1D->Draw("e");
h_LAL_qinv_ss_1D->GetXaxis()->SetRangeUser(0.0,2.9999999);
h_LAL_qinv_ss_1D->GetYaxis()->SetRangeUser(0.6,1.4999999);
h_LAL_qinv_ss_1D->SetMarkerStyle(24);
h_LAL_qinv_ss_1D->SetMarkerSize(1.6);
h_LAL_qinv_ss_1D->SetMarkerColor(1);
h_LAL_qinv_ss_1D->SetLineColor(1);
h_LAL_qinv_ss_1D->SetMarkerColorAlpha(1,1.2);
h_LAL_qinv_ss_1D->GetYaxis()->SetTitle(Form("#scale[1.3]{Single Ratio / %s MeV}",oLAL.str(  ).c_str(  )));
h_LAL_qinv_ss_1D->GetYaxis()->CenterTitle(1);
h_LAL_qinv_ss_1D->GetYaxis()->SetNdivisions(810);
h_LAL_qinv_ss_1D->GetXaxis()->SetTickSize(0.03);
h_LAL_qinv_ss_1D->GetYaxis()->SetTickSize(0.03);
h_LAL_qinv_ss_1D->GetXaxis()->SetLabelSize(0.042);
h_LAL_qinv_ss_1D->GetYaxis()->SetLabelSize(0.042);
h_LAL_qinv_ss_1D->GetYaxis()->SetLabelOffset(0.01);
h_LAL_qinv_ss_1D->GetYaxis()->SetTitleOffset(1.3);
h_LAL_qinv_ss_1D->GetXaxis()->SetTitleOffset(1.3);
h_LAL_qinv_ss_1D->GetXaxis()->SetTitle("#scale[1.3]{q_{inv} [GeV]}");
h_LAL_qinv_ss_1D->GetXaxis()->CenterTitle(1);
std::ostringstream os_LALSR; 
TLatex* texCMS_yLALSR = new TLatex;
texCMS_yLALSR->SetNDC();
texCMS_yLALSR->SetTextFont(42);
texCMS_yLALSR->SetTextSize(0.05);
texCMS_yLALSR->SetTextAlign(12);
texCMS_yLALSR->DrawLatex(0.13,0.93," #font[61]{CMS} #it{Work in Progress}");
texCMS_yLALSR->DrawLatex(0.57,0.93, "#scale[0.9]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS_yLALSR->DrawLatex(0.60,0.8, "#scale[1.2]{#Lambda#bar{#Lambda}}");  
os_LALSR << "#scale[1.0]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_yLALSR->DrawLatex(0.48,0.70, os_LALSR.str(  ).c_str(  )); 
os_LALSR.str(std::string()); 
os_LALSR << "#scale[1.0]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_yLALSR->DrawLatex(0.48,0.64, os_LALSR.str(  ).c_str(  ));  
os_LALSR.str(std::string()); 
TLine* line_yxLALSR = new TLine(0.0, 1, 2.99999, 1);
line_yxLALSR->SetLineStyle(2);
line_yxLALSR->Draw("same");
c_LAL_SR->Update();
c_LAL_SR->Print("results/SR_LAL_"+V0+"_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");

//K0sK0s

TCanvas *c_KsKs_SR = new TCanvas("c_KsKs_SR","",800,600);
gPad->SetTickx();
gPad->SetTicky();
std::ostringstream oKsKsSR; 
oKsKsSR << 1000*h_K0s_qinv_ss_1D->GetBinWidth(1) << std::setprecision(2);
h_K0s_qinv_ss_1D->Draw("e");
h_K0s_qinv_ss_1D->GetXaxis()->SetRangeUser(0.0,2.9999999);
h_K0s_qinv_ss_1D->GetYaxis()->SetRangeUser(0.9,2.6999999);
h_K0s_qinv_ss_1D->SetMarkerStyle(28);
h_K0s_qinv_ss_1D->SetMarkerSize(1.6);
h_K0s_qinv_ss_1D->SetMarkerColor(1);
h_K0s_qinv_ss_1D->SetLineColor(1);
h_K0s_qinv_ss_1D->SetMarkerColorAlpha(1,1.0);
h_K0s_qinv_ss_1D->GetYaxis()->SetTitle(Form("#scale[1.3]{Single Ratio / %s MeV}",oLAL.str(  ).c_str(  )));
h_K0s_qinv_ss_1D->GetYaxis()->CenterTitle(1);
h_K0s_qinv_ss_1D->GetYaxis()->SetNdivisions(810);
h_K0s_qinv_ss_1D->GetXaxis()->SetTickSize(0.03);
h_K0s_qinv_ss_1D->GetYaxis()->SetTickSize(0.03);
h_K0s_qinv_ss_1D->GetXaxis()->SetLabelSize(0.042);
h_K0s_qinv_ss_1D->GetYaxis()->SetLabelSize(0.042);
h_K0s_qinv_ss_1D->GetYaxis()->SetLabelOffset(0.01);
h_K0s_qinv_ss_1D->GetYaxis()->SetTitleOffset(1.3);
h_K0s_qinv_ss_1D->GetXaxis()->SetTitleOffset(1.3);
h_K0s_qinv_ss_1D->GetXaxis()->SetTitle("#scale[1.3]{q_{inv} [GeV]}");
h_K0s_qinv_ss_1D->GetXaxis()->CenterTitle(1);
std::ostringstream os_KsKsSR; 
TLatex* texCMS_yKsKsSR = new TLatex;
texCMS_yKsKsSR->SetNDC();
texCMS_yKsKsSR->SetTextFont(42);
texCMS_yKsKsSR->SetTextSize(0.05);
texCMS_yKsKsSR->SetTextAlign(12);
texCMS_yKsKsSR->DrawLatex(0.13,0.93," #font[61]{CMS} #it{Work in Progress}");
texCMS_yKsKsSR->DrawLatex(0.57,0.93, "#scale[0.9]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");  
texCMS_yKsKsSR->DrawLatex(0.50,0.75, "#scale[1.2]{K^{0}_{S}K^{0}_{S}}");  
os_KsKsSR << "#scale[1.0]{ " << std::setprecision(4) << ntkmin << std::setprecision(4) <<" < N^{offline}_{trk} < "  << ntkmax << " }" << std::setprecision(2); 
texCMS_yKsKsSR->DrawLatex(0.48,0.60, os_KsKsSR.str(  ).c_str(  )); 
os_KsKsSR.str(std::string()); 
os_KsKsSR << "#scale[1.0]{ " << std::setprecision(4) << ktmin << std::setprecision(4) <<" < k_{T} < "  << ktmax << " GeV}" << std::setprecision(2); 
texCMS_yKsKsSR->DrawLatex(0.48,0.54, os_KsKsSR.str(  ).c_str(  ));  
os_KsKsSR.str(std::string()); 
TLine* line_yxKsKsSR = new TLine(0.0, 1, 2.99999, 1);
line_yxKsKsSR->SetLineStyle(2);
line_yxKsKsSR->Draw("same");
texCMS_yKLSR->DrawLatex(0.37,0.28, "#scale[0.9]{f_{2}'(1525)}");  
TArrow *ar_1 = new TArrow(1.15,1.1,1.15,1.2,0.01,"<");
ar_1->SetAngle(40);
ar_1->SetLineWidth(2);
ar_1->Draw();
c_KsKs_SR->Update();
c_KsKs_SR->Print("results/SR_KsKs_"+V0+"_ktmin_"+ktmin+"_ktmax_"+ktmax+"_ntrkmin"+ntkmin+"_ntkmax_"+ntkmax+".pdf");

TFile save("results/histos"+V0+"ktbin_"+ktbin+"_multbin_"+multbin+"_syst_"+syst+".root","RECREATE");
//without purity correction

K0s_hist_data->Write("K0sK0s_CF_nopurity");
LL_hist_data->Write("LLALAL_CF_nopurity");
LAL_hist_data->Write("LAL_CF_nopurity");
K0sLam_hist_data->Write("K0sLamK0sALam_CF_nopurity");

//with purity correction

h_K0s_qinv_ss_1D->Write("K0sK0s_CF_purity");
h_Lam_plus_ALam_qinv_ss_1D->Write("LLALAL_CF_purity");
h_LAL_qinv_ss_1D->Write("LAL_CF_purity");
h_K_plus_LAL_qinv_ss_1D->Write("K0sLamK0sALam_CF_purity");

h_K0s_qinv_ss_graph.Write("K0sK0s_CF_purityG");
h_Lam_qinv_ss_graph.Write("LLALAL_CF_purityG");
h_LAL_qinv_ss_graph.Write("LAL_CF_purityG");
h_K_plus_LAL_qinv_ss_graph.Write("K0sLamK0sALam_CF_purityG");


//qinvhisto after purity correction

h_sig_K0sK0s_befpur_ss->Write("K0sK0s_sig_ss");
h_sig_K0sK0s_befpur_sb->Write("K0sK0s_sig_sb");
h_sig_K0sK0s_befpur_bb->Write("K0sK0s_sig_bb");
h_mix_K0sK0s_befpur_ss->Write("K0sK0s_mix_ss");
h_mix_K0sK0s_befpur_sb->Write("K0sK0s_mix_sb");
h_mix_K0sK0s_befpur_bb->Write("K0sK0s_mix_bb");

h_sig_K0sK0s->Write("K0sK0s_sig_purity");
h_mix_K0sK0s->Write("K0sK0s_mix_purity");

h_sig_LLALAL_befpur_ss->Write("LLALAL_sig_ss");
h_sig_LLALAL_befpur_sb->Write("LLALAL_sig_sb");
h_sig_LLALAL_befpur_bb->Write("LLALAL_sig_bb");
h_mix_LLALAL_befpur_ss->Write("LLALAL_mix_ss");
h_mix_LLALAL_befpur_sb->Write("LLALAL_mix_sb");
h_mix_LLALAL_befpur_bb->Write("LLALAL_mix_bb");

h_sig_LLALAL->Write("LLALAL_sig_purity");
h_mix_LLALAL->Write("LLALAL_mix_purity");

h_sig_K0sLAL_befpur_ss->Write("K0sLAL_sig_ss");
h_sig_K0sLAL_befpur_sb->Write("K0sLAL_sig_sb");
h_sig_K0sLAL_befpur_bb->Write("K0sLAL_sig_bb");
h_mix_K0sLAL_befpur_ss->Write("K0sLAL_mix_ss");
h_mix_K0sLAL_befpur_sb->Write("K0sLAL_mix_sb");
h_mix_K0sLAL_befpur_bb->Write("K0sLAL_mix_bb");

h_sig_K0sLAL->Write("K0sLAL_sig_purity");
h_mix_K0sLAL->Write("K0sLAL_mix_purity");

h_sig_LAL_befpur_ss->Write("LAL_sig_ss");
h_sig_LAL_befpur_sb->Write("LAL_sig_sb");
h_sig_LAL_befpur_bb->Write("LAL_sig_bb");
h_mix_LAL_befpur_ss->Write("LAL_mix_ss");
h_mix_LAL_befpur_sb->Write("LAL_mix_sb");
h_mix_LAL_befpur_bb->Write("LAL_mix_bb");

h_sig_LAL->Write("LAL_sig_purity");
h_mix_LAL->Write("LAL_mix_purity");

save.Close();


} //end macro
