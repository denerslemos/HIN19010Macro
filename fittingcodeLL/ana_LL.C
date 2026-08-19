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
#include "TH1F.h"
#include "TBenchmark.h"
#include "TSystem.h"
//#include "/home/ddesouza/utils/setSPRACEColors.C"
#include "./func.h"
#include "./plots.h"

#include "TStyle.h"

using namespace std;


void ana_LL(){
    
// For systematics, have to update    
const int nsys = 1;

string systematic[nsys] = {
"Nominal"
};

 
//root stuff
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
    
int Ntrkrange = 0;// -> 0 means integrated; -> 1 means MB [0,120]; -> 2 means HM [185,250]; 
    
//qinv range
double qmin = 0.0;
double qmax = 5.999999;

//qinv normalization range
double qmin_Norm = 0.0;
double qmax_Norm = 5.999999;


//kT bin
double ktmin = 0.0;
double ktmax = 1.9999999;

//ntrk bin
double ntkmin;
double ntkmax;
double Rpar_LAL;

double a = 0.263588;
double b = 0.606614;

if(Ntrkrange==0){
ntkmin = 0.0;
ntkmax = 399.9999999;
Rpar_LAL = 2.36619;
//Rpar_LAL = a*pow(100.,1./3.) + b;
}else if(Ntrkrange==1){
ntkmin = 0.0;
ntkmax = 119.9999999;
Rpar_LAL = 2.12899;
//Rpar_LAL = a*pow(70.,1./3.) + b;
}else if(Ntrkrange==2){
ntkmin = 185.0;
ntkmax = 249.9999999;    
Rpar_LAL = 3.24457;
//Rpar_LAL = a*pow(200.,1./3.) + b;    
}

vector<double> Rinv;
vector<double> Rinv_unc;
vector<double> Lamb;
vector<double> Lamb_unc;
vector<double> Ref0;
vector<double> Ref0_unc;
vector<double> Imf0;
vector<double> Imf0_unc;
vector<double> d0ef;
vector<double> d0ef_unc;


for(Int_t i=0; i<nsys; i++){
    
if(i>0 && Ntrkrange>0) break;
    
// if(i!=0) break;

TString syste;

if(i>=42){syste = "Standard";}else{syste=systematic[i];}

//Reading the input file
TString filename; 
filename = Form("histos_V0_ktbin_%i_multbin_%i_syst_%i.root",0,Ntrkrange,i);
TFile* f = new TFile(Form("%s",filename.Data()),"READ");
cout << "File Read" << endl;


if(i==39) qmax_Norm = 4.499999;
if(i==40) qmax_Norm = 2.999999;


//read data and MC
TH1D* h_LAL_qinv_ss_1D=(TH1D*)f->Get("LLALAL_CF_purity");
TH1D* h_LAL_qinv_GEN_1D=(TH1D*)f->Get("hist_CF_LamLam_EPOS");


ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
ROOT::Math::MinimizerOptions::SetDefaultErrorDef(0.95);

cout << "----------------------------------------- Start chi2 fit -----------------------------------------" << endl;
cout << endl;

double N_MC_LamALam;
double B_MC_LamALam;
double invSig_MC_LamALam;
double delta_MC_LamALam;


N_MC_LamALam = 1.0;
B_MC_LamALam = -0.3;
invSig_MC_LamALam = 0.2;
delta_MC_LamALam = 0.0;

//making fits
int iter = 5;
std::pair<double, double>  full_qrange(0.04, 4.999999);
std::pair<double, double>  cluster(0.5, 2.99999);

if(i==33) cluster.first = 0.6;
if(i==34) cluster.first = 0.4;
if(i==35) cluster.second = 3.999999; full_qrange.second = 3.999999;
if(i==36) cluster.second = 4.999999; full_qrange.second = 4.999999;

// if(i==39 || i==40)  h_LAL_qinv_ss_1D->GetXaxis()->SetRangeUser(full_qrange.first,full_qrange.second);


//Start with LamALam
TF1 *f_exp = new TF1("f_exp",fit_strong_LL,full_qrange.first, full_qrange.second, 11);
f_exp->SetParameters(1.0,3.0,0.8,0.0,3.0,-0.2,0.2,0.0,0.5,0,1.0);
f_exp->SetParName(0,"N");
f_exp->SetParName(1,"R_{Exp}");
f_exp->SetParName(2,"Ref0");
f_exp->SetParName(3,"Imf0");
f_exp->SetParName(4,"d0");
f_exp->SetParName(5,"B");
f_exp->SetParName(6,"#sigma");
f_exp->SetParName(7,"#delta");
f_exp->SetParName(8,"#lambda");
f_exp->SetParName(9,"isID");
f_exp->SetParName(10,"NClust");
f_exp->SetLineColor(kRed); 
f_exp->SetLineWidth(3);
f_exp->FixParameter(9,-0.5);//0 for nonidentical and -1 for identical fermions and +1 for bosons

//cluster fit parameters -> start with MC
f_exp->SetParameter(10, N_MC_LamALam);
f_exp->SetParameter(5, B_MC_LamALam);
f_exp->SetParameter(6, invSig_MC_LamALam);
f_exp->SetParameter(7, delta_MC_LamALam);

//fix the parameters in order to get only Non femto bkg
f_exp->FixParameter(0, 1.);
f_exp->FixParameter(1, 1.);
f_exp->FixParameter(2, 1.);
f_exp->FixParameter(3, 0.);
f_exp->FixParameter(4, 0.);
f_exp->FixParameter(8, 0.);
f_exp->FixParameter(9, 0.);

h_LAL_qinv_ss_1D->Fit(f_exp, "Rq0", "",cluster.first,cluster.second);
h_LAL_qinv_ss_1D->Fit(f_exp, "Rq0", "",cluster.first,cluster.second);
h_LAL_qinv_ss_1D->Fit(f_exp, "Rq0", "",cluster.first,cluster.second);

double err;
if(i==37){err = 0.05;}else if(i==38){err = 0.15;}else{err = 0.1;}

f_exp->FixParameter(10, f_exp->GetParameter(10));
f_exp->FixParameter(5, f_exp->GetParameter(5));
f_exp->FixParameter(6, f_exp->GetParameter(6));
f_exp->FixParameter(7, f_exp->GetParameter(7));


TF1 *f_expQS = new TF1("f_expQS",fit_QS_LL,full_qrange.first, full_qrange.second, 7);
f_expQS->SetParameters(1.0,5.0,0.5,-0.3,0.2,0.0,1.0);
f_expQS->SetParName(0,"N");
f_expQS->SetParName(1,"R_{Exp}");
f_expQS->SetParName(2,"#lambda");
f_expQS->SetParName(3,"B");
f_expQS->SetParName(4,"#sigma");
f_expQS->SetParName(5,"#delta");
f_expQS->SetParName(6,"NClust");
f_expQS->SetLineColor(kRed); 
f_expQS->SetLineWidth(2);

f_expQS->FixParameter(6, f_exp->GetParameter(10));
f_expQS->FixParameter(3, f_exp->GetParameter(5));
f_expQS->FixParameter(4, f_exp->GetParameter(6));
f_expQS->FixParameter(5, f_exp->GetParameter(7));

f_expQS->ReleaseParameter(1);

h_LAL_qinv_ss_1D->Fit(f_expQS, "Rq0");


f_exp->FixParameter(9,-0.5);
f_exp->FixParameter(0,1.0);
f_exp->FixParameter(3,0.0);

f_exp->ReleaseParameter(1);
f_exp->SetParameter(1,fabs(f_expQS->GetParameter(1)));
f_exp->SetParLimits(1,fabs(f_expQS->GetParameter(1)),5.0);

f_exp->ReleaseParameter(2);
f_exp->SetParameter(2,1.);
f_exp->SetParLimits(2, 0.2,5.0);

f_exp->ReleaseParameter(4);
f_exp->SetParLimits(4, 0,20);

f_exp->ReleaseParameter(8);
f_exp->SetParameter(8,1.5);
f_exp->SetParLimits(8, 2.2,3.0);
 

h_LAL_qinv_ss_1D->Fit(f_exp, "Rq0", "", full_qrange.first, full_qrange.second);
h_LAL_qinv_ss_1D->Fit(f_exp, "Rq0", "", full_qrange.first, full_qrange.second);
h_LAL_qinv_ss_1D->Fit(f_exp, "Rq0", "", full_qrange.first, full_qrange.second);

for ( int i = 0; i < iter; i++ ) {
   h_LAL_qinv_ss_1D->Fit(f_exp, "Rq0", "", full_qrange.first, 1.0);
}

TFitResultPtr res;
res = h_LAL_qinv_ss_1D->Fit(f_exp, "SRq");
Int_t fitStatus = res;
if(fitStatus==0){
cout << "======================" << endl;
cout << "The Fit Converge :)" << endl; 
cout << "======================" << endl;
}else{
cout << "======================" << endl;
cout << "Fit Problems" << endl; 
cout << "======================" << endl;    
}


TF1 * fit_SI   = 0;
TF1 * fit_CL_LR   = 0;
TF1 * fit_QS   = 0;
fit_SI  = new TF1("fit_SI" , fit_strong_LAL, full_qrange.first, full_qrange.second,11);
fit_CL_LR = new TF1("fit_CL_LR" , fit_cluster, full_qrange.first, full_qrange.second,4);
fit_QS  = new TF1("fit_SI" , fit_QS_LL, full_qrange.first, full_qrange.second,7);

fit_SI->SetParameter(0, 1.0);
fit_SI->SetParameter(1, f_exp->GetParameter(1));
fit_SI->SetParameter(2, f_exp->GetParameter(2));
fit_SI->SetParameter(3, f_exp->GetParameter(3));
fit_SI->SetParameter(4, f_exp->GetParameter(4));
fit_SI->SetParameter(5, 0.0);
fit_SI->SetParameter(6, 0.5);
fit_SI->SetParameter(7, 0.0);
fit_SI->SetParameter(8, 1.0);
fit_SI->SetParameter(9, 0.0);
fit_SI->SetParameter(10, 1.0);

fit_CL_LR->SetParameter(0, f_exp->GetParameter(5));
fit_CL_LR->SetParameter(1, f_exp->GetParameter(6));
fit_CL_LR->SetParameter(2, f_exp->GetParameter(7));
fit_CL_LR->SetParameter(3, f_exp->GetParameter(10));


fit_QS->SetParameter(0, 1.0);
fit_QS->SetParameter(1, f_exp->GetParameter(1));
fit_QS->SetParameter(2, 1.0);
fit_QS->SetParameter(3, 0.0);
fit_QS->SetParameter(4, 1.0);
fit_QS->SetParameter(5, 0.0);
fit_QS->SetParameter(6, 1.0);


//LLplot(h_LAL_qinv_ss_1D, f_exp, fit_QS, fit_SI, fit_CL_LR, ktmin, ktmax, ntkmin, ntkmax, Form("LL_chi2_%s",systematic[i].c_str()),"Preliminary", "#chi^{2}",Form("%s",systematic[i].c_str()));

if(i==0 && Ntrkrange==0)LLplot_sys(h_LAL_qinv_ss_1D, f_exp, fit_QS, fit_SI, fit_CL_LR, ktmin, ktmax, ntkmin, ntkmax, Form("LL_chi2_sys_%s",systematic[i].c_str()),"Preliminary", "#chi^{2}");


Rinv.push_back(f_exp->GetParameter(1));
Rinv_unc.push_back(f_exp->GetParError(1));
Lamb.push_back(f_exp->GetParameter(8));
Lamb_unc.push_back(f_exp->GetParError(1));
Ref0.push_back(f_exp->GetParameter(2));
Ref0_unc.push_back(f_exp->GetParError(2));
Imf0.push_back(f_exp->GetParameter(3));
Imf0_unc.push_back(f_exp->GetParError(3));
d0ef.push_back(f_exp->GetParameter(4));
d0ef_unc.push_back(f_exp->GetParError(4));

}


int npoints=Rinv.size();

std::vector<double> Lnew(npoints);
std::vector<double> Lnew_unc(npoints);
std::vector<double> Rinvnew(npoints);
std::vector<double> Rinvnew_unc(npoints);
std::vector<double> Ref0new(npoints);
std::vector<double> Ref0new_unc(npoints);
std::vector<double> Imf0new(npoints);
std::vector<double> Imf0new_unc(npoints);
std::vector<double> d0new(npoints);
std::vector<double> d0new_unc(npoints);

double sumquad_Lam=0;
double sumquad_R=0;
double sumquad_Ref0=0;
double sumquad_Imf0=0;
double sumquad_d0=0;

for(int j=0; j<npoints; j++){

	Rinvnew[j] = Rinv[j];
	Rinvnew_unc[j] = Rinv_unc[j];
	Lnew[j] = Lamb[j];
	Lnew_unc[j] = Lamb_unc[j];
	Ref0new[j] = Ref0[j];
	Ref0new_unc[j] = Ref0_unc[j];
	Imf0new[j] = Imf0[j];
	Imf0new_unc[j] = Imf0_unc[j];
	d0new[j] = d0ef[j];
	d0new_unc[j] = d0ef_unc[j];
	if(j==0 && npoints > 1){
	    cout << "================== Systematics for LAL ==================" << endl;
	    cout << "Source              lambda              Rinv              Ref0              Imf0              d0" << endl;
	}
	if(j!=0 && npoints > 1){
	    double Lsys = 100.*fabs((Lnew[j] - Lnew[0])/Lnew[0]);
	    double Rsys = 100.*fabs((Rinvnew[j] - Rinvnew[0])/Rinvnew[0]);
	    double Rf0sys = 100.*fabs((Ref0new[j] - Ref0new[0])/Ref0new[0]);
	    double If0sys = 100.*fabs((Imf0new[j] - Imf0new[0])/Imf0new[0]);
	    double d0sys = 100.*fabs((d0new[j] - d0new[0])/d0new[0]);
		cout  <<systematic[j] << "              " << Lsys  <<  "              " << Rsys <<  "              " <<  Rf0sys << "              " << If0sys <<  "              " << d0sys << endl;
		sumquad_Lam+=Lsys*Lsys;
		sumquad_R+=Rsys*Rsys;
		sumquad_Ref0+=Rf0sys*Rf0sys;
		sumquad_Imf0+=If0sys*If0sys;
		sumquad_d0+=d0sys*d0sys;
	}
	if(npoints > 1) cout  << "Total:"<< "              " << sqrt(sumquad_Lam)  <<  "              " << sqrt(sumquad_R) <<  "              " <<  sqrt(sumquad_Ref0) << "              " << sqrt(sumquad_Imf0) <<  "              " << sqrt(sumquad_d0) << endl;

    
}

cout << "Rinv: " << Rinv[0] << " +/- " << Rinv_unc[0] << " (stat) +/- " << fabs(Rinv[0]) * sqrt(sumquad_R) / 100.0 << " (syst)" << endl;
cout << "Lambda: " << Lamb[0] << " +/- " << Lamb_unc[0] << " (stat) +/- " << fabs(Lamb[0]) * sqrt(sumquad_Lam) / 100.0 << " (syst)" << endl;
cout << "Ref0: " << Ref0[0] << " +/- " << Ref0_unc[0] << " (stat) +/- " << fabs(Ref0[0]) * sqrt(sumquad_Ref0) / 100.0 << " (syst)" << endl;
cout << "Imf0: " << Imf0[0] << " +/- " << Imf0_unc[0] << " (stat) +/- " << fabs(Imf0[0]) * sqrt(sumquad_Imf0) / 100.0 << " (syst)" << endl;
cout << "d0: " << d0ef[0] << " +/- " << d0ef_unc[0] << " (stat) +/- " << fabs(d0ef[0]) * sqrt(sumquad_d0) / 100.0 << " (syst)" << endl;



}
