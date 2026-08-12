
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
//functions
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "./Faddeeva.cc"
#include "./Faddeeva.hh"

//ctes

#define k0s_mass ( 0.497614 )
#define la_mass ( 1.1156 )
#define pi_mass ( 0.13957 )
#define eta_mass ( 0.547862 )
#define hbarc ( 0.197326968 )
#define PI ( 4.*atan(1.) )

//===========================================
//Lednicky F1 and F2
//===========================================

double GetLednickyF1(double z){
  double result = Faddeeva::Dawson(z)/(z);
  return result;
}

double GetLednickyF2(double z){
  double result = (1-exp(-z*z))/(z);
  return result;
}

//===========================================
//cluster
//===========================================

double fit_cluster(double* x, double* par){
    double CLR = par[3]*(1 + par[0]*exp(-pow(fabs(par[1]*x[0]/hbarc),2)))*(1+par[2]*x[0]);
    double C1D=0;
    C1D = CLR;
    return C1D;
}

double fit_cluster_mc(double* x, double* par){
    double CLR = par[3]*(1 + par[0]*exp(-pow(fabs(par[1]*x[0]/hbarc),2)))*(1+par[2]*x[0]/hbarc);
    double C1D=0;
    C1D = CLR;
    return C1D;
}

double fit_cluster2(double* x, double* par){
    double CLR = par[3]*(1 + par[0]*exp(-pow(fabs(par[1]*x[0]/hbarc),2)))*(1+par[2]*x[0]);
    double C1D=0;
    C1D = CLR;
    return C1D;
}


double fit_cluster_K0s(double* x, double* par){
    double reject_range_min=0.15;
    double reject_range_max=0.30;
    double C1D=0;
    double CLR = par[3]*(1 + par[0]*exp(-pow(fabs(par[1]*x[0]/hbarc),2)))*(1+par[2]*x[0]);
    if( (reject_range_min<x[0] && x[0]<reject_range_max)){
        TF1::RejectPoint();
        C1D = CLR;
    }else{
        C1D = CLR;
    }
    return C1D;
}

double fit_cluster_mc_K0s(double* x, double* par){
    double reject_range_min=50.15;
    double reject_range_max=50.30;
    double C1D=0;
    double CLR = par[3]*(1 + par[0]*exp(-pow(fabs(par[1]*x[0]/hbarc),2)))*(1+par[2]*x[0]/hbarc);
    if( (reject_range_min<x[0] && x[0]<reject_range_max)){
        TF1::RejectPoint();
        C1D = CLR;
    }else{
        C1D = CLR;
    }
    return C1D;
}


double fit_cluster_K0sch(double* x, double* par){
    double reject_range_min=0.3;
    double reject_range_max=0.65;
    double C1D=0;
    double CLR = par[3]*(1 + par[0]*exp(-pow(fabs(par[1]*x[0]/hbarc),2)))*(1+par[2]*x[0]);
    if( (reject_range_min<x[0] && x[0]<reject_range_max)){
        TF1::RejectPoint();
        C1D = CLR;
    }else{
        C1D = CLR;
    }
    return C1D;
}

// double fit_cluster_LL_withLR(double* x, double* par){
// double C1D=0;
// double CLR = (par[0]+(par[1]-par[0])/pow((1.+pow(x[0]/par[2],par[3])),par[4]))*(1.+par[5]*x[0]);
// C1D = CLR;
// return C1D;
// }

double fit_cluster_LL_sys(double* x, double* par){
double C1D=0;
double CLR = par[5]*(par[0]+(par[1]-par[0])/(1.+pow(x[0]/par[2],par[3])))*(1.+par[4]*x[0]);
C1D = CLR;
return C1D;
}

Double_t fit_cluster_LL(Double_t* x, Double_t* par){
Double_t C1D=0;
Double_t CLR = (par[0]+(par[1]-par[0])/(1.+pow(x[0]/par[2],par[3])))*(1+par[4]*x[0]);
C1D = CLR;
return C1D;
}

Double_t fit_cluster_LL2(Double_t* x, Double_t* par){
Double_t C1D=0;
Double_t CLR = par[5]*(par[0]+(par[1]-par[0])/(1.+pow(x[0]/par[2],par[3])))*(1+par[4]*x[0]);
C1D = CLR;
return C1D;
}

//===========================================
//Strong -> LL and LAL
//===========================================

double fit_QS_LL(double* x, double* par){
 
    double C1D=0;
    
    //define quantities
    double kstar = 0.5 * x[0] / hbarc; //k* -> 1/fm
    double Rinv = par[1]; //fm
    double qR = 2.* kstar * Rinv;     // qinvRinv or 2k*Rinv-> no unit
    double lambda = par[2]; //lambda parameter
    double normalize = par[0]; //normalization
    //cluster parameters bellow
    double Bpar = par[3];
    double sigmainv = par[4];
    double epsilon = par[5];
    double Nclust = par[6];
    
    double CQS = -0.5*exp(-pow(qR,2)); //quantum statistics

    double CLR = Nclust*(1. + Bpar*exp(-pow(fabs(sigmainv*x[0])/hbarc,2)))*(1.+epsilon*x[0]); //cluster
    
    C1D = normalize * (1. + lambda * (CQS) ) * CLR; //total

    return C1D;
    
}

double fit_QS_KK_Gauss(double* x, double* par){
 
    double C1D=0;
    
    double reject_range_min=0.9;
    double reject_range_max=1.4;
    
    //define quantities
    double kstar = 0.5 * x[0] / hbarc; //k* -> 1/fm
    double Rinv = par[1]; //fm
    double qR = 2.* kstar * Rinv;     // qinvRinv or 2k*Rinv-> no unit
    double lambda = par[2]; //lambda parameter
    double normalize = par[0]; //normalization
    //cluster parameters bellow
    double Bpar = par[3];
    double sigmainv = par[4];
    double epsilon = par[5];
    double Nclust = par[6];
    
    double CQS = exp(-pow(qR,2)); //quantum statistics

    double CLR = Nclust*(1. + Bpar*exp(-pow(fabs(sigmainv*x[0])/hbarc,2)))*(1.+epsilon*x[0]); //cluster
    if( (reject_range_min<x[0] && x[0]<reject_range_max)){
        TF1::RejectPoint();    
    C1D = normalize * (1. + lambda * (CQS) ) * CLR; //total
    }else{
    C1D = normalize * (1. + lambda * (CQS) ) * CLR; //total
    }
    return C1D;
    
}

double fit_QS_KK_Exp(double* x, double* par){
 
    double C1D=0;
    
    double reject_range_min=0.9;
    double reject_range_max=1.4;
    
    //define quantities
    double kstar = 0.5 * x[0] / hbarc; //k* -> 1/fm
    double Rinv = par[1]; //fm
    double qR = 2.* kstar * Rinv;     // qinvRinv or 2k*Rinv-> no unit
    double lambda = par[2]; //lambda parameter
    double normalize = par[0]; //normalization
    //cluster parameters bellow
    double Bpar = par[3];
    double sigmainv = par[4];
    double epsilon = par[5];
    double Nclust = par[6];
    
    double CQS = exp(-pow(qR,1)); //quantum statistics

    double CLR = Nclust*(1. + Bpar*exp(-pow(fabs(sigmainv*x[0])/hbarc,2)))*(1.+epsilon*x[0]); //cluster
    if( (reject_range_min<x[0] && x[0]<reject_range_max)){
        TF1::RejectPoint();    
    C1D = normalize * (1. + lambda * (CQS) ) * CLR; //total
    }else{
    C1D = normalize * (1. + lambda * (CQS) ) * CLR; //total
    }
    return C1D;
    
}

double fit_strong_LAL(double* x, double* par){
    
    double C1D=0;
    
    //define quantities
    double kstar = 0.5 * x[0] / hbarc; //k* -> 1/fm
    double Rinv = par[1]; //fm
    double qR = 2.* kstar * Rinv;     // qinvRinv or 2k*Rinv-> no unit
    double Ref0 = par[2]; //real part of f0 in fm
    double Imf0 = par[3]; //imaginary part of f0 in fm
    double d0 = par[4]; //effective range in fm
    double alpha = par[9]; //0 for non identical, -0.5 for identical baryons
    double lambda = par[8]; //lambda parameter
    double normalize = par[0]; //normalization
    //cluster parameters bellow
    double Bpar = par[5];
    double sigmainv = par[6];
    double epsilon = par[7];
    double Nclust = par[10];
    
    const complex<double> i(0.0, 1.0); //imaginary number

    const complex<double> f0(Ref0, fabs(Imf0)); 
    complex<double> f = pow(pow(f0,-1) + 0.5*d0*pow(kstar,2) - i*kstar,-1);
    double fnorm = std::norm(f);
    double ReF = f.real();
    double ImF = f.imag();

    double CSI_p1 = 0.5 * ( fnorm / pow(Rinv,2.) ) * (1. - (d0 / ( 2. * sqrt(PI) * Rinv ) ) );
    double CSI_p2 = (2. * ReF / (sqrt(PI) * Rinv) ) * GetLednickyF1(qR);
    double CSI_p3 = (- ImF / Rinv) * GetLednickyF2(qR);
    double CSI = (1. + alpha) * ( CSI_p1 + CSI_p2 + CSI_p3  ); //strong interaction

    double CQS = alpha * exp(-pow(kstar*Rinv,2)); //quantum statistics

    double CLR = Nclust*(1. + Bpar*exp(-pow(fabs(sigmainv*x[0])/hbarc,2)))*(1.+epsilon*x[0]); //cluster
    
    C1D = normalize * (1. + lambda * (CQS+CSI) ) * CLR; //total

    return C1D;
    
}

double fit_strong_LALX(double* x, double* par){
    
    double C1D=0;
    
    //define quantities
    double kstar = 0.5 * x[0] / hbarc; //k* -> 1/fm
    double Rinv = par[1]; //fm
    double qR = 2.* kstar * Rinv;     // qinvRinv or 2k*Rinv-> no unit
    double Ref0 = par[2]; //real part of f0 in fm
    double Imf0 = par[3]; //imaginary part of f0 in fm
    double d0 = par[4]; //effective range in fm
    double alpha = par[9]; //0 for non identical, -0.5 for identical baryons
    double lambda = par[8]; //lambda parameter
    double normalize = par[0]; //normalization
    //cluster parameters bellow
    double Bpar = par[5];
    double sigmainv = par[6];
    double epsilon = par[7];
    double Nclust = par[10];
    double Aclust = par[11];
    
    const complex<double> i(0.0, 1.0); //imaginary number

    const complex<double> f0(Ref0, fabs(Imf0)); 
    complex<double> f = pow(pow(f0,-1) + 0.5*d0*pow(kstar,2) - i*kstar,-1);
    double fnorm = std::norm(f);
    double ReF = f.real();
    double ImF = f.imag();

    double CSI_p1 = 0.5 * ( fnorm / pow(Rinv,2.) ) * (1. - (d0 / ( 2. * sqrt(PI) * Rinv ) ) );
    double CSI_p2 = (2. * ReF / (sqrt(PI) * Rinv) ) * GetLednickyF1(qR);
    double CSI_p3 = (- ImF / Rinv) * GetLednickyF2(qR);
    double CSI = (1. + alpha)*( CSI_p1 + CSI_p2 + CSI_p3  ); //strong interaction

    double CQS = alpha * exp(-pow(qR,2)); //quantum statistics

//     double CLR = Nclust*(1. + Bpar*exp(-pow(fabs(sigmainv*x[0])/hbarc,2)))*(1.+epsilon*x[0]); //cluster
    double_t CLR = (Nclust+(Bpar-Nclust)/(1.+pow(x[0]/sigmainv,Aclust)))*(1+epsilon*x[0]);
    
    C1D = normalize * (1. + lambda * (CQS+CSI) ) * CLR; //total

    return C1D;
    
}

double fit_strong_LAL_func(double q, double R, double f0real, double f0imag, double d0sc, double alphacorr, double lamb, double norm, double Bcluster, double sinvcluster, double ecluster, double Ncluster){
    
    double C1D=0;
    
    double reject_range_min=0.48;
    double reject_range_max=0.65;
    
    //define quantities
    double kstar = 0.5 * q / hbarc; //k* -> 1/fm
    double Rinv = R; //k* -> 1/fm
    double qR = 2.* kstar * Rinv;     // qinvRinv or 2k*Rinv-> no unit
    double Ref0 = f0real; //real part of f0 in fm
    double Imf0 = f0imag; //imaginary part of f0 in fm
    double d0 = d0sc; //effective range
    double alpha = alphacorr; //0 for non identical, -0.5 for identical baryons
    double lambda = lamb; //lambda parameter
    //cluster parameters bellow
    double Bpar = Bcluster;
    double sigmainv = sinvcluster;
    double epsilon = ecluster;
    double Nclust = Ncluster;
    double normalize = norm; //normalization
    
    const complex<double> i(0.0, 1.0); //imaginary number

    const complex<double> f0(Ref0, fabs(Imf0)); 
    complex<double> f = pow(pow(f0,-1) + 0.5*d0*pow(kstar,2) - i*kstar,-1);
    double fnorm = std::norm(f);
    double ReF = f.real();
    double ImF = f.imag();

    double CSI_p1 = 0.5 * ( fnorm / pow(Rinv,2.) ) * (1. - (d0 / ( 2. * sqrt(PI) * Rinv ) ) );
    double CSI_p2 = (2. * ReF / (sqrt(PI) * Rinv) ) * GetLednickyF1(qR);
    double CSI_p3 = (- ImF / Rinv) * GetLednickyF2(qR);
    double CSI = (1. + alpha)*( CSI_p1 + CSI_p2 + CSI_p3  ); //strong interaction

    double CQS = alpha * exp(- pow(qR,2)); //quantum statistics

    double CLR = Nclust*(1 + Bpar*exp(-pow(fabs(sigmainv*q)/hbarc,2)))*(1+epsilon*q); //cluster
    
    C1D = normalize * (1. + lambda * (CQS+CSI) ) * CLR; //total

    return C1D;
    
}

double fit_strong_LL(double* x, double* par){
    
    double C1D=0;
    
    double reject_range_min=0.25;
    double reject_range_max=0.44;
//    double reject_range_min=0.24;
//    double reject_range_max=0.4; //for bin 60
    
    //define quantities
    double kstar = 0.5 * x[0] / hbarc; //k* -> 1/fm
    double Rinv = par[1]; //fm
    double qR = 2.* kstar * Rinv;     // qinvRinv or 2k*Rinv-> no unit
    double Ref0 = par[2]; //real part of f0 in fm
    double Imf0 = par[3]; //imaginary part of f0 in fm
    double d0 = par[4]; //effective range in fm
    double alpha = par[9]; //0 for non identical, -0.5 for identical baryons
    double lambda = par[8]; //lambda parameter
    double normalize = par[0]; //normalization
    //cluster parameters bellow
    double Bpar = par[5];
    double sigmainv = par[6];
    double epsilon = par[7];
    double Nclust = par[10];
    
    const complex<double> i(0.0, 1.0); //imaginary number

    const complex<double> f0(Ref0, fabs(Imf0)); 
    complex<double> f = pow(pow(f0,-1) + 0.5*d0*pow(kstar,2) - i*kstar,-1);
    double fnorm = std::norm(f);
    double ReF = f.real();
    double ImF = f.imag();

    double CSI_p1 = 0.5 * ( fnorm / pow(Rinv,2.) ) * (1. - (d0 / ( 2. * sqrt(PI) * Rinv ) ) );
    double CSI_p2 = (2. * ReF / (sqrt(PI) * Rinv) ) * GetLednickyF1(qR);
    double CSI_p3 = (- ImF / Rinv) * GetLednickyF2(qR);
    double CSI = (1. + alpha) * ( CSI_p1 + CSI_p2 + CSI_p3  ); //strong interaction

    double CQS = alpha * exp(-pow(qR,2)); //quantum statistics

    double CLR = Nclust*(1. + Bpar*exp(-pow(fabs(sigmainv*x[0])/hbarc,2)))*(1.+epsilon*x[0]); //cluster
    
    if( (reject_range_min<x[0] && x[0]<reject_range_max)){
        TF1::RejectPoint();    
        C1D = normalize * (1. + lambda * (CQS+CSI) ) * CLR; //total
    }else{
        C1D = normalize * (1. + lambda * (CQS+CSI) ) * CLR; //total
    }

    return C1D;
    
}

//===========================================
//Strong -> K0s-Lambda+K0s-AntiLambda
//===========================================

double fit_strong_KL(double* x, double* par){
    
    double C1D=0;
    
    double reject_range_min=0.48;
    double reject_range_max=0.65;
   
    //define quantities
    double kstar = 0.5 * x[0] / hbarc; //k* -> 1/fm
    double Rinv = par[1]; //k* -> 1/fm
    double qR = 2.* kstar * Rinv;     // qinvRinv or 2k*Rinv-> no unit
    double Ref0 = par[2]; //real part of f0 in fm
    double Imf0 = par[3]; //imaginary part of f0 in fm
    double d0 = par[4]; //effective range
    double alpha = par[9]; //0 for non identical, 0.5 for identical baryons
    double lambda = par[8]; //lambda parameter
    double normalize = par[0]; //normalization
    //cluster parameters bellow
    double Bpar = par[5];
    double sigmainv = par[6];
    double epsilon = par[7];
    double Nclust = par[10];
    
    const complex<double> i(0.0, 1.0); //imaginary number

    const complex<double> f0(Ref0, fabs(Imf0)); 
    complex<double> f = pow(pow(f0,-1) + 0.5*par[4]*pow(kstar,2) - i*kstar,-1);
    double fnorm = norm(f);
    double ReF = f.real();
    double ImF = f.imag();

    double CSI_p1 = 0.5 * ( fnorm / pow(Rinv,2.) ) * (1. - (d0 / ( 2. * sqrt(PI) * Rinv ) ) );
    double CSI_p2 = (2. * ReF / (sqrt(PI) * Rinv) ) * GetLednickyF1(qR);
    double CSI_p3 = (- ImF / Rinv) * GetLednickyF2(qR);
    double CSI = (1. + alpha) * ( CSI_p1 + CSI_p2 + CSI_p3  ); //strong interaction

    double CQS = alpha * exp(- pow(qR,2)); //quantum statistics
    
    double CLR = Nclust*(1 + Bpar*exp(-pow(fabs(sigmainv*x[0])/hbarc,2)))*(1+epsilon*x[0]); //cluster
    if( (reject_range_min<x[0] && x[0]<reject_range_max)){
        TF1::RejectPoint();    
        C1D = normalize * (1. + lambda * (CQS+CSI) ) * CLR; //total
    }else{
        C1D = normalize * (1. + lambda * (CQS+CSI) ) * CLR; //total
    }

//     if(Rinv<d0) break;
    
    return C1D;
    
}


//===========================================
//Strong -> K0s-K0s+K0s-K0s
//===========================================

double fit_strong_K0K0(double* x, double* par){
    
    double C1D=0;
    
    complex<double> i(0.0, 1.0); 

    double reject_range_min=0.9;   //f2prime
    double reject_range_max=1.4;  //f2prime

//    double reject_range_min2=0.15;  //Phi
//    double reject_range_max2=0.35; //Phi
    
    //define quantities
    double normalize = par[0]; //normalization
    double Rinv = par[1]; //k* -> 1/fm
    double lambda = par[2]; //lambda parameter

    double kstar = 0.5 * x[0] / hbarc; //k* -> 1/fm
    double qR = 2 *  kstar * Rinv;     // qinvRinv or 2k*Rinv-> no unit

    //set1 -> Antonelli
    //double m_f0 = 0.973, gamma_f0_KK = 2.763, gamma_f0_pipi = 0.5283;
    //double m_a0 = 0.985, gamma_a0_KK = 0.4038, gamma_a0_pieta = 0.3711;

    //set2 -> Achasov1
    //double m_f0 = 0.996, gamma_f0_KK = 1.305, gamma_f0_pipi = 0.2684;
    //double m_a0 = 0.992, gamma_a0_KK = 0.5555, gamma_a0_pieta = 0.4401;

    //set3 -> Achasov2
    //double m_f0 = 0.996, gamma_f0_KK = 1.305, gamma_f0_pipi = 0.2684;
    //double m_a0 = 1.003, gamma_a0_KK = 0.8365, gamma_a0_pieta = 0.4580;

    //set4 -> Martin
    //double m_f0 = 0.978, gamma_f0_KK = 0.792, gamma_f0_pipi = 0.1990;
    //double m_a0 = 0.974, gamma_a0_KK = 0.3330, gamma_a0_pieta = 0.2220;    

    double mass_f0 = par[3];
    double gamma_f0_KK = par[4];
    double gamma_f0_pipi = par[5];

    double mass_a0 = par[6];
    double gamma_a0_KK = par[7];
    double gamma_a0_pieta = par[8];

    //cluster parameters bellow
    double Bpar = par[9];
    double sigmainv = par[10];
    double epsilon = par[11];
    double Nclust = par[12];
    
    double kstarPrimepipi = sqrt((pow(mass_f0,2.)  - pow(pi_mass+pi_mass,2.))*(pow(mass_f0,2.) - pow(pi_mass-pi_mass,2.)))/(2.*mass_f0);
    double kstarPrimepieta = sqrt((pow(mass_a0,2.)  - pow(pi_mass+eta_mass,2.))*(pow(mass_a0,2.) - pow(pi_mass-eta_mass,2.)))/(2.*mass_a0);    

    double s = 4*( pow(k0s_mass,2.) + pow(0.5*x[0],2.) ); //mandelstam s
    
    complex<double> f_f0 = (gamma_f0_KK) * pow( pow(mass_f0,2.) - s - i*(gamma_f0_KK*0.5*x[0] + kstarPrimepipi*gamma_f0_pipi),-1);
    complex<double> f_a0 = (gamma_a0_KK) * pow( pow(mass_a0,2.) - s - i*(gamma_a0_KK*0.5*x[0] + kstarPrimepieta*gamma_a0_pieta),-1);

    complex<double> f = 0.5*(f_f0 + f_a0)*hbarc; //-> GeV^-1 to fm   
    
    double fnorm = norm(f);
    double ReF = f.real();
    double ImF = f.imag();

    double CSI_p1 = 0.5 * ( fnorm / pow(Rinv,2.) );
    double CSI_p2 = (2. * ReF / (sqrt(PI) * Rinv) ) * GetLednickyF1(qR);
    double CSI_p3 = (- ImF / Rinv) * GetLednickyF2(qR);
    double CSI = par[13]*( CSI_p1 + CSI_p2 + CSI_p3  ); //strong interaction

    double CQS = exp(- pow(qR,2)); //quantum statistics

    double CLR = Nclust*(1 + Bpar*exp(-pow(fabs(sigmainv*x[0])/hbarc,2)))*(1+epsilon*x[0]); //cluster
    
    if( (reject_range_min<x[0] && x[0]<reject_range_max)){
        TF1::RejectPoint();    
        C1D = normalize * (1. + lambda * (CQS+CSI) ) * CLR; //total
    }else{
        C1D = normalize * (1. + lambda * (CQS+CSI) ) * CLR; //total
    }

    return C1D;
 
    
}


double fit_strong_K0K0X(double* x, double* par){
    
    double C1D=0;
    
    complex<double> i(0.0, 1.0); 

    double reject_range_min=1.0;   //f2prime
    double reject_range_max=1.3;  //f2prime

//    double reject_range_min2=0.15;  //Phi
//    double reject_range_max2=0.35; //Phi
    
    //define quantities
    double normalize = par[0]; //normalization
    double Rinv = par[1]; //k* -> 1/fm
    double lambda = par[2]; //lambda parameter

    double kstar = 0.5 * x[0] / hbarc; //k* -> 1/fm
    double qR = 2.* kstar * Rinv;     // qinvRinv or 2k*Rinv-> no unit

    //set1 -> Antonelli
    //double m_f0 = 0.973, gamma_f0_KK = 2.763, gamma_f0_pipi = 0.5283;
    //double m_a0 = 0.985, gamma_a0_KK = 0.4038, gamma_a0_pieta = 0.3711;

    //set2 -> Achasov1
    //double m_f0 = 0.996, gamma_f0_KK = 1.305, gamma_f0_pipi = 0.2684;
    //double m_a0 = 0.992, gamma_a0_KK = 0.5555, gamma_a0_pieta = 0.4401;

    //set3 -> Achasov2
    //double m_f0 = 0.996, gamma_f0_KK = 1.305, gamma_f0_pipi = 0.2684;
    //double m_a0 = 1.003, gamma_a0_KK = 0.8365, gamma_a0_pieta = 0.4580;

    //set4 -> Martin
    //double m_f0 = 0.978, gamma_f0_KK = 0.792, gamma_f0_pipi = 0.1990;
    //double m_a0 = 0.974, gamma_a0_KK = 0.3330, gamma_a0_pieta = 0.2220;    

    double mass_f0 = par[3];
    double gamma_f0_KK = par[4];
    double gamma_f0_pipi = par[5];

    double mass_a0 = par[6];
    double gamma_a0_KK = par[7];
    double gamma_a0_pieta = par[8];

    //cluster parameters bellow
    double Bpar = par[9];
    double sigmainv = par[10];
    double epsilon = par[11];
    double Nclust = par[12];
    
    double kstarPrimepipi = sqrt((pow(mass_f0,2.)  - pow(pi_mass+pi_mass,2.))*(pow(mass_f0,2.) - pow(pi_mass-pi_mass,2.)))/(2.*mass_f0);
    double kstarPrimepieta = sqrt((pow(mass_a0,2.)  - pow(pi_mass+eta_mass,2.))*(pow(mass_a0,2.) - pow(pi_mass-eta_mass,2.)))/(2.*mass_a0);    

    double s = 4*( pow(k0s_mass,2.) + pow(0.5*x[0],2.) ); //mandelstam s
    
    complex<double> f_f0 = (gamma_f0_KK) * pow( pow(mass_f0,2.) - s - i*(gamma_f0_KK*0.5*x[0] + kstarPrimepipi*gamma_f0_pipi),-1);
    complex<double> f_a0 = (gamma_a0_KK) * pow( pow(mass_a0,2.) - s - i*(gamma_a0_KK*0.5*x[0] + kstarPrimepieta*gamma_a0_pieta),-1);

    complex<double> f = 0.5*(f_f0 + f_a0)*hbarc; //-> GeV^-1 to fm   
    
    double fnorm = norm(f);
    double ReF = f.real();
    double ImF = f.imag();

    double CSI_p1 = 0.5 * ( fnorm / pow(Rinv,2.) );
    double CSI_p2 = (2. * ReF / (sqrt(PI) * Rinv) ) * GetLednickyF1(qR);
    double CSI_p3 = (- ImF / Rinv) * GetLednickyF2(qR);
    double CSI = ( CSI_p1 + CSI_p2 + CSI_p3  ); //strong interaction

    double CQS = exp(- pow(qR,2)); //quantum statistics

    double CLR = Nclust*(1. + Bpar*exp(-pow(fabs(sigmainv*x[0])/hbarc,2)))*(1.+epsilon*x[0]); //cluster
    
    if( (reject_range_min<x[0] && x[0]<reject_range_max)){
        TF1::RejectPoint();    
        C1D = normalize * (1.+lambda * (CQS+CSI) ) * CLR; //total
    }else{
        C1D = normalize * (1.+lambda * (CQS+CSI) ) * CLR; //total
    }

    return C1D;
 
    
}

void normalize(double qminN, double qmaxN, TH1D* aux_histo_num, TH1D* aux_histo_den){

Int_t bin_for_normInt_min = aux_histo_num->GetXaxis()->FindBin(qminN);
Int_t bin_for_normInt_max = aux_histo_num->GetXaxis()->FindBin(qmaxN);
double int_num_controlRegion = aux_histo_num->Integral(bin_for_normInt_min,bin_for_normInt_max);
double int_den_controlRegion = aux_histo_den->Integral(bin_for_normInt_min,bin_for_normInt_max);
aux_histo_den->Scale(int_num_controlRegion/int_den_controlRegion);

}

void sideband(TH1D* aux_histo_ss, TH1D* aux_histo_sb, TH1D* aux_histo_bb, double fss, double fsb, double fbb){

aux_histo_ss->Scale(1./aux_histo_ss->Integral());
aux_histo_sb->Scale(fsb/aux_histo_sb->Integral());
aux_histo_bb->Scale(fbb/aux_histo_bb->Integral());
aux_histo_ss->Add(aux_histo_sb);
aux_histo_ss->Add(aux_histo_bb);
aux_histo_ss->Scale(1/fss);

}

double lineatone(double* x, double* par){
double C1D=0;
C1D = par[0];
return C1D;
}


//LAL

///signal and bkg qinv distributions

TH1D * hist_sig_global;
TH1D * hist_bkg_global;

TGraphErrors* gr_cluster_sig;
TGraphErrors* gr_cluster_bkg;

TGraphErrors* gr_cluster_sig_global;
TGraphErrors* gr_cluster_bkg_global;

//sideband fit

//Compute Log-Likelihood ratio for cluster
void log_likelihood_ratio_LAL_cluster(int &npar, double *gin, double &f, double *par, int iflag){

   double theLikelihood = 0.0; //-2*ln(L)
   //need to do a loop in the signal and bkg histos -- assuming 100 bins 0to2GeV in qinv
   double middle_of_bin = 0.0;
   
   for (Int_t ibin=0; ibin<75; ibin++){

      middle_of_bin = 0.02 + ibin*0.04; 

      if(middle_of_bin<=0.5)continue; //do not use first bin

      double hsig_bin_content = (double)gr_cluster_sig->Eval(middle_of_bin);
      double hbkg_bin_content = (double)gr_cluster_bkg->Eval(middle_of_bin);

      double Ai = hsig_bin_content;
      double Bi = hbkg_bin_content;
	  double Ci = fit_strong_LAL(&middle_of_bin,par);
     
      //-2*ln(L). See Section 4.3 of ATLAS pPb paper
      theLikelihood = theLikelihood + 2.*( Ai*TMath::Log( ((1.+Ci)*Ai)/(Ci*(Ai+Bi+2.)) ) + (Bi+2.)*TMath::Log( ((1.+Ci)*(Bi+2.))/(Ai+Bi+2.) ) );
      //theLikelihood = theLikelihood - 2.*(Ai - Bi*Ci + (Sai*Sai + Sbi*Sbi*Ci*Ci)*TMath::Log( 1. - (Ai - Bi*Ci)/(Sai*Sai + Sbi*Sbi*Ci*Ci) ));
   }
   f = theLikelihood;   
}

void log_likelihood_ratio_LAL_X(int &npar, double *gin, double &f, double *par, int iflag){

   double theLikelihood = 0.0; //-2*ln(L)
   //need to do a loop in the signal and bkg histos -- assuming 100 bins 0to2GeV in qinv
   double middle_of_bin = 0.0;
   
   for (Int_t ibin=0; ibin<75; ibin++){

      middle_of_bin = 0.02 + ibin*0.04; 

      if(ibin==0)continue; //do not use first bin
      
      double hsig_bin_content = (double)gr_cluster_sig_global->Eval(middle_of_bin);
      double hbkg_bin_content = (double)gr_cluster_bkg_global->Eval(middle_of_bin);

      double Ai = hsig_bin_content;
      double Bi = hbkg_bin_content;
	  double Ci = fit_strong_LAL(&middle_of_bin,par);
     
      //-2*ln(L). See Section 4.3 of ATLAS pPb paper
      theLikelihood = theLikelihood + 2.*( Ai*TMath::Log( ((1.+Ci)*Ai)/(Ci*(Ai+Bi+2.)) ) + (Bi+2.)*TMath::Log( ((1.+Ci)*(Bi+2.))/(Ai+Bi+2.) ) );
      //theLikelihood = theLikelihood - 2.*(Ai - Bi*Ci + (Sai*Sai + Sbi*Sbi*Ci*Ci)*TMath::Log( 1. - (Ai - Bi*Ci)/(Sai*Sai + Sbi*Sbi*Ci*Ci) ));
   }
   f = theLikelihood;   
}


//Compute Log-Likelihood ratio
void log_likelihood_ratio_LAL(int &npar, double *gin, double &f, double *par, int iflag){

   double theLikelihood = 0.0; //-2*ln(L)
   double middle_of_bin = 0.0;

   for (Int_t ibin=0; ibin<75; ibin++){

      middle_of_bin = 0.02 + ibin*0.04; 

      if(ibin==0)continue; //do not use first bin

      Int_t hsig_bin_number = (Int_t) hist_sig_global->FindBin(middle_of_bin);
      double hsig_bin_content = (double) hist_sig_global->GetBinContent(hsig_bin_number);   
      Int_t hbkg_bin_number = (Int_t) hist_bkg_global->FindBin(middle_of_bin);
      double hbkg_bin_content = (double) hist_bkg_global->GetBinContent(hbkg_bin_number);
      double Ai = hsig_bin_content;
      double Bi = hbkg_bin_content;
	  double Ci = fit_strong_LAL(&middle_of_bin,par);
     
      //-2*ln(L). See Section 4.3 of ATLAS pPb paper
      theLikelihood = theLikelihood + 2.*( Ai*TMath::Log( ((1.+Ci)*Ai)/(Ci*(Ai+Bi+2.)) ) + (Bi+2.)*TMath::Log( ((1.+Ci)*(Bi+2.))/(Ai+Bi+2.) ) );
      //theLikelihood = theLikelihood - 2.*(Ai - Bi*Ci + (Sai*Sai + Sbi*Sbi*Ci*Ci)*TMath::Log( 1. - (Ai - Bi*Ci)/(Sai*Sai + Sbi*Sbi*Ci*Ci) ));
       
   }
   f = theLikelihood;   
}


void print_test(int a){

if(a==0){
cout << endl;
cout << "------------------------------------"<< endl;
cout << "      not calculated at all"         << endl;
cout << "------------------------------------"<< endl;
cout << endl;
}else if(a==1){
cout << endl;
cout << "------------------------------------"<< endl;
cout << "  approximation only, not accurate"  << endl;
cout << "------------------------------------"<< endl;
cout << endl;
}else if(a==2){
cout << endl;
cout << "------------------------------------"<< endl;
cout << "full matrix, but forced positive-definite"<< endl;
cout << "------------------------------------"<< endl;
cout << endl;
}else if(a==3){
cout << endl;
cout << "------------------------------------"<< endl;
cout << "full accurate covariance matrix"<< endl;
cout << "------------------------------------"<< endl;
cout << endl;
}
}


void Do_Fitting(TMinuit* minuit, bool dominus, bool showcorr){

  // Run the fit procedure

  cout<<"Running Fit"<<endl;
  assert(minuit); // minuit needs to have been created by now
 
  double arglist[5] = {0,0,0,0,0}; //Arguments that can be passed with Minuit commands
  Int_t errFlag = 0;
  double tolerance = 10e-4;
  double MaxCalls  = 10000;//5000
  double verbosity = 0; //Print option: {-1,0,1,2}
  double strategy = 1; //Print option: {0,1,2}
  Int_t LLFlag = 2; //llist
  Int_t sigma = 1; // error in sigma

  //Define the error (sigmas)
  minuit->SetErrorDef(sigma); 

  // Set how verbose the output is (from no output at -1, to max at 3)
  arglist[0] = verbosity;
  minuit->mnexcm("SET PRINT", arglist, LLFlag, errFlag);

  // Set strategy
  arglist[0] = strategy;
  minuit->mnexcm("SET STR", arglist , LLFlag, errFlag);

  // Run Migrad with a specified max number of calls
  arglist[0] = MaxCalls; 
  arglist[1] = tolerance; 
  minuit->mnexcm("MINIMIZE", arglist, LLFlag, errFlag);
  if(showcorr){minuit->mnexcm("SHOw CORrelations", arglist, LLFlag, errFlag);minuit->mnexcm("SHOw COVariance", arglist, LLFlag, errFlag);}
  minuit->mnexcm("HESSE", arglist, LLFlag, errFlag);
  if(showcorr){minuit->mnexcm("SHOw CORrelations", arglist, LLFlag, errFlag);minuit->mnexcm("SHOw COVariance", arglist, LLFlag, errFlag);}

  if(dominus){minuit->mnexcm("MINOS", arglist, LLFlag, errFlag);}
  if(showcorr){minuit->mnexcm("SHOw CORrelations", arglist, LLFlag, errFlag);minuit->mnexcm("SHOw COVariance", arglist, LLFlag, errFlag);}

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
