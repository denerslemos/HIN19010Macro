#include "TCanvas.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
#include <sstream>

// Adjust relative path if needed for color settings
#include "./setSPRACEColors.C"

void LALMain(double Rf0, double Rf0_stat, double Rf0_sys, double If0, double If0_stat, double If0_sys, double d0, double d0_stat, double d0_sys) {

    const int LALpoints = 1;

    // =========================================================================
    // External Data Points & Models
    // =========================================================================

    // ALICE Lambda-antiLambda
    double LAL_ALICE_Rf0[LALpoints]      = {-0.90};
    double LAL_ALICE_Rf0_stat[LALpoints] = {0.04};
    double LAL_ALICE_Rf0_sys[LALpoints]  = {0.16};

    double LAL_ALICE_If0[LALpoints]      = {0.40};
    double LAL_ALICE_If0_stat[LALpoints] = {0.06};
    double LAL_ALICE_If0_sys[LALpoints]  = {0.18};

    double LAL_ALICE_d0[LALpoints]      = {2.76};
    double LAL_ALICE_d0_stat[LALpoints] = {0.29};
    double LAL_ALICE_d0_sys[LALpoints]  = {0.73};

    // ALICE p-antiLambda
    double pL_ALICE_Rf0[LALpoints]      = {-1.15};
    double pL_ALICE_Rf0_stat[LALpoints] = {0.05};
    double pL_ALICE_Rf0_sys[LALpoints]  = {0.23};

    double pL_ALICE_If0[LALpoints]      = {0.53};
    double pL_ALICE_If0_stat[LALpoints] = {0.04};
    double pL_ALICE_If0_sys[LALpoints]  = {0.15};

    double pL_ALICE_d0[LALpoints]      = {3.06};
    double pL_ALICE_d0_stat[LALpoints] = {0.14};
    double pL_ALICE_d0_sys[LALpoints]  = {0.98};

    // Models: p-pbar
    double ppbar_Grach_Rf0[LALpoints]      = {-1.2};
    double ppbar_Grach_Rf0_stat[LALpoints] = {0.0};
    double ppbar_Grach_If0[LALpoints]      = {0.9};
    double ppbar_Grach_If0_stat[LALpoints] = {0.0};
    double ppbar_Grach_d0[LALpoints]       = {1.0};
    double ppbar_Grach_d0_stat[LALpoints]  = {0.0};

    double ppbar_Pirner_Rf0[LALpoints] = {-1.0992908};
    double ppbar_Pirner_If0[LALpoints] = {0.40271493212669673};

    double ppbar_Klempt_Rf0[LALpoints]      = {-0.87949203};
    double ppbar_Klempt_Rf0_stat[LALpoints] = {0.0};
    double ppbar_Klempt_If0[LALpoints]      = {0.6402714932126695};
    double ppbar_Klempt_If0_stat[LALpoints] = {0.6945701357466062 - 0.6402714932126695};

    // CMS Lambda-antiLambda (from arguments)
    double LAL_CMS_Rf0[LALpoints]      = {Rf0};
    double LAL_CMS_Rf0_stat[LALpoints] = {Rf0_stat};
    double LAL_CMS_Rf0_sys[LALpoints]  = {Rf0_sys};

    double LAL_CMS_If0[LALpoints]      = {If0};
    double LAL_CMS_If0_stat[LALpoints] = {If0_stat};
    double LAL_CMS_If0_sys[LALpoints]  = {If0_sys};

    double LAL_CMS_d0[LALpoints]      = {d0};
    double LAL_CMS_d0_stat[LALpoints] = {d0_stat};
    double LAL_CMS_d0_sys[LALpoints]  = {d0_sys};

    // =========================================================================
    // Canvas 1: Im(f0) vs Re(f0)
    // =========================================================================

    // ALICE L-Lbar
    auto gr_LAL_ALICE_RvsI = new TGraphErrors(LALpoints, LAL_ALICE_Rf0, LAL_ALICE_If0, LAL_ALICE_Rf0_stat, LAL_ALICE_If0_stat);
    gr_LAL_ALICE_RvsI->SetMarkerStyle(21);
    gr_LAL_ALICE_RvsI->SetMarkerColor(kBlue+1);
    gr_LAL_ALICE_RvsI->SetMarkerSize(1.5);
    gr_LAL_ALICE_RvsI->SetLineColor(kBlue+1);
    
    auto gr_LAL_ALICE_RvsI_sys = new TGraphErrors(LALpoints, LAL_ALICE_Rf0, LAL_ALICE_If0, LAL_ALICE_Rf0_sys, LAL_ALICE_If0_sys);
    gr_LAL_ALICE_RvsI_sys->SetLineColor(kBlue+1);
    gr_LAL_ALICE_RvsI_sys->SetMarkerColor(0);
    gr_LAL_ALICE_RvsI_sys->SetFillColor(0);
    gr_LAL_ALICE_RvsI_sys->SetFillStyle(1);

    // ALICE p-Lbar
    auto gr_pL_ALICE_RvsI = new TGraphErrors(LALpoints, pL_ALICE_Rf0, pL_ALICE_If0, pL_ALICE_Rf0_stat, pL_ALICE_If0_stat);
    gr_pL_ALICE_RvsI->SetMarkerStyle(29);
    gr_pL_ALICE_RvsI->SetMarkerColor(kBlack);
    gr_pL_ALICE_RvsI->SetMarkerSize(2.0);
    gr_pL_ALICE_RvsI->SetLineColor(kBlack);

    auto gr_pL_ALICE_RvsI_sys = new TGraphErrors(LALpoints, pL_ALICE_Rf0, pL_ALICE_If0, pL_ALICE_Rf0_sys, pL_ALICE_If0_sys);
    gr_pL_ALICE_RvsI_sys->SetLineColor(kBlack);
    gr_pL_ALICE_RvsI_sys->SetMarkerColor(0);
    gr_pL_ALICE_RvsI_sys->SetFillColor(0);
    gr_pL_ALICE_RvsI_sys->SetFillStyle(1);

    // Theoretical Models
    auto gr_ppbar_Klempt_RvsI = new TGraphErrors(LALpoints, ppbar_Klempt_Rf0, ppbar_Klempt_If0, ppbar_Klempt_Rf0_stat, ppbar_Klempt_If0_stat);
    gr_ppbar_Klempt_RvsI->SetMarkerStyle(26);
    gr_ppbar_Klempt_RvsI->SetMarkerColor(kGreen+2);
    gr_ppbar_Klempt_RvsI->SetMarkerSize(1.5);
    gr_ppbar_Klempt_RvsI->SetLineColor(kGreen+2);

    auto gr_ppbar_Grach_RvsI = new TGraphErrors(LALpoints, ppbar_Grach_Rf0, ppbar_Grach_If0, ppbar_Grach_Rf0_stat, ppbar_Grach_If0_stat);
    gr_ppbar_Grach_RvsI->SetMarkerStyle(28);
    gr_ppbar_Grach_RvsI->SetMarkerColor(kMagenta);
    gr_ppbar_Grach_RvsI->SetMarkerSize(1.5);
    gr_ppbar_Grach_RvsI->SetLineColor(kMagenta);

    auto gr_ppbar_Pirner_RvsI = new TGraphErrors(LALpoints, ppbar_Pirner_Rf0, ppbar_Pirner_If0, 0, 0);
    gr_ppbar_Pirner_RvsI->SetMarkerStyle(44);
    gr_ppbar_Pirner_RvsI->SetMarkerColor(kGray+1);
    gr_ppbar_Pirner_RvsI->SetMarkerSize(1.5);
    gr_ppbar_Pirner_RvsI->SetLineColor(kGray+1);
    gr_ppbar_Pirner_RvsI->SetLineWidth(0);
    // CMS pPb
    auto gr_CMS = new TGraphErrors(LALpoints, LAL_CMS_Rf0, LAL_CMS_If0, LAL_CMS_Rf0_stat, LAL_CMS_If0_stat);
    gr_CMS->SetMarkerStyle(20);
    gr_CMS->SetMarkerColor(SPcolors[SPdred]);
    gr_CMS->SetMarkerSize(1.5);
    gr_CMS->SetLineColor(SPcolors[SPdred]);

    auto gr_CMS_sys = new TGraphErrors(LALpoints, LAL_CMS_Rf0, LAL_CMS_If0, LAL_CMS_Rf0_sys, LAL_CMS_If0_sys);
    gr_CMS_sys->SetMarkerColor(0);
    gr_CMS_sys->SetFillColor(0);
    gr_CMS_sys->SetFillStyle(1);
    gr_CMS_sys->SetLineColor(SPcolors[SPdred]);

    // Multigraph Assembly & Drawing (Canvas 1)
    TMultiGraph *mg1 = new TMultiGraph();
    mg1->Add(gr_CMS, "P");
    mg1->Add(gr_LAL_ALICE_RvsI, "P");
    mg1->Add(gr_pL_ALICE_RvsI, "P");
    mg1->Add(gr_ppbar_Klempt_RvsI, "P");
    mg1->Add(gr_ppbar_Grach_RvsI, "P");
    mg1->Add(gr_ppbar_Pirner_RvsI, "P");


    TCanvas *c1 = new TCanvas("c1", "", 1000, 700);
    mg1->Draw("a5");
    gr_pL_ALICE_RvsI_sys->Draw("5SAME");
    gr_LAL_ALICE_RvsI_sys->Draw("5SAME");
    gr_CMS->Draw("PSAME");
    gr_CMS_sys->Draw("5SAME");

    // Axis Formatting
    mg1->GetYaxis()->SetRangeUser(-0.5, 2.0);
    mg1->GetXaxis()->SetLimits(-1.5, -0.6);
    mg1->GetYaxis()->SetTitle("#scale[1.6]{#Jgothic f_{0} [fm]}");
    mg1->GetXaxis()->SetTitle("#scale[1.6]{#Rgothic f_{0} [fm]}");
    mg1->GetXaxis()->SetNdivisions(506);
    mg1->GetYaxis()->SetNdivisions(506);
    mg1->GetXaxis()->CenterTitle(true);
    mg1->GetYaxis()->CenterTitle(true);
    mg1->GetXaxis()->SetTickSize(0.03);
    mg1->GetYaxis()->SetTickSize(0.03);
    mg1->GetXaxis()->SetLabelSize(0.05);
    mg1->GetYaxis()->SetLabelSize(0.05);
    mg1->GetYaxis()->SetLabelOffset(0.0025);
    mg1->GetXaxis()->SetLabelOffset(0.0025);
    mg1->GetYaxis()->SetTitleOffset(1.3);
    mg1->GetXaxis()->SetTitleOffset(1.3);

    // Legend 1
    TLegend *legend1 = new TLegend(0.4, 0.55, 0.88, 0.88);
    legend1->SetBorderSize(0);
    legend1->AddEntry(gr_CMS, "#Lambda#bar{#Lambda}: CMS, pPb, #sqrt{s_{NN}} = 8.16 TeV", "P");
    legend1->AddEntry(gr_LAL_ALICE_RvsI, "#Lambda#bar{#Lambda}: ALICE, PbPb, #sqrt{s_{NN}} = 2.76, 5.02 TeV", "P");
    legend1->AddEntry(gr_pL_ALICE_RvsI, "p#bar{#Lambda} #oplus #bar{p}#Lambda: ALICE, PbPb, #sqrt{s_{NN}} = 2.76, 5.02 TeV", "P");
    legend1->AddEntry(gr_ppbar_Grach_RvsI, "p#bar{p} (singlet): PLB 208, 309", "P");
    legend1->AddEntry(gr_ppbar_Klempt_RvsI, "p#bar{p} (singlet): PR 368, 119", "P");
    legend1->AddEntry(gr_ppbar_Pirner_RvsI, "p#bar{p} (singlet): ZPA 338, 111" "");
//    legend1->AddEntry(gr_LAL_ALICE_RvsI_sys, "Systematic Uncertainties", "F");
//    legend1->AddEntry(gr_LAL_ALICE_RvsI, "Statistical Uncertainties", "L");
    legend1->Draw();

    // Labels
    TLatex texCMS1;
    texCMS1.SetNDC();
    texCMS1.SetTextFont(42);
    texCMS1.SetTextSize(0.05);
    texCMS1.SetTextAlign(12);
    texCMS1.DrawLatex(0.54, 0.17, "#font[61]{baryon-antibaryon}");
    texCMS1.DrawLatex(0.10, 0.93, "#font[61]{CMS} #it{Work in Progress}");
    texCMS1.DrawLatex(0.60, 0.931, "#scale[0.95]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");

    c1->Update();
    c1->Print("Rf0vsIf0LAL.pdf");

    // =========================================================================
    // Canvas 2: d0 vs Re(f0)
    // =========================================================================

    // ALICE L-Lbar
    auto gr_LAL_ALICE_Rvsd0 = new TGraphErrors(LALpoints, LAL_ALICE_Rf0, LAL_ALICE_d0, LAL_ALICE_Rf0_stat, LAL_ALICE_d0_stat);
    gr_LAL_ALICE_Rvsd0->SetMarkerStyle(21);
    gr_LAL_ALICE_Rvsd0->SetMarkerColor(kBlue+1);
    gr_LAL_ALICE_Rvsd0->SetMarkerSize(1.5);
    gr_LAL_ALICE_Rvsd0->SetLineColor(kBlue+1);
    
    auto gr_LAL_ALICE_Rvsd0_sys = new TGraphErrors(LALpoints, LAL_ALICE_Rf0, LAL_ALICE_d0, LAL_ALICE_Rf0_sys, LAL_ALICE_d0_sys);
    gr_LAL_ALICE_Rvsd0_sys->SetMarkerColor(0);
    gr_LAL_ALICE_Rvsd0_sys->SetFillColor(0);
    gr_LAL_ALICE_Rvsd0_sys->SetFillStyle(1);
    gr_LAL_ALICE_Rvsd0_sys->SetLineColor(kBlue+1);

    // ALICE p-Lbar
    auto gr_pL_ALICE_Rvsd0 = new TGraphErrors(LALpoints, pL_ALICE_Rf0, pL_ALICE_d0, pL_ALICE_Rf0_stat, pL_ALICE_d0_stat);
    gr_pL_ALICE_Rvsd0->SetMarkerStyle(29);
    gr_pL_ALICE_Rvsd0->SetMarkerColor(kBlack);
    gr_pL_ALICE_Rvsd0->SetMarkerSize(2.0);
    gr_pL_ALICE_Rvsd0->SetLineColor(kBlack);

    auto gr_pL_ALICE_Rvsd0_sys = new TGraphErrors(LALpoints, pL_ALICE_Rf0, pL_ALICE_d0, pL_ALICE_Rf0_sys, pL_ALICE_d0_sys);
    gr_pL_ALICE_Rvsd0_sys->SetMarkerColor(0);
    gr_pL_ALICE_Rvsd0_sys->SetFillColor(0);
    gr_pL_ALICE_Rvsd0_sys->SetFillStyle(1);
    gr_pL_ALICE_Rvsd0_sys->SetLineColor(SPcolors[SPgreen]);

    // Theoretical Models
    auto gr_ppbar_Grach_Rvsd0 = new TGraphErrors(LALpoints, ppbar_Grach_Rf0, ppbar_Grach_d0, ppbar_Grach_Rf0_stat, ppbar_Grach_d0_stat);
    gr_ppbar_Grach_Rvsd0->SetMarkerStyle(28);
    gr_ppbar_Grach_Rvsd0->SetMarkerColor(kMagenta);
    gr_ppbar_Grach_Rvsd0->SetMarkerSize(1.5);
    gr_ppbar_Grach_Rvsd0->SetLineColor(kMagenta);
    
    // CMS pPb
    auto gr_CMS_d0 = new TGraphErrors(LALpoints, LAL_CMS_Rf0, LAL_CMS_d0, LAL_CMS_Rf0_stat, LAL_CMS_d0_stat);
    gr_CMS_d0->SetMarkerStyle(20);
    gr_CMS_d0->SetMarkerColor(SPcolors[SPdred]);
    gr_CMS_d0->SetMarkerSize(1.5);
    gr_CMS_d0->SetLineColor(SPcolors[SPdred]);

    auto gr_CMS_d0_sys = new TGraphErrors(LALpoints, LAL_CMS_Rf0, LAL_CMS_d0, LAL_CMS_Rf0_sys, LAL_CMS_d0_sys);
    gr_CMS_d0_sys->SetMarkerColor(0);
    gr_CMS_d0_sys->SetFillColor(0);
    gr_CMS_d0_sys->SetFillStyle(1);
    gr_CMS_d0_sys->SetLineColor(SPcolors[SPred]);

    // Multigraph Assembly & Drawing (Canvas 2)
    TMultiGraph *mg2 = new TMultiGraph();
    mg2->Add(gr_LAL_ALICE_Rvsd0, "P");
    mg2->Add(gr_pL_ALICE_Rvsd0, "P");
    mg2->Add(gr_ppbar_Grach_Rvsd0, "P");
    mg2->Add(gr_CMS_d0, "P");

    TCanvas *c2 = new TCanvas("c2", "", 1000, 700);
    mg2->Draw("a");
    gr_LAL_ALICE_Rvsd0_sys->Draw("5SAME");
    gr_pL_ALICE_Rvsd0_sys->Draw("5SAME");
    gr_CMS_d0_sys->Draw("5SAME");
    gr_CMS_d0->Draw("PSAME");

    // Axis Formatting
    mg2->GetYaxis()->SetRangeUser(-0.05, 6.0);
    mg2->GetXaxis()->SetLimits(-1.5, -0.6);
    mg2->GetYaxis()->SetTitle("#scale[1.3]{d_{0} [fm]}");
    mg2->GetXaxis()->SetTitle("#scale[1.3]{#Rgothic f_{0} [fm]}");
    mg2->GetXaxis()->SetNdivisions(506);
    mg2->GetYaxis()->SetNdivisions(506);
    mg2->GetXaxis()->CenterTitle(true);
    mg2->GetYaxis()->CenterTitle(true);
    mg2->GetXaxis()->SetTickSize(0.03);
    mg2->GetYaxis()->SetTickSize(0.03);
    mg2->GetXaxis()->SetLabelSize(0.05);
    mg2->GetYaxis()->SetLabelSize(0.05);
    mg2->GetYaxis()->SetLabelOffset(0.0025);
    mg2->GetXaxis()->SetLabelOffset(0.0025);
    mg2->GetYaxis()->SetTitleOffset(1.3);
    mg2->GetXaxis()->SetTitleOffset(1.3);

    // Legend 2
    TLegend *legend2 = new TLegend(0.15, 0.7, 0.8, 0.88);
    legend2->SetBorderSize(0);
    legend2->AddEntry(gr_CMS_d0, "#Lambda#bar{#Lambda}: CMS, pPb, #sqrt{s_{NN}} = 8.16 TeV", "P");
    legend2->AddEntry(gr_LAL_ALICE_Rvsd0, "#Lambda#bar{#Lambda}: ALICE, PbPb, #sqrt{s_{NN}} = 2.76, 5.02 TeV", "P");
    legend2->AddEntry(gr_pL_ALICE_Rvsd0, "p#bar{#Lambda} #oplus #bar{p}#Lambda: ALICE, PbPb, #sqrt{s_{NN}} = 2.76, 5.02 TeV", "P");
    legend2->AddEntry(gr_ppbar_Grach_Rvsd0, "p#bar{p} (singlet): PLB 208, 309", "P");
//    legend2->AddEntry(gr_LAL_ALICE_Rvsd0_sys, "Systematic Uncertainties", "F");
//    legend2->AddEntry(gr_LAL_ALICE_Rvsd0, "Statistical Uncertainties", "L");
    legend2->Draw();

    // Labels
    TLatex texCMS2;
    texCMS2.SetNDC();
    texCMS2.SetTextFont(42);
    texCMS2.SetTextSize(0.05);
    texCMS2.SetTextAlign(12);
    texCMS2.DrawLatex(0.10, 0.93, "#font[61]{CMS} #it{Work in Progress}");
    texCMS2.DrawLatex(0.54, 0.17, "#font[61]{baryon-antibaryon}");
    texCMS2.DrawLatex(0.60, 0.931, "#scale[0.95]{pPb, #sqrt{s_{NN}} = 8.16 TeV}");

    c2->Update();
    c2->Print("Rf0vsd0LAL.pdf");
}