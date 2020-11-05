#include<TApplication.h>
#include<TFile.h>
#include<TMath.h>
#include<TMinuit.h>
#include<TROOT.h>
#include<TSystem.h>
#include<TTree.h>
#include<TVector2.h>

#include<TCanvas.h>
#include<TF1.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include<TLegend.h>
#include<TLine.h>
#include<TH2F.h>
#include<TPaveText.h>
#include<TStyle.h>

#include<algorithm>
#include<chrono>
#include<ctime>
#include<fstream>
#include<iostream>
#include<mutex>
#include<string>
#include<thread>
#include<vector>

#include<TEfficiency.h>
#include<TProfile.h>
#include<TPaveStats.h>
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "RooRealProxy.h"

//#include "./CruijffPdf.h"
//#include "./DoubleCBPdf.h"

using namespace std;
using namespace RooFit;

//DEFINE HISTOGRAMS

   //Total
   TH1F* h_Energy_EB_old;
   TH1F* h_Energy_EE_old;
   TH1F* h_Energy_EB_new;
   TH1F* h_Energy_EE_new;
   TH1F* h_Et_EB_old;
   TH1F* h_Et_EE_old;
   TH1F* h_Et_EB_new;
   TH1F* h_Et_EE_new;
   TH1F* h_Eta_old;
   TH1F* h_Eta_new;
   TH1F* h_Phi_EB_old;
   TH1F* h_Phi_EE_old;
   TH1F* h_Phi_EB_new;
   TH1F* h_Phi_EE_new;
   TH1F* h_EtaWidth_EB_old;
   TH1F* h_EtaWidth_EE_old;
   TH1F* h_EtaWidth_EB_new;
   TH1F* h_EtaWidth_EE_new;
   TH1F* h_PhiWidth_EB_old;
   TH1F* h_PhiWidth_EE_old;
   TH1F* h_PhiWidth_EB_new;
   TH1F* h_PhiWidth_EE_new;
   TH1F* h_nPFClusters_old;
   TH1F* h_nPFClusters_new;
   TH1F* h_nPFClusters_EB_old;
   TH1F* h_nPFClusters_EE_old;
   TH1F* h_nPFClusters_EB_new;
   TH1F* h_nPFClusters_EE_new;
   TH1F* h_R9_EB_old;
   TH1F* h_R9_EE_old;
   TH1F* h_R9_EB_new;
   TH1F* h_R9_EE_new;
   TH1F* h_full5x5_R9_EB_old;
   TH1F* h_full5x5_R9_EE_old;
   TH1F* h_full5x5_R9_EB_new;
   TH1F* h_full5x5_R9_EE_new;
   TH1F* h_sigmaIetaIeta_EB_old;
   TH1F* h_sigmaIetaIeta_EE_old;
   TH1F* h_sigmaIetaIeta_EB_new;
   TH1F* h_sigmaIetaIeta_EE_new;
   TH1F* h_full5x5_sigmaIetaIeta_EB_old;
   TH1F* h_full5x5_sigmaIetaIeta_EE_old;
   TH1F* h_full5x5_sigmaIetaIeta_EB_new;
   TH1F* h_full5x5_sigmaIetaIeta_EE_new;
   TH1F* h_sigmaIetaIphi_EB_old;
   TH1F* h_sigmaIetaIphi_EE_old;
   TH1F* h_sigmaIetaIphi_EB_new;
   TH1F* h_sigmaIetaIphi_EE_new;
   TH1F* h_full5x5_sigmaIetaIphi_EB_old;
   TH1F* h_full5x5_sigmaIetaIphi_EE_old;
   TH1F* h_full5x5_sigmaIetaIphi_EB_new;
   TH1F* h_full5x5_sigmaIetaIphi_EE_new;
   TH1F* h_sigmaIphiIphi_EB_old;
   TH1F* h_sigmaIphiIphi_EE_old;
   TH1F* h_sigmaIphiIphi_EB_new;
   TH1F* h_sigmaIphiIphi_EE_new;
   TH1F* h_full5x5_sigmaIphiIphi_EB_old;
   TH1F* h_full5x5_sigmaIphiIphi_EE_old;
   TH1F* h_full5x5_sigmaIphiIphi_EB_new;
   TH1F* h_full5x5_sigmaIphiIphi_EE_new;

   //seed matched plots
   TH1F* h_Energy_EB_seedMatched_old;
   TH1F* h_Energy_EE_seedMatched_old;
   TH1F* h_Energy_EB_seedMatched_new;
   TH1F* h_Energy_EE_seedMatched_new;
   TH1F* h_Et_EB_seedMatched_old;
   TH1F* h_Et_EE_seedMatched_old;
   TH1F* h_Et_EB_seedMatched_new;
   TH1F* h_Et_EE_seedMatched_new;
   TH1F* h_Eta_seedMatched_old;
   TH1F* h_Eta_seedMatched_new;
   TH1F* h_Phi_EB_seedMatched_old;
   TH1F* h_Phi_EE_seedMatched_old;
   TH1F* h_Phi_EB_seedMatched_new;
   TH1F* h_Phi_EE_seedMatched_new;
   TH1F* h_EtaWidth_EB_seedMatched_old;
   TH1F* h_EtaWidth_EE_seedMatched_old;
   TH1F* h_EtaWidth_EB_seedMatched_new;
   TH1F* h_EtaWidth_EE_seedMatched_new;
   TH1F* h_PhiWidth_EB_seedMatched_old;
   TH1F* h_PhiWidth_EE_seedMatched_old;
   TH1F* h_PhiWidth_EB_seedMatched_new;
   TH1F* h_PhiWidth_EE_seedMatched_new;
   TH1F* h_nPFClusters_seedMatched_old;
   TH1F* h_nPFClusters_seedMatched_new;
   TH1F* h_nPFClusters_EB_seedMatched_old;
   TH1F* h_nPFClusters_EE_seedMatched_old;
   TH1F* h_nPFClusters_EB_seedMatched_new;
   TH1F* h_nPFClusters_EE_seedMatched_new;
   TH1F* h_R9_EB_seedMatched_old;
   TH1F* h_R9_EE_seedMatched_old;
   TH1F* h_R9_EB_seedMatched_new;
   TH1F* h_R9_EE_seedMatched_new;
   TH1F* h_full5x5_R9_EB_seedMatched_old;
   TH1F* h_full5x5_R9_EE_seedMatched_old;
   TH1F* h_full5x5_R9_EB_seedMatched_new;
   TH1F* h_full5x5_R9_EE_seedMatched_new;
   TH1F* h_sigmaIetaIeta_EB_seedMatched_old;
   TH1F* h_sigmaIetaIeta_EE_seedMatched_old;
   TH1F* h_sigmaIetaIeta_EB_seedMatched_new;
   TH1F* h_sigmaIetaIeta_EE_seedMatched_new;
   TH1F* h_full5x5_sigmaIetaIeta_EB_seedMatched_old;
   TH1F* h_full5x5_sigmaIetaIeta_EE_seedMatched_old;
   TH1F* h_full5x5_sigmaIetaIeta_EB_seedMatched_new;
   TH1F* h_full5x5_sigmaIetaIeta_EE_seedMatched_new;
   TH1F* h_sigmaIetaIphi_EB_seedMatched_old;
   TH1F* h_sigmaIetaIphi_EE_seedMatched_old;
   TH1F* h_sigmaIetaIphi_EB_seedMatched_new;
   TH1F* h_sigmaIetaIphi_EE_seedMatched_new;
   TH1F* h_full5x5_sigmaIetaIphi_EB_seedMatched_old;
   TH1F* h_full5x5_sigmaIetaIphi_EE_seedMatched_old;
   TH1F* h_full5x5_sigmaIetaIphi_EB_seedMatched_new;
   TH1F* h_full5x5_sigmaIetaIphi_EE_seedMatched_new;
   TH1F* h_sigmaIphiIphi_EB_seedMatched_old;
   TH1F* h_sigmaIphiIphi_EE_seedMatched_old;
   TH1F* h_sigmaIphiIphi_EB_seedMatched_new;
   TH1F* h_sigmaIphiIphi_EE_seedMatched_new;
   TH1F* h_full5x5_sigmaIphiIphi_EB_seedMatched_old;
   TH1F* h_full5x5_sigmaIphiIphi_EE_seedMatched_old;
   TH1F* h_full5x5_sigmaIphiIphi_EB_seedMatched_new;
   TH1F* h_full5x5_sigmaIphiIphi_EE_seedMatched_new;

   //caloMatched
   TH1F* h_Energy_EB_caloMatched_old;
   TH1F* h_Energy_EE_caloMatched_old;
   TH1F* h_Energy_EB_caloMatched_new;
   TH1F* h_Energy_EE_caloMatched_new;
   TH1F* h_EoEtrue_EB_old;
   TH1F* h_EoEtrue_EE_old;
   TH1F* h_EoEtrue_EB_new;
   TH1F* h_EoEtrue_EE_new;
   TH1F* h_EoEgen_EB_old;
   TH1F* h_EoEgen_EE_old;
   TH1F* h_EoEgen_EB_new;
   TH1F* h_EoEgen_EE_new;
   TH1F* h_Et_EB_caloMatched_old;
   TH1F* h_Et_EE_caloMatched_old;
   TH1F* h_Et_EB_caloMatched_new;
   TH1F* h_Et_EE_caloMatched_new;
   TH1F* h_Eta_caloMatched_old;
   TH1F* h_Eta_caloMatched_new;
   TH1F* h_Phi_EB_caloMatched_old;
   TH1F* h_Phi_EE_caloMatched_old;
   TH1F* h_Phi_EB_caloMatched_new;
   TH1F* h_Phi_EE_caloMatched_new;
   TH1F* h_EtaWidth_EB_caloMatched_old;
   TH1F* h_EtaWidth_EE_caloMatched_old;
   TH1F* h_EtaWidth_EB_caloMatched_new;
   TH1F* h_EtaWidth_EE_caloMatched_new;
   TH1F* h_PhiWidth_EB_caloMatched_old;
   TH1F* h_PhiWidth_EE_caloMatched_old;
   TH1F* h_PhiWidth_EB_caloMatched_new;
   TH1F* h_PhiWidth_EE_caloMatched_new;
   TH1F* h_nPFClusters_caloMatched_old;
   TH1F* h_nPFClusters_caloMatched_new;
   TH1F* h_nPFClusters_EB_caloMatched_old;
   TH1F* h_nPFClusters_EE_caloMatched_old;
   TH1F* h_nPFClusters_EB_caloMatched_new;
   TH1F* h_nPFClusters_EE_caloMatched_new;
   TH1F* h_R9_EB_caloMatched_old;
   TH1F* h_R9_EE_caloMatched_old;
   TH1F* h_R9_EB_caloMatched_new;
   TH1F* h_R9_EE_caloMatched_new;
   TH1F* h_full5x5_R9_EB_caloMatched_old;
   TH1F* h_full5x5_R9_EE_caloMatched_old;
   TH1F* h_full5x5_R9_EB_caloMatched_new;
   TH1F* h_full5x5_R9_EE_caloMatched_new;
   TH1F* h_sigmaIetaIeta_EB_caloMatched_old;
   TH1F* h_sigmaIetaIeta_EE_caloMatched_old;
   TH1F* h_sigmaIetaIeta_EB_caloMatched_new;
   TH1F* h_sigmaIetaIeta_EE_caloMatched_new;
   TH1F* h_full5x5_sigmaIetaIeta_EB_caloMatched_old;
   TH1F* h_full5x5_sigmaIetaIeta_EE_caloMatched_old;
   TH1F* h_full5x5_sigmaIetaIeta_EB_caloMatched_new;
   TH1F* h_full5x5_sigmaIetaIeta_EE_caloMatched_new;
   TH1F* h_sigmaIetaIphi_EB_caloMatched_old;
   TH1F* h_sigmaIetaIphi_EE_caloMatched_old;
   TH1F* h_sigmaIetaIphi_EB_caloMatched_new;
   TH1F* h_sigmaIetaIphi_EE_caloMatched_new;
   TH1F* h_full5x5_sigmaIetaIphi_EB_caloMatched_old;
   TH1F* h_full5x5_sigmaIetaIphi_EE_caloMatched_old;
   TH1F* h_full5x5_sigmaIetaIphi_EB_caloMatched_new;
   TH1F* h_full5x5_sigmaIetaIphi_EE_caloMatched_new;
   TH1F* h_sigmaIphiIphi_EB_caloMatched_old;
   TH1F* h_sigmaIphiIphi_EE_caloMatched_old;
   TH1F* h_sigmaIphiIphi_EB_caloMatched_new;
   TH1F* h_sigmaIphiIphi_EE_caloMatched_new;
   TH1F* h_full5x5_sigmaIphiIphi_EB_caloMatched_old;
   TH1F* h_full5x5_sigmaIphiIphi_EE_caloMatched_old;
   TH1F* h_full5x5_sigmaIphiIphi_EB_caloMatched_new;
   TH1F* h_full5x5_sigmaIphiIphi_EE_caloMatched_new;

   //caloMatched && seedMatched
   TH1F* h_EoEtrue_EB_seedMatched_old;
   TH1F* h_EoEtrue_EE_seedMatched_old;
   TH1F* h_EoEtrue_EB_seedMatched_new;
   TH1F* h_EoEtrue_EE_seedMatched_new;
   TH1F* h_EoEgen_EB_seedMatched_old;
   TH1F* h_EoEgen_EE_seedMatched_old;
   TH1F* h_EoEgen_EB_seedMatched_new;
   TH1F* h_EoEgen_EE_seedMatched_new;

   //caloUnmatched
   TH1F* h_Energy_EB_caloUnmatched_old;
   TH1F* h_Energy_EE_caloUnmatched_old;
   TH1F* h_Energy_EB_caloUnmatched_new;
   TH1F* h_Energy_EE_caloUnmatched_new;
   TH1F* h_Et_EB_caloUnmatched_old;
   TH1F* h_Et_EE_caloUnmatched_old;
   TH1F* h_Et_EB_caloUnmatched_new;
   TH1F* h_Et_EE_caloUnmatched_new;
   TH1F* h_Eta_caloUnmatched_old;
   TH1F* h_Eta_caloUnmatched_new;
   TH1F* h_Phi_EB_caloUnmatched_old;
   TH1F* h_Phi_EE_caloUnmatched_old;
   TH1F* h_Phi_EB_caloUnmatched_new;
   TH1F* h_Phi_EE_caloUnmatched_new;
   TH1F* h_EtaWidth_EB_caloUnmatched_old;
   TH1F* h_EtaWidth_EE_caloUnmatched_old;
   TH1F* h_EtaWidth_EB_caloUnmatched_new;
   TH1F* h_EtaWidth_EE_caloUnmatched_new;
   TH1F* h_PhiWidth_EB_caloUnmatched_old;
   TH1F* h_PhiWidth_EE_caloUnmatched_old;
   TH1F* h_PhiWidth_EB_caloUnmatched_new;
   TH1F* h_PhiWidth_EE_caloUnmatched_new;
   TH1F* h_nPFClusters_caloUnmatched_old;
   TH1F* h_nPFClusters_caloUnmatched_new;
   TH1F* h_nPFClusters_EB_caloUnmatched_old;
   TH1F* h_nPFClusters_EE_caloUnmatched_old;
   TH1F* h_nPFClusters_EB_caloUnmatched_new;
   TH1F* h_nPFClusters_EE_caloUnmatched_new;
   TH1F* h_R9_EB_caloUnmatched_old;
   TH1F* h_R9_EE_caloUnmatched_old;
   TH1F* h_R9_EB_caloUnmatched_new;
   TH1F* h_R9_EE_caloUnmatched_new;
   TH1F* h_full5x5_R9_EB_caloUnmatched_old;
   TH1F* h_full5x5_R9_EE_caloUnmatched_old;
   TH1F* h_full5x5_R9_EB_caloUnmatched_new;
   TH1F* h_full5x5_R9_EE_caloUnmatched_new;
   TH1F* h_sigmaIetaIeta_EB_caloUnmatched_old;
   TH1F* h_sigmaIetaIeta_EE_caloUnmatched_old;
   TH1F* h_sigmaIetaIeta_EB_caloUnmatched_new;
   TH1F* h_sigmaIetaIeta_EE_caloUnmatched_new;
   TH1F* h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old;
   TH1F* h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old;
   TH1F* h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new;
   TH1F* h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new;
   TH1F* h_sigmaIetaIphi_EB_caloUnmatched_old;
   TH1F* h_sigmaIetaIphi_EE_caloUnmatched_old;
   TH1F* h_sigmaIetaIphi_EB_caloUnmatched_new;
   TH1F* h_sigmaIetaIphi_EE_caloUnmatched_new;
   TH1F* h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old;
   TH1F* h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old;
   TH1F* h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new;
   TH1F* h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new;
   TH1F* h_sigmaIphiIphi_EB_caloUnmatched_old;
   TH1F* h_sigmaIphiIphi_EE_caloUnmatched_old;
   TH1F* h_sigmaIphiIphi_EB_caloUnmatched_new;
   TH1F* h_sigmaIphiIphi_EE_caloUnmatched_new;
   TH1F* h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old;
   TH1F* h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old;
   TH1F* h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new;
   TH1F* h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new;

//efficiencies
   TH1F* h_Eta_Calo_Denum_old;
   TH1F* h_Eta_Calo_Denum_new;
   TH1F* h_Eta_Calo_old;
   TH1F* h_Eta_Calo_new;
   TH1F* h_Eta_Gen_Denum_old;
   TH1F* h_Eta_Gen_Denum_new;
   TH1F* h_Eta_Gen_old;
   TH1F* h_Eta_Gen_new;
   TH1F* h_Eta_Seed_Denum_old;
   TH1F* h_Eta_Seed_Denum_new;
   TH1F* h_Eta_Seed_old;
   TH1F* h_Eta_Seed_new;
   TH1F* h_Et_Calo_EB_Denum_old;
   TH1F* h_Et_Calo_EB_Denum_new;
   TH1F* h_Et_Calo_EB_old;
   TH1F* h_Et_Calo_EB_new;
   TH1F* h_Et_Calo_EE_Denum_old;
   TH1F* h_Et_Calo_EE_Denum_new;
   TH1F* h_Et_Calo_EE_old;
   TH1F* h_Et_Calo_EE_new;
   TH1F* h_Et_Gen_EB_Denum_old;
   TH1F* h_Et_Gen_EB_Denum_new;
   TH1F* h_Et_Gen_EB_old;
   TH1F* h_Et_Gen_EB_new;
   TH1F* h_Et_Gen_EE_Denum_old;
   TH1F* h_Et_Gen_EE_Denum_new;
   TH1F* h_Et_Gen_EE_old;
   TH1F* h_Et_Gen_EE_new;
   TH1F* h_Et_Seed_EB_Denum_old;
   TH1F* h_Et_Seed_EB_Denum_new;
   TH1F* h_Et_Seed_EB_old;
   TH1F* h_Et_Seed_EB_new;
   TH1F* h_Et_Seed_EE_Denum_old;
   TH1F* h_Et_Seed_EE_Denum_new;
   TH1F* h_Et_Seed_EE_old;
   TH1F* h_Et_Seed_EE_new;
   TEfficiency* eff_SuperCluster_vs_EtaCalo;
   TEfficiency* eff_DeepSuperCluster_vs_EtaCalo;    
   TEfficiency* eff_SuperCluster_vs_EtaSeed;
   TEfficiency* eff_DeepSuperCluster_vs_EtaSeed;    
   TEfficiency* eff_SuperCluster_vs_EtCalo_EB;
   TEfficiency* eff_DeepSuperCluster_vs_EtCalo_EB;
   TEfficiency* eff_SuperCluster_vs_EtSeed_EB;
   TEfficiency* eff_DeepSuperCluster_vs_EtSeed_EB;
   TEfficiency* eff_SuperCluster_vs_EtCalo_EE;
   TEfficiency* eff_DeepSuperCluster_vs_EtCalo_EE;
   TEfficiency* eff_SuperCluster_vs_EtSeed_EE;
   TEfficiency* eff_DeepSuperCluster_vs_EtSeed_EE;

    TH1F* h_calo_pfClusters_eta_new;
    TH1F* h_calo_pfClusters_eta_old;
    TH1F* h_sc_pfClusters_eta_new;
    TH1F* h_sc_pfClusters_eta_old;
    TH1F* h_sc_pfClusters_caloInSC_new;
    TH1F* h_sc_pfClusters_caloInSC_old;

    TEfficiency* eff_Reco_SC_calo;
    TEfficiency* eff_Reco_DeepSC_calo;
    TEfficiency* eff_Reco_SC_gen;
    TEfficiency* eff_Reco_DeepSC_gen;
    TEfficiency* eff_Reco_SC_calo_dr;
    TEfficiency* eff_Reco_DeepSC_calo_dr;
    TEfficiency* eff_Reco_SC_gen_dr;
    TEfficiency* eff_Reco_DeepSC_gen_dr;
    TH1F* h_Eta_SC_old;
    TH1F* h_Eta_SC_new;
    TH1F* h_Eta_SCGen_old;
    TH1F* h_Eta_SCGen_new;
    TH1F* h_Eta_SC_old_dr;
    TH1F* h_Eta_SC_new_dr;
    TH1F* h_Eta_SCGen_old_dr;
    TH1F* h_Eta_SCGen_new_dr;
    TH1F* h_Eta_Calo_Denum_2;
    TH2F* h_caloPT_vs_Eta_new;
    TH2F* h_caloPT_vs_Eta_old;

//DEFINE PROFILES
   TProfile* prof_EoEtrue_vs_Eta_Calo_old;
   TProfile* prof_EoEtrue_vs_Eta_Calo_new;
   TProfile* prof_EoEtrue_vs_Eta_Seed_old;
   TProfile* prof_EoEtrue_vs_Eta_Seed_new;
   TProfile* prof_EoEtrue_vs_Et_Calo_EB_old;
   TProfile* prof_EoEtrue_vs_Et_Calo_EB_new;
   TProfile* prof_EoEtrue_vs_Et_Calo_EE_old;
   TProfile* prof_EoEtrue_vs_Et_Calo_EE_new;
   TProfile* prof_EoEtrue_vs_Et_Seed_EB_old;
   TProfile* prof_EoEtrue_vs_Et_Seed_EB_new;
   TProfile* prof_EoEtrue_vs_Et_Seed_EE_old;
   TProfile* prof_EoEtrue_vs_Et_Seed_EE_new;
   TProfile* prof_EoEtrue_vs_Energy_Calo_EB_old;
   TProfile* prof_EoEtrue_vs_Energy_Calo_EB_new;
   TProfile* prof_EoEtrue_vs_Energy_Calo_EE_old;
   TProfile* prof_EoEtrue_vs_Energy_Calo_EE_new;
   TProfile* prof_EoEtrue_vs_nVtx_EB_old;
   TProfile* prof_EoEtrue_vs_nVtx_EB_new;
   TProfile* prof_EoEtrue_vs_nVtx_EE_old;
   TProfile* prof_EoEtrue_vs_nVtx_EE_new;
   TProfile* prof_EoEtrue_vs_Rho_EB_old;
   TProfile* prof_EoEtrue_vs_Rho_EB_new;
   TProfile* prof_EoEtrue_vs_Rho_EE_old;
   TProfile* prof_EoEtrue_vs_Rho_EE_new;
   TProfile* prof_EoEgen_vs_Eta_Gen_old;
   TProfile* prof_EoEgen_vs_Eta_Gen_new;
   TProfile* prof_EoEgen_vs_Eta_Seed_old;
   TProfile* prof_EoEgen_vs_Eta_Seed_new;
   TProfile* prof_EoEgen_vs_Et_Gen_EB_old;
   TProfile* prof_EoEgen_vs_Et_Gen_EB_new;
   TProfile* prof_EoEgen_vs_Et_Gen_EE_old;
   TProfile* prof_EoEgen_vs_Et_Gen_EE_new;
   TProfile* prof_EoEgen_vs_Et_Seed_EB_old;
   TProfile* prof_EoEgen_vs_Et_Seed_EB_new;
   TProfile* prof_EoEgen_vs_Et_Seed_EE_old;
   TProfile* prof_EoEgen_vs_Et_Seed_EE_new;
   TProfile* prof_EoEgen_vs_Energy_Gen_EB_old;
   TProfile* prof_EoEgen_vs_Energy_Gen_EB_new;
   TProfile* prof_EoEgen_vs_Energy_Gen_EE_old;
   TProfile* prof_EoEgen_vs_Energy_Gen_EE_new;
   TProfile* prof_EoEgen_vs_nVtx_EB_old;
   TProfile* prof_EoEgen_vs_nVtx_EB_new;
   TProfile* prof_EoEgen_vs_nVtx_EE_old;
   TProfile* prof_EoEgen_vs_nVtx_EE_new;
   TProfile* prof_EoEgen_vs_Rho_EB_old;
   TProfile* prof_EoEgen_vs_Rho_EB_new;
   TProfile* prof_EoEgen_vs_Rho_EE_old;
   TProfile* prof_EoEgen_vs_Rho_EE_new;

//DEFINE HISTOGRAM VECTORS 
    std::vector<TH1F*> EoEtrue_vs_EtaWidth_old;
    std::vector<TH1F*> EoEtrue_vs_EtaWidth_new;
    std::vector<TH1F*> EoEtrue_vs_PhiWidth_old;
    std::vector<TH1F*> EoEtrue_vs_PhiWidth_new;

    std::vector<TH1F*> EoEtrue_vs_EtaWidth_EB_old;
    std::vector<TH1F*> EoEtrue_vs_EtaWidth_EB_new;
    std::vector<TH1F*> EoEtrue_vs_PhiWidth_EB_old;
    std::vector<TH1F*> EoEtrue_vs_PhiWidth_EB_new;

    std::vector<TH1F*> EoEtrue_vs_EtaWidth_EE_old;
    std::vector<TH1F*> EoEtrue_vs_EtaWidth_EE_new;
    std::vector<TH1F*> EoEtrue_vs_PhiWidth_EE_old;
    std::vector<TH1F*> EoEtrue_vs_PhiWidth_EE_new;

   std::vector<TH1F*> EoEtrue_vs_Eta_Calo_old;
   std::vector<TH1F*> EoEtrue_vs_Eta_Calo_new;
   std::vector<TH1F*> EoEtrue_vs_Eta_Seed_old;
   std::vector<TH1F*> EoEtrue_vs_Eta_Seed_new;
   std::vector<TH1F*> EoEtrue_vs_Et_Calo_EB_old;
   std::vector<TH1F*> EoEtrue_vs_Et_Calo_EB_new;
   std::vector<TH1F*> EoEtrue_vs_Et_Seed_EB_old;
   std::vector<TH1F*> EoEtrue_vs_Et_Seed_EB_new;
   std::vector<TH1F*> EoEtrue_vs_Energy_Calo_EB_old;
   std::vector<TH1F*> EoEtrue_vs_Energy_Calo_EB_new;
   std::vector<TH1F*> EoEtrue_vs_nVtx_EB_old;
   std::vector<TH1F*> EoEtrue_vs_nVtx_EB_new;
   std::vector<TH1F*> EoEtrue_vs_Rho_EB_old;
   std::vector<TH1F*> EoEtrue_vs_Rho_EB_new;
   std::vector<TH1F*> EoEtrue_vs_Et_Calo_EE_old;
   std::vector<TH1F*> EoEtrue_vs_Et_Calo_EE_new;
   std::vector<TH1F*> EoEtrue_vs_Et_Seed_EE_old;
   std::vector<TH1F*> EoEtrue_vs_Et_Seed_EE_new;
   std::vector<TH1F*> EoEtrue_vs_Energy_Calo_EE_old;
   std::vector<TH1F*> EoEtrue_vs_Energy_Calo_EE_new;
   std::vector<TH1F*> EoEtrue_vs_nVtx_EE_old;
   std::vector<TH1F*> EoEtrue_vs_nVtx_EE_new;
   std::vector<TH1F*> EoEtrue_vs_Rho_EE_old;
   std::vector<TH1F*> EoEtrue_vs_Rho_EE_new;
   std::vector<TH1F*> EoEgen_vs_Eta_Gen_old;
   std::vector<TH1F*> EoEgen_vs_Eta_Gen_new;
   std::vector<TH1F*> EoEgen_vs_Eta_Seed_old;
   std::vector<TH1F*> EoEgen_vs_Eta_Seed_new;
   std::vector<TH1F*> EoEgen_vs_Et_Gen_EB_old;
   std::vector<TH1F*> EoEgen_vs_Et_Gen_EB_new;
   std::vector<TH1F*> EoEgen_vs_Et_Seed_EB_old;
   std::vector<TH1F*> EoEgen_vs_Et_Seed_EB_new;
   std::vector<TH1F*> EoEgen_vs_Energy_Gen_EB_old;
   std::vector<TH1F*> EoEgen_vs_Energy_Gen_EB_new;
   std::vector<TH1F*> EoEgen_vs_nVtx_EB_old;
   std::vector<TH1F*> EoEgen_vs_nVtx_EB_new;
   std::vector<TH1F*> EoEgen_vs_Rho_EB_old;
   std::vector<TH1F*> EoEgen_vs_Rho_EB_new;
   std::vector<TH1F*> EoEgen_vs_Et_Gen_EE_old;
   std::vector<TH1F*> EoEgen_vs_Et_Gen_EE_new;
   std::vector<TH1F*> EoEgen_vs_Et_Seed_EE_old;
   std::vector<TH1F*> EoEgen_vs_Et_Seed_EE_new;
   std::vector<TH1F*> EoEgen_vs_Energy_Gen_EE_old;
   std::vector<TH1F*> EoEgen_vs_Energy_Gen_EE_new;
   std::vector<TH1F*> EoEgen_vs_nVtx_EE_old;
   std::vector<TH1F*> EoEgen_vs_nVtx_EE_new;
   std::vector<TH1F*> EoEgen_vs_Rho_EE_old;
   std::vector<TH1F*> EoEgen_vs_Rho_EE_new;
   std::vector<std::vector<TH1F*>> EoEtrue_vs_seedEt_seedEta_old;
   std::vector<std::vector<TH1F*>> EoEtrue_vs_seedEt_seedEta_new;

//DEFINE RESOLUTION MAPS 
   TH2F* h2_EoEtrue_Mean_old;
   TH2F* h2_EoEtrue_Mean_new;
   TH2F* h2_EoEtrue_Mean_Effective_old;
   TH2F* h2_EoEtrue_Mean_Effective_new;
   TH2F* h2_EoEtrue_Resolution_old;
   TH2F* h2_EoEtrue_Resolution_new;
   TH2F* h2_EoEtrue_Resolution_Effective_old;
   TH2F* h2_EoEtrue_Resolution_Effective_new;



//

class DoubleCBPdf;
class CruijffPdf;

class CruijffPdf : public RooAbsPdf 
{
  public:
    explicit CruijffPdf(const char *name, const char *title, RooAbsReal& _m,
		RooAbsReal& _m0, 
		RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
		RooAbsReal& _alphaL, RooAbsReal& _alphaR);

  CruijffPdf(const CruijffPdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { 
    return new CruijffPdf(*this,newname); } 
  virtual ~CruijffPdf() { }

  protected:
    RooRealProxy m;
    RooRealProxy m0;
    RooRealProxy sigmaL;
    RooRealProxy sigmaR;
    RooRealProxy alphaL;
    RooRealProxy alphaR;

  double evaluate() const;

};

CruijffPdf::CruijffPdf(const char *name, const char *title,
	     RooAbsReal& _m, RooAbsReal& _m0, 
	     RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
	     RooAbsReal& _alphaL, RooAbsReal& _alphaR)
  :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigmaL("sigmaL", "SigmaL", this, _sigmaL),
  sigmaR("sigmaR", "SigmaR", this, _sigmaR),
  alphaL("alphaL", "AlphaL", this, _alphaL),
  alphaR("alphaR", "AlphaR", this, _alphaR)
{
}

CruijffPdf::CruijffPdf(const CruijffPdf& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigmaL("sigmaL", this, other.sigmaL), sigmaR("sigmaR", this, other.sigmaR), 
  alphaL("alphaL", this, other.alphaL), alphaR("alphaR", this, other.alphaR)
{
}


double CruijffPdf::evaluate() const 
{
  double dx = (m-m0) ;
  double sigma = dx<0 ? sigmaL: sigmaR ;
  double alpha = dx<0 ? alphaL: alphaR ;
  double f = 2*sigma*sigma + alpha*dx*dx ;
  return exp(-dx*dx/f) ;
}

class DoubleCBPdf : public RooAbsPdf 
{
  public:
    explicit DoubleCBPdf(const char *name, const char *title,
                RooAbsReal& _x,
                RooAbsReal& _mean,
                RooAbsReal& _width,
                RooAbsReal& _alpha1,
                RooAbsReal& _n1,
                RooAbsReal& _alpha2,
                RooAbsReal& _n2);
  DoubleCBPdf(const DoubleCBPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { 
    return new DoubleCBPdf(*this,newname); }
  virtual ~DoubleCBPdf() { }
  
  protected:
    RooRealProxy x ;
    RooRealProxy mean;
    RooRealProxy width;
    RooRealProxy alpha1;
    RooRealProxy n1;
    RooRealProxy alpha2;
    RooRealProxy n2;
  
  double evaluate() const ;

};

DoubleCBPdf::DoubleCBPdf(const char *name, const char *title, 
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha1,
              RooAbsReal& _n1,
              RooAbsReal& _alpha2,
              RooAbsReal& _n2) 
  :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  mean("mean","mean",this,_mean),
  width("width","width",this,_width),
  alpha1("alpha1","alpha1",this,_alpha1),
  n1("n1","n1",this,_n1),
  alpha2("alpha2","alpha2",this,_alpha2),
  n2("n2","n2",this,_n2)
{ 
} 


DoubleCBPdf::DoubleCBPdf(const DoubleCBPdf& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  mean("mean",this,other.mean),
  width("width",this,other.width),
  alpha1("alpha1",this,other.alpha1),
  n1("n1",this,other.n1),
  alpha2("alpha2",this,other.alpha2),
  n2("n2",this,other.n2)
{ 
} 
 
double DoubleCBPdf::evaluate() const 
{ 
  double t = (x-mean)/width;
  double val = -99.;
  if(t>-alpha1 && t<alpha2){
     val = exp(-0.5*t*t);
  }else if(t<=-alpha1){
     double alpha1invn1 = alpha1/n1;
     val = exp(-0.5*alpha1*alpha1)*pow(1. - alpha1invn1*(alpha1+t), -n1);
  }else if(t>=alpha2){
     double alpha2invn2 = alpha2/n2;
     val = exp(-0.5*alpha2*alpha2)*pow(1. - alpha2invn2*(alpha2-t), -n2);        
  }
  //if(!std::isnormal(val)) {
  //   printf("bad val: x = %5f, t = %5f, mean = %5f, sigma = %5f, alpha1 = %5f, n1 = %5f, alpha2 = %5f, n2 = %5f\n",double(x), t, double(mean),double(width),double(alpha1),double(n1),double(alpha2), double(n2));
  //   printf("val = %5f\n",val);
  //}  
  return val;
} 


void ReadInfile(string SuperClusterRef, string SuperClusterVal, string SCRef_File, string SCVal_File){
    cout<<"Reading input file"<<endl;

    std::vector<float> etCuts {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
    std::vector<float> etaCuts {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.479,1.75,2.0,2.25,3.0};

    TFile *hist_infile_SCRef = TFile::Open(SCRef_File.c_str());;
    TFile *hist_infile_SCVal = TFile::Open(SCVal_File.c_str());


    string ref_1 = "", ref_2 = "", val_1 = "", val_2 = "";
    if(SuperClusterRef == "SingleSet"){
        ref_1 = "DeepSC";
        ref_2 = "new";
    }
    else if(SuperClusterRef == "Mustache"){
        ref_1 = "SC";
        ref_2 = "old";
    }

    if(SuperClusterVal == "SingleSet"){
        val_1 = "DeepSC";
        val_2 = "new";
    }
    else if(SuperClusterVal == "LocalSets"){
        val_1 = "DeepSC";
        val_2 = "new";
    }

    //Super Cluster Ref

    hist_infile_SCRef->GetObject(("h_Eta_"+ref_1).c_str(), h_Eta_old);     h_Eta_old->SetName("h_Eta_old");
    //neeed new
    hist_infile_SCRef ->GetObject("h_Eta_Calo_Denum", h_Eta_Calo_Denum_old);
    h_Eta_Calo_Denum_old->SetName("h_Eta_Calo_Denum_old");

    hist_infile_SCRef ->GetObject("h_Et_Calo_EB_Denum", h_Et_Calo_EB_Denum_old);
    h_Et_Calo_EB_Denum_old->SetName("h_Et_Calo_EB_Denum_old");
    hist_infile_SCRef ->GetObject("h_Et_Calo_EE_Denum", h_Et_Calo_EE_Denum_old);
    h_Et_Calo_EE_Denum_old->SetName("h_Et_Calo_EE_Denum_old");


    //single Set Params -> delimited by old

    hist_infile_SCRef->GetObject(("h_Eta_"+ref_1).c_str(), h_Eta_old);     h_Eta_old->SetName("h_Eta_old");
        //Eta
    hist_infile_SCRef->GetObject(("h_Eta_"+ref_1).c_str(), h_Eta_old); h_Eta_old->SetName("h_Eta_old");
    hist_infile_SCRef->GetObject(("h_Eta_Calo_"+ref_1).c_str(), h_Eta_Calo_old); h_Eta_Calo_old->SetName("h_Eta_Calo_old");
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_Eta_Calo_"+ref_1).c_str(), prof_EoEtrue_vs_Eta_Calo_old); prof_EoEtrue_vs_Eta_Calo_old->SetName("prof_EoEtrue_vs_Eta_Calo_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_Eta_Gen_"+ref_1).c_str(), prof_EoEgen_vs_Eta_Gen_old); prof_EoEgen_vs_Eta_Gen_old->SetName("prof_EoEgen_vs_Eta_Gen_old");
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_Eta_Seed_"+ref_1).c_str(), prof_EoEtrue_vs_Eta_Seed_old); prof_EoEtrue_vs_Eta_Seed_old->SetName("prof_EoEtrue_vs_Eta_Seed_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_Eta_Seed_"+ref_1).c_str(), prof_EoEgen_vs_Eta_Seed_old); prof_EoEgen_vs_Eta_Seed_old->SetName("prof_EoEgen_vs_Eta_Seed_old");
    hist_infile_SCRef->GetObject(("h_Eta_seedMatched_"+ref_1).c_str(), h_Eta_seedMatched_old); h_Eta_seedMatched_old->SetName("h_Eta_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_Eta_caloMatched_"+ref_1).c_str(), h_Eta_caloMatched_old); h_Eta_caloMatched_old->SetName("h_Eta_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_Eta_caloUnmatched_"+ref_1).c_str(), h_Eta_caloUnmatched_old); h_Eta_caloUnmatched_old->SetName("h_Eta_caloUnmatched_old");
        //nPFClusters
    hist_infile_SCRef->GetObject(("h_nPFClusters_"+ref_1).c_str(), h_nPFClusters_old); h_nPFClusters_old->SetName("h_nPFClusters_old");
    hist_infile_SCRef->GetObject(("h_nPFClusters_seedMatched_"+ref_1).c_str(), h_nPFClusters_seedMatched_old); h_nPFClusters_seedMatched_old->SetName("h_nPFClusters_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_nPFClusters_caloMatched_"+ref_1).c_str(), h_nPFClusters_caloMatched_old); h_nPFClusters_caloMatched_old->SetName("h_nPFClusters_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_nPFClusters_caloUnmatched_"+ref_1).c_str(), h_nPFClusters_caloUnmatched_old); h_nPFClusters_caloUnmatched_old->SetName("h_nPFClusters_caloUnmatched_old");
    
    hist_infile_SCRef->GetObject(("h_nPFClusters_EB_"+ref_1).c_str(), h_nPFClusters_EB_old); h_nPFClusters_EB_old->SetName("h_nPFClusters_EB_old");
    hist_infile_SCRef->GetObject(("h_nPFClusters_EB_seedMatched_"+ref_1).c_str(), h_nPFClusters_EB_seedMatched_old); h_nPFClusters_EB_seedMatched_old->SetName("h_nPFClusters_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_nPFClusters_EB_caloMatched_"+ref_1).c_str(), h_nPFClusters_EB_caloMatched_old); h_nPFClusters_EB_caloMatched_old->SetName("h_nPFClusters_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_nPFClusters_EB_caloUnmatched_"+ref_1).c_str(), h_nPFClusters_EB_caloUnmatched_old); h_nPFClusters_EB_caloUnmatched_old->SetName("h_nPFClusters_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_nPFClusters_EE_"+ref_1).c_str(), h_nPFClusters_EE_old); h_nPFClusters_EE_old->SetName("h_nPFClusters_EE_old");
    hist_infile_SCRef->GetObject(("h_nPFClusters_EE_seedMatched_"+ref_1).c_str(), h_nPFClusters_EE_seedMatched_old); h_nPFClusters_EE_seedMatched_old->SetName("h_nPFClusters_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_nPFClusters_EE_caloMatched_"+ref_1).c_str(), h_nPFClusters_EE_caloMatched_old); h_nPFClusters_EE_caloMatched_old->SetName("h_nPFClusters_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_nPFClusters_EE_caloUnmatched_"+ref_1).c_str(), h_nPFClusters_EE_caloUnmatched_old); h_nPFClusters_EE_caloUnmatched_old->SetName("h_nPFClusters_EE_caloUnmatched_old");
         
        //nVtx
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_nVtx_EB_"+ref_1).c_str(), prof_EoEtrue_vs_nVtx_EB_old); prof_EoEtrue_vs_nVtx_EB_old->SetName("prof_EoEtrue_vs_nVtx_EB_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_nVtx_EB_"+ref_1).c_str(), prof_EoEgen_vs_nVtx_EB_old); prof_EoEgen_vs_nVtx_EB_old->SetName("prof_EoEgen_vs_nVtx_EB_old");
             
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_nVtx_EE_"+ref_1).c_str(), prof_EoEtrue_vs_nVtx_EE_old); prof_EoEtrue_vs_nVtx_EE_old->SetName("prof_EoEtrue_vs_nVtx_EE_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_nVtx_EE_"+ref_1).c_str(), prof_EoEgen_vs_nVtx_EE_old); prof_EoEgen_vs_nVtx_EE_old->SetName("prof_EoEgen_vs_nVtx_EE_old");
             
        //Rho
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_Rho_EB_"+ref_1).c_str(), prof_EoEtrue_vs_Rho_EB_old); prof_EoEtrue_vs_Rho_EB_old->SetName("prof_EoEtrue_vs_Rho_EB_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_Rho_EB_"+ref_1).c_str(), prof_EoEgen_vs_Rho_EB_old); prof_EoEgen_vs_Rho_EB_old->SetName("prof_EoEgen_vs_Rho_EB_old");
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_Rho_EE_"+ref_1).c_str(), prof_EoEtrue_vs_Rho_EE_old); prof_EoEtrue_vs_Rho_EE_old->SetName("prof_EoEtrue_vs_Rho_EE_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_Rho_EE_"+ref_1).c_str(), prof_EoEgen_vs_Rho_EE_old); prof_EoEgen_vs_Rho_EE_old->SetName("prof_EoEgen_vs_Rho_EE_old");

        //Energy
    hist_infile_SCRef->GetObject(("h_Energy_EB_"+ref_1).c_str(), h_Energy_EB_old); h_Energy_EB_old->SetName("h_Energy_EB_old");
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_Energy_Calo_EB_"+ref_1).c_str(), prof_EoEtrue_vs_Energy_Calo_EB_old); prof_EoEtrue_vs_Energy_Calo_EB_old->SetName("prof_EoEtrue_vs_Energy_Calo_EB_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_Energy_Gen_EB_"+ref_1).c_str(), prof_EoEgen_vs_Energy_Gen_EB_old); prof_EoEgen_vs_Energy_Gen_EB_old->SetName("prof_EoEgen_vs_Energy_Gen_EB_old");
    hist_infile_SCRef->GetObject(("h_Energy_EB_seedMatched_"+ref_1).c_str(), h_Energy_EB_seedMatched_old); h_Energy_EB_seedMatched_old->SetName("h_Energy_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_Energy_EB_caloMatched_"+ref_1).c_str(), h_Energy_EB_caloMatched_old); h_Energy_EB_caloMatched_old->SetName("h_Energy_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_Energy_EB_caloUnmatched_"+ref_1).c_str(), h_Energy_EB_caloUnmatched_old); h_Energy_EB_caloUnmatched_old->SetName("h_Energy_EB_caloUnmatched_old");

    hist_infile_SCRef->GetObject(("h_Energy_EE_"+ref_1).c_str(), h_Energy_EE_old); h_Energy_EE_old->SetName("h_Energy_EE_old");
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_Energy_Calo_EE_"+ref_1).c_str(), prof_EoEtrue_vs_Energy_Calo_EE_old); prof_EoEtrue_vs_Energy_Calo_EE_old->SetName("prof_EoEtrue_vs_Energy_Calo_EE_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_Energy_Gen_EE_"+ref_1).c_str(), prof_EoEgen_vs_Energy_Gen_EE_old); prof_EoEgen_vs_Energy_Gen_EE_old->SetName("prof_EoEgen_vs_Energy_Gen_EE_old");
    hist_infile_SCRef->GetObject(("h_Energy_EE_seedMatched_"+ref_1).c_str(), h_Energy_EE_seedMatched_old); h_Energy_EE_seedMatched_old->SetName("h_Energy_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_Energy_EE_caloMatched_"+ref_1).c_str(), h_Energy_EE_caloMatched_old); h_Energy_EE_caloMatched_old->SetName("h_Energy_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_Energy_EE_caloUnmatched_"+ref_1).c_str(), h_Energy_EE_caloUnmatched_old); h_Energy_EE_caloUnmatched_old->SetName("h_Energy_EE_caloUnmatched_old");
                        
        //EoEtrue
    hist_infile_SCRef->GetObject(("h_EoEtrue_EB_"+ref_1).c_str(), h_EoEtrue_EB_old); h_EoEtrue_EB_old->SetName("h_EoEtrue_EB_old");
    hist_infile_SCRef->GetObject(("h_EoEtrue_EB_seedMatched_"+ref_1).c_str(), h_EoEtrue_EB_seedMatched_old); h_EoEtrue_EB_seedMatched_old->SetName("h_EoEtrue_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_EoEtrue_EE_"+ref_1).c_str(), h_EoEtrue_EE_old); h_EoEtrue_EE_old->SetName("h_EoEtrue_EE_old");
    hist_infile_SCRef->GetObject(("h_EoEtrue_EE_seedMatched_"+ref_1).c_str(), h_EoEtrue_EE_seedMatched_old); h_EoEtrue_EE_seedMatched_old->SetName("h_EoEtrue_EE_seedMatched_old");
             
        //EoEgen
    hist_infile_SCRef->GetObject(("h_EoEgen_EB_"+ref_1).c_str(), h_EoEgen_EB_old); h_EoEgen_EB_old->SetName("h_EoEgen_EB_old");
    hist_infile_SCRef->GetObject(("h_EoEgen_EB_seedMatched_"+ref_1).c_str(), h_EoEgen_EB_seedMatched_old); h_EoEgen_EB_seedMatched_old->SetName("h_EoEgen_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_EoEgen_EE_"+ref_1).c_str(), h_EoEgen_EE_old); h_EoEgen_EE_old->SetName("h_EoEgen_EE_old");
    hist_infile_SCRef->GetObject(("h_EoEgen_EE_seedMatched_"+ref_1).c_str(), h_EoEgen_EE_seedMatched_old); h_EoEgen_EE_seedMatched_old->SetName("h_EoEgen_EE_seedMatched_old");
             
        //Et
    hist_infile_SCRef->GetObject(("h_Et_EB_"+ref_1).c_str(), h_Et_EB_old); h_Et_EB_old->SetName("h_Et_EB_old");
    hist_infile_SCRef->GetObject(("h_Et_Calo_EB_"+ref_1).c_str(), h_Et_Calo_EB_old); h_Et_Calo_EB_old->SetName("h_Et_Calo_EB_old");
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_Et_Calo_EB_"+ref_1).c_str(), prof_EoEtrue_vs_Et_Calo_EB_old); prof_EoEtrue_vs_Et_Calo_EB_old->SetName("prof_EoEtrue_vs_Et_Calo_EB_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_Et_Gen_EB_"+ref_1).c_str(), prof_EoEgen_vs_Et_Gen_EB_old); prof_EoEgen_vs_Et_Gen_EB_old->SetName("prof_EoEgen_vs_Et_Gen_EB_old");
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_Et_Seed_EB_"+ref_1).c_str(), prof_EoEtrue_vs_Et_Seed_EB_old); prof_EoEtrue_vs_Et_Seed_EB_old->SetName("prof_EoEtrue_vs_Et_Seed_EB_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_Et_Seed_EB_"+ref_1).c_str(), prof_EoEgen_vs_Et_Seed_EB_old); prof_EoEgen_vs_Et_Seed_EB_old->SetName("prof_EoEgen_vs_Et_Seed_EB_old");
    hist_infile_SCRef->GetObject(("h_Et_EB_seedMatched_"+ref_1).c_str(), h_Et_EB_seedMatched_old); h_Et_EB_seedMatched_old->SetName("h_Et_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_Et_EB_caloMatched_"+ref_1).c_str(), h_Et_EB_caloMatched_old); h_Et_EB_caloMatched_old->SetName("h_Et_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_Et_EB_caloUnmatched_"+ref_1).c_str(), h_Et_EB_caloUnmatched_old); h_Et_EB_caloUnmatched_old->SetName("h_Et_EB_caloUnmatched_old");

    hist_infile_SCRef->GetObject(("h_Et_EE_"+ref_1).c_str(), h_Et_EE_old); h_Et_EE_old->SetName("h_Et_EE_old");
    hist_infile_SCRef->GetObject(("h_Et_Calo_EE_"+ref_1).c_str(), h_Et_Calo_EE_old); h_Et_Calo_EE_old->SetName("h_Et_Calo_EE_old");
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_Et_Calo_EE_"+ref_1).c_str(), prof_EoEtrue_vs_Et_Calo_EE_old); prof_EoEtrue_vs_Et_Calo_EE_old->SetName("prof_EoEtrue_vs_Et_Calo_EE_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_Et_Gen_EE_"+ref_1).c_str(), prof_EoEgen_vs_Et_Gen_EE_old); prof_EoEgen_vs_Et_Gen_EE_old->SetName("prof_EoEgen_vs_Et_Gen_EE_old");
    hist_infile_SCRef->GetObject(("prof_EoEtrue_vs_Et_Seed_EE_"+ref_1).c_str(), prof_EoEtrue_vs_Et_Seed_EE_old); prof_EoEtrue_vs_Et_Seed_EE_old->SetName("prof_EoEtrue_vs_Et_Seed_EE_old");
    hist_infile_SCRef->GetObject(("prof_EoEgen_vs_Et_Seed_EE_"+ref_1).c_str(), prof_EoEgen_vs_Et_Seed_EE_old); prof_EoEgen_vs_Et_Seed_EE_old->SetName("prof_EoEgen_vs_Et_Seed_EE_old");
    hist_infile_SCRef->GetObject(("h_Et_EE_seedMatched_"+ref_1).c_str(), h_Et_EE_seedMatched_old); h_Et_EE_seedMatched_old->SetName("h_Et_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_Et_EE_caloMatched_"+ref_1).c_str(), h_Et_EE_caloMatched_old); h_Et_EE_caloMatched_old->SetName("h_Et_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_Et_EE_caloUnmatched_"+ref_1).c_str(), h_Et_EE_caloUnmatched_old); h_Et_EE_caloUnmatched_old->SetName("h_Et_EE_caloUnmatched_old");


        //Phi
    hist_infile_SCRef->GetObject(("h_Phi_EB_"+ref_1).c_str(), h_Phi_EB_old); h_Phi_EB_old->SetName("h_Phi_EB_old");
    hist_infile_SCRef->GetObject(("h_Phi_EB_seedMatched_"+ref_1).c_str(), h_Phi_EB_seedMatched_old); h_Phi_EB_seedMatched_old->SetName("h_Phi_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_Phi_EB_caloMatched_"+ref_1).c_str(), h_Phi_EB_caloMatched_old); h_Phi_EB_caloMatched_old->SetName("h_Phi_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_Phi_EB_caloUnmatched_"+ref_1).c_str(), h_Phi_EB_caloUnmatched_old); h_Phi_EB_caloUnmatched_old->SetName("h_Phi_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_Phi_EE_"+ref_1).c_str(), h_Phi_EE_old); h_Phi_EE_old->SetName("h_Phi_EE_old");
    hist_infile_SCRef->GetObject(("h_Phi_EE_seedMatched_"+ref_1).c_str(), h_Phi_EE_seedMatched_old); h_Phi_EE_seedMatched_old->SetName("h_Phi_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_Phi_EE_caloMatched_"+ref_1).c_str(), h_Phi_EE_caloMatched_old); h_Phi_EE_caloMatched_old->SetName("h_Phi_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_Phi_EE_caloUnmatched_"+ref_1).c_str(), h_Phi_EE_caloUnmatched_old); h_Phi_EE_caloUnmatched_old->SetName("h_Phi_EE_caloUnmatched_old");
             
        //EtaWidth
    hist_infile_SCRef->GetObject(("h_EtaWidth_EB_"+ref_1).c_str(), h_EtaWidth_EB_old); h_EtaWidth_EB_old->SetName("h_EtaWidth_EB_old");
    hist_infile_SCRef->GetObject(("h_EtaWidth_EB_seedMatched_"+ref_1).c_str(), h_EtaWidth_EB_seedMatched_old); h_EtaWidth_EB_seedMatched_old->SetName("h_EtaWidth_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_EtaWidth_EB_caloMatched_"+ref_1).c_str(), h_EtaWidth_EB_caloMatched_old); h_EtaWidth_EB_caloMatched_old->SetName("h_EtaWidth_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_EtaWidth_EB_caloUnmatched_"+ref_1).c_str(), h_EtaWidth_EB_caloUnmatched_old); h_EtaWidth_EB_caloUnmatched_old->SetName("h_EtaWidth_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_EtaWidth_EE_"+ref_1).c_str(), h_EtaWidth_EE_old); h_EtaWidth_EE_old->SetName("h_EtaWidth_EE_old");
    hist_infile_SCRef->GetObject(("h_EtaWidth_EE_seedMatched_"+ref_1).c_str(), h_EtaWidth_EE_seedMatched_old); h_EtaWidth_EE_seedMatched_old->SetName("h_EtaWidth_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_EtaWidth_EE_caloMatched_"+ref_1).c_str(), h_EtaWidth_EE_caloMatched_old); h_EtaWidth_EE_caloMatched_old->SetName("h_EtaWidth_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_EtaWidth_EE_caloUnmatched_"+ref_1).c_str(), h_EtaWidth_EE_caloUnmatched_old); h_EtaWidth_EE_caloUnmatched_old->SetName("h_EtaWidth_EE_caloUnmatched_old");
             
        //PhiWidth
    hist_infile_SCRef->GetObject(("h_PhiWidth_EB_"+ref_1).c_str(), h_PhiWidth_EB_old); h_PhiWidth_EB_old->SetName("h_PhiWidth_EB_old");
    hist_infile_SCRef->GetObject(("h_PhiWidth_EB_seedMatched_"+ref_1).c_str(), h_PhiWidth_EB_seedMatched_old); h_PhiWidth_EB_seedMatched_old->SetName("h_PhiWidth_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_PhiWidth_EB_caloMatched_"+ref_1).c_str(), h_PhiWidth_EB_caloMatched_old); h_PhiWidth_EB_caloMatched_old->SetName("h_PhiWidth_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_PhiWidth_EB_caloUnmatched_"+ref_1).c_str(), h_PhiWidth_EB_caloUnmatched_old); h_PhiWidth_EB_caloUnmatched_old->SetName("h_PhiWidth_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_PhiWidth_EE_"+ref_1).c_str(), h_PhiWidth_EE_old); h_PhiWidth_EE_old->SetName("h_PhiWidth_EE_old");
    hist_infile_SCRef->GetObject(("h_PhiWidth_EE_seedMatched_"+ref_1).c_str(), h_PhiWidth_EE_seedMatched_old); h_PhiWidth_EE_seedMatched_old->SetName("h_PhiWidth_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_PhiWidth_EE_caloMatched_"+ref_1).c_str(), h_PhiWidth_EE_caloMatched_old); h_PhiWidth_EE_caloMatched_old->SetName("h_PhiWidth_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_PhiWidth_EE_caloUnmatched_"+ref_1).c_str(), h_PhiWidth_EE_caloUnmatched_old); h_PhiWidth_EE_caloUnmatched_old->SetName("h_PhiWidth_EE_caloUnmatched_old");
             
        //full5x5_R9
    hist_infile_SCRef->GetObject(("h_full5x5_R9_EB_"+ref_1).c_str(), h_full5x5_R9_EB_old); h_full5x5_R9_EB_old->SetName("h_full5x5_R9_EB_old");
    hist_infile_SCRef->GetObject(("h_full5x5_R9_EB_seedMatched_"+ref_1).c_str(), h_full5x5_R9_EB_seedMatched_old); h_full5x5_R9_EB_seedMatched_old->SetName("h_full5x5_R9_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_R9_EB_caloMatched_"+ref_1).c_str(), h_full5x5_R9_EB_caloMatched_old); h_full5x5_R9_EB_caloMatched_old->SetName("h_full5x5_R9_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_R9_EB_caloUnmatched_"+ref_1).c_str(), h_full5x5_R9_EB_caloUnmatched_old); h_full5x5_R9_EB_caloUnmatched_old->SetName("h_full5x5_R9_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_full5x5_R9_EE_"+ref_1).c_str(), h_full5x5_R9_EE_old); h_full5x5_R9_EE_old->SetName("h_full5x5_R9_EE_old");
    hist_infile_SCRef->GetObject(("h_full5x5_R9_EE_seedMatched_"+ref_1).c_str(), h_full5x5_R9_EE_seedMatched_old); h_full5x5_R9_EE_seedMatched_old->SetName("h_full5x5_R9_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_R9_EE_caloMatched_"+ref_1).c_str(), h_full5x5_R9_EE_caloMatched_old); h_full5x5_R9_EE_caloMatched_old->SetName("h_full5x5_R9_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_R9_EE_caloUnmatched_"+ref_1).c_str(), h_full5x5_R9_EE_caloUnmatched_old); h_full5x5_R9_EE_caloUnmatched_old->SetName("h_full5x5_R9_EE_caloUnmatched_old");
             
        //R9
    hist_infile_SCRef->GetObject(("h_R9_EB_"+ref_1).c_str(), h_R9_EB_old); h_R9_EB_old->SetName("h_R9_EB_old");
    hist_infile_SCRef->GetObject(("h_R9_EB_seedMatched_"+ref_1).c_str(), h_R9_EB_seedMatched_old); h_R9_EB_seedMatched_old->SetName("h_R9_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_R9_EB_caloMatched_"+ref_1).c_str(), h_R9_EB_caloMatched_old); h_R9_EB_caloMatched_old->SetName("h_R9_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_R9_EB_caloUnmatched_"+ref_1).c_str(), h_R9_EB_caloUnmatched_old); h_R9_EB_caloUnmatched_old->SetName("h_R9_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_R9_EE_"+ref_1).c_str(), h_R9_EE_old); h_R9_EE_old->SetName("h_R9_EE_old");
    hist_infile_SCRef->GetObject(("h_R9_EE_seedMatched_"+ref_1).c_str(), h_R9_EE_seedMatched_old); h_R9_EE_seedMatched_old->SetName("h_R9_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_R9_EE_caloMatched_"+ref_1).c_str(), h_R9_EE_caloMatched_old); h_R9_EE_caloMatched_old->SetName("h_R9_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_R9_EE_caloUnmatched_"+ref_1).c_str(), h_R9_EE_caloUnmatched_old); h_R9_EE_caloUnmatched_old->SetName("h_R9_EE_caloUnmatched_old");
             
        //full5x5_SigmaIetaIeta
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIeta_EB_"+ref_1).c_str(), h_full5x5_sigmaIetaIeta_EB_old); h_full5x5_sigmaIetaIeta_EB_old->SetName("h_full5x5_SigmaIetaIeta_EB_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIeta_EB_seedMatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIeta_EB_seedMatched_old); h_full5x5_sigmaIetaIeta_EB_seedMatched_old->SetName("h_full5x5_SigmaIetaIeta_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIeta_EB_caloMatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIeta_EB_caloMatched_old); h_full5x5_sigmaIetaIeta_EB_caloMatched_old->SetName("h_full5x5_SigmaIetaIeta_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIeta_EB_caloUnmatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old); h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old->SetName("h_full5x5_SigmaIetaIeta_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIeta_EE_"+ref_1).c_str(), h_full5x5_sigmaIetaIeta_EE_old); h_full5x5_sigmaIetaIeta_EE_old->SetName("h_full5x5_SigmaIetaIeta_EE_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIeta_EE_seedMatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIeta_EE_seedMatched_old); h_full5x5_sigmaIetaIeta_EE_seedMatched_old->SetName("h_full5x5_SigmaIetaIeta_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIeta_EE_caloMatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIeta_EE_caloMatched_old); h_full5x5_sigmaIetaIeta_EE_caloMatched_old->SetName("h_full5x5_SigmaIetaIeta_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIeta_EE_caloUnmatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old); h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old->SetName("h_full5x5_SigmaIetaIeta_EE_caloUnmatched_old");
             
        //full5x5_SigmaIetaIphi
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIphi_EB_"+ref_1).c_str(), h_full5x5_sigmaIetaIphi_EB_old); h_full5x5_sigmaIetaIphi_EB_old->SetName("h_full5x5_SigmaIetaIphi_EB_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIphi_EB_seedMatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIphi_EB_seedMatched_old); h_full5x5_sigmaIetaIphi_EB_seedMatched_old->SetName("h_full5x5_SigmaIetaIphi_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIphi_EB_caloMatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIphi_EB_caloMatched_old); h_full5x5_sigmaIetaIphi_EB_caloMatched_old->SetName("h_full5x5_SigmaIetaIphi_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIphi_EB_caloUnmatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old); h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old->SetName("h_full5x5_SigmaIetaIphi_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIphi_EE_"+ref_1).c_str(), h_full5x5_sigmaIetaIphi_EE_old); h_full5x5_sigmaIetaIphi_EE_old->SetName("h_full5x5_SigmaIetaIphi_EE_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIphi_EE_seedMatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIphi_EE_seedMatched_old); h_full5x5_sigmaIetaIphi_EE_seedMatched_old->SetName("h_full5x5_SigmaIetaIphi_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIphi_EE_caloMatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIphi_EE_caloMatched_old); h_full5x5_sigmaIetaIphi_EE_caloMatched_old->SetName("h_full5x5_SigmaIetaIphi_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIetaIphi_EE_caloUnmatched_"+ref_1).c_str(), h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old); h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old->SetName("h_full5x5_SigmaIetaIphi_EE_caloUnmatched_old");
             
        //SigmaIetaIeta
    hist_infile_SCRef->GetObject(("h_SigmaIetaIeta_EB_"+ref_1).c_str(), h_sigmaIetaIeta_EB_old); h_sigmaIetaIeta_EB_old->SetName("h_SigmaIetaIeta_EB_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIeta_EB_seedMatched_"+ref_1).c_str(), h_sigmaIetaIeta_EB_seedMatched_old); h_sigmaIetaIeta_EB_seedMatched_old->SetName("h_SigmaIetaIeta_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIeta_EB_caloMatched_"+ref_1).c_str(), h_sigmaIetaIeta_EB_caloMatched_old); h_sigmaIetaIeta_EB_caloMatched_old->SetName("h_SigmaIetaIeta_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIeta_EB_caloUnmatched_"+ref_1).c_str(), h_sigmaIetaIeta_EB_caloUnmatched_old); h_sigmaIetaIeta_EB_caloUnmatched_old->SetName("h_SigmaIetaIeta_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_SigmaIetaIeta_EE_"+ref_1).c_str(), h_sigmaIetaIeta_EE_old); h_sigmaIetaIeta_EE_old->SetName("h_SigmaIetaIeta_EE_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIeta_EE_seedMatched_"+ref_1).c_str(), h_sigmaIetaIeta_EE_seedMatched_old); h_sigmaIetaIeta_EE_seedMatched_old->SetName("h_SigmaIetaIeta_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIeta_EE_caloMatched_"+ref_1).c_str(), h_sigmaIetaIeta_EE_caloMatched_old); h_sigmaIetaIeta_EE_caloMatched_old->SetName("h_SigmaIetaIeta_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIeta_EE_caloUnmatched_"+ref_1).c_str(), h_sigmaIetaIeta_EE_caloUnmatched_old); h_sigmaIetaIeta_EE_caloUnmatched_old->SetName("h_SigmaIetaIeta_EE_caloUnmatched_old");
             
        //SigmaIetaIphi
    hist_infile_SCRef->GetObject(("h_SigmaIetaIphi_EB_"+ref_1).c_str(), h_sigmaIetaIphi_EB_old); h_sigmaIetaIphi_EB_old->SetName("h_SigmaIetaIphi_EB_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIphi_EB_seedMatched_"+ref_1).c_str(), h_sigmaIetaIphi_EB_seedMatched_old); h_sigmaIetaIphi_EB_seedMatched_old->SetName("h_SigmaIetaIphi_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIphi_EB_caloMatched_"+ref_1).c_str(), h_sigmaIetaIphi_EB_caloMatched_old); h_sigmaIetaIphi_EB_caloMatched_old->SetName("h_SigmaIetaIphi_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIphi_EB_caloUnmatched_"+ref_1).c_str(), h_sigmaIetaIphi_EB_caloUnmatched_old); h_sigmaIetaIphi_EB_caloUnmatched_old->SetName("h_SigmaIetaIphi_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_SigmaIetaIphi_EE_"+ref_1).c_str(), h_sigmaIetaIphi_EE_old); h_sigmaIetaIphi_EE_old->SetName("h_SigmaIetaIphi_EE_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIphi_EE_seedMatched_"+ref_1).c_str(), h_sigmaIetaIphi_EE_seedMatched_old); h_sigmaIetaIphi_EE_seedMatched_old->SetName("h_SigmaIetaIphi_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIphi_EE_caloMatched_"+ref_1).c_str(), h_sigmaIetaIphi_EE_caloMatched_old); h_sigmaIetaIphi_EE_caloMatched_old->SetName("h_SigmaIetaIphi_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIetaIphi_EE_caloUnmatched_"+ref_1).c_str(), h_sigmaIetaIphi_EE_caloUnmatched_old); h_sigmaIetaIphi_EE_caloUnmatched_old->SetName("h_SigmaIetaIphi_EE_caloUnmatched_old");
             
        //full5x5_SigmaIphiIphi
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIphiIphi_EB_"+ref_1).c_str(), h_full5x5_sigmaIphiIphi_EB_old); h_full5x5_sigmaIphiIphi_EB_old->SetName("h_full5x5_SigmaIphiIphi_EB_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIphiIphi_EB_seedMatched_"+ref_1).c_str(), h_full5x5_sigmaIphiIphi_EB_seedMatched_old); h_full5x5_sigmaIphiIphi_EB_seedMatched_old->SetName("h_full5x5_SigmaIphiIphi_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIphiIphi_EB_caloMatched_"+ref_1).c_str(), h_full5x5_sigmaIphiIphi_EB_caloMatched_old); h_full5x5_sigmaIphiIphi_EB_caloMatched_old->SetName("h_full5x5_SigmaIphiIphi_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIphiIphi_EB_caloUnmatched_"+ref_1).c_str(), h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old); h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old->SetName("h_full5x5_SigmaIphiIphi_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIphiIphi_EE_"+ref_1).c_str(), h_full5x5_sigmaIphiIphi_EE_old); h_full5x5_sigmaIphiIphi_EE_old->SetName("h_full5x5_SigmaIphiIphi_EE_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIphiIphi_EE_seedMatched_"+ref_1).c_str(), h_full5x5_sigmaIphiIphi_EE_seedMatched_old); h_full5x5_sigmaIphiIphi_EE_seedMatched_old->SetName("h_full5x5_SigmaIphiIphi_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIphiIphi_EE_caloMatched_"+ref_1).c_str(), h_full5x5_sigmaIphiIphi_EE_caloMatched_old); h_full5x5_sigmaIphiIphi_EE_caloMatched_old->SetName("h_full5x5_SigmaIphiIphi_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_full5x5_SigmaIphiIphi_EE_caloUnmatched_"+ref_1).c_str(), h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old); h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old->SetName("h_full5x5_SigmaIphiIphi_EE_caloUnmatched_old");
             
        //SigmaIphiIphi
    hist_infile_SCRef->GetObject(("h_SigmaIphiIphi_EB_"+ref_1).c_str(), h_sigmaIphiIphi_EB_old); h_sigmaIphiIphi_EB_old->SetName("h_SigmaIphiIphi_EB_old");
    hist_infile_SCRef->GetObject(("h_SigmaIphiIphi_EB_seedMatched_"+ref_1).c_str(), h_sigmaIphiIphi_EB_seedMatched_old); h_sigmaIphiIphi_EB_seedMatched_old->SetName("h_SigmaIphiIphi_EB_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIphiIphi_EB_caloMatched_"+ref_1).c_str(), h_sigmaIphiIphi_EB_caloMatched_old); h_sigmaIphiIphi_EB_caloMatched_old->SetName("h_SigmaIphiIphi_EB_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIphiIphi_EB_caloUnmatched_"+ref_1).c_str(), h_sigmaIphiIphi_EB_caloUnmatched_old); h_sigmaIphiIphi_EB_caloUnmatched_old->SetName("h_SigmaIphiIphi_EB_caloUnmatched_old");
          
    hist_infile_SCRef->GetObject(("h_SigmaIphiIphi_EE_"+ref_1).c_str(), h_sigmaIphiIphi_EE_old); h_sigmaIphiIphi_EE_old->SetName("h_SigmaIphiIphi_EE_old");
    hist_infile_SCRef->GetObject(("h_SigmaIphiIphi_EE_seedMatched_"+ref_1).c_str(), h_sigmaIphiIphi_EE_seedMatched_old); h_sigmaIphiIphi_EE_seedMatched_old->SetName("h_SigmaIphiIphi_EE_seedMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIphiIphi_EE_caloMatched_"+ref_1).c_str(), h_sigmaIphiIphi_EE_caloMatched_old); h_sigmaIphiIphi_EE_caloMatched_old->SetName("h_SigmaIphiIphi_EE_caloMatched_old");
    hist_infile_SCRef->GetObject(("h_SigmaIphiIphi_EE_caloUnmatched_"+ref_1).c_str(), h_sigmaIphiIphi_EE_caloUnmatched_old); h_sigmaIphiIphi_EE_caloUnmatched_old->SetName("h_SigmaIphiIphi_EE_caloUnmatched_old");


    int nBins_Eta = 100;//binOpts[findOption(std::string("EtaBins"),binOpts)].second[0];
    EoEtrue_vs_Eta_Calo_old.resize(nBins_Eta);
    for(int iBin=0; iBin<nBins_Eta; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_Eta_Calo_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Eta_Calo_old[iBin]); 
        EoEtrue_vs_Eta_Calo_old[iBin]->SetName(("EoEtrue_vs_Eta_Calo_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_Eta_Gen_old.resize(nBins_Eta);
    for(int iBin=0; iBin<nBins_Eta; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_Eta_Gen_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Eta_Gen_old[iBin]); 
        EoEgen_vs_Eta_Gen_old[iBin]->SetName(("EoEgen_vs_Eta_Gen_old_"+to_string(iBin)).c_str());
    }EoEtrue_vs_Eta_Seed_old.resize(nBins_Eta);
    for(int iBin=0; iBin<nBins_Eta; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_Eta_Seed_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Eta_Seed_old[iBin]); 
        EoEtrue_vs_Eta_Seed_old[iBin]->SetName(("EoEtrue_vs_Eta_Seed_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_Eta_Seed_old.resize(nBins_Eta);
    for(int iBin=0; iBin<nBins_Eta; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_Eta_Seed_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Eta_Seed_old[iBin]); 
        EoEgen_vs_Eta_Seed_old[iBin]->SetName(("EoEgen_vs_Eta_Seed_old_"+to_string(iBin)).c_str());
      } 
    int nBins_Et_EB = 50;//binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[0];
    int nBins_Et_EE = 50;//binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[0];
    EoEtrue_vs_Et_Calo_EB_old.resize(nBins_Et_EB);
    for(int iBin=0; iBin<nBins_Et_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_Et_Calo_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Et_Calo_EB_old[iBin]); 
        EoEtrue_vs_Et_Calo_EB_old[iBin]->SetName(("EoEtrue_vs_Et_Calo_EB_old_"+to_string(iBin)).c_str());
    }EoEtrue_vs_Et_Calo_EE_old.resize(nBins_Et_EE);
    for(int iBin=0; iBin<nBins_Et_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_Et_Calo_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Et_Calo_EE_old[iBin]); 
        EoEtrue_vs_Et_Calo_EE_old[iBin]->SetName(("EoEtrue_vs_Et_Calo_EE_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_Et_Gen_EB_old.resize(nBins_Et_EB);
    for(int iBin=0; iBin<nBins_Et_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_Et_Gen_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Et_Gen_EB_old[iBin]); 
        EoEgen_vs_Et_Gen_EB_old[iBin]->SetName(("EoEgen_vs_Et_Gen_EB_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_Et_Gen_EE_old.resize(nBins_Et_EE);
    for(int iBin=0; iBin<nBins_Et_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_Et_Gen_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Et_Gen_EE_old[iBin]); 
        EoEgen_vs_Et_Gen_EE_old[iBin]->SetName(("EoEgen_vs_Et_Gen_EE_old_"+to_string(iBin)).c_str());
    }EoEtrue_vs_Et_Seed_EB_old.resize(nBins_Et_EB);
    for(int iBin=0; iBin<nBins_Et_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_Et_Seed_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Et_Seed_EB_old[iBin]); 
        EoEtrue_vs_Et_Seed_EB_old[iBin]->SetName(("EoEtrue_vs_Et_Seed_EB_old_"+to_string(iBin)).c_str());
    }EoEtrue_vs_Et_Seed_EE_old.resize(nBins_Et_EE);
    for(int iBin=0; iBin<nBins_Et_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_Et_Seed_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Et_Seed_EE_old[iBin]); 
        EoEtrue_vs_Et_Seed_EE_old[iBin]->SetName(("EoEtrue_vs_Et_Seed_EE_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_Et_Seed_EB_old.resize(nBins_Et_EB);
    for(int iBin=0; iBin<nBins_Et_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_Et_Seed_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Et_Seed_EB_old[iBin]); 
        EoEgen_vs_Et_Seed_EB_old[iBin]->SetName(("EoEgen_vs_Et_Seed_EB_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_Et_Seed_EE_old.resize(nBins_Et_EE);
    for(int iBin=0; iBin<nBins_Et_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_Et_Seed_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Et_Seed_EE_old[iBin]); 
        EoEgen_vs_Et_Seed_EE_old[iBin]->SetName(("EoEgen_vs_Et_Seed_EE_old_"+to_string(iBin)).c_str());
       } 
    int nBins_Energy_EB = 150;//binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[0];
    int nBins_Energy_EE = 100;//binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[0];
    EoEtrue_vs_Energy_Calo_EB_old.resize(nBins_Energy_EB);
    for(int iBin=0; iBin<nBins_Energy_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_Energy_Calo_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Energy_Calo_EB_old[iBin]); 
        EoEtrue_vs_Energy_Calo_EB_old[iBin]->SetName(("EoEtrue_vs_Energy_Calo_EB_old_"+to_string(iBin)).c_str());
    }EoEtrue_vs_Energy_Calo_EE_old.resize(nBins_Energy_EE);
    for(int iBin=0; iBin<nBins_Energy_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_Energy_Calo_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Energy_Calo_EE_old[iBin]); 
        EoEtrue_vs_Energy_Calo_EE_old[iBin]->SetName(("EoEtrue_vs_Energy_Calo_EE_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_Energy_Gen_EB_old.resize(nBins_Energy_EB);
    for(int iBin=0; iBin<nBins_Energy_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_Energy_Calo_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Energy_Gen_EB_old[iBin]); 
        EoEgen_vs_Energy_Gen_EB_old[iBin]->SetName(("EoEgen_vs_Energy_Calo_EB_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_Energy_Gen_EE_old.resize(nBins_Energy_EE);
    for(int iBin=0; iBin<nBins_Energy_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_Energy_Calo_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Energy_Gen_EE_old[iBin]); 
        EoEgen_vs_Energy_Gen_EE_old[iBin]->SetName(("EoEgen_vs_Energy_Calo_EE_old_"+to_string(iBin)).c_str());
    }
    int nBins_nVtx_EB = 60;//binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[0];
    int nBins_nVtx_EE = 60;//binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[0];
    EoEtrue_vs_nVtx_EB_old.resize(nBins_nVtx_EB);
    for(int iBin=0; iBin<nBins_nVtx_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_nVtx_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_nVtx_EB_old[iBin]); 
        EoEtrue_vs_nVtx_EB_old[iBin]->SetName(("EoEtrue_vs_nVtx_EB_old_"+to_string(iBin)).c_str());
    }EoEtrue_vs_nVtx_EE_old.resize(nBins_nVtx_EE);
    for(int iBin=0; iBin<nBins_nVtx_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_nVtx_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_nVtx_EE_old[iBin]); 
        EoEtrue_vs_nVtx_EE_old[iBin]->SetName(("EoEtrue_vs_nVtx_EE_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_nVtx_EB_old.resize(nBins_nVtx_EB);
    for(int iBin=0; iBin<nBins_nVtx_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_nVtx_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_nVtx_EB_old[iBin]); 
        EoEgen_vs_nVtx_EB_old[iBin]->SetName(("EoEgen_vs_nVtx_EB_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_nVtx_EE_old.resize(nBins_nVtx_EE);
    for(int iBin=0; iBin<nBins_nVtx_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_nVtx_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_nVtx_EE_old[iBin]); 
        EoEgen_vs_nVtx_EE_old[iBin]->SetName(("EoEgen_vs_nVtx_EE_old_"+to_string(iBin)).c_str());
    }
    int nBins_Rho_EB = 40;//binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[0]; 
    int nBins_Rho_EE = 40;//binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[0];
    EoEtrue_vs_Rho_EB_old.resize(nBins_Rho_EB);
    for(int iBin=0; iBin<nBins_Rho_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_Rho_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Rho_EB_old[iBin]); 
        EoEtrue_vs_Rho_EB_old[iBin]->SetName(("EoEtrue_vs_Rho_EB_old_"+to_string(iBin)).c_str());
    }EoEtrue_vs_Rho_EE_old.resize(nBins_Rho_EE);
    for(int iBin=0; iBin<nBins_Rho_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_Rho_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Rho_EE_old[iBin]); 
        EoEtrue_vs_Rho_EE_old[iBin]->SetName(("EoEtrue_vs_Rho_EE_old_"+to_string(iBin)).c_str());
    }EoEgen_vs_Rho_EB_old.resize(nBins_Rho_EB);
    for(int iBin=0; iBin<nBins_Rho_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_Rho_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Rho_EB_old[iBin]); 
        EoEgen_vs_Rho_EB_old[iBin]->SetName(("EoEgen_vs_Rho_EB_old_"+to_string(iBin)).c_str());
   } EoEgen_vs_Rho_EE_old.resize(nBins_Rho_EE);
    for(int iBin=0; iBin<nBins_Rho_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEgen_vs_Rho_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Rho_EE_old[iBin]); 
        EoEgen_vs_Rho_EE_old[iBin]->SetName(("EoEgen_vs_Rho_EE_old_"+to_string(iBin)).c_str());
    
    }EoEtrue_vs_seedEt_seedEta_old.resize(etCuts.size()-1);
    for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
    {
        EoEtrue_vs_seedEt_seedEta_old[iBin].resize(etaCuts.size()-1); 
        for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
        {
            hist_infile_SCRef->GetObject(("EoEtrue_vs_seedEt_"+to_string(etCuts.at(iBin))+"_"+to_string(etCuts.at(iBin+1))+"_seedEta_"+to_string(etaCuts.at(jBin))+"_"+to_string(etaCuts.at(jBin+1))+"_new").c_str(), EoEtrue_vs_seedEt_seedEta_old[iBin][jBin]); 
            EoEtrue_vs_seedEt_seedEta_old[iBin][jBin]->SetName(("EoEtrue_vs_seedEt_"+to_string(etCuts.at(iBin))+"_"+to_string(etCuts.at(iBin+1))+"_seedEta_"+to_string(etaCuts.at(jBin))+"_"+to_string(etaCuts.at(jBin+1))+"_old").c_str());
            }
    }


    int nBins_EtaWidth_EB = 100;//binOpts[findOption(std::string("EtaWidthBins_Barrel"),binOpts)].second[0]; 
    int nBins_EtaWidth_EE = 100;//binOpts[findOption(std::string("EtaWidthBins_Endcap"),binOpts)].second[0];
    
    EoEtrue_vs_EtaWidth_EB_old.resize(nBins_EtaWidth_EB);
    for(int iBin=0; iBin<nBins_EtaWidth_EB; iBin++){
        cout<<"A "<<iBin<<endl;
        hist_infile_SCRef->GetObject(("EoEtrue_vs_EtaWidth_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_EtaWidth_EB_old[iBin]); 
        cout<<"B "<<iBin<<endl;
        EoEtrue_vs_EtaWidth_EB_old[iBin]->SetName(("EoEtrue_vs_EtaWidth_EB_old_"+to_string(iBin)).c_str());
    }
    cout<<"C"<<endl;
    EoEtrue_vs_EtaWidth_EE_old.resize(nBins_EtaWidth_EE);
    for(int iBin=0; iBin<nBins_EtaWidth_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_EtaWidth_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_EtaWidth_EE_old[iBin]); 
        EoEtrue_vs_EtaWidth_EE_old[iBin]->SetName(("EoEtrue_vs_EtaWidth_EE_old_"+to_string(iBin)).c_str());
    }
    int nBins_PhiWidth_EB = 100;//binOpts[findOption(std::string("PhiWidthBins_Barrel"),binOpts)].second[0]; 
    int nBins_PhiWidth_EE = 100;//binOpts[findOption(std::string("PhiWidthBins_Endcap"),binOpts)].second[0];
    
    EoEtrue_vs_PhiWidth_EB_old.resize(nBins_PhiWidth_EB);
    for(int iBin=0; iBin<nBins_PhiWidth_EB; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_PhiWidth_EB_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_PhiWidth_EB_old[iBin]); 
        EoEtrue_vs_PhiWidth_EB_old[iBin]->SetName(("EoEtrue_vs_PhiWidth_EB_old_"+to_string(iBin)).c_str());
    }
    EoEtrue_vs_PhiWidth_EE_old.resize(nBins_PhiWidth_EE);
    for(int iBin=0; iBin<nBins_PhiWidth_EE; iBin++){
        hist_infile_SCRef->GetObject(("EoEtrue_vs_PhiWidth_EE_"+ref_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_PhiWidth_EE_old[iBin]); 
        EoEtrue_vs_PhiWidth_EE_old[iBin]->SetName(("EoEtrue_vs_PhiWidth_EE_old_"+to_string(iBin)).c_str());
    }
    
    hist_infile_SCRef->GetObject(("h2_EoEtrue_Mean_"+ref_2).c_str(), h2_EoEtrue_Mean_old); h2_EoEtrue_Mean_old->SetName("h2_EoEtrue_Mean_old");
    hist_infile_SCRef->GetObject(("h2_EoEtrue_Mean_Effective_"+ref_2).c_str(), h2_EoEtrue_Mean_Effective_old); h2_EoEtrue_Mean_Effective_old->SetName("h2_EoEtrue_Mean_Effective_old");
    hist_infile_SCRef->GetObject(("h2_EoEtrue_Resolution_"+ref_2).c_str(), h2_EoEtrue_Resolution_old); h2_EoEtrue_Resolution_old->SetName("h2_EoEtrue_Resolution_old");
    hist_infile_SCRef->GetObject(("h2_EoEtrue_Resolution_Effective_"+ref_2).c_str(), h2_EoEtrue_Resolution_Effective_old); h2_EoEtrue_Resolution_Effective_old->SetName("h2_EoEtrue_Resolution_Effective_old");

    //local Set Params -> delimited by new
    hist_infile_SCVal->GetObject(("h_Eta_"+val_1).c_str(), h_Eta_new);
    hist_infile_SCVal->GetObject(("h_Eta_Calo_"+val_1).c_str(), h_Eta_Calo_new);

        //Eta
    hist_infile_SCVal->GetObject(("h_Eta_"+val_1).c_str(), h_Eta_new);
    hist_infile_SCVal->GetObject(("h_Eta_Calo_"+val_1).c_str(), h_Eta_Calo_new);
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_Eta_Calo_"+val_1).c_str(), prof_EoEtrue_vs_Eta_Calo_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_Eta_Gen_"+val_1).c_str(), prof_EoEgen_vs_Eta_Gen_new);
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_Eta_Seed_"+val_1).c_str(), prof_EoEtrue_vs_Eta_Seed_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_Eta_Seed_"+val_1).c_str(), prof_EoEgen_vs_Eta_Seed_new);
    hist_infile_SCVal->GetObject(("h_Eta_seedMatched_"+val_1).c_str(), h_Eta_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_Eta_caloMatched_"+val_1).c_str(), h_Eta_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_Eta_caloUnmatched_"+val_1).c_str(), h_Eta_caloUnmatched_new);
        //nPFClusters
    hist_infile_SCVal->GetObject(("h_nPFClusters_"+val_1).c_str(), h_nPFClusters_new);
    hist_infile_SCVal->GetObject(("h_nPFClusters_seedMatched_"+val_1).c_str(), h_nPFClusters_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_nPFClusters_caloMatched_"+val_1).c_str(), h_nPFClusters_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_nPFClusters_caloUnmatched_"+val_1).c_str(), h_nPFClusters_caloUnmatched_new);
    
    hist_infile_SCVal->GetObject(("h_nPFClusters_EB_"+val_1).c_str(), h_nPFClusters_EB_new);
    hist_infile_SCVal->GetObject(("h_nPFClusters_EB_seedMatched_"+val_1).c_str(), h_nPFClusters_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_nPFClusters_EB_caloMatched_"+val_1).c_str(), h_nPFClusters_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_nPFClusters_EB_caloUnmatched_"+val_1).c_str(), h_nPFClusters_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_nPFClusters_EE_"+val_1).c_str(), h_nPFClusters_EE_new);
    hist_infile_SCVal->GetObject(("h_nPFClusters_EE_seedMatched_"+val_1).c_str(), h_nPFClusters_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_nPFClusters_EE_caloMatched_"+val_1).c_str(), h_nPFClusters_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_nPFClusters_EE_caloUnmatched_"+val_1).c_str(), h_nPFClusters_EE_caloUnmatched_new);
         
        //nVtx
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_nVtx_EB_"+val_1).c_str(), prof_EoEtrue_vs_nVtx_EB_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_nVtx_EB_"+val_1).c_str(), prof_EoEgen_vs_nVtx_EB_new);
             
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_nVtx_EE_"+val_1).c_str(), prof_EoEtrue_vs_nVtx_EE_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_nVtx_EE_"+val_1).c_str(), prof_EoEgen_vs_nVtx_EE_new);
             
        //Rho
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_Rho_EB_"+val_1).c_str(), prof_EoEtrue_vs_Rho_EB_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_Rho_EB_"+val_1).c_str(), prof_EoEgen_vs_Rho_EB_new);
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_Rho_EE_"+val_1).c_str(), prof_EoEtrue_vs_Rho_EE_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_Rho_EE_"+val_1).c_str(), prof_EoEgen_vs_Rho_EE_new);

        //Energy
    hist_infile_SCVal->GetObject(("h_Energy_EB_"+val_1).c_str(), h_Energy_EB_new);
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_Energy_Calo_EB_"+val_1).c_str(), prof_EoEtrue_vs_Energy_Calo_EB_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_Energy_Gen_EB_"+val_1).c_str(), prof_EoEgen_vs_Energy_Gen_EB_new);
    hist_infile_SCVal->GetObject(("h_Energy_EB_seedMatched_"+val_1).c_str(), h_Energy_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_Energy_EB_caloMatched_"+val_1).c_str(), h_Energy_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_Energy_EB_caloUnmatched_"+val_1).c_str(), h_Energy_EB_caloUnmatched_new);

    hist_infile_SCVal->GetObject(("h_Energy_EE_"+val_1).c_str(), h_Energy_EE_new);
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_Energy_Calo_EE_"+val_1).c_str(), prof_EoEtrue_vs_Energy_Calo_EE_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_Energy_Gen_EE_"+val_1).c_str(), prof_EoEgen_vs_Energy_Gen_EE_new);
    hist_infile_SCVal->GetObject(("h_Energy_EE_seedMatched_"+val_1).c_str(), h_Energy_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_Energy_EE_caloMatched_"+val_1).c_str(), h_Energy_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_Energy_EE_caloUnmatched_"+val_1).c_str(), h_Energy_EE_caloUnmatched_new);
                        
        //EoEtrue
    hist_infile_SCVal->GetObject(("h_EoEtrue_EB_"+val_1).c_str(), h_EoEtrue_EB_new);
    hist_infile_SCVal->GetObject(("h_EoEtrue_EB_seedMatched_"+val_1).c_str(), h_EoEtrue_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_EoEtrue_EE_"+val_1).c_str(), h_EoEtrue_EE_new);
    hist_infile_SCVal->GetObject(("h_EoEtrue_EE_seedMatched_"+val_1).c_str(), h_EoEtrue_EE_seedMatched_new);
             
        //EoEgen
    hist_infile_SCVal->GetObject(("h_EoEgen_EB_"+val_1).c_str(), h_EoEgen_EB_new);
    hist_infile_SCVal->GetObject(("h_EoEgen_EB_seedMatched_"+val_1).c_str(), h_EoEgen_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_EoEgen_EE_"+val_1).c_str(), h_EoEgen_EE_new);
    hist_infile_SCVal->GetObject(("h_EoEgen_EE_seedMatched_"+val_1).c_str(), h_EoEgen_EE_seedMatched_new);
             
        //Et
    hist_infile_SCVal->GetObject(("h_Et_EB_"+val_1).c_str(), h_Et_EB_new);
    hist_infile_SCVal->GetObject(("h_Et_Calo_EB_"+val_1).c_str(), h_Et_Calo_EB_new);
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_Et_Calo_EB_"+val_1).c_str(), prof_EoEtrue_vs_Et_Calo_EB_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_Et_Gen_EB_"+val_1).c_str(), prof_EoEgen_vs_Et_Gen_EB_new);
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_Et_Seed_EB_"+val_1).c_str(), prof_EoEtrue_vs_Et_Seed_EB_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_Et_Seed_EB_"+val_1).c_str(), prof_EoEgen_vs_Et_Seed_EB_new);
    hist_infile_SCVal->GetObject(("h_Et_EB_seedMatched_"+val_1).c_str(), h_Et_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_Et_EB_caloMatched_"+val_1).c_str(), h_Et_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_Et_EB_caloUnmatched_"+val_1).c_str(), h_Et_EB_caloUnmatched_new);

    hist_infile_SCVal->GetObject(("h_Et_EE_"+val_1).c_str(), h_Et_EE_new);
    hist_infile_SCVal->GetObject(("h_Et_Calo_EE_"+val_1).c_str(), h_Et_Calo_EE_new);
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_Et_Calo_EE_"+val_1).c_str(), prof_EoEtrue_vs_Et_Calo_EE_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_Et_Gen_EE_"+val_1).c_str(), prof_EoEgen_vs_Et_Gen_EE_new);
    hist_infile_SCVal->GetObject(("prof_EoEtrue_vs_Et_Seed_EE_"+val_1).c_str(), prof_EoEtrue_vs_Et_Seed_EE_new);
    hist_infile_SCVal->GetObject(("prof_EoEgen_vs_Et_Seed_EE_"+val_1).c_str(), prof_EoEgen_vs_Et_Seed_EE_new);
    hist_infile_SCVal->GetObject(("h_Et_EE_seedMatched_"+val_1).c_str(), h_Et_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_Et_EE_caloMatched_"+val_1).c_str(), h_Et_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_Et_EE_caloUnmatched_"+val_1).c_str(), h_Et_EE_caloUnmatched_new);


        //Phi
    hist_infile_SCVal->GetObject(("h_Phi_EB_"+val_1).c_str(), h_Phi_EB_new);
    hist_infile_SCVal->GetObject(("h_Phi_EB_seedMatched_"+val_1).c_str(), h_Phi_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_Phi_EB_caloMatched_"+val_1).c_str(), h_Phi_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_Phi_EB_caloUnmatched_"+val_1).c_str(), h_Phi_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_Phi_EE_"+val_1).c_str(), h_Phi_EE_new);
    hist_infile_SCVal->GetObject(("h_Phi_EE_seedMatched_"+val_1).c_str(), h_Phi_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_Phi_EE_caloMatched_"+val_1).c_str(), h_Phi_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_Phi_EE_caloUnmatched_"+val_1).c_str(), h_Phi_EE_caloUnmatched_new);
             
        //EtaWidth
    hist_infile_SCVal->GetObject(("h_EtaWidth_EB_"+val_1).c_str(), h_EtaWidth_EB_new);
    hist_infile_SCVal->GetObject(("h_EtaWidth_EB_seedMatched_"+val_1).c_str(), h_EtaWidth_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_EtaWidth_EB_caloMatched_"+val_1).c_str(), h_EtaWidth_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_EtaWidth_EB_caloUnmatched_"+val_1).c_str(), h_EtaWidth_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_EtaWidth_EE_"+val_1).c_str(), h_EtaWidth_EE_new);
    hist_infile_SCVal->GetObject(("h_EtaWidth_EE_seedMatched_"+val_1).c_str(), h_EtaWidth_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_EtaWidth_EE_caloMatched_"+val_1).c_str(), h_EtaWidth_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_EtaWidth_EE_caloUnmatched_"+val_1).c_str(), h_EtaWidth_EE_caloUnmatched_new);
             
        //PhiWidth
    hist_infile_SCVal->GetObject(("h_PhiWidth_EB_"+val_1).c_str(), h_PhiWidth_EB_new);
    hist_infile_SCVal->GetObject(("h_PhiWidth_EB_seedMatched_"+val_1).c_str(), h_PhiWidth_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_PhiWidth_EB_caloMatched_"+val_1).c_str(), h_PhiWidth_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_PhiWidth_EB_caloUnmatched_"+val_1).c_str(), h_PhiWidth_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_PhiWidth_EE_"+val_1).c_str(), h_PhiWidth_EE_new);
    hist_infile_SCVal->GetObject(("h_PhiWidth_EE_seedMatched_"+val_1).c_str(), h_PhiWidth_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_PhiWidth_EE_caloMatched_"+val_1).c_str(), h_PhiWidth_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_PhiWidth_EE_caloUnmatched_"+val_1).c_str(), h_PhiWidth_EE_caloUnmatched_new);
             
        //full5x5_R9
    hist_infile_SCVal->GetObject(("h_full5x5_R9_EB_"+val_1).c_str(), h_full5x5_R9_EB_new);
    hist_infile_SCVal->GetObject(("h_full5x5_R9_EB_seedMatched_"+val_1).c_str(), h_full5x5_R9_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_R9_EB_caloMatched_"+val_1).c_str(), h_full5x5_R9_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_R9_EB_caloUnmatched_"+val_1).c_str(), h_full5x5_R9_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_full5x5_R9_EE_"+val_1).c_str(), h_full5x5_R9_EE_new);
    hist_infile_SCVal->GetObject(("h_full5x5_R9_EE_seedMatched_"+val_1).c_str(), h_full5x5_R9_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_R9_EE_caloMatched_"+val_1).c_str(), h_full5x5_R9_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_R9_EE_caloUnmatched_"+val_1).c_str(), h_full5x5_R9_EE_caloUnmatched_new);
             
        //R9
    hist_infile_SCVal->GetObject(("h_R9_EB_"+val_1).c_str(), h_R9_EB_new);
    hist_infile_SCVal->GetObject(("h_R9_EB_seedMatched_"+val_1).c_str(), h_R9_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_R9_EB_caloMatched_"+val_1).c_str(), h_R9_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_R9_EB_caloUnmatched_"+val_1).c_str(), h_R9_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_R9_EE_"+val_1).c_str(), h_R9_EE_new);
    hist_infile_SCVal->GetObject(("h_R9_EE_seedMatched_"+val_1).c_str(), h_R9_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_R9_EE_caloMatched_"+val_1).c_str(), h_R9_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_R9_EE_caloUnmatched_"+val_1).c_str(), h_R9_EE_caloUnmatched_new);
             
        //full5x5_SigmaIetaIeta
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIeta_EB_"+val_1).c_str(), h_full5x5_sigmaIetaIeta_EB_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIeta_EB_seedMatched_"+val_1).c_str(), h_full5x5_sigmaIetaIeta_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIeta_EB_caloMatched_"+val_1).c_str(), h_full5x5_sigmaIetaIeta_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIeta_EB_caloUnmatched_"+val_1).c_str(), h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIeta_EE_"+val_1).c_str(), h_full5x5_sigmaIetaIeta_EE_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIeta_EE_seedMatched_"+val_1).c_str(), h_full5x5_sigmaIetaIeta_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIeta_EE_caloMatched_"+val_1).c_str(), h_full5x5_sigmaIetaIeta_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIeta_EE_caloUnmatched_"+val_1).c_str(), h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new);
             
        //full5x5_SigmaIetaIphi
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIphi_EB_"+val_1).c_str(), h_full5x5_sigmaIetaIphi_EB_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIphi_EB_seedMatched_"+val_1).c_str(), h_full5x5_sigmaIetaIphi_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIphi_EB_caloMatched_"+val_1).c_str(), h_full5x5_sigmaIetaIphi_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIphi_EB_caloUnmatched_"+val_1).c_str(), h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIphi_EE_"+val_1).c_str(), h_full5x5_sigmaIetaIphi_EE_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIphi_EE_seedMatched_"+val_1).c_str(), h_full5x5_sigmaIetaIphi_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIphi_EE_caloMatched_"+val_1).c_str(), h_full5x5_sigmaIetaIphi_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIetaIphi_EE_caloUnmatched_"+val_1).c_str(), h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new);
             
        //SigmaIetaIeta
    hist_infile_SCVal->GetObject(("h_SigmaIetaIeta_EB_"+val_1).c_str(), h_sigmaIetaIeta_EB_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIeta_EB_seedMatched_"+val_1).c_str(), h_sigmaIetaIeta_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIeta_EB_caloMatched_"+val_1).c_str(), h_sigmaIetaIeta_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIeta_EB_caloUnmatched_"+val_1).c_str(), h_sigmaIetaIeta_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_SigmaIetaIeta_EE_"+val_1).c_str(), h_sigmaIetaIeta_EE_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIeta_EE_seedMatched_"+val_1).c_str(), h_sigmaIetaIeta_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIeta_EE_caloMatched_"+val_1).c_str(), h_sigmaIetaIeta_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIeta_EE_caloUnmatched_"+val_1).c_str(), h_sigmaIetaIeta_EE_caloUnmatched_new);
             
        //SigmaIetaIphi
    hist_infile_SCVal->GetObject(("h_SigmaIetaIphi_EB_"+val_1).c_str(), h_sigmaIetaIphi_EB_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIphi_EB_seedMatched_"+val_1).c_str(), h_sigmaIetaIphi_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIphi_EB_caloMatched_"+val_1).c_str(), h_sigmaIetaIphi_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIphi_EB_caloUnmatched_"+val_1).c_str(), h_sigmaIetaIphi_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_SigmaIetaIphi_EE_"+val_1).c_str(), h_sigmaIetaIphi_EE_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIphi_EE_seedMatched_"+val_1).c_str(), h_sigmaIetaIphi_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIphi_EE_caloMatched_"+val_1).c_str(), h_sigmaIetaIphi_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIetaIphi_EE_caloUnmatched_"+val_1).c_str(), h_sigmaIetaIphi_EE_caloUnmatched_new);
             
        //full5x5_SigmaIphiIphi
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIphiIphi_EB_"+val_1).c_str(), h_full5x5_sigmaIphiIphi_EB_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIphiIphi_EB_seedMatched_"+val_1).c_str(), h_full5x5_sigmaIphiIphi_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIphiIphi_EB_caloMatched_"+val_1).c_str(), h_full5x5_sigmaIphiIphi_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIphiIphi_EB_caloUnmatched_"+val_1).c_str(), h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIphiIphi_EE_"+val_1).c_str(), h_full5x5_sigmaIphiIphi_EE_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIphiIphi_EE_seedMatched_"+val_1).c_str(), h_full5x5_sigmaIphiIphi_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIphiIphi_EE_caloMatched_"+val_1).c_str(), h_full5x5_sigmaIphiIphi_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_full5x5_SigmaIphiIphi_EE_caloUnmatched_"+val_1).c_str(), h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new);
             
        //SigmaIphiIphi
    hist_infile_SCVal->GetObject(("h_SigmaIphiIphi_EB_"+val_1).c_str(), h_sigmaIphiIphi_EB_new);
    hist_infile_SCVal->GetObject(("h_SigmaIphiIphi_EB_seedMatched_"+val_1).c_str(), h_sigmaIphiIphi_EB_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIphiIphi_EB_caloMatched_"+val_1).c_str(), h_sigmaIphiIphi_EB_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIphiIphi_EB_caloUnmatched_"+val_1).c_str(), h_sigmaIphiIphi_EB_caloUnmatched_new);
          
    hist_infile_SCVal->GetObject(("h_SigmaIphiIphi_EE_"+val_1).c_str(), h_sigmaIphiIphi_EE_new);
    hist_infile_SCVal->GetObject(("h_SigmaIphiIphi_EE_seedMatched_"+val_1).c_str(), h_sigmaIphiIphi_EE_seedMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIphiIphi_EE_caloMatched_"+val_1).c_str(), h_sigmaIphiIphi_EE_caloMatched_new);
    hist_infile_SCVal->GetObject(("h_SigmaIphiIphi_EE_caloUnmatched_"+val_1).c_str(), h_sigmaIphiIphi_EE_caloUnmatched_new);


    EoEtrue_vs_Eta_Calo_new.resize(nBins_Eta);
    for(int iBin=0; iBin<nBins_Eta; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_Eta_Calo_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Eta_Calo_new[iBin]);
    EoEgen_vs_Eta_Gen_new.resize(nBins_Eta);
    for(int iBin=0; iBin<nBins_Eta; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_Eta_Gen_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Eta_Gen_new[iBin]);
    EoEtrue_vs_Eta_Seed_new.resize(nBins_Eta);
    for(int iBin=0; iBin<nBins_Eta; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_Eta_Seed_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Eta_Seed_new[iBin]);
    EoEgen_vs_Eta_Seed_new.resize(nBins_Eta);
    for(int iBin=0; iBin<nBins_Eta; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_Eta_Seed_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Eta_Seed_new[iBin]);
       
    EoEtrue_vs_Et_Calo_EB_new.resize(nBins_Et_EB);
    for(int iBin=0; iBin<nBins_Et_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_Et_Calo_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Et_Calo_EB_new[iBin]);
    EoEtrue_vs_Et_Calo_EE_new.resize(nBins_Et_EE);
    for(int iBin=0; iBin<nBins_Et_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_Et_Calo_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Et_Calo_EE_new[iBin]);
    EoEgen_vs_Et_Gen_EB_new.resize(nBins_Et_EB);
    for(int iBin=0; iBin<nBins_Et_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_Et_Gen_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Et_Gen_EB_new[iBin]);
    EoEgen_vs_Et_Gen_EE_new.resize(nBins_Et_EE);
    for(int iBin=0; iBin<nBins_Et_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_Et_Gen_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Et_Gen_EE_new[iBin]);
    EoEtrue_vs_Et_Seed_EB_new.resize(nBins_Et_EB);
    for(int iBin=0; iBin<nBins_Et_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_Et_Seed_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Et_Seed_EB_new[iBin]);
    EoEtrue_vs_Et_Seed_EE_new.resize(nBins_Et_EE);
    for(int iBin=0; iBin<nBins_Et_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_Et_Seed_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Et_Seed_EE_new[iBin]);
    EoEgen_vs_Et_Seed_EB_new.resize(nBins_Et_EB);
    for(int iBin=0; iBin<nBins_Et_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_Et_Seed_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Et_Seed_EB_new[iBin]);
    EoEgen_vs_Et_Seed_EE_new.resize(nBins_Et_EE);
    for(int iBin=0; iBin<nBins_Et_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_Et_Seed_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Et_Seed_EE_new[iBin]);
        
    EoEtrue_vs_Energy_Calo_EB_new.resize(nBins_Energy_EB);
    for(int iBin=0; iBin<nBins_Energy_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_Energy_Calo_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Energy_Calo_EB_new[iBin]);
    EoEtrue_vs_Energy_Calo_EE_new.resize(nBins_Energy_EE);
    for(int iBin=0; iBin<nBins_Energy_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_Energy_Calo_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Energy_Calo_EE_new[iBin]);
    EoEgen_vs_Energy_Gen_EB_new.resize(nBins_Energy_EB);
    for(int iBin=0; iBin<nBins_Energy_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_Energy_Calo_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Energy_Gen_EB_new[iBin]);
    EoEgen_vs_Energy_Gen_EE_new.resize(nBins_Energy_EE);
    for(int iBin=0; iBin<nBins_Energy_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_Energy_Calo_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Energy_Gen_EE_new[iBin]);
    
    EoEtrue_vs_nVtx_EB_new.resize(nBins_nVtx_EB);
    for(int iBin=0; iBin<nBins_nVtx_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_nVtx_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_nVtx_EB_new[iBin]);
    EoEtrue_vs_nVtx_EE_new.resize(nBins_nVtx_EE);
    for(int iBin=0; iBin<nBins_nVtx_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_nVtx_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_nVtx_EE_new[iBin]);
    EoEgen_vs_nVtx_EB_new.resize(nBins_nVtx_EB);
    for(int iBin=0; iBin<nBins_nVtx_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_nVtx_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_nVtx_EB_new[iBin]);
    EoEgen_vs_nVtx_EE_new.resize(nBins_nVtx_EE);
    for(int iBin=0; iBin<nBins_nVtx_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_nVtx_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_nVtx_EE_new[iBin]);
    
    EoEtrue_vs_Rho_EB_new.resize(nBins_Rho_EB);
    for(int iBin=0; iBin<nBins_Rho_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_Rho_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Rho_EB_new[iBin]);
    EoEtrue_vs_Rho_EE_new.resize(nBins_Rho_EE);
    for(int iBin=0; iBin<nBins_Rho_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_Rho_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_Rho_EE_new[iBin]);
    EoEgen_vs_Rho_EB_new.resize(nBins_Rho_EB);
    for(int iBin=0; iBin<nBins_Rho_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_Rho_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Rho_EB_new[iBin]);
    EoEgen_vs_Rho_EE_new.resize(nBins_Rho_EE);
    for(int iBin=0; iBin<nBins_Rho_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEgen_vs_Rho_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEgen_vs_Rho_EE_new[iBin]);
    
    EoEtrue_vs_seedEt_seedEta_new.resize(etCuts.size()-1);
    for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
    {
        EoEtrue_vs_seedEt_seedEta_new[iBin].resize(etaCuts.size()-1); 
        for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
        {
            hist_infile_SCVal->GetObject(("EoEtrue_vs_seedEt_"+to_string(etCuts.at(iBin))+"_"+to_string(etCuts.at(iBin+1))+"_seedEta_"+to_string(etaCuts.at(jBin))+"_"+to_string(etaCuts.at(jBin+1))+"_new").c_str(), EoEtrue_vs_seedEt_seedEta_new[iBin][jBin]);
            }
    }

    EoEtrue_vs_EtaWidth_EB_new.resize(nBins_EtaWidth_EB);
    for(int iBin=0; iBin<nBins_EtaWidth_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_EtaWidth_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_EtaWidth_EB_new[iBin]);
    EoEtrue_vs_EtaWidth_EE_new.resize(nBins_EtaWidth_EE);
    for(int iBin=0; iBin<nBins_EtaWidth_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_EtaWidth_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_EtaWidth_EE_new[iBin]);
    
    
    EoEtrue_vs_PhiWidth_EB_new.resize(nBins_PhiWidth_EB);
    for(int iBin=0; iBin<nBins_PhiWidth_EB; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_PhiWidth_EB_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_PhiWidth_EB_new[iBin]);
    EoEtrue_vs_PhiWidth_EE_new.resize(nBins_PhiWidth_EE);
    for(int iBin=0; iBin<nBins_PhiWidth_EE; iBin++)
        hist_infile_SCVal->GetObject(("EoEtrue_vs_PhiWidth_EE_"+val_2+"_"+to_string(iBin)).c_str(), EoEtrue_vs_PhiWidth_EE_new[iBin]);
   
    hist_infile_SCVal->GetObject(("h2_EoEtrue_Mean_"+val_2).c_str(), h2_EoEtrue_Mean_new);
    hist_infile_SCVal->GetObject(("h2_EoEtrue_Mean_Effective_"+val_2).c_str(), h2_EoEtrue_Mean_Effective_new);
    hist_infile_SCVal->GetObject(("h2_EoEtrue_Resolution_"+val_2).c_str(), h2_EoEtrue_Resolution_new);
    hist_infile_SCVal->GetObject(("h2_EoEtrue_Resolution_Effective_"+val_2).c_str(), h2_EoEtrue_Resolution_Effective_new);
        
    //neeed new
    hist_infile_SCVal->GetObject("h_Eta_Calo_Denum", h_Eta_Calo_Denum_new);
    hist_infile_SCVal->GetObject("h_Et_Calo_EB_Denum", h_Et_Calo_EB_Denum_new);
    hist_infile_SCVal->GetObject("h_Et_Calo_EE_Denum", h_Et_Calo_EE_Denum_new);

    hist_infile_SCRef ->GetObject("h_Eta_SCGen_new", h_Eta_SCGen_old);
    h_Eta_SCGen_old->SetName("h_Eta_SCGen_old");
    hist_infile_SCVal->GetObject("h_Eta_SCGen_new", h_Eta_SCGen_new);

    hist_infile_SCRef ->GetObject("h_Eta_SCGen_new_dr", h_Eta_SCGen_old_dr);
    h_Eta_SCGen_old_dr->SetName("h_Eta_SCGen_old_dr");
    hist_infile_SCVal->GetObject("h_Eta_SCGen_new_dr", h_Eta_SCGen_new_dr);

    hist_infile_SCRef ->GetObject("h_Eta_Gen_Denum", h_Eta_Gen_Denum_old);
    h_Eta_Gen_Denum_old->SetName("h_Eta_Gen_Denum_old");
    hist_infile_SCVal->GetObject("h_Eta_Gen_Denum", h_Eta_Gen_Denum_new);

    cout<<"Closed input file"<<endl;
}

//set efficiencies
void setEfficiencies()
{
   eff_SuperCluster_vs_EtaCalo = new TEfficiency(*h_Eta_Calo_old,*h_Eta_Calo_Denum_old);
   eff_DeepSuperCluster_vs_EtaCalo = new TEfficiency(*h_Eta_Calo_new,*h_Eta_Calo_Denum_new);   
   eff_SuperCluster_vs_EtCalo_EB = new TEfficiency(*h_Et_Calo_EB_old,*h_Et_Calo_EB_Denum_old);
   eff_DeepSuperCluster_vs_EtCalo_EB = new TEfficiency(*h_Et_Calo_EB_new,*h_Et_Calo_EB_Denum_new);
   eff_SuperCluster_vs_EtCalo_EE = new TEfficiency(*h_Et_Calo_EE_old,*h_Et_Calo_EE_Denum_old);
   eff_DeepSuperCluster_vs_EtCalo_EE = new TEfficiency(*h_Et_Calo_EE_new,*h_Et_Calo_EE_Denum_new);

}

bool float_equals(float a, float b, float epsilon = 0.00001)
{
    return std::fabs(a - b) < epsilon;
}

double my2sideCrystalBall(double* x, double* par) {

  //a priori we allow for different shape of right and left tail, thus two values of alpha and n 

  Double_t xcur = x[0];
  Double_t alphaL = par[0];
  Double_t nL = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t alphaR = par[5];
  Double_t nR = par[6];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlphaL = fabs((Double_t)alphaL);
  Double_t invAbsAlphaL = 1./absAlphaL;
  Double_t absAlphaR = fabs((Double_t)alphaR);
  Double_t invAbsAlphaR = 1./absAlphaR;

  
  if ( t<-absAlphaL ) {
    //cout<<"checkpoint dscb left"<<endl;
    Double_t AL = TMath::Power(nL*invAbsAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    Double_t BL = nL*invAbsAlphaL - absAlphaL;
    return N*AL*TMath::Power(BL-t,-nL);
  } else if ( t <= absAlphaR )  {
    //cout<<"checkpoint dscb gaussian"<<endl;
    return N*exp(-0.5*t*t);
  } else {
    //cout<<"checkpoint dscb right"<<endl;
    Double_t AR = TMath::Power(nR*invAbsAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    Double_t BR = nR*invAbsAlphaR - absAlphaR;
    return N*AR*TMath::Power(BR+t,-nR);
  }

}

double mycruijff(double* x, double* par) 
{
  Double_t m = x[0];
  Double_t m0 = par[0];
  Double_t sigmaL = par[1];
  Double_t sigmaR = par[2];
  Double_t alphaL = par[3];
  Double_t alphaR = par[4];
  Double_t N = par[5];
  Double_t dx =  (m-m0) ;
  Double_t sigma = dx<0 ? sigmaL: sigmaR ;
  Double_t alpha = dx<0 ? alphaL: alphaR ;
  Double_t f = 2*sigma*sigma + alpha*dx*dx ; 
  
  return N*exp(-dx*dx/f) ;
}

TF1* makeCruijffFit(TH1* hist, float xmin, float xmax, float meanSet=1., float sigmaLSet=0.1, float sigmaRSet=0.1, float alphaLSet=0.1, float alphaRSet=0.1)
{
  //hist->Scale(1./hist->GetEntries()); 
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;
  float nrEntries = hist->Integral();
  RooRealVar  nsig("N_{S}", "#signal events", nrEntries, nrEntries*0.1, nrEntries*10.);
  RooRealVar mean( "#DeltaE", "mean_{cb}", meanSet ,meanSet-3.*hist->GetMean(),meanSet+3.*hist->GetMean(),"");
  RooRealVar sigmaL("#sigma_{L}","#sigma_{L}", sigmaLSet, 0., 0.5);
  RooRealVar sigmaR("#sigma_{R}","#sigma_{R}", sigmaRSet, 0., 0.5);
  RooRealVar alphaL( "alpha_{L}", "alpha_{L}", alphaLSet, 0., 20.);
  RooRealVar alphaR( "alpha_{R}", "alpha_{R}", alphaRSet, 0., 20.);
  CruijffPdf cruijff("cruijff","cruijff",res,mean,sigmaL,sigmaR,alphaL,alphaR);

  RooDataHist data("res","E^{reco}/E^{gen}",res,hist);
   
  RooAddPdf model("model", "model", RooArgList(cruijff), RooArgList(nsig));
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1),RooFit::PrintEvalErrors(-1),RooFit::SumW2Error(true));
  auto fitRes = model.fitTo(data,RooFit::Optimize(1),RooFit::Timer(0),RooFit::PrintEvalErrors(-1),RooFit::Save(1));
  if( (sigmaL.getValV()<=0.001 || sigmaR.getValV()<=0.001) || fitRes->edm()>10){
    std::cout <<" fit status "<<fitRes->status()<<" sigma "<<sigmaL.getValV()<<" "<<sigmaR.getValV()<<" nrsig "<<nsig.getValV()<<" edm "<<fitRes->edm()<<std::endl;
    double maxSigma = std::max(sigmaL.getValV(),sigmaR.getValV());
    sigmaL.setVal(maxSigma);
    sigmaR.setVal(maxSigma);
    fitRes = model.fitTo(data,RooFit::Optimize(1),RooFit::Timer(0),RooFit::PrintEvalErrors(-1),RooFit::Save(1));
    std::cout <<"trying again "<<fitRes->status()<<" sigma "<<sigmaL.getValV()<<" "<<sigmaR.getValV()<<" nrsig "<<nsig.getValV()<<" edm "<<fitRes->edm()<<std::endl;
  }
  
  double errors[6] ={(mean.getAsymErrorHi()+mean.getAsymErrorLo())/2., (sigmaL.getAsymErrorHi()+sigmaL.getAsymErrorLo())/2., (sigmaR.getAsymErrorHi()+sigmaR.getAsymErrorLo())/2., (alphaL.getAsymErrorHi()+alphaL.getAsymErrorLo())/2., (alphaR.getAsymErrorHi()+alphaR.getAsymErrorLo())/2., 0.};
  TF1* Cruijff = new TF1("Cruijff",&mycruijff,xmin,xmax,6);
  Cruijff->SetParameters(mean.getVal(), sigmaL.getVal(), sigmaR.getVal(), alphaL.getVal(), alphaR.getVal(), hist->GetMaximum());
  Cruijff->SetParErrors(errors);

  return Cruijff;
}

TF1* makeDoubleCBFit(TH1* hist,float xmin,float xmax, float meanSet=1., float sigmaSet=0.1, float alpha1Set=0.1, float n1Set=2., float alpha2Set=0.1, float n2Set=2.)
{
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;

  float nrEntries = hist->Integral();
  RooRealVar nsig("N_{S}", "#signal events", nrEntries,nrEntries*0.5,nrEntries*2.);
  RooRealVar mean( "#DeltaE", "mean_{cb}", meanSet ,meanSet-3.*hist->GetMean(),meanSet+3.*hist->GetMean(),"");
  RooRealVar cbSigma("#sigma_{CB}","CB Width", sigmaSet, 0., 0.5);
  RooRealVar alpha1( "alpha_{1}", "alpha_{1}", alpha1Set, 0., 20.);
  RooRealVar alpha2( "alpha_{2}", "alpha_{2}", alpha2Set,0.,20.);
  RooRealVar n1( "n_{1}", "n_{1}", n1Set,0.,5000.);
  RooRealVar n2( "n_{2}", "n_{2}", n2Set,0.,5000.);

  DoubleCBPdf doubleCB("doubleCB","doubleCB",res,mean,cbSigma,alpha1,n1,alpha2,n2);
  RooAddPdf model("model", "model", RooArgList(doubleCB), RooArgList(nsig));

  RooDataHist data("res","E^{reco}/E^{gen}",res,hist);
 
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));
  model.fitTo(data,RooFit::Optimize(1),RooFit::Timer(0),RooFit::PrintEvalErrors(-1),RooFit::Save(1)); 
  
  double errors[7] ={(mean.getAsymErrorHi()+mean.getAsymErrorLo())/2., (cbSigma.getAsymErrorHi()+cbSigma.getAsymErrorLo())/2., (alpha1.getAsymErrorHi()+alpha1.getAsymErrorLo())/2., (n1.getAsymErrorHi()+n1.getAsymErrorLo())/2., (alpha2.getAsymErrorHi()+alpha2.getAsymErrorLo())/2., (n2.getAsymErrorHi()+n2.getAsymErrorLo())/2., 0.};
   
  hist->GetXaxis()->SetRangeUser(0.2,1.8); 
  double maximum = hist->GetMaximum();
  hist->GetXaxis()->SetRangeUser(0.,2.); 
  TF1* DoubleCB = new TF1("DoubleCB",&my2sideCrystalBall,xmin,xmax,7);
  DoubleCB->SetParameters(mean.getVal(), cbSigma.getVal(), alpha1.getVal(), n1.getVal(), alpha2.getVal(), n2.getVal(), maximum);
  DoubleCB->SetParErrors(errors);

  return DoubleCB;
}

TF1* fitHisto(TH1* hist, std::string fitFunction_="cruijff")
{
   hist->GetXaxis()->SetRangeUser(0.2,1.8); 
   int binmax = hist->GetMaximumBin(); 
   hist->GetXaxis()->SetRangeUser(0.,2.); 
   double xMAX = hist->GetXaxis()->GetBinCenter(binmax);  

   TF1* fit1;
   TF1* fit2;
   TF1* fit3; 
   std::cout << "fitHisto: " << hist->GetName() << std::endl;
   if(fitFunction_=="doubleCB"){ 
      fit1 = makeDoubleCBFit(hist,xMAX-1.5*hist->GetRMS(),xMAX+1.5*hist->GetRMS());
      fit2 = makeDoubleCBFit(hist, fit1->GetParameter(0)-5.*fit1->GetParameter(1), fit1->GetParameter(0)+5.*fit1->GetParameter(1), fit1->GetParameter(0), fit1->GetParameter(1), fit1->GetParameter(2), fit1->GetParameter(3), fit1->GetParameter(4), fit1->GetParameter(5));
      fit3 = makeDoubleCBFit(hist, fit2->GetParameter(0)-5.*fit2->GetParameter(1), fit2->GetParameter(0)+5.*fit2->GetParameter(1), fit2->GetParameter(0), fit2->GetParameter(1), fit2->GetParameter(2), fit2->GetParameter(3), fit2->GetParameter(4), fit2->GetParameter(5));
      return fit3;
   }else{
      float sigma = 0;
      fit1 = makeCruijffFit(hist,xMAX-1.5*hist->GetRMS(),xMAX+1.5*hist->GetRMS(), 1., hist->GetRMS()/2., hist->GetRMS()/2., 0.1, 0.1);
      sigma = (fit1->GetParameter(1)+fit1->GetParameter(2))/2.;
      fit2 = makeCruijffFit(hist, fit1->GetParameter(0)-5.*sigma, fit1->GetParameter(0)+5.*sigma, fit1->GetParameter(0), fit1->GetParameter(1), fit1->GetParameter(2), fit1->GetParameter(3), fit1->GetParameter(4)); 
      sigma = (fit2->GetParameter(1)+fit2->GetParameter(2))/2.;
      fit3 = makeCruijffFit(hist, fit2->GetParameter(0)-3.*sigma, fit2->GetParameter(0)+3.*sigma, fit2->GetParameter(0), fit2->GetParameter(1), fit2->GetParameter(2), fit2->GetParameter(3), fit2->GetParameter(4));  
      return fit3;
   }
}

void drawHistFunc(TH1F* hist, TF1* func, std::string x_label, std::string Name)
{
   gStyle->SetOptStat(0000); 
   hist->SetMaximum(hist->GetMaximum()*1.05);
   hist->SetLineColor(kRed+1);
   hist->SetMarkerColor(kRed+1);
   hist->SetLineWidth(2);
   hist->GetXaxis()->SetTitle(x_label.c_str());
   hist->GetXaxis()->SetRangeUser(hist->GetMean()-1.,hist->GetMean()+1.);   

   func->SetLineColor(kBlue+1);

   TLegend* legend = new TLegend(0.57, 0.77, 0.72, 0.89);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);
   legend -> AddEntry(hist,std::string("#mu = "+to_string(func->GetParameter(0))+" +/- "+to_string(func->GetParError(0))).c_str(),"");
   legend -> AddEntry(hist,std::string("#sigma = "+to_string((func->GetParameter(2)+func->GetParameter(1))/2.)+" +/- "+to_string(0.5*sqrt(func->GetParError(2)*func->GetParError(2) +func->GetParError(1)*func->GetParError(1)))).c_str(),"");
     
   TCanvas* c = new TCanvas();
   hist->Draw("E");    
   func->Draw("L,same");
   legend -> Draw("same");  
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf"); 
   gStyle->SetOptStat(1110); 
   
}

std::pair<float,float> computeRange(TGraphErrors* graph)
{
   std::vector<double> y_points;
   double x,y;
   for(int i = 0; i<graph->GetN(); i++){ 
       graph->GetPoint(i,x,y);   
       if(y>=0.) y_points.push_back(y); 
   } 
   std::sort(y_points.begin(),y_points.end());
   if(y_points.size() == 0) return make_pair(0.,2.); 
   else return make_pair(y_points.at(0),y_points.at(y_points.size()-1)); 
}

double computeMean(TH1 * hist, int imin, int imax)
{
   if(imin<1) imin = 1;
   if(imax>hist->GetNbinsX()) imax = hist->GetNbinsX();
   
   double val = 0.;
   double total = 0.;
   for(int ibin=imin; ibin<imax+1; ibin++){
       val+=hist->GetXaxis()->GetBinCenter(ibin) *hist->GetBinContent(ibin);
       total+=hist->GetBinContent(ibin); 
   } 

   if(total==0) return -1.;
   else return val/total;   
}

std::pair<double,double> computeEffectiveSigma(TH1 * hist)
{
    TAxis *xaxis = hist->GetXaxis();
    int nb = xaxis->GetNbins();
    if(nb < 10) {
       cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
       return std::make_pair(-1.,-1.);
    }

    double bwid = xaxis->GetBinWidth(1);
    if(bwid == 0) {
       cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
       return std::make_pair(-1.,-1.);
    }

    //double xmax = xaxis->GetXmax();
    double xmin = xaxis->GetXmin();
    double ave = hist->GetMean();
    double rms = hist->GetRMS();

    double total=0.;
    for(int i=0; i<nb+2; i++) {
        total+=hist->GetBinContent(i);
    }
    
    int ierr=0;
    int ismin=999;

    double rlim=0.68269*total;
    int nrms=rms/(bwid);    // Set scan size to +/- rms
    if(nrms > nb/10) nrms=nb/10; // Could be tuned...

    double widmin=9999999.;
    int jbin = nb;
    int kbin = 1;  
     

    for(int iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre

        int ibm=(ave-xmin)/bwid+1+iscan;
        double x=(ibm-0.5)*bwid+xmin;
        double xj=x;
        double xk=x;
        int jbm=ibm;
        int kbm=ibm;
        double bin=hist->GetBinContent(ibm);
        total=bin;

        for(int j=1;j<nb;j++){
            if(jbm < nb) {
               jbm++;
               xj+=bwid;
               bin=hist->GetBinContent(jbm);
               total+=bin;
               jbin = jbm; 
               if(total > rlim) break;
            }else ierr=1;

            if(kbm > 0) {
               kbm--;
               xk-=bwid;
               bin=hist->GetBinContent(kbm);
               total+=bin;
               kbin = kbm;  
               if(total > rlim) break;
            }else ierr=1;
        }

        double dxf=(total-rlim)*bwid/bin;
        double wid=(xj-xk+bwid-dxf)*0.5;

        if(wid < widmin){
           widmin=wid;
           ismin=iscan;
        }
    }

    if(ismin == nrms || ismin == -nrms) ierr=3;
    if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;

    //std::cout << "EFFECTIVE SIGMA: " << kbin << " - " << jbin << " - " << nb << " - " << (float)hist->Integral(kbin,jbin)/(float)hist->Integral() << std::endl;
    return std::make_pair(computeMean(hist,kbin,jbin),widmin);
}

std::pair<TGraphErrors*,TGraphErrors*> makeFitProfile(std::vector<TH1F*>* vecHist,double min, double max, std::string xTitle, std::string fitFunction_="cruijff", bool doEffective = false)
{
   TF1* doubleCB;
   int nBins = vecHist->size();
   double x[nBins], x_error[nBins], y_mean[nBins], y_meanError[nBins], y_sigma[nBins], y_sigmaError[nBins];
   float delta = fabs(max-min)/nBins; 

   for(unsigned int iBin=0; iBin<vecHist->size(); iBin++)
   {
       x[iBin] = min + (iBin+0.5)*delta; 
       x_error[iBin] = 0.;
       
       if(vecHist->at(iBin)->Integral()>20){
          if(!doEffective){
             doubleCB = fitHisto(vecHist->at(iBin), fitFunction_);
             y_mean[iBin] = doubleCB->GetParameter(0);
             y_meanError[iBin] = doubleCB->GetParError(0);   
             if(y_meanError[iBin]>0.1){
                y_mean[iBin] = -1.;
                y_meanError[iBin] = 0.;     
             }  
             y_sigma[iBin] = (doubleCB->GetParameter(1)+doubleCB->GetParameter(2))/2.;
             y_sigmaError[iBin] = sqrt(doubleCB->GetParError(1)*doubleCB->GetParError(1)+doubleCB->GetParError(2)*doubleCB->GetParError(2))/2.;
             if(y_sigmaError[iBin]>0.03){
                y_sigma[iBin] = -1.;
                y_sigmaError[iBin] = 0.;     
             }    
             if(y_meanError[iBin]<0.1 && y_sigmaError[iBin]<0.03) drawHistFunc(vecHist->at(iBin),doubleCB, std::string(""), std::string(vecHist->at(iBin)->GetName())); 
          }else{
            std::pair<double,double> effective = computeEffectiveSigma(vecHist->at(iBin));
            y_mean[iBin] = effective.first;
            y_meanError[iBin] = 0.;
            y_sigma[iBin] =  effective.second;
            y_sigmaError[iBin] = 0.; 
            if(effective.first>10.) y_mean[iBin] = -1.; 
            if(effective.second>10.) y_sigma[iBin] = -1.; 
          }   
       }else{
          y_mean[iBin] = -1.;
          y_meanError[iBin] = 0.;
          y_sigma[iBin] = -1.;
          y_sigmaError[iBin] = 0.;
       } 
   }
         
   TGraphErrors* gr_Mean = new TGraphErrors(nBins,x,y_mean,x_error,y_meanError);
   TGraphErrors* gr_Sigma = new TGraphErrors(nBins,x,y_sigma,x_error,y_sigmaError);;

   //delete doubleCB;
   return std::make_pair(gr_Mean,gr_Sigma);
   
}

void drawHisto(TH1F* h_old, TH1F* h_new, std::string x_label, std::string drawType, std::string Name, bool log, bool fit=false, std::string fitFunc_="cruijff", std::string refLegend="Mustache", std::string valLegend="DeepSC", double rat_min=0.7, double rat_max=1.3)
{

   //h_old->Scale(1./h_old->GetEntries());
   //h_new->Scale(1./h_new->GetEntries());

   TFitResultPtr frp_old;
   TFitResultPtr frp_new;
   TF1* doubleCB_old = new TF1();
   TF1* doubleCB_new = new TF1();

   if(fit){
      TH1F* h_old_clone = (TH1F*)h_old->Clone();  
      TH1F* h_new_clone = (TH1F*)h_new->Clone();    

      doubleCB_old = fitHisto(h_old_clone, fitFunc_);
      doubleCB_old->SetLineColor(kRed+1);

      doubleCB_new = fitHisto(h_new_clone, fitFunc_);
      doubleCB_new->SetLineColor(kBlue+1);
   }

   //h_old->SetTitle("");
   //h_new->SetTitle("");

   std::vector<float> maxima;
   maxima.resize(4);
   maxima[0] = h_old->GetMaximum(); 
   maxima[1] = h_new->GetMaximum();
   if(fit) maxima[2] = doubleCB_old->GetMaximum();
   if(fit) maxima[3] = doubleCB_new->GetMaximum(); 
   std::sort(maxima.begin(),maxima.end());  
   
   h_old->SetMaximum(maxima.at(maxima.size()-1)*1.05);
   h_old->SetLineColor(kRed+1);
   h_old->SetMarkerColor(kRed+1);
   h_old->SetLineWidth(2);
   h_old->GetXaxis()->SetTitle(x_label.c_str());

   h_new->SetLineColor(kBlue+1);
   h_new->SetMarkerColor(kBlue+1);
   h_new->SetLineWidth(2);

   /*float maximum = -1.;
   if( h_old -> GetMaximum() > h_new -> GetMaximum()) maximum=h_old -> GetMaximum();
   else maximum=h_new -> GetMaximum();  
   h_old -> SetMaximum( 1.1*maximum );*/
    
   TLegend* legend = new TLegend(0.60, 0.82, 0.75, 0.94);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);

   TPaveStats* st_old = new TPaveStats();
   TPaveStats* st_new = new TPaveStats();
   TPaveStats* st_ratio = new TPaveStats();
   
   TCanvas* c = new TCanvas();
 
   TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.36,1.00,1.00);
   TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.36);
   
   cUp->SetBottomMargin(0.01); 
   cDown->SetTopMargin(0.01); 
   cDown->SetBottomMargin(0.2); 
    
   cUp->Draw();
   if(log) cUp->SetLogy();
   cDown->Draw();
     
   cUp->cd();
   if(!fit){
      TLegend* legend = new TLegend(0.60, 0.82, 0.75, 0.94);
      legend -> SetFillColor(kWhite);
      legend -> SetFillStyle(1000);
      legend -> SetLineWidth(0);
      legend -> SetLineColor(kWhite);
      legend -> SetTextFont(42);  
      legend -> SetTextSize(0.04); 
      legend -> AddEntry(h_old,refLegend.c_str(),"L");
      legend -> AddEntry(h_new,valLegend.c_str(),"L");
      h_old->Draw(drawType.c_str());
      gPad -> Update();
      st_old= (TPaveStats*)(h_old->GetListOfFunctions()->FindObject("stats"));
      st_old->SetX1NDC(0.82); //new x start position
      st_old->SetX2NDC(0.99); //new x end position
      st_old->SetY1NDC(0.82); //new y start position
      st_old->SetY2NDC(0.94); //new y end position
      st_old->SetTextColor(kRed+1);
      st_old->Draw("sames");
      h_new->Draw(std::string(drawType+",sames").c_str());
      gPad -> Update();
      st_new= (TPaveStats*)(h_new->GetListOfFunctions()->FindObject("stats"));
      st_new->SetX1NDC(0.82); //new x start position
      st_new->SetX2NDC(0.99); //new x end position
      st_new->SetY1NDC(0.68); //new y start position
      st_new->SetY2NDC(0.80); //new y end position
      st_new->SetTextColor(kBlue+1);
      st_new->Draw("sames");
      legend -> Draw("same");
   }else{
      TLegend* legend = new TLegend(0.57, 0.77, 0.72, 0.89);
      legend -> SetFillColor(kWhite);
      legend -> SetFillStyle(1000);
      legend -> SetLineWidth(0);
      legend -> SetLineColor(kWhite);
      legend -> SetTextFont(42);  
      legend -> SetTextSize(0.05);
      legend -> AddEntry(h_old,std::string(refLegend+": "+to_string(doubleCB_old->GetParameter(0))+" +/- "+to_string((doubleCB_old->GetParameter(2)+doubleCB_old->GetParameter(1))/2.)).c_str(),"L");
      legend -> AddEntry(h_new,std::string(valLegend+": "+to_string(doubleCB_new->GetParameter(0))+" +/- "+to_string((doubleCB_new->GetParameter(2)+doubleCB_new->GetParameter(1))/2.)).c_str(),"L");
      h_old->Draw(std::string(drawType).c_str());    
      h_new->Draw(std::string(drawType+",sames").c_str());
      doubleCB_old->Draw("L,same");
      doubleCB_new->Draw("L,same");
      legend -> Draw("same");   
      gPad -> Update();
      st_old= (TPaveStats*)(h_old->GetListOfFunctions()->FindObject("stats"));
      st_old->SetX1NDC(0.); //new x start position
      st_old->SetX2NDC(0.); //new x end position
      st_old->SetY1NDC(0.); //new y start position
      st_old->SetY2NDC(0.); //new y end position
      st_old->SetTextColor(kBlack);
      st_old->Draw("sames");
      gPad -> Update();
      st_new= (TPaveStats*)(h_new->GetListOfFunctions()->FindObject("stats"));
      st_new->SetX1NDC(0.); //new x start position
      st_new->SetX2NDC(0.); //new x end position
      st_new->SetY1NDC(0.); //new y start position
      st_new->SetY2NDC(0.); //new y end position
      st_new->SetTextColor(kBlack);
      st_new->Draw("sames");
   }
   cDown->cd();
    
   TH1F* histo_ratio=(TH1F*)h_new->Clone("histo_ratio");
   if(histo_ratio->GetSumw2N()<=0) histo_ratio->Sumw2();
   histo_ratio->Divide(h_old);
    
   histo_ratio -> GetXaxis() -> SetTitle(x_label.c_str());
   histo_ratio -> GetYaxis() -> SetTitle(std::string(valLegend+"/"+refLegend).c_str());
   histo_ratio -> SetMaximum(rat_max);
   histo_ratio -> SetMinimum(rat_min);
   histo_ratio -> SetMarkerColor(kBlack);
   histo_ratio -> SetMarkerSize(0.5);
   histo_ratio ->SetTitle("");
   histo_ratio -> GetXaxis() -> SetLabelSize(0.07);
   histo_ratio -> GetYaxis() -> SetLabelSize(0.07);
   histo_ratio -> GetXaxis() -> SetTitleSize(0.07);
   histo_ratio -> GetYaxis() -> SetTitleSize(0.07);
   histo_ratio -> GetYaxis() -> SetTitleOffset(0.7);
   histo_ratio -> Draw("e");
   TF1* f_const = new TF1("f_1", "[0]",histo_ratio->GetBinCenter(1)-histo_ratio->GetBinWidth(1)/2, histo_ratio->GetBinCenter(histo_ratio->GetNbinsX())+histo_ratio->GetBinWidth(histo_ratio->GetNbinsX())/2);
   f_const -> FixParameter(0,1);
   f_const -> SetLineColor(kRed);
   f_const -> SetLineWidth(2);
   f_const -> Draw("same");
    
   gPad -> Update();
   st_ratio= (TPaveStats*)(histo_ratio->GetListOfFunctions()->FindObject("stats"));
   st_ratio->SetX1NDC(0.); //new x start position
   st_ratio->SetX2NDC(0.); //new x end position
   st_ratio->SetY1NDC(0.); //new y start position
   st_ratio->SetY2NDC(0.); //new y end position
   st_ratio->SetTextColor(kBlack);
   st_ratio->Draw("sames");

   if(!log){
      c->SaveAs(std::string(Name+".png").c_str(),"png");
      c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	
   }else{
      c->SaveAs(std::string(Name+"_logY.png").c_str(),"png");
      c->SaveAs(std::string(Name+"_logY.pdf").c_str(),"pdf");
   }

   delete histo_ratio;
}

void drawH2(TH2F* h2, std::string xtitle, std::string ytitle, string ztitle, std::string Name, bool log = false, float z_min=-1, float z_max=-1)
{
   gStyle->SetOptStat(0000); 

   h2->GetXaxis()->SetTitle(xtitle.c_str()); 
   h2->GetYaxis()->SetTitle(ytitle.c_str()); 
   h2->GetZaxis()->SetTitle(ztitle.c_str());
   if(z_min!=-1. && z_max!=-1.) h2->GetZaxis()->SetRangeUser(z_min,z_max);

   TCanvas* c = new TCanvas();
   if(log) c->SetLogz(); 
   h2->Draw("COLZ");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");

   gStyle->SetOptStat(1110); 
}

TGraphErrors* makeRatioGraph(TGraphErrors* gr_SuperCluster, TGraphErrors* gr_DeepSuperCluster)
{
   TGraphErrors* gr_ratio = new TGraphErrors();
   double x,y_ref,y_val,yErr_ref,yErr_val;
   for(int i = 0; i<gr_SuperCluster->GetN(); i++){ 
       gr_SuperCluster->GetPoint(i,x,y_ref); 
       yErr_ref = gr_SuperCluster->GetErrorY(i); 
       gr_DeepSuperCluster->GetPoint(i,x,y_val);  
       yErr_val = gr_DeepSuperCluster->GetErrorY(i);     
       if(y_ref>0. && y_val>0.){
          gr_ratio->SetPoint(i,x,y_val/y_ref);
          gr_ratio->SetPointError(i,0.,y_val/y_ref*sqrt((yErr_ref/y_ref)*(yErr_ref/y_ref)+(yErr_val/y_val)*(yErr_val/y_val)));
       }else{
          gr_ratio->SetPoint(i,x,-1.);
          gr_ratio->SetPointError(i,0.,0.);
       }   
   } 
   return gr_ratio; 
}

void drawGraph(TGraphErrors* gr_SuperCluster, TGraphErrors* gr_DeepSuperCluster, std::string xtitle, std::string ytitle, std::string Name, std::string refLegend="Mustache", std::string valLegend="DeepSC",float rat_min=0.5, float rat_max=1.5,float y_min=-1., float y_max=-1.)
{ 
   gStyle->SetOptStat(0000);  
   gr_SuperCluster->SetTitle(Name.c_str());
   gr_SuperCluster->SetLineColor(kRed+1);
   gr_SuperCluster->SetMarkerStyle(20);
   gr_SuperCluster->SetMarkerSize(0.5);
   gr_SuperCluster->SetMarkerColor(kRed+1);
   gr_SuperCluster->SetLineWidth(2);
   gr_SuperCluster->GetXaxis()->SetTitle(xtitle.c_str()); 
   gr_SuperCluster->GetYaxis()->SetTitle(ytitle.c_str()); 
   
   gr_DeepSuperCluster->SetLineColor(kBlue+1);
   gr_DeepSuperCluster->SetMarkerStyle(20);
   gr_DeepSuperCluster->SetMarkerSize(0.5);
   gr_DeepSuperCluster->SetMarkerColor(kBlue+1);
   gr_DeepSuperCluster->SetLineWidth(2);

   float min = y_min;
   float max = y_max;
   if(y_min<0. || y_max<0.){
      min = 0.9*computeRange(gr_SuperCluster).first; 
      max = 1.1*computeRange(gr_SuperCluster).second;
   }
   gr_SuperCluster->GetYaxis()->SetRangeUser(min,max);  

   TLegend* legend = new TLegend(0.799, 0.77, 0.999, 0.95);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.03);
   legend -> AddEntry(gr_SuperCluster,refLegend.c_str(),"L");
   legend -> AddEntry(gr_DeepSuperCluster,valLegend.c_str(),"L");

   TCanvas* c = new TCanvas();

   TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.36,1.00,1.00);
   TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.36);
   
   cUp->SetBottomMargin(0.01); 
   cDown->SetTopMargin(0.01); 
   cDown->SetBottomMargin(0.2); 
    
   cUp->Draw();
   cDown->Draw();
     
   cUp->cd();
   gr_SuperCluster->Draw("AP");
   gr_DeepSuperCluster->Draw("P, same");
   legend -> Draw("same");

   cDown->cd();
    
   TGraphErrors* gr_ratio = makeRatioGraph(gr_SuperCluster,gr_DeepSuperCluster); 
   std::pair<float,float> ratioRange = computeRange(gr_ratio);
   min = rat_min;//0.65; 
   max = rat_max;//1.35; 
   //if(ratioRange.first*0.9> 0.65) min = ratioRange.first*0.9; 
   //if(ratioRange.second*1.1<1.35) max = ratioRange.second*1.1; 
   gr_ratio->GetYaxis() -> SetRangeUser(min,max);   
   gr_ratio->GetXaxis()->SetTitle(xtitle.c_str()); 
   gr_ratio->GetYaxis() -> SetTitle(std::string(valLegend+"/"+refLegend).c_str());
   gr_ratio->SetMarkerColor(kBlack);
   gr_ratio->SetMarkerStyle(20);
   gr_ratio->SetMarkerSize(0.5);
   gr_ratio->SetTitle("");
   gr_ratio->GetXaxis() -> SetLabelSize(0.07);
   gr_ratio->GetYaxis() -> SetLabelSize(0.07);
   gr_ratio->GetXaxis() -> SetTitleSize(0.07);
   gr_ratio->GetYaxis() -> SetTitleSize(0.07);
   gr_ratio->GetYaxis() -> SetTitleOffset(0.7);
   gr_ratio->Draw("AP");
   TF1* f_const = new TF1("f_1", "[0]",gr_ratio->GetXaxis()->GetBinCenter(1)-gr_ratio->GetXaxis()->GetBinWidth(1)/2, gr_ratio->GetXaxis()->GetBinCenter(gr_ratio->GetXaxis()->GetNbins())+gr_ratio->GetXaxis()->GetBinWidth(gr_ratio->GetXaxis()->GetNbins())/2);
   f_const -> FixParameter(0,1);
   f_const -> SetLineColor(kRed);
   //f_const -> SetLineWidth(2);
   f_const -> Draw("same");

   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	

   gStyle->SetOptStat(1110);
}

void drawEfficiency(TEfficiency* eff_SuperCluster, TEfficiency* eff_DeepSuperCluster, std::string xtitle, std::string Name, std::string refLegend="Mustache", std::string valLegend="DeepSC")
{
   eff_SuperCluster->SetLineColor(kRed+1);
   eff_SuperCluster->SetLineWidth(2);
   eff_SuperCluster->SetTitle(std::string(Name+"; "+xtitle+" ; Efficiency").c_str()); 
   
   eff_DeepSuperCluster->SetLineColor(kBlue+1);
   eff_DeepSuperCluster->SetLineWidth(2);

   TLegend* legend = new TLegend(0.365, 0.12, 0.635, 0.34);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);
   legend -> AddEntry(eff_SuperCluster,refLegend.c_str(),"L");
   legend -> AddEntry(eff_DeepSuperCluster,valLegend.c_str(),"L");

   TCanvas* c = new TCanvas();
   eff_SuperCluster->Draw("APL");
   eff_DeepSuperCluster->Draw("PL, same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	
}

void drawProfile(TProfile* prof_SuperCluster, TProfile* prof_DeepSuperCluster, std::string xtitle, std::string ytitle, std::string Name, std::string refLegend="Mustache", std::string valLegend="DeepSC")
{ 
   gStyle->SetOptStat(0000);  
   prof_SuperCluster->SetLineColor(kRed+1);
   prof_SuperCluster->SetLineWidth(2);
   prof_SuperCluster->GetXaxis()->SetTitle(xtitle.c_str()); 
   prof_SuperCluster->GetYaxis()->SetTitle(ytitle.c_str());   
   
   prof_DeepSuperCluster->SetLineColor(kBlue+1);
   prof_DeepSuperCluster->SetLineWidth(2);

   float min = prof_SuperCluster->GetMinimum();
   if(prof_DeepSuperCluster->GetMinimum()<min) min = prof_DeepSuperCluster->GetMinimum(); 
   float max = prof_SuperCluster->GetMaximum();
   if(prof_DeepSuperCluster->GetMaximum()>max) max = prof_DeepSuperCluster->GetMaximum(); 
   min = min*0.9;
   max = max*1.1;

   prof_SuperCluster->GetYaxis()->SetRangeUser(min,max);

   TLegend* legend = new TLegend(0.365, 0.72, 0.635, 0.90);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.03);
   legend -> AddEntry(prof_SuperCluster,refLegend.c_str(),"L");
   legend -> AddEntry(prof_DeepSuperCluster,valLegend.c_str(),"L");

   TCanvas* c = new TCanvas();
   prof_SuperCluster->Draw("PL");
   prof_DeepSuperCluster->Draw("PL, same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	

   gStyle->SetOptStat(1110);
}

void drawTH2Fits(std::vector<std::vector<TH1F*>> EoEtrue_vs_seedEt_seedEta_old, std::vector<std::vector<TH1F*>> EoEtrue_vs_seedEt_seedEta_new, std::string superClusterRef_, string superClusterVal_, string fitFunction_,  bool doEffective=false)
{
   TF1* doubleCB;
   for(unsigned int iBin=0; iBin<EoEtrue_vs_seedEt_seedEta_old.size(); iBin++){
       for(unsigned int jBin=0; jBin<EoEtrue_vs_seedEt_seedEta_old[iBin].size(); jBin++) 
       {
           if(!doEffective){
              doubleCB = fitHisto(EoEtrue_vs_seedEt_seedEta_old[iBin][jBin], fitFunction_);
              float y_mean = doubleCB->GetParameter(0);
              float y_meanError = doubleCB->GetParError(0);   
              if(y_meanError>0.1){
                 y_mean = -1.;
                 y_meanError = 0.;     
              }  
              float y_sigma = (doubleCB->GetParameter(1)+doubleCB->GetParameter(2))/2.;
              float y_sigmaError = sqrt(doubleCB->GetParError(1)*doubleCB->GetParError(1)+doubleCB->GetParError(2)*doubleCB->GetParError(2))/2.;
              if(y_sigmaError>0.03){
                 y_sigma = -1.;
                 y_sigmaError = 0.;     
              }    
              if(y_meanError<0.1 && y_sigmaError<0.03) drawHistFunc(EoEtrue_vs_seedEt_seedEta_old[iBin][jBin],doubleCB, std::string(""), std::string(EoEtrue_vs_seedEt_seedEta_old[iBin][jBin]->GetName())); 
              if(y_meanError<0.1 && y_sigmaError<0.03 && y_mean<5. && y_sigma<5. && y_mean>0. && y_sigma>0.){
                 h2_EoEtrue_Mean_old->SetBinContent(iBin+1,jBin+1,y_mean); 
                 h2_EoEtrue_Resolution_old->SetBinContent(iBin+1,jBin+1,y_sigma); 
              }
           }else{
              std::pair<double,double> effective = computeEffectiveSigma(EoEtrue_vs_seedEt_seedEta_old[iBin][jBin]);
              float y_mean = effective.first;
              float y_sigma =  effective.second;
              if(effective.first>10.) y_mean = -1.; 
              if(effective.second>10.) y_sigma = -1.; 
              if(y_mean>0. && y_sigma>0.  && y_mean<5. && y_sigma>5.){
                 h2_EoEtrue_Mean_Effective_old->SetBinContent(iBin+1,jBin+1,y_mean); 
                 h2_EoEtrue_Resolution_Effective_old->SetBinContent(iBin+1,jBin+1,y_sigma); 
              } 
           } 
       }
   }

   for(unsigned int iBin=0; iBin<EoEtrue_vs_seedEt_seedEta_new.size(); iBin++){
       for(unsigned int jBin=0; jBin<EoEtrue_vs_seedEt_seedEta_new[iBin].size(); jBin++) 
       {
           if(!doEffective){
              doubleCB = fitHisto(EoEtrue_vs_seedEt_seedEta_new[iBin][jBin], fitFunction_);
              float y_mean = doubleCB->GetParameter(0);
              float y_meanError = doubleCB->GetParError(0);   
              if(y_meanError>0.1){
                 y_mean = -1.;
                 y_meanError = 0.;     
              }  
              float y_sigma = (doubleCB->GetParameter(1)+doubleCB->GetParameter(2))/2.;
              float y_sigmaError = sqrt(doubleCB->GetParError(1)*doubleCB->GetParError(1)+doubleCB->GetParError(2)*doubleCB->GetParError(2))/2.;
              if(y_sigmaError>0.03){
                 y_sigma = -1.;
                 y_sigmaError = 0.;     
              }    
              if(y_meanError<0.1 && y_sigmaError<0.03) drawHistFunc(EoEtrue_vs_seedEt_seedEta_new[iBin][jBin],doubleCB, std::string(""), std::string(EoEtrue_vs_seedEt_seedEta_new[iBin][jBin]->GetName())); 
              if(y_meanError<0.1 && y_sigmaError<0.03 && y_mean<5. && y_sigma<5. && y_mean>0. && y_sigma>0.){
                 h2_EoEtrue_Mean_new->SetBinContent(iBin+1,jBin+1,y_mean); 
                 h2_EoEtrue_Resolution_new->SetBinContent(iBin+1,jBin+1,y_sigma); 
              }
           }else{
              std::pair<double,double> effective = computeEffectiveSigma(EoEtrue_vs_seedEt_seedEta_new[iBin][jBin]);
              float y_mean = effective.first;
              float y_sigma =  effective.second;
              if(effective.first>10.) y_mean = -1.; 
              if(effective.second>10.) y_sigma = -1.; 
              if(y_mean>0. && y_sigma>0. && y_mean<5. && y_sigma>5.){
                 h2_EoEtrue_Mean_Effective_new->SetBinContent(iBin+1,jBin+1,y_mean); 
                 h2_EoEtrue_Resolution_Effective_new->SetBinContent(iBin+1,jBin+1,y_sigma); 
              } 
           } 
       }
   }

   std::cout << "SONO QUI 1" << std::endl;
   TH2F* h2_EoEtrue_Mean_Ratio = (TH2F*)h2_EoEtrue_Mean_old->Clone("h2_EoEtrue_Mean_Ratio"); 
   TH2F* h2_EoEtrue_Mean_Effective_Ratio = (TH2F*)h2_EoEtrue_Mean_old->Clone("h2_EoEtrue_Mean_Effective_Ratio"); 
   TH2F* h2_EoEtrue_Resolution_Ratio = (TH2F*)h2_EoEtrue_Mean_old->Clone("h2_EoEtrue_Resolution_Ratio"); 
   TH2F* h2_EoEtrue_Resolution_Effective_Ratio = (TH2F*)h2_EoEtrue_Mean_old->Clone("h2_EoEtrue_Resolution_Effective_Ratio"); 
  
   h2_EoEtrue_Mean_Ratio->SetTitle("h2_EoEtrue_Mean_Ratio");
   h2_EoEtrue_Mean_Effective_Ratio->SetTitle("h2_EoEtrue_Mean_Effective_Ratio");
   h2_EoEtrue_Resolution_Ratio->SetTitle("h2_EoEtrue_Resolution_Ratio");
   h2_EoEtrue_Resolution_Effective_Ratio->SetTitle("h2_EoEtrue_Resolution_Effective_Ratio");
   h2_EoEtrue_Mean_Ratio->Divide(h2_EoEtrue_Mean_new,h2_EoEtrue_Mean_old);
   h2_EoEtrue_Mean_Effective_Ratio->Divide(h2_EoEtrue_Mean_Effective_new,h2_EoEtrue_Mean_Effective_old);
   h2_EoEtrue_Resolution_Ratio->Divide(h2_EoEtrue_Resolution_new,h2_EoEtrue_Resolution_old);
   h2_EoEtrue_Resolution_Effective_Ratio->Divide(h2_EoEtrue_Resolution_Effective_new,h2_EoEtrue_Resolution_Effective_old);

   std::cout << "SONO QUI 2" << std::endl;
   drawH2(h2_EoEtrue_Mean_old, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Mean_"+superClusterRef_));
   drawH2(h2_EoEtrue_Mean_new, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Mean_"+superClusterVal_));
   drawH2(h2_EoEtrue_Mean_Effective_old, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Mean_Effective_"+superClusterRef_));
   drawH2(h2_EoEtrue_Mean_Effective_new, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Mean_Effective_"+superClusterVal_));
   drawH2(h2_EoEtrue_Resolution_old, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Resolution_"+superClusterRef_));
   drawH2(h2_EoEtrue_Resolution_new, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Resolution_"+superClusterVal_));
   drawH2(h2_EoEtrue_Resolution_Effective_old, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Resolution_Effective_"+superClusterRef_)); 
   drawH2(h2_EoEtrue_Resolution_Effective_new, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Resolution_Effective_"+superClusterVal_));
   std::cout << "SONO QUI 3" << std::endl;
   drawH2(h2_EoEtrue_Mean_Ratio, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Mean_Ratio")); 
   drawH2(h2_EoEtrue_Mean_Effective_Ratio, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Mean_Effective_Ratio")); 
   drawH2(h2_EoEtrue_Resolution_Ratio, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Resolution_Ratio")); 
   drawH2(h2_EoEtrue_Resolution_Effective_Ratio, std::string("seedEt (GeV)"), std::string("seedEta"), string(""), std::string("h2_EoEtrue_Resolution_Effective_Ratio"));  
   std::cout << "SONO QUI 4" << std::endl;
   
}

//draw plots
void drawPlots(std::string fitFunction_, string superClusterRef_, string superClusterVal_)
{
   
   drawEfficiency(eff_SuperCluster_vs_EtaCalo, eff_DeepSuperCluster_vs_EtaCalo, std::string("caloParticle_#eta"), std::string("Efficiency_vs_CaloEta"), superClusterRef_, superClusterVal_); 
   drawEfficiency(eff_SuperCluster_vs_EtCalo_EB, eff_DeepSuperCluster_vs_EtCalo_EB, std::string("caloParticle_Et (GeV)"), std::string("Efficiency_vs_CaloEt_EB"), superClusterRef_, superClusterVal_); 
   drawEfficiency(eff_SuperCluster_vs_EtCalo_EE, eff_DeepSuperCluster_vs_EtCalo_EE, std::string("caloParticle_Et (GeV)"), std::string("Efficiency_vs_CaloEt_EE"), superClusterRef_, superClusterVal_); 

   drawTH2Fits(EoEtrue_vs_seedEt_seedEta_old, EoEtrue_vs_seedEt_seedEta_new, superClusterRef_, superClusterVal_, fitFunction_, true);
   drawTH2Fits(EoEtrue_vs_seedEt_seedEta_old, EoEtrue_vs_seedEt_seedEta_new, superClusterRef_, superClusterVal_, fitFunction_, false);
   
   drawProfile(prof_EoEtrue_vs_Eta_Calo_old, prof_EoEtrue_vs_Eta_Calo_new, std::string("caloParticle_#eta"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEta"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Eta_Seed_old, prof_EoEtrue_vs_Eta_Seed_new, std::string("seed_#eta"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_SeedEta"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Et_Calo_EB_old, prof_EoEtrue_vs_Et_Calo_EB_new, std::string("caloParticle_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEt_EB"), superClusterRef_, superClusterVal_);  
   drawProfile(prof_EoEtrue_vs_Et_Seed_EB_old, prof_EoEtrue_vs_Et_Seed_EB_new, std::string("seed_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_SeedEt_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Energy_Calo_EB_old, prof_EoEtrue_vs_Energy_Calo_EB_new, std::string("caloParticle_Energy (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEnergy_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_nVtx_EB_old, prof_EoEtrue_vs_nVtx_EB_new, std::string("nVtx"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_nVtx_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Rho_EB_old, prof_EoEtrue_vs_Rho_EB_new, std::string("#rho"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_Rho_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Et_Calo_EE_old, prof_EoEtrue_vs_Et_Calo_EE_new, std::string("caloParticle_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEt_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Et_Seed_EE_old, prof_EoEtrue_vs_Et_Seed_EE_new, std::string("seed_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_SeedEt_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Energy_Calo_EE_old, prof_EoEtrue_vs_Energy_Calo_EE_new, std::string("caloParticle_Energy (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEnergy_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_nVtx_EE_old, prof_EoEtrue_vs_nVtx_EE_new, std::string("nVtx"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_nVtx_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Rho_EE_old, prof_EoEtrue_vs_Rho_EE_new, std::string("#rho"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_Rho_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Eta_Gen_old, prof_EoEgen_vs_Eta_Gen_new, std::string("genParticle_#eta"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_GenEta"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Eta_Seed_old, prof_EoEgen_vs_Eta_Seed_new, std::string("seed_#eta"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_SeedEta"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Et_Gen_EB_old, prof_EoEgen_vs_Et_Gen_EB_new, std::string("genParticle_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_GenEt_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Et_Seed_EB_old, prof_EoEgen_vs_Et_Seed_EB_new, std::string("seed_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_SeedEt_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Energy_Gen_EB_old, prof_EoEgen_vs_Energy_Gen_EB_new, std::string("genParticle_Energy (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_GenEnergy_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_nVtx_EB_old, prof_EoEgen_vs_nVtx_EB_new, std::string("nVtx"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_nVtx_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Rho_EB_old, prof_EoEgen_vs_Rho_EB_new, std::string("#rho"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_Rho_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Et_Gen_EE_old, prof_EoEgen_vs_Et_Gen_EE_new, std::string("genParticle_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_GenEt_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Et_Seed_EE_old, prof_EoEgen_vs_Et_Seed_EE_new, std::string("seed_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_SeedEt_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Energy_Gen_EE_old, prof_EoEgen_vs_Energy_Gen_EE_new, std::string("genParticle_Energy (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_GenEnergy_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_nVtx_EE_old, prof_EoEgen_vs_nVtx_EE_new, std::string("nVtx"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_nVtx_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Rho_EE_old, prof_EoEgen_vs_Rho_EE_new, std::string("#rho"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_Rho_EE"), superClusterRef_, superClusterVal_);
   
   
    
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_EtaWidth_EB_old= makeFitProfile(&EoEtrue_vs_EtaWidth_EB_old,-3.,3.,std::string("eta_width"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_EtaWidth_EB_new= makeFitProfile(&EoEtrue_vs_EtaWidth_EB_new,-3.,3.,std::string("eta_width"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_EtaWidth_EE_old= makeFitProfile(&EoEtrue_vs_EtaWidth_EE_old,-3.,3.,std::string("eta_width"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_EtaWidth_EE_new= makeFitProfile(&EoEtrue_vs_EtaWidth_EE_new,-3.,3.,std::string("eta_width"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_PhiWidth_EB_old= makeFitProfile(&EoEtrue_vs_PhiWidth_EB_old,-3.,3.,std::string("phi_width"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_PhiWidth_EB_new= makeFitProfile(&EoEtrue_vs_PhiWidth_EB_new,-3.,3.,std::string("phi_width"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_PhiWidth_EE_old= makeFitProfile(&EoEtrue_vs_PhiWidth_EE_old,-3.,3.,std::string("phi_width"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_PhiWidth_EE_new= makeFitProfile(&EoEtrue_vs_PhiWidth_EE_new,-3.,3.,std::string("phi_width"),fitFunction_);
    
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Eta_old= makeFitProfile(&EoEtrue_vs_Eta_Calo_old,-3.,3.,std::string("caloParticle_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Eta_new= makeFitProfile(&EoEtrue_vs_Eta_Calo_new,-3.,3.,std::string("caloParticle_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEta_old= makeFitProfile(&EoEtrue_vs_Eta_Seed_old,-3.,3.,std::string("seed_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEta_new= makeFitProfile(&EoEtrue_vs_Eta_Seed_new,-3.,3.,std::string("seed_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EB_old= makeFitProfile(&EoEtrue_vs_Et_Calo_EB_old,0.,100.,std::string("caloParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EB_new= makeFitProfile(&EoEtrue_vs_Et_Calo_EB_new,0.,100.,std::string("caloParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EB_old= makeFitProfile(&EoEtrue_vs_Et_Seed_EB_old,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EB_new= makeFitProfile(&EoEtrue_vs_Et_Seed_EB_new,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EB_old= makeFitProfile(&EoEtrue_vs_Energy_Calo_EB_old,0.,250.,std::string("caloParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EB_new= makeFitProfile(&EoEtrue_vs_Energy_Calo_EB_new,0.,250.,std::string("caloParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EB_old= makeFitProfile(&EoEtrue_vs_nVtx_EB_old,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EB_new= makeFitProfile(&EoEtrue_vs_nVtx_EB_new,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EB_old= makeFitProfile(&EoEtrue_vs_Rho_EB_old,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EB_new= makeFitProfile(&EoEtrue_vs_Rho_EB_new,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EE_old= makeFitProfile(&EoEtrue_vs_Et_Calo_EE_old,0.,100.,std::string("caloParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EE_new= makeFitProfile(&EoEtrue_vs_Et_Calo_EE_new,0.,100.,std::string("caloParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EE_old= makeFitProfile(&EoEtrue_vs_Et_Seed_EE_old,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EE_new= makeFitProfile(&EoEtrue_vs_Et_Seed_EE_new,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EE_old= makeFitProfile(&EoEtrue_vs_Energy_Calo_EE_old,0.,1000.,std::string("caloParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EE_new= makeFitProfile(&EoEtrue_vs_Energy_Calo_EE_new,0.,1000.,std::string("caloParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EE_old= makeFitProfile(&EoEtrue_vs_nVtx_EE_old,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EE_new= makeFitProfile(&EoEtrue_vs_nVtx_EE_new,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EE_old= makeFitProfile(&EoEtrue_vs_Rho_EE_old,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EE_new= makeFitProfile(&EoEtrue_vs_Rho_EE_new,0.,80.,std::string("#rho (GeV)"),fitFunction_); 
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Eta_old= makeFitProfile(&EoEgen_vs_Eta_Gen_old,-3.,3.,std::string("genParticle_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Eta_new= makeFitProfile(&EoEgen_vs_Eta_Gen_new,-3.,3.,std::string("genParticle_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEta_old= makeFitProfile(&EoEgen_vs_Eta_Seed_old,-3.,3.,std::string("seed_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEta_new= makeFitProfile(&EoEgen_vs_Eta_Seed_new,-3.,3.,std::string("seed_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EB_old= makeFitProfile(&EoEgen_vs_Et_Gen_EB_old,0.,100.,std::string("genParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EB_new= makeFitProfile(&EoEgen_vs_Et_Gen_EB_new,0.,100.,std::string("genParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EB_old= makeFitProfile(&EoEgen_vs_Et_Seed_EB_old,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EB_new= makeFitProfile(&EoEgen_vs_Et_Seed_EB_new,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EB_old= makeFitProfile(&EoEgen_vs_Energy_Gen_EB_old,0.,250.,std::string("genParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EB_new= makeFitProfile(&EoEgen_vs_Energy_Gen_EB_new,0.,250.,std::string("genParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EB_old= makeFitProfile(&EoEgen_vs_nVtx_EB_old,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EB_new= makeFitProfile(&EoEgen_vs_nVtx_EB_new,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EB_old= makeFitProfile(&EoEgen_vs_Rho_EB_old,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EB_new= makeFitProfile(&EoEgen_vs_Rho_EB_new,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EE_old= makeFitProfile(&EoEgen_vs_Et_Gen_EE_old,0.,100.,std::string("genParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EE_new= makeFitProfile(&EoEgen_vs_Et_Gen_EE_new,0.,100.,std::string("genParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EE_old= makeFitProfile(&EoEgen_vs_Et_Seed_EE_old,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EE_new= makeFitProfile(&EoEgen_vs_Et_Seed_EE_new,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EE_old= makeFitProfile(&EoEgen_vs_Energy_Gen_EE_old,0.,1000.,std::string("genParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EE_new= makeFitProfile(&EoEgen_vs_Energy_Gen_EE_new,0.,1000.,std::string("genParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EE_old= makeFitProfile(&EoEgen_vs_nVtx_EE_old,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EE_new= makeFitProfile(&EoEgen_vs_nVtx_EE_new,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EE_old= makeFitProfile(&EoEgen_vs_Rho_EE_old,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EE_new= makeFitProfile(&EoEgen_vs_Rho_EE_new,0.,80.,std::string("#rho (GeV)"),fitFunction_);
    
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Eta_old_eff = makeFitProfile(&EoEtrue_vs_Eta_Calo_old,-3.,3.,std::string("caloParticle_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Eta_new_eff = makeFitProfile(&EoEtrue_vs_Eta_Calo_new,-3.,3.,std::string("caloParticle_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEta_old_eff = makeFitProfile(&EoEtrue_vs_Eta_Seed_old,-3.,3.,std::string("seed_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEta_new_eff = makeFitProfile(&EoEtrue_vs_Eta_Seed_new,-3.,3.,std::string("seed_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EB_old_eff= makeFitProfile(&EoEtrue_vs_Et_Calo_EB_old,0.,100.,std::string("caloParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EB_new_eff= makeFitProfile(&EoEtrue_vs_Et_Calo_EB_new,0.,100.,std::string("caloParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EB_old_eff= makeFitProfile(&EoEtrue_vs_Et_Seed_EB_old,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EB_new_eff= makeFitProfile(&EoEtrue_vs_Et_Seed_EB_new,0.,100.,std::string("seed_Et"),fitFunction_,true); 
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EB_old_eff= makeFitProfile(&EoEtrue_vs_Energy_Calo_EB_old,0.,250.,std::string("caloParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EB_new_eff= makeFitProfile(&EoEtrue_vs_Energy_Calo_EB_new,0.,250.,std::string("caloParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EB_old_eff= makeFitProfile(&EoEtrue_vs_nVtx_EB_old,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EB_new_eff= makeFitProfile(&EoEtrue_vs_nVtx_EB_new,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EB_old_eff= makeFitProfile(&EoEtrue_vs_Rho_EB_old,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EB_new_eff= makeFitProfile(&EoEtrue_vs_Rho_EB_new,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EE_old_eff= makeFitProfile(&EoEtrue_vs_Et_Calo_EE_old,0.,100.,std::string("caloParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EE_new_eff= makeFitProfile(&EoEtrue_vs_Et_Calo_EE_new,0.,100.,std::string("caloParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EE_old_eff= makeFitProfile(&EoEtrue_vs_Et_Seed_EE_old,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EE_new_eff= makeFitProfile(&EoEtrue_vs_Et_Seed_EE_new,0.,100.,std::string("seed_Et"),fitFunction_,true); 
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EE_old_eff= makeFitProfile(&EoEtrue_vs_Energy_Calo_EE_old,0.,1000.,std::string("caloParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EE_new_eff= makeFitProfile(&EoEtrue_vs_Energy_Calo_EE_new,0.,1000.,std::string("caloParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EE_old_eff= makeFitProfile(&EoEtrue_vs_nVtx_EE_old,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EE_new_eff= makeFitProfile(&EoEtrue_vs_nVtx_EE_new,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EE_old_eff= makeFitProfile(&EoEtrue_vs_Rho_EE_old,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EE_new_eff= makeFitProfile(&EoEtrue_vs_Rho_EE_new,0.,80.,std::string("#rho (GeV)"),fitFunction_,true); 
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Eta_old_eff= makeFitProfile(&EoEgen_vs_Eta_Gen_old,-3.,3.,std::string("genParticle_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Eta_new_eff= makeFitProfile(&EoEgen_vs_Eta_Gen_new,-3.,3.,std::string("genParticle_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEta_old_eff= makeFitProfile(&EoEgen_vs_Eta_Seed_old,-3.,3.,std::string("seed_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEta_new_eff= makeFitProfile(&EoEgen_vs_Eta_Seed_new,-3.,3.,std::string("seed_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EB_old_eff= makeFitProfile(&EoEgen_vs_Et_Gen_EB_old,0.,100.,std::string("genParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EB_new_eff= makeFitProfile(&EoEgen_vs_Et_Gen_EB_new,0.,100.,std::string("genParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EB_old_eff= makeFitProfile(&EoEgen_vs_Et_Seed_EB_old,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EB_new_eff= makeFitProfile(&EoEgen_vs_Et_Seed_EB_new,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EB_old_eff= makeFitProfile(&EoEgen_vs_Energy_Gen_EB_old,0.,250.,std::string("genParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EB_new_eff= makeFitProfile(&EoEgen_vs_Energy_Gen_EB_new,0.,250.,std::string("genParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EB_old_eff= makeFitProfile(&EoEgen_vs_nVtx_EB_old,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EB_new_eff= makeFitProfile(&EoEgen_vs_nVtx_EB_new,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EB_old_eff= makeFitProfile(&EoEgen_vs_Rho_EB_old,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EB_new_eff= makeFitProfile(&EoEgen_vs_Rho_EB_new,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EE_old_eff= makeFitProfile(&EoEgen_vs_Et_Gen_EE_old,0.,100.,std::string("genParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EE_new_eff= makeFitProfile(&EoEgen_vs_Et_Gen_EE_new,0.,100.,std::string("genParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EE_old_eff= makeFitProfile(&EoEgen_vs_Et_Seed_EE_old,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EE_new_eff= makeFitProfile(&EoEgen_vs_Et_Seed_EE_new,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EE_old_eff= makeFitProfile(&EoEgen_vs_Energy_Gen_EE_old,0.,1000.,std::string("genParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EE_new_eff= makeFitProfile(&EoEgen_vs_Energy_Gen_EE_new,0.,1000.,std::string("genParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EE_old_eff= makeFitProfile(&EoEgen_vs_nVtx_EE_old,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EE_new_eff= makeFitProfile(&EoEgen_vs_nVtx_EE_new,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EE_old_eff= makeFitProfile(&EoEgen_vs_Rho_EE_old,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EE_new_eff= makeFitProfile(&EoEgen_vs_Rho_EE_new,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   
   
   drawGraph(gr_EoEtrue_vs_Eta_old.first, gr_EoEtrue_vs_Eta_new.first, std::string("caloParticle_#eta"), std::string("#mu"), std::string("EoEtrue_vs_caloEta_Mean"), superClusterRef_, superClusterVal_,0.97,1.03);
   drawGraph(gr_EoEtrue_vs_SeedEta_old.first, gr_EoEtrue_vs_SeedEta_new.first, std::string("seed_#eta"), std::string("#mu"), std::string("EoEtrue_vs_seedEta_Mean"), superClusterRef_, superClusterVal_,0.97,1.03);
   drawGraph(gr_EoEtrue_vs_Eta_old.second, gr_EoEtrue_vs_Eta_new.second, std::string("caloParticle_#eta"), std::string("#sigma"), std::string("EoEtrue_vs_caloEta_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEta_old.second, gr_EoEtrue_vs_Eta_new.second, std::string("seed_#eta"), std::string("#sigma"), std::string("EoEtrue_vs_seedEta_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_Et_EB_old.first, gr_EoEtrue_vs_Et_EB_new.first, std::string("caloParticle_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEt_EB_Mean"), superClusterRef_, superClusterVal_,0.99,1.01);
   drawGraph(gr_EoEtrue_vs_SeedEt_EB_old.first, gr_EoEtrue_vs_SeedEt_EB_new.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_seedEt_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Et_EB_old.second, gr_EoEtrue_vs_Et_EB_new.second, std::string("caloParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEt_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EB_old.second, gr_EoEtrue_vs_SeedEt_EB_new.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_seedEt_EB_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_Et_EE_old.first, gr_EoEtrue_vs_Et_EE_new.first, std::string("caloParticle_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEt_EE_Mean"), superClusterRef_, superClusterVal_,0.99,1.01);
   drawGraph(gr_EoEtrue_vs_SeedEt_EE_old.first, gr_EoEtrue_vs_SeedEt_EE_new.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_seedEt_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Et_EE_old.second, gr_EoEtrue_vs_Et_EE_new.second, std::string("caloParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEt_EE_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EE_old.second, gr_EoEtrue_vs_SeedEt_EE_new.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_seedEt_EE_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Energy_EB_old.first, gr_EoEtrue_vs_Energy_EB_new.first, std::string("caloParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEnergy_EB_Mean"), superClusterRef_, superClusterVal_,0.99,1.01);
   drawGraph(gr_EoEtrue_vs_Energy_EB_old.second, gr_EoEtrue_vs_Energy_EB_new.second, std::string("caloParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEnergy_EB_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_Energy_EE_old.first, gr_EoEtrue_vs_Energy_EE_new.first, std::string("caloParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEnergy_EE_Mean"), superClusterRef_, superClusterVal_,0.99,1.01);
   drawGraph(gr_EoEtrue_vs_Energy_EE_old.second, gr_EoEtrue_vs_Energy_EE_new.second, std::string("caloParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEnergy_EE_Resolution"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEtrue_vs_nVtx_EB_old.first, gr_EoEtrue_vs_nVtx_EB_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EB_old.first, gr_EoEtrue_vs_nVtx_EB_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EB_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.97,1.03);
   drawGraph(gr_EoEtrue_vs_nVtx_EB_old.second, gr_EoEtrue_vs_nVtx_EB_new.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEtrue_vs_nVtx_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EE_old.first, gr_EoEtrue_vs_nVtx_EE_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EE_old.first, gr_EoEtrue_vs_nVtx_EE_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EE_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.95,1.03);
   drawGraph(gr_EoEtrue_vs_nVtx_EE_old.second, gr_EoEtrue_vs_nVtx_EE_new.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEtrue_vs_nVtx_EE_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EB_old.first, gr_EoEtrue_vs_Rho_EB_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EB_old.first, gr_EoEtrue_vs_Rho_EB_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EB_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.97,1.03);
   drawGraph(gr_EoEtrue_vs_Rho_EB_old.second, gr_EoEtrue_vs_Rho_EB_new.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_Rho_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EE_old.first, gr_EoEtrue_vs_Rho_EE_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EE_old.first, gr_EoEtrue_vs_Rho_EE_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EE_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.95,1.03);
   drawGraph(gr_EoEtrue_vs_Rho_EE_old.second, gr_EoEtrue_vs_Rho_EE_new.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_Rho_EE_Resolution"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEgen_vs_Eta_old.first, gr_EoEgen_vs_Eta_new.first, std::string("genParticle_#eta"), std::string("#mu"), std::string("EoEgen_vs_genEta_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEta_old.first, gr_EoEgen_vs_SeedEta_new.first, std::string("seed_#eta"), std::string("#mu"), std::string("EoEgen_vs_seedEta_Mean"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Eta_old.second, gr_EoEgen_vs_Eta_new.second, std::string("genParticle_#eta"), std::string("#sigma"), std::string("EoEgen_vs_genEta_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEta_old.second, gr_EoEgen_vs_SeedEta_new.second, std::string("seed_#eta"), std::string("#sigma"), std::string("EoEgen_vs_seedEta_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Et_EB_old.first, gr_EoEgen_vs_Et_EB_new.first, std::string("genParticle_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEt_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EB_old.first, gr_EoEgen_vs_Et_EB_new.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_seedEt_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Et_EB_old.second, gr_EoEgen_vs_Et_EB_new.second, std::string("genParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEt_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EB_old.second, gr_EoEgen_vs_SeedEt_EB_new.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_seedEt_EB_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Et_EE_old.first, gr_EoEgen_vs_Et_EE_new.first, std::string("genParticle_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEt_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EE_old.first, gr_EoEgen_vs_SeedEt_EE_new.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_seedEt_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Et_EE_old.second, gr_EoEgen_vs_Et_EE_new.second, std::string("genParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEt_EE_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EE_old.second, gr_EoEgen_vs_SeedEt_EE_new.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_seedEt_EE_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Energy_EB_old.first, gr_EoEgen_vs_Energy_EB_new.first, std::string("genParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEnergy_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Energy_EB_old.second, gr_EoEgen_vs_Energy_EB_new.second, std::string("genParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEnergy_EB_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Energy_EE_old.first, gr_EoEgen_vs_Energy_EE_new.first, std::string("genParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEnergy_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Energy_EE_old.second, gr_EoEgen_vs_Energy_EE_new.second, std::string("genParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEnergy_EE_Resolution"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEgen_vs_nVtx_EB_old.first, gr_EoEgen_vs_nVtx_EB_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EB_old.first, gr_EoEgen_vs_nVtx_EB_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EB_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.97,1.03);
   drawGraph(gr_EoEgen_vs_nVtx_EB_old.second, gr_EoEgen_vs_nVtx_EB_new.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEgen_vs_nVtx_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EE_old.first, gr_EoEgen_vs_nVtx_EE_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EE_old.first, gr_EoEgen_vs_nVtx_EE_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EE_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.93,1.);
   drawGraph(gr_EoEgen_vs_nVtx_EE_old.second, gr_EoEgen_vs_nVtx_EE_new.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEgen_vs_nVtx_EE_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EB_old.first, gr_EoEgen_vs_Rho_EB_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EB_old.first, gr_EoEgen_vs_Rho_EB_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EB_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.97,1.03);
   drawGraph(gr_EoEgen_vs_Rho_EB_old.second, gr_EoEgen_vs_Rho_EB_new.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_Rho_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EE_old.first, gr_EoEgen_vs_Rho_EE_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EE_old.first, gr_EoEgen_vs_Rho_EE_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EE_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.93,1.0);
   drawGraph(gr_EoEgen_vs_Rho_EE_old.second, gr_EoEgen_vs_Rho_EE_new.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_Rho_EE_Resolution"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEtrue_vs_Eta_old_eff.first, gr_EoEtrue_vs_Eta_new_eff.first, std::string("caloParticle_#eta"), std::string("#mu"), std::string("EoEtrue_vs_caloEta_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEta_old_eff.first, gr_EoEtrue_vs_SeedEta_new_eff.first, std::string("seed_#eta"), std::string("#mu"), std::string("EoEtrue_vs_seedEta_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Eta_old_eff.second, gr_EoEtrue_vs_Eta_new_eff.second, std::string("caloParticle_#eta"), std::string("#sigma"), std::string("EoEtrue_vs_caloEta_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEta_old_eff.second, gr_EoEtrue_vs_SeedEta_new_eff.second, std::string("seed_#eta"), std::string("#sigma"), std::string("EoEtrue_vs_seedEta_Resolution_Effective"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_Et_EB_old_eff.first, gr_EoEtrue_vs_Et_EB_new_eff.first, std::string("caloParticle_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEt_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EB_old_eff.first, gr_EoEtrue_vs_SeedEt_EB_new_eff.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_seedEt_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Et_EB_old_eff.second, gr_EoEtrue_vs_Et_EB_new_eff.second, std::string("caloParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEt_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EB_old_eff.second, gr_EoEtrue_vs_SeedEt_EB_new_eff.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_seedEt_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEtrue_vs_Et_EE_old_eff.first, gr_EoEtrue_vs_Et_EE_new_eff.first, std::string("caloParticle_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEt_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EE_old_eff.first, gr_EoEtrue_vs_SeedEt_EE_new_eff.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_seedEt_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Et_EE_old_eff.second, gr_EoEtrue_vs_Et_EE_new_eff.second, std::string("caloParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEt_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EE_old_eff.second, gr_EoEtrue_vs_SeedEt_EE_new_eff.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_seedEt_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Energy_EB_old_eff.first, gr_EoEtrue_vs_Energy_EB_new_eff.first, std::string("caloParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEnergy_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Energy_EB_old_eff.second, gr_EoEtrue_vs_Energy_EB_new_eff.second, std::string("caloParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEnergy_EB_Resolution_Effective"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_Energy_EE_old_eff.first, gr_EoEtrue_vs_Energy_EE_new_eff.first, std::string("caloParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEnergy_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Energy_EE_old_eff.second, gr_EoEtrue_vs_Energy_EE_new_eff.second, std::string("caloParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEnergy_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEtrue_vs_nVtx_EB_old_eff.first, gr_EoEtrue_vs_nVtx_EB_new_eff.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EB_old_eff.second, gr_EoEtrue_vs_nVtx_EB_new_eff.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEtrue_vs_nVtx_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EE_old_eff.first, gr_EoEtrue_vs_nVtx_EE_new_eff.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EE_old_eff.second, gr_EoEtrue_vs_nVtx_EE_new_eff.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEtrue_vs_nVtx_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EB_old_eff.first, gr_EoEtrue_vs_Rho_EB_new_eff.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EB_old_eff.second, gr_EoEtrue_vs_Rho_EB_new_eff.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_Rho_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EE_old_eff.first, gr_EoEtrue_vs_Rho_EE_new_eff.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EE_old_eff.second, gr_EoEtrue_vs_Rho_EE_new_eff.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_Rho_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);  

   drawGraph(gr_EoEgen_vs_Eta_old_eff.first, gr_EoEgen_vs_Eta_new_eff.first, std::string("genParticle_#eta"), std::string("#mu"), std::string("EoEgen_vs_genEta_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEta_old_eff.first, gr_EoEgen_vs_SeedEta_new_eff.first, std::string("seed_#eta"), std::string("#mu"), std::string("EoEgen_vs_seedEta_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Eta_old_eff.second, gr_EoEgen_vs_Eta_new_eff.second, std::string("genParticle_#eta"), std::string("#sigma"), std::string("EoEgen_vs_genEta_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEta_old_eff.second, gr_EoEgen_vs_SeedEta_new_eff.second, std::string("seed_#eta"), std::string("#sigma"), std::string("EoEgen_vs_seedEta_Resolution_Effective"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Et_EB_old_eff.first, gr_EoEgen_vs_Et_EB_new_eff.first, std::string("genParticle_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEt_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EB_old_eff.first, gr_EoEgen_vs_SeedEt_EB_new_eff.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_seedEt_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Et_EB_old_eff.second, gr_EoEgen_vs_Et_EB_new_eff.second, std::string("genParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEt_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EB_old_eff.second, gr_EoEgen_vs_SeedEt_EB_new_eff.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_seedEt_EB_Resolution_Effective"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Et_EE_old_eff.first, gr_EoEgen_vs_Et_EE_new_eff.first, std::string("genParticle_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEt_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EE_old_eff.first, gr_EoEgen_vs_SeedEt_EE_new_eff.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_seedEt_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Et_EE_old_eff.second, gr_EoEgen_vs_Et_EE_new_eff.second, std::string("genParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEt_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EE_old_eff.second, gr_EoEgen_vs_SeedEt_EE_new_eff.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_seedEt_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Energy_EB_old_eff.first, gr_EoEgen_vs_Energy_EB_new_eff.first, std::string("genParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEnergy_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Energy_EB_old_eff.second, gr_EoEgen_vs_Energy_EB_new_eff.second, std::string("genParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEnergy_EB_Resolution_Effective"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Energy_EE_old_eff.first, gr_EoEgen_vs_Energy_EE_new_eff.first, std::string("genParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEnergy_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Energy_EE_old_eff.second, gr_EoEgen_vs_Energy_EE_new_eff.second, std::string("genParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEnergy_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEgen_vs_nVtx_EB_old_eff.first, gr_EoEgen_vs_nVtx_EB_new_eff.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EB_old_eff.second, gr_EoEgen_vs_nVtx_EB_new_eff.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEgen_vs_nVtx_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EE_old_eff.first, gr_EoEgen_vs_nVtx_EE_new_eff.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EE_old_eff.second, gr_EoEgen_vs_nVtx_EE_new_eff.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEgen_vs_nVtx_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EB_old_eff.first, gr_EoEgen_vs_Rho_EB_new_eff.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EB_old_eff.second, gr_EoEgen_vs_Rho_EB_new_eff.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_Rho_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EE_old_eff.first, gr_EoEgen_vs_Rho_EE_new_eff.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EE_old_eff.second, gr_EoEgen_vs_Rho_EE_new_eff.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_Rho_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);

   drawHisto(h_EoEtrue_EB_old, h_EoEtrue_EB_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EB"), 0, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEtrue_EB_old, h_EoEtrue_EB_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EB"), 1, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEtrue_EE_old, h_EoEtrue_EE_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EE"), 0, true, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EoEtrue_EE_old, h_EoEtrue_EE_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EE"), 1, true, fitFunction_);
   drawHisto(h_EoEtrue_EB_seedMatched_old, h_EoEtrue_EB_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EB_seedMatched"), 0, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEtrue_EB_seedMatched_old, h_EoEtrue_EB_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EB_seedMatched"), 1, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEtrue_EE_seedMatched_old, h_EoEtrue_EE_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EE_seedMatched"), 0, true, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EoEtrue_EE_seedMatched_old, h_EoEtrue_EE_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EE_seedMatched"), 1, true, fitFunction_);
   drawHisto(h_EoEgen_EB_old, h_EoEgen_EB_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EB"), 0, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEgen_EB_old, h_EoEgen_EB_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EB"), 1, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEgen_EE_old, h_EoEgen_EE_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EE"), 0, true, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EoEgen_EE_old, h_EoEgen_EE_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EE"), 1, true, fitFunction_);
   drawHisto(h_EoEgen_EB_seedMatched_old, h_EoEgen_EB_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EB_seedMatched"), 0, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEgen_EB_seedMatched_old, h_EoEgen_EB_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EB_seedMatched"), 1, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEgen_EE_seedMatched_old, h_EoEgen_EE_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EE_seedMatched"), 0, true, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EoEgen_EE_seedMatched_old, h_EoEgen_EE_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EE_seedMatched"), 1, true, fitFunction_);
   drawHisto(h_Energy_EB_old, h_Energy_EB_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EB_old, h_Energy_EB_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_old, h_Energy_EE_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_old, h_Energy_EE_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_old, h_Eta_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_old, h_Eta_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EB_old, h_Phi_EB_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_Phi_EB_old, h_Phi_EB_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EE_old, h_Phi_EE_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_Phi_EE_old, h_Phi_EE_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EB_old, h_EtaWidth_EB_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EtaWidth_EB_old, h_EtaWidth_EB_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EtaWidth_EE_old, h_EtaWidth_EE_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EE_old, h_EtaWidth_EE_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EB_old, h_PhiWidth_EB_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_PhiWidth_EB_old, h_PhiWidth_EB_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_PhiWidth_EE_old, h_PhiWidth_EE_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EE_old, h_PhiWidth_EE_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_nPFClusters_old, h_nPFClusters_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_old, h_nPFClusters_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_old, h_nPFClusters_EB_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_old, h_nPFClusters_EB_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_old, h_nPFClusters_EE_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_old, h_nPFClusters_EE_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_old, h_R9_EB_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_old, h_R9_EB_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_old, h_R9_EE_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_old, h_R9_EE_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_old, h_full5x5_R9_EB_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_old, h_full5x5_R9_EB_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_old, h_full5x5_R9_EE_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_old, h_full5x5_R9_EE_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_old, h_sigmaIetaIeta_EB_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_old, h_sigmaIetaIeta_EB_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_old, h_sigmaIetaIeta_EE_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_old, h_sigmaIetaIeta_EE_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   
   drawHisto(h_sigmaIetaIphi_EB_old, h_sigmaIetaIphi_EB_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_old, h_sigmaIetaIphi_EB_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_old, h_sigmaIetaIphi_EE_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_old, h_sigmaIetaIphi_EE_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_old, h_sigmaIphiIphi_EB_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_old, h_sigmaIphiIphi_EB_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_old, h_sigmaIphiIphi_EE_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_old, h_sigmaIphiIphi_EE_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   
   //make seedMatched plots
   drawHisto(h_Energy_EB_seedMatched_old, h_Energy_EB_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EB_seedMatched_old, h_Energy_EB_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_seedMatched_old, h_Energy_EE_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_seedMatched_old, h_Energy_EE_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_seedMatched_old, h_Eta_seedMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_seedMatched_old, h_Eta_seedMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EB_seedMatched_old, h_Phi_EB_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_Phi_EB_seedMatched_old, h_Phi_EB_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EE_seedMatched_old, h_Phi_EE_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_Phi_EE_seedMatched_old, h_Phi_EE_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EB_seedMatched_old, h_EtaWidth_EB_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EtaWidth_EB_seedMatched_old, h_EtaWidth_EB_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EtaWidth_EE_seedMatched_old, h_EtaWidth_EE_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EE_seedMatched_old, h_EtaWidth_EE_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EB_seedMatched_old, h_PhiWidth_EB_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_PhiWidth_EB_seedMatched_old, h_PhiWidth_EB_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_PhiWidth_EE_seedMatched_old, h_PhiWidth_EE_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EE_seedMatched_old, h_PhiWidth_EE_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_nPFClusters_seedMatched_old, h_nPFClusters_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_seedMatched_old, h_nPFClusters_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_seedMatched_old, h_nPFClusters_EB_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_seedMatched_old, h_nPFClusters_EB_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_seedMatched_old, h_nPFClusters_EE_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_seedMatched_old, h_nPFClusters_EE_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_seedMatched_old, h_R9_EB_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_seedMatched_old, h_R9_EB_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_seedMatched_old, h_R9_EE_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_seedMatched_old, h_R9_EE_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_seedMatched_old, h_full5x5_R9_EB_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_seedMatched_old, h_full5x5_R9_EB_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_seedMatched_old, h_full5x5_R9_EE_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_seedMatched_old, h_full5x5_R9_EE_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_seedMatched_old, h_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_seedMatched_old, h_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_seedMatched_old, h_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_seedMatched_old, h_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_seedMatched_old, h_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_seedMatched_old, h_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_seedMatched_old, h_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_seedMatched_old, h_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_seedMatched_old, h_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_seedMatched_old, h_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_seedMatched_old, h_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_seedMatched_old, h_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   
   //make caloMatched plots
   drawHisto(h_Energy_EB_caloMatched_old, h_Energy_EB_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EB_caloMatched_old, h_Energy_EB_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_caloMatched_old, h_Energy_EE_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_caloMatched_old, h_Energy_EE_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_caloMatched_old, h_Eta_caloMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_caloMatched_old, h_Eta_caloMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EB_caloMatched_old, h_Phi_EB_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_Phi_EB_caloMatched_old, h_Phi_EB_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EE_caloMatched_old, h_Phi_EE_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_Phi_EE_caloMatched_old, h_Phi_EE_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EB_caloMatched_old, h_EtaWidth_EB_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_,0.7,2.5);
   drawHisto(h_EtaWidth_EB_caloMatched_old, h_EtaWidth_EB_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_,0.7,2.5); 
   drawHisto(h_EtaWidth_EE_caloMatched_old, h_EtaWidth_EE_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EE_caloMatched_old, h_EtaWidth_EE_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EB_caloMatched_old, h_PhiWidth_EB_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_,0.7,2.5);
   drawHisto(h_PhiWidth_EB_caloMatched_old, h_PhiWidth_EB_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_,0.7,2.5); 
   drawHisto(h_PhiWidth_EE_caloMatched_old, h_PhiWidth_EE_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EE_caloMatched_old, h_PhiWidth_EE_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_nPFClusters_caloMatched_old, h_nPFClusters_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_caloMatched_old, h_nPFClusters_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_caloMatched_old, h_nPFClusters_EB_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_caloMatched_old, h_nPFClusters_EB_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_caloMatched_old, h_nPFClusters_EE_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_caloMatched_old, h_nPFClusters_EE_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_caloMatched_old, h_R9_EB_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_caloMatched_old, h_R9_EB_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_caloMatched_old, h_R9_EE_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_caloMatched_old, h_R9_EE_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_caloMatched_old, h_full5x5_R9_EB_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_caloMatched_old, h_full5x5_R9_EB_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_caloMatched_old, h_full5x5_R9_EE_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_caloMatched_old, h_full5x5_R9_EE_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_caloMatched_old, h_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_caloMatched_old, h_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_caloMatched_old, h_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_caloMatched_old, h_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_caloMatched_old, h_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_caloMatched_old, h_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_caloMatched_old, h_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_caloMatched_old, h_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_caloMatched_old, h_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_caloMatched_old, h_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_caloMatched_old, h_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_caloMatched_old, h_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   
   //make caloUnmatched plots
   drawHisto(h_Energy_EB_caloUnmatched_old, h_Energy_EB_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EB_caloUnmatched_old, h_Energy_EB_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_caloUnmatched_old, h_Energy_EE_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_caloUnmatched_old, h_Energy_EE_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_caloUnmatched_old, h_Eta_caloUnmatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_caloUnmatched_old, h_Eta_caloUnmatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EB_caloUnmatched_old, h_Phi_EB_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_Phi_EB_caloUnmatched_old, h_Phi_EB_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EE_caloUnmatched_old, h_Phi_EE_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_Phi_EE_caloUnmatched_old, h_Phi_EE_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EB_caloUnmatched_old, h_EtaWidth_EB_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EtaWidth_EB_caloUnmatched_old, h_EtaWidth_EB_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EtaWidth_EE_caloUnmatched_old, h_EtaWidth_EE_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EE_caloUnmatched_old, h_EtaWidth_EE_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EB_caloUnmatched_old, h_PhiWidth_EB_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_PhiWidth_EB_caloUnmatched_old, h_PhiWidth_EB_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_PhiWidth_EE_caloUnmatched_old, h_PhiWidth_EE_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EE_caloUnmatched_old, h_PhiWidth_EE_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_nPFClusters_caloUnmatched_old, h_nPFClusters_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_caloUnmatched_old, h_nPFClusters_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_caloUnmatched_old, h_nPFClusters_EB_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_caloUnmatched_old, h_nPFClusters_EB_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_caloUnmatched_old, h_nPFClusters_EE_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_caloUnmatched_old, h_nPFClusters_EE_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_caloUnmatched_old, h_R9_EB_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_caloUnmatched_old, h_R9_EB_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_caloUnmatched_old, h_R9_EE_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_caloUnmatched_old, h_R9_EE_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_caloUnmatched_old, h_full5x5_R9_EB_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_caloUnmatched_old, h_full5x5_R9_EB_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_caloUnmatched_old, h_full5x5_R9_EE_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_caloUnmatched_old, h_full5x5_R9_EE_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_caloUnmatched_old, h_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_caloUnmatched_old, h_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_caloUnmatched_old, h_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_caloUnmatched_old, h_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_caloUnmatched_old, h_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_caloUnmatched_old, h_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_caloUnmatched_old, h_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_caloUnmatched_old, h_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_caloUnmatched_old, h_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_caloUnmatched_old, h_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_caloUnmatched_old, h_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_caloUnmatched_old, h_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   
   
   drawGraph(gr_EoEtrue_vs_EtaWidth_EB_old.first, gr_EoEtrue_vs_EtaWidth_EB_new.first, std::string("eta_width"), std::string("#mu"), std::string("EoEtrue_vs_EtaWidth_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_EtaWidth_EB_old.second, gr_EoEtrue_vs_EtaWidth_EB_new.second, std::string("eta_width"), std::string("#sigma"), std::string("EoEtrue_vs_EtaWidth_EB_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_EtaWidth_EE_old.first, gr_EoEtrue_vs_EtaWidth_EE_new.first, std::string("eta_width"), std::string("#mu"), std::string("EoEtrue_vs_EtaWidth_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_EtaWidth_EE_old.second, gr_EoEtrue_vs_EtaWidth_EE_new.second, std::string("eta_width"), std::string("#sigma"), std::string("EoEtrue_vs_EtaWidth_EE_Resolution"), superClusterRef_, superClusterVal_); 

   drawGraph(gr_EoEtrue_vs_PhiWidth_EB_old.first, gr_EoEtrue_vs_PhiWidth_EB_new.first, std::string("phi_width"), std::string("#mu"), std::string("EoEtrue_vs_PhiWidth_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_PhiWidth_EB_old.second, gr_EoEtrue_vs_PhiWidth_EB_new.second, std::string("phi_width"), std::string("#sigma"), std::string("EoEtrue_vs_PhiWidth_EB_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_PhiWidth_EE_old.first, gr_EoEtrue_vs_PhiWidth_EE_new.first, std::string("phi_width"), std::string("#mu"), std::string("EoEtrue_vs_PhiWidth_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_PhiWidth_EE_old.second, gr_EoEtrue_vs_PhiWidth_EE_new.second, std::string("phi_width"), std::string("#sigma"), std::string("EoEtrue_vs_PhiWidth_EE_Resolution"), superClusterRef_, superClusterVal_);
    

   drawHisto(h_full5x5_sigmaIetaIeta_EB_old, h_full5x5_sigmaIetaIeta_EB_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_old, h_full5x5_sigmaIetaIeta_EB_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_old, h_full5x5_sigmaIetaIeta_EE_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_old, h_full5x5_sigmaIetaIeta_EE_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_old, h_full5x5_sigmaIetaIphi_EB_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_old, h_full5x5_sigmaIetaIphi_EB_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_old, h_full5x5_sigmaIetaIphi_EE_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_old, h_full5x5_sigmaIetaIphi_EE_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_old, h_full5x5_sigmaIphiIphi_EB_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_old, h_full5x5_sigmaIphiIphi_EB_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_old, h_full5x5_sigmaIphiIphi_EE_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_old, h_full5x5_sigmaIphiIphi_EE_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);
    drawHisto(h_full5x5_sigmaIetaIeta_EB_seedMatched_old, h_full5x5_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_seedMatched_old, h_full5x5_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_seedMatched_old, h_full5x5_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_seedMatched_old, h_full5x5_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_seedMatched_old, h_full5x5_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_seedMatched_old, h_full5x5_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_seedMatched_old, h_full5x5_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_seedMatched_old, h_full5x5_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_seedMatched_old, h_full5x5_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_seedMatched_old, h_full5x5_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_seedMatched_old, h_full5x5_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_seedMatched_old, h_full5x5_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);
    drawHisto(h_full5x5_sigmaIetaIeta_EB_caloMatched_old, h_full5x5_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_caloMatched_old, h_full5x5_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloMatched_old, h_full5x5_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloMatched_old, h_full5x5_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloMatched_old, h_full5x5_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloMatched_old, h_full5x5_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloMatched_old, h_full5x5_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloMatched_old, h_full5x5_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloMatched_old, h_full5x5_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloMatched_old, h_full5x5_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloMatched_old, h_full5x5_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloMatched_old, h_full5x5_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);
    drawHisto(h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
    
   //eff ratios
   h_Eta_Calo_old->Divide(h_Eta_Calo_Denum_old);
   h_Eta_Calo_new->Divide(h_Eta_Calo_Denum_new);
   h_Et_Calo_EB_old->Divide(h_Et_Calo_EB_Denum_old);
   h_Et_Calo_EB_new->Divide(h_Et_Calo_EB_Denum_new);
   h_Et_Calo_EE_old->Divide(h_Et_Calo_EE_Denum_old);
   h_Et_Calo_EE_new->Divide(h_Et_Calo_EE_Denum_new);

   drawHisto(h_Eta_Calo_old, h_Eta_Calo_new, std::string("caloParticle_#eta"), std::string("hist"), std::string("Efficiency_vs_CaloEta_ratio"), 0, false, fitFunction_, superClusterRef_, superClusterVal_,0.90,1.02); 
   drawHisto(h_Et_Calo_EB_old, h_Et_Calo_EB_new, std::string("caloParticle_Et (GeV)"), std::string("hist"), std::string("Efficiency_vs_CaloEt_EB_ratio"), 0, false, fitFunction_, superClusterRef_, superClusterVal_,0.90,1.02); 
   drawHisto(h_Et_Calo_EE_old, h_Et_Calo_EE_new, std::string("caloParticle_Et (GeV)"), std::string("hist"), std::string("Efficiency_vs_CaloEt_EE_ratio"), 0, false, fitFunction_, superClusterRef_, superClusterVal_,0.90,1.02); 

   

}

void SCAlgoValidation_localSetVsSingleSet(){

    //Single vs Mustache
    /*string superClusterRef_label = "Mustache";
    string superClusterRef_collection = "Mustache";
    string SCRef_File = "./SCAlgoPlots_SingleSetParams.root";

    string superClusterVal_label = "retunedMustache";
    string superClusterVal_collection = "SingleSet";
    string SCVal_File = "./SCAlgoPlots_SingleSetParams.root";
    */

    //Single vs Local
    string superClusterRef_label = "SingleSet";
    string superClusterRef_collection = "SingleSet";
    string SCRef_File = "./SCAlgoPlots_SingleSetParams.root";

    string superClusterVal_label = "LocalSets";
    string superClusterVal_collection = "LocalSets";
    string SCVal_File = "./SCAlgoPlots_LocalSetParams.root";
    

    //Local vs Mustache
    /*string superClusterRef_label = "Mustache";
    string superClusterRef_collection = "Mustache";
    string SCRef_File = "./SCAlgoPlots_LocalSetParams.root";

    string superClusterVal_label = "retunedMustache";
    string superClusterVal_collection = "LocalSets";
    string SCVal_File = "./SCAlgoPlots_LocalSetParams.root";
    */

    ReadInfile(superClusterRef_collection, superClusterVal_collection, SCRef_File, SCVal_File);
    setEfficiencies();
    drawPlots("cruijff",superClusterRef_label, superClusterVal_label);
}
