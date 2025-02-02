#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "PhysicsTools/Utilities/macros/setTDRStyle.C"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoSimStudies/Dumpers/interface/SCAlgoValidation.h"
#include "RecoSimStudies/Dumpers/interface/CruijffPdf.h"

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include <algorithm> 
#include <iostream>
#include <utility>

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

using namespace std;

int main(int argc, char** argv)
{
    const edm::ParameterSet &process         = edm::readPSetsFrom( argv[1] )->getParameter<edm::ParameterSet>( "process" );
    const edm::ParameterSet &filesOpt        = process.getParameter<edm::ParameterSet>( "ioFilesOpt" );
    const edm::ParameterSet &histOpt        = process.getParameter<edm::ParameterSet>( "histOpt" );  

    //Read config
    string inputFiles_ = filesOpt.getParameter<string>( "inputFiles" );
    vector<string> inputFiles = split(inputFiles_,',');

    string outputDir_ = filesOpt.getParameter<string>( "outputDir" );
    if( outputDir_ == "" ) outputDir_ = "output/"; 

    int maxEvents_ = filesOpt.getUntrackedParameter<int>( "maxEvents" );
    string fitFunction_ = filesOpt.getParameter<string>( "fitFunction" );
    string superClusterRef_ = filesOpt.getParameter<string>( "superClusterRef" );
    string superClusterVal_ = filesOpt.getParameter<string>( "superClusterVal" );

    std::vector<float> etCuts {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
    std::vector<float> etaCuts {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.479,1.75,2.0,2.25,3.0};
    std::vector<std::pair<std::string,std::vector<double>>> binOpts = getBinOpts(histOpt); 
    setHistograms(binOpts, etCuts, etaCuts);
    binOpts = getBinOpts(histOpt); //reset binOpts 

    //useful vars 
    vector<float> simEnergy;
    vector<int> superCluster_sim_fraction_MatchedIndex;
    vector<int> deepSuperCluster_sim_fraction_MatchedIndex;
    vector<int> caloParticle_superCluster_sim_fraction_MatchedIndex;
    vector<int> caloParticle_deepSuperCluster_sim_fraction_MatchedIndex;
    std::map<int,float> SuperCluster_RecoEnergy_sim_fraction;
    std::map<int,float> DeepSuperCluster_RecoEnergy_sim_fraction;
    std::map<int,std::map<int,double>> Calo_SuperCluster_sim_fraction;
    std::map<int,std::map<int,double>> SuperCluster_Calo_sim_fraction;
    std::map<int,std::map<int,double>> Calo_DeepSuperCluster_sim_fraction;
    std::map<int,std::map<int,double>> DeepSuperCluster_Calo_sim_fraction; 
    std::vector<int> superCluster_seeds;
    std::vector<int> deepSuperCluster_seeds;
    std::vector<int> common_seeds; 

    bool gotoMain = false;
    for(unsigned int iFile=0; iFile<inputFiles.size() && !gotoMain; ++iFile)
    {
       TFile* inFile = TFile::Open(inputFiles[iFile].c_str());
       TTree* inTree = (TTree*)inFile->Get("recosimdumper/caloTree");
       if(inFile)
       {
          cout << "\n--- Reading file: " << inputFiles[iFile].c_str() << endl;

          //setBranches
          setTreeBranches(inTree,superClusterRef_,superClusterVal_);

          // Loop over all entries of the TTree
          for(int entry = 0; entry < inTree->GetEntries(); entry++){
     
             if(entry>maxEvents_ && maxEvents_>0) continue;
             //if(entry<59259) continue;
             if(entry%10000==0) std::cout << "--- Reading tree = " << entry << std::endl;
             inTree->GetEntry(entry);

             simEnergy.clear();
             superCluster_sim_fraction_MatchedIndex.clear();
             deepSuperCluster_sim_fraction_MatchedIndex.clear();
             caloParticle_superCluster_sim_fraction_MatchedIndex.clear();
             caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.clear(); 
             SuperCluster_RecoEnergy_sim_fraction.clear(); 
             DeepSuperCluster_RecoEnergy_sim_fraction.clear();
             superCluster_seeds.clear();
             deepSuperCluster_seeds.clear();
             common_seeds.clear();

             //get SC with highest score
             Calo_SuperCluster_sim_fraction.clear();
             SuperCluster_Calo_sim_fraction.clear();
             for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++)
                 for(unsigned int iSC=0; iSC<superCluster_rawEnergy->size(); iSC++){
                     Calo_SuperCluster_sim_fraction[iCalo][iSC]=0.;
                     SuperCluster_Calo_sim_fraction[iSC][iCalo]=0.;
                 } 
       
             Calo_DeepSuperCluster_sim_fraction.clear();
             DeepSuperCluster_Calo_sim_fraction.clear();
             for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++)
                 for(unsigned int iSC=0; iSC<deepSuperCluster_rawEnergy->size(); iSC++){
                     Calo_DeepSuperCluster_sim_fraction[iCalo][iSC]=0.;
                     DeepSuperCluster_Calo_sim_fraction[iSC][iCalo]=0.;
                 } 
           
             for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++){
                 for(unsigned int iSC=0; iSC<superCluster_rawEnergy->size(); iSC++){
                     for(unsigned int iPF=0; iPF<superCluster_pfClustersIndex->at(iSC).size(); iPF++)
                     {  
                         int pfCluster_index = superCluster_pfClustersIndex->at(iSC).at(iPF);    
                         if(pfCluster_sim_fraction->at(pfCluster_index).at(iCalo) > 0.){
                            Calo_SuperCluster_sim_fraction[iCalo][iSC]+=pfCluster_sim_fraction->at(pfCluster_index).at(iCalo);
                            SuperCluster_Calo_sim_fraction[iSC][iCalo]+=pfCluster_sim_fraction->at(pfCluster_index).at(iCalo);
                         }else{
                            Calo_SuperCluster_sim_fraction[iCalo][iSC]+=0.;
                            SuperCluster_Calo_sim_fraction[iSC][iCalo]+=0.;            
                         }  
                     }
                 }
                 for(unsigned int iSC=0; iSC<deepSuperCluster_rawEnergy->size(); iSC++){
                     for(unsigned int iPF=0; iPF<deepSuperCluster_pfClustersIndex->at(iSC).size(); iPF++)
                     {  
                         int pfCluster_index = deepSuperCluster_pfClustersIndex->at(iSC).at(iPF);    
                         if(pfCluster_sim_fraction->at(pfCluster_index).at(iCalo) > 0.){
                            Calo_DeepSuperCluster_sim_fraction[iCalo][iSC]+=pfCluster_sim_fraction->at(pfCluster_index).at(iCalo);
                            DeepSuperCluster_Calo_sim_fraction[iSC][iCalo]+=pfCluster_sim_fraction->at(pfCluster_index).at(iCalo);
                         }else{
                            Calo_DeepSuperCluster_sim_fraction[iCalo][iSC]+=0.;
                            DeepSuperCluster_Calo_sim_fraction[iSC][iCalo]+=0.;
                         }
                     }
                 }
             }
   
             std::pair<std::vector<int>,std::vector<int>> matchSuperCluster = matchParticles(&Calo_SuperCluster_sim_fraction, &SuperCluster_Calo_sim_fraction);
             caloParticle_superCluster_sim_fraction_MatchedIndex = matchSuperCluster.first;
             superCluster_sim_fraction_MatchedIndex = matchSuperCluster.second;

             if(superCluster_sim_fraction_MatchedIndex.size()!=superCluster_rawEnergy->size()) std::cout << "SIZE SC: " << superCluster_sim_fraction_MatchedIndex.size() << " - " << superCluster_rawEnergy->size() << std::endl; 

             std::pair<std::vector<int>,std::vector<int>> matchDeepSuperCluster = matchParticles(&Calo_DeepSuperCluster_sim_fraction, &DeepSuperCluster_Calo_sim_fraction);
             caloParticle_deepSuperCluster_sim_fraction_MatchedIndex = matchDeepSuperCluster.first;
             deepSuperCluster_sim_fraction_MatchedIndex = matchDeepSuperCluster.second; 

             if(deepSuperCluster_sim_fraction_MatchedIndex.size()!=deepSuperCluster_rawEnergy->size()) std::cout << "SIZE deepSC: " << deepSuperCluster_sim_fraction_MatchedIndex.size() << " - " << deepSuperCluster_rawEnergy->size() << std::endl;

             //match SC and DeepSC seeds         
             superCluster_seeds.clear();
             for(unsigned int iPF=0; iPF<superCluster_rawEnergy->size(); iPF++)
                 superCluster_seeds.push_back(superCluster_seedIndex->at(iPF));      
             deepSuperCluster_seeds.clear();
             for(unsigned int iPF=0; iPF<deepSuperCluster_rawEnergy->size(); iPF++)
                 deepSuperCluster_seeds.push_back(deepSuperCluster_seedIndex->at(iPF));
             common_seeds.clear(); 
             for(unsigned int iPF=0; iPF<superCluster_seeds.size(); iPF++)
                 for(unsigned int iDeepPF=0; iDeepPF<deepSuperCluster_seeds.size(); iDeepPF++)
                     if(deepSuperCluster_seeds.at(iDeepPF)==superCluster_seeds.at(iPF)) common_seeds.push_back(superCluster_seeds.at(iPF));
            
             //fill total and seeMatched histograms
             for(unsigned int iPF=0; iPF<superCluster_sim_fraction_MatchedIndex.size(); iPF++){                
                if(superCluster_sim_fraction_MatchedIndex.at(iPF)>=0)   
                   SuperCluster_RecoEnergy_sim_fraction[superCluster_sim_fraction_MatchedIndex.at(iPF)]+=superCluster_rawEnergy->at(iPF); 

                h_Eta_old->Fill(superCluster_eta->at(iPF));
                h_nPFClusters_old->Fill(superCluster_nPFClusters->at(iPF)); 
                if(fabs(superCluster_eta->at(iPF))<1.4442){ 
                   h_EtaWidth_EB_old->Fill(superCluster_etaWidth->at(iPF)); 
                   h_Phi_EB_old->Fill(superCluster_phi->at(iPF)); 
                   h_PhiWidth_EB_old->Fill(superCluster_phiWidth->at(iPF));   
             	   h_Energy_EB_old->Fill(superCluster_rawEnergy->at(iPF));  
                   h_nPFClusters_EB_old->Fill(superCluster_nPFClusters->at(iPF));   
                   h_R9_EB_old->Fill(superCluster_r9->at(iPF));    
                   h_full5x5_R9_EB_old->Fill(superCluster_full5x5_r9->at(iPF));   
                   h_sigmaIetaIeta_EB_old->Fill(superCluster_sigmaIetaIeta->at(iPF));   
                   h_full5x5_sigmaIetaIeta_EB_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));      
                   h_sigmaIetaIphi_EB_old->Fill(superCluster_sigmaIetaIphi->at(iPF));   
                   h_full5x5_sigmaIetaIphi_EB_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                   h_sigmaIphiIphi_EB_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                   h_full5x5_sigmaIphiIphi_EB_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));         
                }else if(fabs(superCluster_eta->at(iPF))>1.566){
                   h_EtaWidth_EE_old->Fill(superCluster_etaWidth->at(iPF));
                   h_Phi_EE_old->Fill(superCluster_phi->at(iPF));
                   h_PhiWidth_EE_old->Fill(superCluster_phiWidth->at(iPF)); 
                   h_Energy_EE_old->Fill(superCluster_rawEnergy->at(iPF)); 
                   h_nPFClusters_EE_old->Fill(superCluster_nPFClusters->at(iPF));  
                   h_R9_EE_old->Fill(superCluster_r9->at(iPF));             
                   h_full5x5_R9_EE_old->Fill(superCluster_full5x5_r9->at(iPF));
                   h_sigmaIetaIeta_EE_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                   h_full5x5_sigmaIetaIeta_EE_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                   h_sigmaIetaIphi_EE_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                   h_full5x5_sigmaIetaIphi_EE_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                   h_sigmaIphiIphi_EE_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                   h_full5x5_sigmaIphiIphi_EE_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));                             
                } 
         
                bool sameSeed=false;  
                for(unsigned int iSeed=0; iSeed<common_seeds.size(); iSeed++)
                    if(superCluster_seedIndex->at(iPF) == common_seeds.at(iSeed)){
                       sameSeed=true;
                       break;
                    }    
          
                if(sameSeed==false) continue;
         
                h_Eta_seedMatched_old->Fill(superCluster_eta->at(iPF));
                h_nPFClusters_seedMatched_old->Fill(superCluster_nPFClusters->at(iPF));  
                if(fabs(superCluster_eta->at(iPF))<1.4442){ 
                   h_EtaWidth_EB_seedMatched_old->Fill(superCluster_etaWidth->at(iPF));  
                   h_Phi_EB_seedMatched_old->Fill(superCluster_phi->at(iPF)); 
                   h_PhiWidth_EB_seedMatched_old->Fill(superCluster_phiWidth->at(iPF));  
                   h_Energy_EB_seedMatched_old->Fill(superCluster_rawEnergy->at(iPF));  
                   h_nPFClusters_EB_seedMatched_old->Fill(superCluster_nPFClusters->at(iPF));   
                   h_R9_EB_seedMatched_old->Fill(superCluster_r9->at(iPF));    
                   h_full5x5_R9_EB_seedMatched_old->Fill(superCluster_full5x5_r9->at(iPF));   
                   h_sigmaIetaIeta_EB_seedMatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                   h_full5x5_sigmaIetaIeta_EB_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                   h_sigmaIetaIphi_EB_seedMatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                   h_full5x5_sigmaIetaIphi_EB_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                   h_sigmaIphiIphi_EB_seedMatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                   h_full5x5_sigmaIphiIphi_EB_seedMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));              
                }else if(fabs(superCluster_eta->at(iPF))>1.566){
                   h_EtaWidth_EE_seedMatched_old->Fill(superCluster_etaWidth->at(iPF));
                   h_Phi_EE_seedMatched_old->Fill(superCluster_phi->at(iPF));
                   h_PhiWidth_EE_seedMatched_old->Fill(superCluster_phiWidth->at(iPF)); 
                   h_Energy_EE_seedMatched_old->Fill(superCluster_rawEnergy->at(iPF)); 
                   h_nPFClusters_EE_seedMatched_old->Fill(superCluster_nPFClusters->at(iPF));  
                   h_R9_EE_seedMatched_old->Fill(superCluster_r9->at(iPF));             
                   h_full5x5_R9_EE_seedMatched_old->Fill(superCluster_full5x5_r9->at(iPF));
                   h_sigmaIetaIeta_EE_seedMatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                   h_full5x5_sigmaIetaIeta_EE_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                   h_sigmaIetaIphi_EE_seedMatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                   h_full5x5_sigmaIetaIphi_EE_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                   h_sigmaIphiIphi_EE_seedMatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                   h_full5x5_sigmaIphiIphi_EE_seedMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));                             
                }  
             } 

             for(unsigned int iPF=0; iPF<deepSuperCluster_sim_fraction_MatchedIndex.size(); iPF++){
                 if(deepSuperCluster_sim_fraction_MatchedIndex.at(iPF)>=0)   
                 DeepSuperCluster_RecoEnergy_sim_fraction[deepSuperCluster_sim_fraction_MatchedIndex.at(iPF)]+=deepSuperCluster_rawEnergy->at(iPF); 
          
                 h_Eta_new->Fill(deepSuperCluster_eta->at(iPF));
                 h_nPFClusters_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                 if(fabs(deepSuperCluster_eta->at(iPF))<1.4442){ 
                    h_EtaWidth_EB_new->Fill(deepSuperCluster_etaWidth->at(iPF));  
                    h_Phi_EB_new->Fill(deepSuperCluster_phi->at(iPF)); 
                    h_PhiWidth_EB_new->Fill(deepSuperCluster_phiWidth->at(iPF));  
                    h_Energy_EB_new->Fill(deepSuperCluster_rawEnergy->at(iPF));  
                    h_nPFClusters_EB_new->Fill(deepSuperCluster_nPFClusters->at(iPF));   
                    h_R9_EB_new->Fill(deepSuperCluster_r9->at(iPF));    
                    h_full5x5_R9_EB_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));   
                    h_sigmaIetaIeta_EB_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                    h_full5x5_sigmaIetaIeta_EB_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                    h_sigmaIetaIphi_EB_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                    h_full5x5_sigmaIetaIphi_EB_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                    h_sigmaIphiIphi_EB_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                    h_full5x5_sigmaIphiIphi_EB_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));              
                 }else if(fabs(deepSuperCluster_eta->at(iPF))>1.566){
                    h_EtaWidth_EE_new->Fill(deepSuperCluster_etaWidth->at(iPF));
                    h_Phi_EE_new->Fill(deepSuperCluster_phi->at(iPF));
                    h_PhiWidth_EE_new->Fill(deepSuperCluster_phiWidth->at(iPF)); 
                    h_Energy_EE_new->Fill(deepSuperCluster_rawEnergy->at(iPF)); 
                    h_nPFClusters_EE_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                    h_R9_EE_new->Fill(deepSuperCluster_r9->at(iPF));             
                    h_full5x5_R9_EE_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));
                    h_sigmaIetaIeta_EE_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                    h_full5x5_sigmaIetaIeta_EE_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                    h_sigmaIetaIphi_EE_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                    h_full5x5_sigmaIetaIphi_EE_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                    h_sigmaIphiIphi_EE_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                    h_full5x5_sigmaIphiIphi_EE_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));                             
                 }  
        
                 bool sameSeed=false;  
                 for(unsigned int iSeed=0; iSeed<common_seeds.size(); iSeed++)
                     if(deepSuperCluster_seedIndex->at(iPF) == common_seeds.at(iSeed)){
                        sameSeed=true;
                        break;
                     }    
          
                 if(sameSeed==false) continue;

                 h_Eta_seedMatched_new->Fill(deepSuperCluster_eta->at(iPF));
                 h_nPFClusters_seedMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                 if(fabs(deepSuperCluster_eta->at(iPF))<1.4442){ 
                    h_EtaWidth_EB_seedMatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));  
                    h_Phi_EB_seedMatched_new->Fill(deepSuperCluster_phi->at(iPF)); 
                    h_PhiWidth_EB_seedMatched_new->Fill(deepSuperCluster_phiWidth->at(iPF));  
                    h_Energy_EB_seedMatched_new->Fill(deepSuperCluster_rawEnergy->at(iPF));  
                    h_nPFClusters_EB_seedMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));   
                    h_R9_EB_seedMatched_new->Fill(deepSuperCluster_r9->at(iPF));    
                    h_full5x5_R9_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));   
                    h_sigmaIetaIeta_EB_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                    h_full5x5_sigmaIetaIeta_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                    h_sigmaIetaIphi_EB_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                    h_full5x5_sigmaIetaIphi_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                    h_sigmaIphiIphi_EB_seedMatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                    h_full5x5_sigmaIphiIphi_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));              
                 }else if(fabs(deepSuperCluster_eta->at(iPF))>1.566){
                    h_EtaWidth_EE_seedMatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));
                    h_Phi_EE_seedMatched_new->Fill(deepSuperCluster_phi->at(iPF));
                    h_PhiWidth_EE_seedMatched_new->Fill(deepSuperCluster_phiWidth->at(iPF)); 
                    h_Energy_EE_seedMatched_new->Fill(deepSuperCluster_rawEnergy->at(iPF)); 
                    h_nPFClusters_EE_seedMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                    h_R9_EE_seedMatched_new->Fill(deepSuperCluster_r9->at(iPF));             
                    h_full5x5_R9_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));
                    h_sigmaIetaIeta_EE_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                    h_full5x5_sigmaIetaIeta_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                    h_sigmaIetaIphi_EE_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                    h_full5x5_sigmaIetaIphi_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                    h_sigmaIphiIphi_EE_seedMatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                    h_full5x5_sigmaIphiIphi_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));                             
                 }  
              }

              //fill caloMatched histograms
              for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++){

                  if(caloParticle_superCluster_sim_fraction_MatchedIndex.size()!=caloParticle_simEnergy->size()) continue;

                  h_Eta_Calo_Denum->Fill(caloParticle_simEta->at(iCalo));
                  if(fabs(caloParticle_genEta->at(iCalo))<1.479) h_Et_Calo_EB_Denum->Fill(caloParticle_simEt->at(iCalo));
                  else h_Et_Calo_EE_Denum->Fill(caloParticle_simEt->at(iCalo)); 

                  int iPF=caloParticle_superCluster_sim_fraction_MatchedIndex.at(iCalo); 
                  if(iPF==-1) continue; 

                  h_Eta_Calo_old->Fill(caloParticle_simEta->at(iCalo)); 
                  if(fabs(caloParticle_genEta->at(iCalo))<1.479) h_Et_Calo_EB_old->Fill(caloParticle_simEt->at(iCalo)); 
                  else h_Et_Calo_EE_old->Fill(caloParticle_simEt->at(iCalo)); 

                  h_Eta_caloMatched_old->Fill(superCluster_eta->at(iPF));
                  h_nPFClusters_caloMatched_old->Fill(superCluster_nPFClusters->at(iPF));  

                  if(fabs(caloParticle_genEta->at(iCalo))<1.479){
                     h_EoEtrue_EB_old->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EoEgen_EB_old->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                     h_EtaWidth_EB_caloMatched_old->Fill(superCluster_etaWidth->at(iPF));  
                     h_Phi_EB_caloMatched_old->Fill(superCluster_phi->at(iPF)); 
                     h_PhiWidth_EB_caloMatched_old->Fill(superCluster_phiWidth->at(iPF));  
                     h_Energy_EB_caloMatched_old->Fill(superCluster_rawEnergy->at(iPF));  
                     h_nPFClusters_EB_caloMatched_old->Fill(superCluster_nPFClusters->at(iPF));   
                     h_R9_EB_caloMatched_old->Fill(superCluster_r9->at(iPF));    
                     h_full5x5_R9_EB_caloMatched_old->Fill(superCluster_full5x5_r9->at(iPF));   
                     h_sigmaIetaIeta_EB_caloMatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                     h_full5x5_sigmaIetaIeta_EB_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                     h_sigmaIetaIphi_EB_caloMatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                     h_full5x5_sigmaIetaIphi_EB_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                     h_sigmaIphiIphi_EB_caloMatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                     h_full5x5_sigmaIphiIphi_EB_caloMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF)); 
                  }else if(fabs(caloParticle_genEta->at(iCalo))>=1.479){
                     h_EoEtrue_EE_old->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EoEgen_EE_old->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                     h_EtaWidth_EE_caloMatched_old->Fill(superCluster_etaWidth->at(iPF));
                     h_Phi_EE_caloMatched_old->Fill(superCluster_phi->at(iPF));
                     h_PhiWidth_EE_caloMatched_old->Fill(superCluster_phiWidth->at(iPF)); 
                     h_Energy_EE_caloMatched_old->Fill(superCluster_rawEnergy->at(iPF)); 
                     h_nPFClusters_EE_caloMatched_old->Fill(superCluster_nPFClusters->at(iPF));  
                     h_R9_EE_caloMatched_old->Fill(superCluster_r9->at(iPF));             
                     h_full5x5_R9_EE_caloMatched_old->Fill(superCluster_full5x5_r9->at(iPF));
                     h_sigmaIetaIeta_EE_caloMatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                     h_full5x5_sigmaIetaIeta_EE_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                     h_sigmaIetaIphi_EE_caloMatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                     h_full5x5_sigmaIetaIphi_EE_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                     h_sigmaIphiIphi_EE_caloMatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                     h_full5x5_sigmaIphiIphi_EE_caloMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));          
                  }

                  bool sameSeed=false;  
                  for(unsigned int iSeed=0; iSeed<common_seeds.size(); iSeed++)
                     if(superCluster_seedIndex->at(iPF) == common_seeds.at(iSeed)){
                        sameSeed=true;
                        break;
                     }    
         
                  if(sameSeed==false) continue;

                  int iSeed = superCluster_seedIndex->at(iPF);
                  float seed_et = pfCluster_energy->at(iSeed)/TMath::CosH(pfCluster_eta->at(iSeed));
                  float seed_eta = pfCluster_eta->at(iSeed);
   
                  prof_EoEtrue_vs_Eta_Calo_old->Fill(caloParticle_simEta->at(iCalo),superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  prof_EoEtrue_vs_Eta_Seed_old->Fill(seed_eta,superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  prof_EoEgen_vs_Eta_Gen_old->Fill(caloParticle_genEta->at(iCalo),superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                  prof_EoEgen_vs_Eta_Seed_old->Fill(seed_eta,superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 

                  int nBins_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[0]; 
                  float min_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[1]; 
                  float max_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[2]; 
                  float delta = fabs(min_Eta-max_Eta)/nBins_Eta;
                  int iBin = int((caloParticle_simEta->at(iCalo)-min_Eta)/delta); 
                  if(caloParticle_simEta->at(iCalo)>=min_Eta && caloParticle_simEta->at(iCalo)<max_Eta) EoEtrue_vs_Eta_Calo_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  iBin = int((caloParticle_genEta->at(iCalo)-min_Eta)/delta);  
                  if(caloParticle_genEta->at(iCalo)>=min_Eta && caloParticle_genEta->at(iCalo)<max_Eta) EoEgen_vs_Eta_Gen_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  iBin = int((seed_eta-min_Eta)/delta);  
                  if(seed_eta>=min_Eta && seed_eta<max_Eta){ 
                     EoEtrue_vs_Eta_Seed_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_Eta_Seed_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  }       
                  int nBins_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[0]; 
                  float min_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[1]; 
                  float max_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[2]; 
                  delta = fabs(min_Et-max_Et)/nBins_Et;
                  iBin = int((caloParticle_simEt->at(iCalo)-min_Et)/delta);             
                  if(caloParticle_simEt->at(iCalo)>=min_Et && caloParticle_simEt->at(iCalo)<max_Et && fabs(caloParticle_genEta->at(iCalo))<1.479)
                     EoEtrue_vs_Et_Calo_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                  iBin = int((caloParticle_genEt->at(iCalo)-min_Et)/delta);     
                  if(caloParticle_genEt->at(iCalo)>=min_Et && caloParticle_genEt->at(iCalo)<max_Et && fabs(caloParticle_genEta->at(iCalo))<1.479)  
         	     EoEgen_vs_Et_Gen_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                  iBin = int((seed_et-min_Et)/delta);  
                  if(seed_et>=min_Et && seed_et<max_Et && fabs(seed_eta)<1.479){ 
                     EoEtrue_vs_Et_Seed_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_Et_Seed_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  }        
                  nBins_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[0]; 
                  min_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[1]; 
                  max_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[2]; 
                  delta = fabs(min_Et-max_Et)/nBins_Et;
                  iBin = int((caloParticle_simEt->at(iCalo)-min_Et)/delta);
                  if(caloParticle_simEt->at(iCalo)>=min_Et && caloParticle_simEt->at(iCalo)<max_Et && fabs(caloParticle_genEta->at(iCalo))>1.479)
                     EoEtrue_vs_Et_Calo_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));   
                  iBin = int((caloParticle_genEt->at(iCalo)-min_Et)/delta);
                  if(caloParticle_genEt->at(iCalo)>=min_Et && caloParticle_genEt->at(iCalo)<max_Et && fabs(caloParticle_genEta->at(iCalo))>1.479)              
         	     EoEgen_vs_Et_Gen_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));     
                  iBin = int((seed_et-min_Et)/delta);  
                  if(seed_et>=min_Et && seed_et<max_Et && fabs(seed_eta)>1.479){ 
                     EoEtrue_vs_Et_Seed_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_Et_Seed_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  }  

                  int seed_etBin = h2_EoEtrue_Mean_old->GetXaxis()->FindBin(seed_et)-1;  
                  int seed_etaBin = h2_EoEtrue_Mean_old->GetYaxis()->FindBin(fabs(seed_eta))-1; 
 
                  if(seed_et>=100) seed_etBin = h2_EoEtrue_Mean_old->GetNbinsX()-1;  
                  if(fabs(seed_eta)>=3.) seed_etaBin = h2_EoEtrue_Mean_old->GetNbinsY()-1; 
                  //std::cout << "Old Et - Eta - " << seed_etBin << " - " << seed_et << " - " << fabs(seed_eta) << " - " << seed_etaBin << std::endl;
                  EoEtrue_vs_seedEt_seedEta_old[seed_etBin][seed_etaBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 

                  int nBins_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[0]; 
                  float min_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[1]; 
                  float max_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_Energy-max_Energy)/nBins_Energy;
                  iBin = int((caloParticle_simEnergy->at(iCalo)-min_Energy)/delta);
                  if(caloParticle_simEnergy->at(iCalo)>=min_Energy && caloParticle_simEnergy->at(iCalo)<max_Energy && fabs(caloParticle_genEta->at(iCalo))<1.479)
                     EoEtrue_vs_Energy_Calo_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                  iBin = int((caloParticle_genEnergy->at(iCalo)-min_Energy)/delta);    
                  if(caloParticle_genEnergy->at(iCalo)>=min_Energy && caloParticle_genEnergy->at(iCalo)<max_Energy && fabs(caloParticle_genEta->at(iCalo))<1.479)    
                     EoEgen_vs_Energy_Gen_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));                  
                  nBins_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[0]; 
                  min_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[1]; 
                  max_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_Energy-max_Energy)/nBins_Energy;
                  iBin = int((caloParticle_simEnergy->at(iCalo)-min_Energy)/delta);
                  if(caloParticle_simEnergy->at(iCalo)>=min_Energy && caloParticle_simEnergy->at(iCalo)<max_Energy && fabs(caloParticle_genEta->at(iCalo))>1.479)
                     EoEtrue_vs_Energy_Calo_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  iBin = int((caloParticle_genEnergy->at(iCalo)-min_Energy)/delta);
                  if(caloParticle_genEnergy->at(iCalo)>=min_Energy && caloParticle_genEnergy->at(iCalo)<max_Energy && fabs(caloParticle_genEta->at(iCalo))>1.479)  
 	             EoEgen_vs_Energy_Gen_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
 
                  int nBins_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[0]; 
                  float min_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[1]; 
                  float max_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_nVtx-max_nVtx)/nBins_nVtx;
                  iBin = int((nVtx-min_nVtx)/delta);
                  if(nVtx>=min_nVtx && nVtx<max_nVtx && fabs(caloParticle_genEta->at(iCalo))<1.479){
                     EoEtrue_vs_nVtx_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_nVtx_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));     
                  }   
                  nBins_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[0]; 
                  min_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[1]; 
                  max_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_nVtx-max_nVtx)/nBins_nVtx;
                  iBin = int((nVtx-min_nVtx)/delta);
                  if(nVtx>=min_nVtx && nVtx<max_nVtx && fabs(caloParticle_genEta->at(iCalo))>1.479){
                     EoEtrue_vs_nVtx_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_nVtx_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                  }    

                  int nBins_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[0]; 
                  float min_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[1]; 
                  float max_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_Rho-max_Rho)/nBins_Rho;
                  iBin = int((rho-min_Rho)/delta);
                  if(rho>=min_Rho && rho<max_Rho && fabs(caloParticle_genEta->at(iCalo))<1.479){
                     EoEtrue_vs_Rho_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                     EoEgen_vs_Rho_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));        
                  }  
                  nBins_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[0]; 
                  min_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[1]; 
                  max_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_Rho-max_Rho)/nBins_Rho;
                  iBin = int((rho-min_Rho)/delta);
                  if(rho>=min_Rho && rho<max_Rho && fabs(caloParticle_genEta->at(iCalo))>1.479){
                     EoEtrue_vs_Rho_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_Rho_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));    
                  } 


                  int nBins_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Barrel"),binOpts)].second[0]; 
                  float min_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Barrel"),binOpts)].second[1]; 
                  float max_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_EtaWidth-max_EtaWidth)/nBins_EtaWidth;
                  iBin = int((superCluster_etaWidth->at(iPF)-min_EtaWidth)/delta);
                  if(superCluster_etaWidth->at(iPF)>=min_EtaWidth && superCluster_etaWidth->at(iPF)<max_EtaWidth && fabs(caloParticle_genEta->at(iCalo))<1.479){
                     EoEtrue_vs_EtaWidth_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                     //EoEgen_vs_EtaWidth_EB_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));        
                  } 
                  nBins_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Endcap"),binOpts)].second[0]; 
                  min_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Endcap"),binOpts)].second[1]; 
                  max_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_EtaWidth-max_EtaWidth)/nBins_EtaWidth;
                  iBin = int((superCluster_etaWidth->at(iPF)-min_EtaWidth)/delta);
                  if(superCluster_etaWidth->at(iPF)>=min_EtaWidth && superCluster_etaWidth->at(iPF)<max_EtaWidth && fabs(caloParticle_genEta->at(iCalo))>1.479){
                     EoEtrue_vs_EtaWidth_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     //EoEgen_vs_EtaWidth_EE_old[iBin]->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));    
                  } 


                  if(fabs(caloParticle_genEta->at(iCalo))<1.479){
                     h_EoEtrue_EB_seedMatched_old->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EoEgen_EB_seedMatched_old->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                     prof_EoEtrue_vs_Et_Calo_EB_old->Fill(caloParticle_simEt->at(iCalo),superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEtrue_vs_Energy_Calo_EB_old->Fill(caloParticle_simEnergy->at(iCalo),superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));  
                     prof_EoEtrue_vs_nVtx_EB_old->Fill(nVtx,superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEtrue_vs_Rho_EB_old->Fill(rho,superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEgen_vs_Et_Gen_EB_old->Fill(caloParticle_genEt->at(iCalo),superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                     prof_EoEgen_vs_Energy_Gen_EB_old->Fill(caloParticle_genEnergy->at(iCalo),superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));  
                     prof_EoEgen_vs_nVtx_EB_old->Fill(nVtx,superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                     prof_EoEgen_vs_Rho_EB_old->Fill(rho,superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  }else if(fabs(caloParticle_genEta->at(iCalo))>=1.479){
                     h_EoEtrue_EE_seedMatched_old->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EoEgen_EE_seedMatched_old->Fill(superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                     prof_EoEtrue_vs_Et_Calo_EE_old->Fill(caloParticle_simEt->at(iCalo),superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEtrue_vs_Energy_Calo_EE_old->Fill(caloParticle_simEnergy->at(iCalo),superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));  
                     prof_EoEtrue_vs_nVtx_EE_old->Fill(nVtx,superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEtrue_vs_Rho_EE_old->Fill(rho,superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEgen_vs_Et_Gen_EE_old->Fill(caloParticle_genEt->at(iCalo),superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                     prof_EoEgen_vs_Energy_Gen_EE_old->Fill(caloParticle_genEnergy->at(iCalo),superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));  
                     prof_EoEgen_vs_nVtx_EE_old->Fill(nVtx,superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                     prof_EoEgen_vs_Rho_EE_old->Fill(rho,superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  }  

                  if(fabs(seed_eta)<1.479){
                     prof_EoEtrue_vs_Et_Seed_EB_old->Fill(seed_et,superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEgen_vs_Et_Seed_EB_old->Fill(seed_et,superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));  
                  }else if(fabs(seed_eta)>1.479){
                     prof_EoEtrue_vs_Et_Seed_EE_old->Fill(seed_et,superCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEgen_vs_Et_Seed_EE_old->Fill(seed_et,superCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  } 
               }

               for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++){

                  if(caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.size()!=caloParticle_simEnergy->size()) continue;

                  int iPF=caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.at(iCalo); 
                  if(iPF==-1) continue; 

                  h_Eta_Calo_new->Fill(caloParticle_simEta->at(iCalo)); 
                  if(fabs(caloParticle_genEta->at(iCalo))<1.479) h_Et_Calo_EB_new->Fill(caloParticle_simEt->at(iCalo)); 
                  else h_Et_Calo_EE_new->Fill(caloParticle_simEt->at(iCalo)); 

                  h_Eta_caloMatched_new->Fill(deepSuperCluster_eta->at(iPF));
                  h_nPFClusters_caloMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));
                  
                  if(fabs(caloParticle_genEta->at(iCalo))<1.479){
                     h_EoEtrue_EB_new->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EoEgen_EB_new->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                     h_EtaWidth_EB_caloMatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));  
                     h_Phi_EB_caloMatched_new->Fill(deepSuperCluster_phi->at(iPF)); 
                     h_PhiWidth_EB_caloMatched_new->Fill(deepSuperCluster_phiWidth->at(iPF));  
                     h_Energy_EB_caloMatched_new->Fill(deepSuperCluster_rawEnergy->at(iPF));  
                     h_nPFClusters_EB_caloMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));   
                     h_R9_EB_caloMatched_new->Fill(deepSuperCluster_r9->at(iPF));    
                     h_full5x5_R9_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));   
                     h_sigmaIetaIeta_EB_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                     h_full5x5_sigmaIetaIeta_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                     h_sigmaIetaIphi_EB_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                     h_full5x5_sigmaIetaIphi_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                     h_sigmaIphiIphi_EB_caloMatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                     h_full5x5_sigmaIphiIphi_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF)); 
                  }else if(fabs(caloParticle_genEta->at(iCalo))>=1.479){
                     h_EoEtrue_EE_new->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EoEgen_EE_new->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                     h_EtaWidth_EE_caloMatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));
                     h_Phi_EE_caloMatched_new->Fill(deepSuperCluster_phi->at(iPF));
                     h_PhiWidth_EE_caloMatched_new->Fill(deepSuperCluster_phiWidth->at(iPF)); 
                     h_Energy_EE_caloMatched_new->Fill(deepSuperCluster_rawEnergy->at(iPF)); 
                     h_nPFClusters_EE_caloMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                     h_R9_EE_caloMatched_new->Fill(deepSuperCluster_r9->at(iPF));             
                     h_full5x5_R9_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));
                     h_sigmaIetaIeta_EE_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                     h_full5x5_sigmaIetaIeta_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                     h_sigmaIetaIphi_EE_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                     h_full5x5_sigmaIetaIphi_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                     h_sigmaIphiIphi_EE_caloMatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                     h_full5x5_sigmaIphiIphi_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));   
                  }

                  bool sameSeed=false;  
                  for(unsigned int iSeed=0; iSeed<common_seeds.size(); iSeed++)
                     if(deepSuperCluster_seedIndex->at(iPF) == common_seeds.at(iSeed)){
                        sameSeed=true;
                        break;
                     } 

                  if(sameSeed==false) continue;

                  int iSeed = deepSuperCluster_seedIndex->at(iPF);
                  float seed_et = pfCluster_energy->at(iSeed)/TMath::CosH(pfCluster_eta->at(iSeed));
                  float seed_eta = pfCluster_eta->at(iSeed);
                 
                  prof_EoEtrue_vs_Eta_Calo_new->Fill(caloParticle_simEta->at(iCalo),deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  prof_EoEtrue_vs_Eta_Seed_new->Fill(seed_eta,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  prof_EoEgen_vs_Eta_Gen_new->Fill(caloParticle_genEta->at(iCalo),deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                  prof_EoEgen_vs_Eta_Seed_new->Fill(seed_eta,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
   
                  int nBins_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[0]; 
                  float min_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[1]; 
                  float max_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[2]; 
                  float delta = fabs(min_Eta-max_Eta)/nBins_Eta;
                  int iBin = int((caloParticle_simEta->at(iCalo)-min_Eta)/delta); 
                  if(caloParticle_simEta->at(iCalo)>=min_Eta && caloParticle_simEta->at(iCalo)<max_Eta) EoEtrue_vs_Eta_Calo_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  iBin = int((caloParticle_genEta->at(iCalo)-min_Eta)/delta);  
                  if(caloParticle_genEta->at(iCalo)>=min_Eta && caloParticle_genEta->at(iCalo)<max_Eta) EoEgen_vs_Eta_Gen_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  iBin = int((seed_eta-min_Eta)/delta);  
                  if(seed_eta>=min_Eta && seed_eta<max_Eta){ 
                     EoEtrue_vs_Eta_Seed_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_Eta_Seed_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  }  
                  
                  int nBins_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[0]; 
                  float min_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[1]; 
                  float max_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[2]; 
                  delta = fabs(min_Et-max_Et)/nBins_Et;
                  iBin = int((caloParticle_simEt->at(iCalo)-min_Et)/delta);             
                  if(caloParticle_simEt->at(iCalo)>=min_Et && caloParticle_simEt->at(iCalo)<max_Et && fabs(caloParticle_genEta->at(iCalo))<1.479)
                     EoEtrue_vs_Et_Calo_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                  iBin = int((caloParticle_genEt->at(iCalo)-min_Et)/delta);     
                  if(caloParticle_genEt->at(iCalo)>=min_Et && caloParticle_genEt->at(iCalo)<max_Et && fabs(caloParticle_genEta->at(iCalo))<1.479)  
         	     EoEgen_vs_Et_Gen_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  iBin = int((seed_et-min_Et)/delta);   
                  if(seed_et>=min_Et && seed_et<max_Et && fabs(seed_eta)<1.479){ 
                     EoEtrue_vs_Et_Seed_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_Et_Seed_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  }        
                  nBins_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[0]; 
                  min_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[1]; 
                  max_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[2]; 
                  delta = fabs(min_Et-max_Et)/nBins_Et;
                  iBin = int((caloParticle_simEt->at(iCalo)-min_Et)/delta);
                  if(caloParticle_simEt->at(iCalo)>=min_Et && caloParticle_simEt->at(iCalo)<max_Et && fabs(caloParticle_genEta->at(iCalo))>1.479)
                     EoEtrue_vs_Et_Calo_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));   
                  iBin = int((caloParticle_genEt->at(iCalo)-min_Et)/delta);
                  if(caloParticle_genEt->at(iCalo)>=min_Et && caloParticle_genEt->at(iCalo)<max_Et && fabs(caloParticle_genEta->at(iCalo))>1.479)              
         	     EoEgen_vs_Et_Gen_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));       
                  iBin = int((seed_et-min_Et)/delta);   
                  if(seed_et>=min_Et && seed_et<max_Et && fabs(seed_eta)>1.479){ 
                     EoEtrue_vs_Et_Seed_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_Et_Seed_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  }  
 
                  int seed_etBin = h2_EoEtrue_Mean_new->GetXaxis()->FindBin(seed_et)-1;  
                  int seed_etaBin = h2_EoEtrue_Mean_new->GetYaxis()->FindBin(fabs(seed_eta))-1;  
                  if(seed_et>=100) seed_etBin = h2_EoEtrue_Mean_new->GetNbinsX()-1;  
                  if(fabs(seed_eta)>=3.) seed_etaBin = h2_EoEtrue_Mean_new->GetNbinsY()-1;  
                  //std::cout << "New Et - Eta - " << seed_etBin << " - " << seed_et << " - " << fabs(seed_eta) << " - " << seed_etaBin << std::endl;
                  EoEtrue_vs_seedEt_seedEta_new[seed_etBin][seed_etaBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
   
                  int nBins_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[0]; 
                  float min_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[1]; 
                  float max_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_Energy-max_Energy)/nBins_Energy;
                  iBin = int((caloParticle_simEnergy->at(iCalo)-min_Energy)/delta);
                  if(caloParticle_simEnergy->at(iCalo)>=min_Energy && caloParticle_simEnergy->at(iCalo)<max_Energy && fabs(caloParticle_genEta->at(iCalo))<1.479)
                     EoEtrue_vs_Energy_Calo_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                  iBin = int((caloParticle_genEnergy->at(iCalo)-min_Energy)/delta);    
                  if(caloParticle_genEnergy->at(iCalo)>=min_Energy && caloParticle_genEnergy->at(iCalo)<max_Energy && fabs(caloParticle_genEta->at(iCalo))<1.479)    
                     EoEgen_vs_Energy_Gen_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));                  
                  nBins_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[0]; 
                  min_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[1]; 
                  max_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_Energy-max_Energy)/nBins_Energy;
                  iBin = int((caloParticle_simEnergy->at(iCalo)-min_Energy)/delta);
                  if(caloParticle_simEnergy->at(iCalo)>=min_Energy && caloParticle_simEnergy->at(iCalo)<max_Energy && fabs(caloParticle_genEta->at(iCalo))>1.479)
                     EoEtrue_vs_Energy_Calo_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  iBin = int((caloParticle_genEnergy->at(iCalo)-min_Energy)/delta);
                  if(caloParticle_genEnergy->at(iCalo)>=min_Energy && caloParticle_genEnergy->at(iCalo)<max_Energy && fabs(caloParticle_genEta->at(iCalo))>1.479)  
 	             EoEgen_vs_Energy_Gen_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
 
                  int nBins_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[0]; 
                  float min_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[1]; 
                  float max_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_nVtx-max_nVtx)/nBins_nVtx;
                  iBin = int((nVtx-min_nVtx)/delta);
                  if(nVtx>=min_nVtx && nVtx<max_nVtx && fabs(caloParticle_genEta->at(iCalo))<1.479){
                     EoEtrue_vs_nVtx_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_nVtx_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));     
                  }   
                  nBins_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[0]; 
                  min_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[1]; 
                  max_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_nVtx-max_nVtx)/nBins_nVtx;
                  iBin = int((nVtx-min_nVtx)/delta);
                  if(nVtx>=min_nVtx && nVtx<max_nVtx && fabs(caloParticle_genEta->at(iCalo))>1.479){
                     EoEtrue_vs_nVtx_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_nVtx_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                  }    

                  int nBins_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[0]; 
                  float min_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[1]; 
                  float max_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_Rho-max_Rho)/nBins_Rho;
                  iBin = int((rho-min_Rho)/delta);
                  if(rho>=min_Rho && rho<max_Rho && fabs(caloParticle_genEta->at(iCalo))<1.479){
                     EoEtrue_vs_Rho_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                     EoEgen_vs_Rho_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));        
                  }  
                  nBins_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[0]; 
                  min_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[1]; 
                  max_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_Rho-max_Rho)/nBins_Rho;
                  iBin = int((rho-min_Rho)/delta);
                  if(rho>=min_Rho && rho<max_Rho && fabs(caloParticle_genEta->at(iCalo))>1.479){
                     EoEtrue_vs_Rho_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     EoEgen_vs_Rho_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));    
                  } 



                  int nBins_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Barrel"),binOpts)].second[0]; 
                  float min_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Barrel"),binOpts)].second[1]; 
                  float max_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_EtaWidth-max_EtaWidth)/nBins_EtaWidth;
                  iBin = int((deepSuperCluster_etaWidth->at(iPF)-min_EtaWidth)/delta);
                  if(deepSuperCluster_etaWidth->at(iPF)>=min_EtaWidth && deepSuperCluster_etaWidth->at(iPF)<max_EtaWidth && fabs(caloParticle_genEta->at(iCalo))<1.479){
                     EoEtrue_vs_EtaWidth_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                     //EoEgen_vs_EtaWidth_EB_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));        
                  } 
                  nBins_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Endcap"),binOpts)].second[0]; 
                  min_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Endcap"),binOpts)].second[1]; 
                  max_EtaWidth = binOpts[findOption(std::string("EtaWidthBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_EtaWidth-max_EtaWidth)/nBins_EtaWidth;
                  iBin = int((deepSuperCluster_etaWidth->at(iPF)-min_EtaWidth)/delta);
                  if(deepSuperCluster_etaWidth->at(iPF)>=min_EtaWidth && deepSuperCluster_etaWidth->at(iPF)<max_EtaWidth && fabs(caloParticle_genEta->at(iCalo))>1.479){
                     EoEtrue_vs_EtaWidth_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     //EoEgen_vs_EtaWidth_EE_new[iBin]->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));    
                  } 

                  if(fabs(caloParticle_genEta->at(iCalo))<1.479){
                     h_EoEtrue_EB_seedMatched_new->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EoEgen_EB_seedMatched_new->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                     prof_EoEtrue_vs_Et_Calo_EB_new->Fill(caloParticle_simEt->at(iCalo),deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEtrue_vs_Energy_Calo_EB_new->Fill(caloParticle_simEnergy->at(iCalo),deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));  
                     prof_EoEtrue_vs_nVtx_EB_new->Fill(nVtx,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEtrue_vs_Rho_EB_new->Fill(rho,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEgen_vs_Et_Gen_EB_new->Fill(caloParticle_genEt->at(iCalo),deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                     prof_EoEgen_vs_Energy_Gen_EB_new->Fill(caloParticle_genEnergy->at(iCalo),deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));  
                     prof_EoEgen_vs_nVtx_EB_new->Fill(nVtx,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                     prof_EoEgen_vs_Rho_EB_new->Fill(rho,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  }else if(fabs(caloParticle_genEta->at(iCalo))>=1.479){
                     h_EoEtrue_EE_seedMatched_new->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EoEgen_EE_seedMatched_new->Fill(deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                     prof_EoEtrue_vs_Et_Calo_EE_new->Fill(caloParticle_simEt->at(iCalo),deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEtrue_vs_Energy_Calo_EE_new->Fill(caloParticle_simEnergy->at(iCalo),deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));  
                     prof_EoEtrue_vs_nVtx_EE_new->Fill(nVtx,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEtrue_vs_Rho_EE_new->Fill(rho,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEgen_vs_Et_Gen_EE_new->Fill(caloParticle_genEt->at(iCalo),deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo)); 
                     prof_EoEgen_vs_Energy_Gen_EE_new->Fill(caloParticle_genEnergy->at(iCalo),deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));  
                     prof_EoEgen_vs_nVtx_EE_new->Fill(nVtx,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                     prof_EoEgen_vs_Rho_EE_new->Fill(rho,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  } 

                  if(fabs(seed_eta)<1.479){
                     prof_EoEtrue_vs_Et_Seed_EB_new->Fill(seed_et,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEgen_vs_Et_Seed_EB_new->Fill(seed_et,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));  
                  }else if(fabs(seed_eta)>1.479){
                     prof_EoEtrue_vs_Et_Seed_EE_new->Fill(seed_et,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEgen_vs_Et_Seed_EE_new->Fill(seed_et,deepSuperCluster_rawEnergy->at(iPF)/caloParticle_genEnergy->at(iCalo));
                  } 
               }

               //fill caloUnmatched histos
               for(unsigned int iPF=0; iPF<superCluster_rawEnergy->size(); iPF++){

                   std::vector<int>::iterator it = std::find(caloParticle_superCluster_sim_fraction_MatchedIndex.begin(), caloParticle_superCluster_sim_fraction_MatchedIndex.end(), iPF); 
                   if(it != caloParticle_superCluster_sim_fraction_MatchedIndex.end()) continue;   
 
                   h_Eta_caloUnmatched_old->Fill(superCluster_eta->at(iPF));
                   h_nPFClusters_caloUnmatched_old->Fill(superCluster_nPFClusters->at(iPF));  
                   if(fabs(superCluster_eta->at(iPF))<1.4442){ 
                      h_EtaWidth_EB_caloUnmatched_old->Fill(superCluster_etaWidth->at(iPF));  
                      h_Phi_EB_caloUnmatched_old->Fill(superCluster_phi->at(iPF)); 
                      h_PhiWidth_EB_caloUnmatched_old->Fill(superCluster_phiWidth->at(iPF));  
                      h_Energy_EB_caloUnmatched_old->Fill(superCluster_rawEnergy->at(iPF));  
                      h_nPFClusters_EB_caloUnmatched_old->Fill(superCluster_nPFClusters->at(iPF));   
                      h_R9_EB_caloUnmatched_old->Fill(superCluster_r9->at(iPF));    
                      h_full5x5_R9_EB_caloUnmatched_old->Fill(superCluster_full5x5_r9->at(iPF));   
                      h_sigmaIetaIeta_EB_caloUnmatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                      h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                      h_sigmaIetaIphi_EB_caloUnmatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                      h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                      h_sigmaIphiIphi_EB_caloUnmatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                      h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));              
                   }else if(fabs(superCluster_eta->at(iPF))>1.566){
                      h_EtaWidth_EE_caloUnmatched_old->Fill(superCluster_etaWidth->at(iPF));
             	      h_Phi_EE_caloUnmatched_old->Fill(superCluster_phi->at(iPF));
                      h_PhiWidth_EE_caloUnmatched_old->Fill(superCluster_phiWidth->at(iPF)); 
                      h_Energy_EE_caloUnmatched_old->Fill(superCluster_rawEnergy->at(iPF)); 
                      h_nPFClusters_EE_caloUnmatched_old->Fill(superCluster_nPFClusters->at(iPF));  
                      h_R9_EE_caloUnmatched_old->Fill(superCluster_r9->at(iPF));             
                      h_full5x5_R9_EE_caloUnmatched_old->Fill(superCluster_full5x5_r9->at(iPF));
                      h_sigmaIetaIeta_EE_caloUnmatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                      h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                      h_sigmaIetaIphi_EE_caloUnmatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                      h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                      h_sigmaIphiIphi_EE_caloUnmatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                      h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));                             
                   } 
               }

               for(unsigned int iPF=0; iPF<deepSuperCluster_rawEnergy->size(); iPF++){

                   std::vector<int>::iterator it = std::find(caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.begin(), caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.end(), iPF); 
                   if(it != caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.end()) continue;  

                   h_Eta_caloUnmatched_new->Fill(deepSuperCluster_eta->at(iPF));
                   h_nPFClusters_caloUnmatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  

                   if(fabs(deepSuperCluster_eta->at(iPF))<1.4442){ 
                      h_EtaWidth_EB_caloUnmatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));  
                      h_Phi_EB_caloUnmatched_new->Fill(deepSuperCluster_phi->at(iPF)); 
                      h_PhiWidth_EB_caloUnmatched_new->Fill(deepSuperCluster_phiWidth->at(iPF));  
                      h_Energy_EB_caloUnmatched_new->Fill(deepSuperCluster_rawEnergy->at(iPF));  
                      h_nPFClusters_EB_caloUnmatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));   
                      h_R9_EB_caloUnmatched_new->Fill(deepSuperCluster_r9->at(iPF));    
                      h_full5x5_R9_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));   
                      h_sigmaIetaIeta_EB_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                      h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                      h_sigmaIetaIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                      h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                      h_sigmaIphiIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                      h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));              
                   }else if(fabs(deepSuperCluster_eta->at(iPF))>1.566){
                      h_EtaWidth_EE_caloUnmatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));
                      h_Phi_EE_caloUnmatched_new->Fill(deepSuperCluster_phi->at(iPF));
                      h_PhiWidth_EE_caloUnmatched_new->Fill(deepSuperCluster_phiWidth->at(iPF)); 
                      h_Energy_EE_caloUnmatched_new->Fill(deepSuperCluster_rawEnergy->at(iPF)); 
                      h_nPFClusters_EE_caloUnmatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                      h_R9_EE_caloUnmatched_new->Fill(deepSuperCluster_r9->at(iPF));             
                      h_full5x5_R9_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));
                      h_sigmaIetaIeta_EE_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                      h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                      h_sigmaIetaIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                      h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                      h_sigmaIphiIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                      h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));  
                   }           
                }       
          }
       }
    }
  
    setEfficiencies();
    
    if(superClusterRef_ == "superCluster") superClusterRef_ = "Mustache";
    if(superClusterRef_ == "retunedSuperCluster") superClusterRef_ = "retunedMustache";
    if(superClusterRef_ == "deepSuperCluster") superClusterRef_ = "DeepSC";
    if(superClusterRef_ == "deepSuperClusterLWP") superClusterRef_ = "deepSC LWP";
    if(superClusterRef_ == "deepSuperClusterTWP") superClusterRef_ = "deepSC TWP";
    if(superClusterVal_ == "superCluster") superClusterRef_ = "Mustache";
    if(superClusterVal_ == "retunedSuperCluster") superClusterVal_ = "retunedMustache";
    if(superClusterVal_ == "deepSuperCluster") superClusterVal_ = "DeepSC";
    if(superClusterVal_ == "deepSuperClusterLWP") superClusterVal_ = "deepSC LWP";
    if(superClusterVal_ == "deepSuperClusterTWP") superClusterVal_ = "deepSC TWP";
    
    //drawPlots(fitFunction_,superClusterRef_,superClusterVal_);
    savePlots();
   
}
