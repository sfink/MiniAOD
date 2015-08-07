#include "MiniAOD/MiniAODHelper/interface/TopTagger.h"

TopTagger::TopTagger(TopTag::Mode mode_, TopTag::SubjetAssign subjetAssign_, std::string filePath_):helper(new MiniAODHelper()),mode(mode_),subjetAssign(subjetAssign_),btagger("pfCombinedInclusiveSecondaryVertexV2BJetTags"){
  
  char* CMSSWPath = getenv("CMSSW_BASE");
  std::string filePath = CMSSWPath;
  filePath = filePath+"/src/MiniAOD/MiniAODHelper/data/toptagger/"+filePath_;
  
  switch(mode){
    case TopTag::HEP:
      break;
      
    case TopTag::Likelihood:
      {
        file = new TFile(filePath.c_str());

        mtop_top_histo=(TH1F*)file->Get("TopJet_Top_M_True");
        mtop_nottop_histo=(TH1F*)file->Get("TopJet_Top_M_False");

        if(subjetAssign == TopTag::Wmass){
          mratio_top_histo=(TH1F*)file->Get("TopJet_MRatio_W_Top_True");
          mratio_nottop_histo=(TH1F*)file->Get("TopJet_MRatio_W_Top_False");
          atan_top_histo=(TH1F*)file->Get("TopJet_MRatio_W_Top_True");
          atan_nottop_histo=(TH1F*)file->Get("TopJet_MRatio_W_Top_False");
        }
        else if(subjetAssign == TopTag::CSV){
          mratio_top_histo=(TH1F*)file->Get("TopJet_MRatio_Wbtag_Top_True");
          mratio_nottop_histo=(TH1F*)file->Get("TopJet_MRatio_Wbtag_Top_False");
          atan_top_histo=(TH1F*)file->Get("TopJet_MRatio_Wbtag_Top_True");
          atan_nottop_histo=(TH1F*)file->Get("TopJet_MRatio_Wbtag_Top_False");
        }
        else{
          std::cout << "Error! No matching Subjet Assignment found!" << std:: endl;
        }
      }
      break;
    
    case TopTag::TMVA:
      {  
        // Get absolute path for TMVA weights
        std::string weightsPath = filePath;

        // Initialize TMVA input variables
        GetTMVAVars(weightsPath);

        // Setup TMVA Reader
        TMVAReader = new TMVA::Reader("Silent");
        
        for(std::vector<string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
          TMVAReader->AddVariable(*itVarName,&TMVAVars[itVarName-TMVAVarNames.begin()]);
        }

        TMVAReader->BookMVA("TMVATopTagger",(weightsPath).c_str());
      }      
      break;

    default:
      std::cout << "Error! No matching Top Tagger Mode found!" << std:: endl;
  } 
}


TopTagger::~TopTagger(){
}


void TopTagger::TopSubjetCSVDef(std::vector<pat::Jet> &subjets){
  
  subjets = helper->GetSortedByCSV(subjets);
  
  pat::Jet Bsubjet = subjets[0];
  subjets.erase(subjets.begin());
  
  subjets = helper->GetSortedByPt(subjets);
  pat::Jet W1subjet = subjets[0];
  pat::Jet W2subjet = subjets[1];
  
  subjets.clear();
  subjets.push_back(Bsubjet);
  subjets.push_back(W1subjet);
  subjets.push_back(W2subjet);
}


void TopTagger::GetTMVAVarNames(std::string filePath_, bool verbose){
  
  TMVAVarNames.clear();
  
  std::string line;
  ifstream inFile (filePath_);
  
  if(inFile.is_open()){
  
    while(getline(inFile,line)){
      if(line.find("</Variables>")!=std::string::npos) break;
      
      if(line.find("<Variable VarIndex")!=std::string::npos){
        size_t varStart = line.find("\"",line.find("Title"))+1;
        size_t varEnd = line.find("\"",varStart);

        TMVAVarNames.push_back(line.substr(varStart,varEnd-varStart));
      }
    }
    
    if(verbose){
      for(std::vector<string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
        std::cout << "Variable " << itVarName-TMVAVarNames.begin() << ": " << *itVarName << std::endl;
      }
    }
    
    inFile.close();
  }
}


void TopTagger::GetTMVAVars(std::string filePath_, bool verbose){
  
  GetTMVAVarNames(filePath_,verbose);
  
  TMVAVars = new float[50];
  
  for(int i=0;i<50;i++){
    TMVAVars[i] = -999;
  }
  
  if(verbose){
    for(std::vector<string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
      std::cout << "Variable " << itVarName-TMVAVarNames.begin() << ": " << *itVarName << ": " << TMVAVars[itVarName-TMVAVarNames.begin()] << std::endl;
    }
  }  
}


void TopTagger::ResetTMVAVars(){
  for(int i=0;i<50;i++){
    TMVAVars[i] = -999;
  }
}


float TopTagger::GetTopTaggerOutput(const boosted::HTTTopJet& topjet, bool verbose){
  
  if(topjet.nonW.pt()==0.) return -1.1;
  if(topjet.W1.pt()==0.) return -1.1;
  if(topjet.W2.pt()==0.) return -1.1;
  
  std::vector<pat::Jet> subjets;
  subjets.push_back(topjet.nonW);
  subjets.push_back(topjet.W1);
  subjets.push_back(topjet.W2);
  
  if(subjetAssign == TopTag::Pt){
    subjets = helper->GetSortedByPt(subjets);
  }
  else if(subjetAssign == TopTag::CSV){
    TopSubjetCSVDef(subjets);
  }
  
  switch(mode){
    case TopTag::HEP:
      {
        float m12 = -999;
        float m23 = -999;
        float m13 = -999;
        float m123 = -999;

        m12 = (subjets[0].p4()+subjets[1].p4()).M();
        m13 = (subjets[0].p4()+subjets[2].p4()).M();
        m23 = (subjets[1].p4()+subjets[2].p4()).M();
        m123 = (subjets[0].p4()+subjets[1].p4()+subjets[2].p4()).M();

        //if(m123<mTopMin) return false;

        float mT = 172.3;
        float mW = 80.4;

        if(subjetAssign == TopTag::Pt){
          float y1 = 1+pow(m13/m12,2);
          float y2 = 1+pow(m12/m13,2);
          float x = 1-pow(m23/m123,2);
          
          float fw1 = fabs((m23/m123*mT/mW)-1);
          float fw2 = fabs((sqrt(x/y1)*mT/mW)-1);
          float fw3 = fabs((sqrt(x/y2)*mT/mW)-1);
            
          float fw = 999;
          
          if((0.2<atan(m23/m123)) && (atan(m23/m123)<1.3)) fw = fmin(fw,fw1);
          if(m23/m123>0.35){
            fw = fmin(fw,fw2);
            fw = fmin(fw,fw3);
          }
          
          if(fw<999) return 1-fw;
          
        }
        else if(subjetAssign == TopTag::CSV){
          float fw = fabs((m23/m123*mT/mW)-1);
          
          if(0.2<m12/m13 && 0.2<m13/m12) return 1-fw;
        }

        return -1.1;
      }
      break;

    case TopTag::Likelihood:   
      {  
        float mW = -999;
        float mBW1 = -999;
        float mBW2 = -999;
        float mTop = -999;

        mBW1 = (subjets[0].p4()+subjets[1].p4()).M();
        mBW2 = (subjets[0].p4()+subjets[2].p4()).M();
        mW = (subjets[1].p4()+subjets[2].p4()).M();
        mTop = (subjets[0].p4()+subjets[1].p4()+subjets[2].p4()).M();

        float mratio=mW/mTop;
        float atan=TMath::ATan(mBW1/mBW2);

        float ptt = mtop_top_histo->Interpolate(mTop);
        float ptf = mtop_nottop_histo->Interpolate(mTop);
        float pat = atan_top_histo->Interpolate(atan);
        float paf = atan_nottop_histo->Interpolate(atan);
        float pmt = mratio_top_histo->Interpolate(mratio);
        float pmf = mratio_nottop_histo->Interpolate(mratio);
        float lr = ptt*pat*pmt/(ptt*pat*pmt+ptf*paf*pmf);

        return 2*lr-1;
      }
      break;

    case TopTag::TMVA:
      {
        ResetTMVAVars();

        for(std::vector<string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
          int iVar = itVarName-TMVAVarNames.begin();
        
          if(*itVarName=="TopJet_Top_M")                TMVAVars[iVar] = topjet.topjet.mass();                                                                                                         
          else if(*itVarName=="TopJet_PrunedMass")      TMVAVars[iVar] = topjet.prunedMass;                                                                                                            
          else if(*itVarName=="TopJet_UnfilteredMass")  TMVAVars[iVar] = topjet.unfilteredMass;                                                                                                        
          else if(*itVarName=="TopJet_fRec")            TMVAVars[iVar] = topjet.fRec;                                                                                                                  
          else if(*itVarName=="TopJet_DRoptRoptCalc")   TMVAVars[iVar] = topjet.Ropt-topjet.RoptCalc;                                                                                                  
          else if(*itVarName=="TopJet_Tau21Filtered")   TMVAVars[iVar] = topjet.tau2Filtered/topjet.tau1Filtered;                                                                                      
          else if(*itVarName=="TopJet_Tau32Filtered")   TMVAVars[iVar] = topjet.tau3Filtered/topjet.tau2Filtered;                                                                                      
          else if(*itVarName=="TopJet_WM" || *itVarName=="TopJet_Wbtag_M")                    TMVAVars[iVar] = (subjets[1].p4()+subjets[2].p4()).M();                                              
          else if(*itVarName=="TopJet_BW1M" || *itVarName=="TopJet_BW1btag_M")                TMVAVars[iVar] = (subjets[0].p4()+subjets[1].p4()).M();                                              
          else if(*itVarName=="TopJet_BW2M" || *itVarName=="TopJet_BW2btag_M")                TMVAVars[iVar] = (subjets[0].p4()+subjets[2].p4()).M();                                              
          else if(*itVarName=="TopJet_BCSV" || *itVarName=="TopJet_Bbtag_CSV")                TMVAVars[iVar] = fmax(subjets[0].bDiscriminator(btagger),-.1);                                       
          else if(*itVarName=="TopJet_W1CSV" || *itVarName=="TopJet_W1btag_CSV")              TMVAVars[iVar] = fmax(subjets[1].bDiscriminator(btagger),-.1);                                       
          else if(*itVarName=="TopJet_W2CSV" || *itVarName=="TopJet_W2btag_CSV")              TMVAVars[iVar] = fmax(subjets[2].bDiscriminator(btagger),-.1);                                       
          else if(*itVarName=="TopJet_MRatio_WTop" || *itVarName=="TopJet_MRatio_Wbtag_Top")  TMVAVars[iVar] = (subjets[1].p4()+subjets[2].p4()).M()/topjet.topjet.mass();                         
          else if(*itVarName=="TopJet_Atan_BW1W2btag" || *itVarName=="TopJet_Atan_BW1W2btag") TMVAVars[iVar] = atan((subjets[0].p4()+subjets[1].p4()).M()/(subjets[0].p4()+subjets[2].p4()).M());
          else std::cout << "Error! No matching Top Tagger Input Variable found!" << std:: endl;
        }

        if(verbose){
          std::cout << "Top Tagger Variables:" << std::endl;
          for(std::vector<string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
            std::cout << "Variable " << itVarName-TMVAVarNames.begin() << ": " << *itVarName << ": " << TMVAVars[itVarName-TMVAVarNames.begin()] << std::endl;
          }
        }

        float TMVAoutput = TMVAReader->EvaluateMVA("TMVATopTagger");

        if(verbose) std::cout << "TMVAOutput: " << TMVAoutput << std::endl;

        return TMVAoutput;
      }  
      break;

    default:
    
      std::cout << "Error! No matching Top Tagger Mode found!" << std:: endl;
      return -999;
  }
}


boosted::HTTTopJetCollection TopTagger::GetSortedByTopTaggerOutput(const boosted::HTTTopJetCollection& topJets, bool verbose){
  
  boosted::HTTTopJetCollection result = topJets;
  
  TopTaggerOutputComparison topTagComp(this);
  std::sort(result.begin(), result.end(),topTagComp);
  
  return result;
}


float TopTagger::GetTopHad(boosted::HTTTopJetCollection& topJets, boosted::HTTTopJet& topHadCand, bool verbose){
  
  float maxTopTag=-9999;
  
  for(boosted::HTTTopJetCollection::iterator itTopJet = topJets.begin() ; itTopJet != topJets.end(); ++itTopJet){
    
    float topTag = GetTopTaggerOutput(*itTopJet,verbose);
    
    if(topTag > maxTopTag) {
      maxTopTag = topTag;
      topHadCand = *itTopJet;
    }
  } 
  
  return maxTopTag;
}
