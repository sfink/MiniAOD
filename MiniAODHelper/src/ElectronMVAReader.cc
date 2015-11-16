#include "MiniAOD/MiniAODHelper/interface/ElectronMVAReader.h"

ElectronMVAReader::ElectronMVAReader(std::string WeightFilePath_){
  
  char* CMSSWPath = getenv("CMSSW_BASE");
  std::string filePath = CMSSWPath;
  weightsPath = filePath+"/src/"+WeightFilePath_;
  
  // Initialize TMVA input variables
  GetTMVAVars(weightsPath);
  // Setup TMVA Reader
  TMVAReader = new TMVA::Reader("Silent");    
        for(std::vector<std::string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
          TMVAReader->AddVariable(*itVarName,&TMVAVars[itVarName-TMVAVarNames.begin()]);
        }
  TMVAReader->BookMVA("ElectronMVA",(weightsPath).c_str());
}


ElectronMVAReader::~ElectronMVAReader(){
}


void ElectronMVAReader::GetTMVAVarNames(std::string filePath_, bool verbose){
  
  TMVAVarNames.clear();
  
  std::string line;
  std::ifstream inFile (filePath_);
  
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
      for(std::vector<std::string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
        std::cout << "Variable " << itVarName-TMVAVarNames.begin() << ": " << *itVarName << std::endl;
      }
    }
    
    inFile.close();
  }
}


void ElectronMVAReader::GetTMVAVars(std::string filePath_, bool verbose){
  
  GetTMVAVarNames(filePath_,verbose);
  
  TMVAVars = new float[50];
  
  for(int i=0;i<50;i++){
    TMVAVars[i] = -999;
  }
  
  if(verbose){
    for(std::vector<std::string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
      std::cout << "Variable " << itVarName-TMVAVarNames.begin() << ": " << *itVarName << ": " << TMVAVars[itVarName-TMVAVarNames.begin()] << std::endl;
    }
  }  
}


void ElectronMVAReader::ResetTMVAVars(){
  for(int i=0;i<50;i++){
    TMVAVars[i] = -999;
  }
}


float ElectronMVAReader::GetElectronMVAReaderOutput(const pat::Electron& iElectron, const edm::Handle< reco::ConversionCollection >& conversions, const edm::Handle< reco::BeamSpot >& theBeamSpot, bool verbose){

        ResetTMVAVars();


        bool validKF= false; 
        reco::TrackRef myTrackRef = iElectron.closestCtfTrackRef();
        validKF = (myTrackRef.isAvailable() && (myTrackRef.isNonnull()) ); 

	reco::ConversionRef conv_ref = ::ConversionTools::matchedConversion(iElectron, conversions, theBeamSpot.product()->position());
	double vertexFitProbability = -1.; 
	if(!conv_ref.isNull()) {
	  const reco::Vertex &vtx = conv_ref.get()->conversionVertex(); if (vtx.isValid()) {
	    vertexFitProbability = TMath::Prob( vtx.chi2(), vtx.ndof());
	  } 
	}


        for(std::vector<std::string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
          int iVar = itVarName-TMVAVarNames.begin();

          if(*itVarName=="ele_oldsigmaietaieta")                  	TMVAVars[iVar] = iElectron.full5x5_sigmaIetaIeta();
          else if(*itVarName=="ele_oldsigmaiphiiphi")                   TMVAVars[iVar] = iElectron.full5x5_sigmaIphiIphi();
          else if(*itVarName=="ele_oldcircularity")                     TMVAVars[iVar] = std::min(std::max(1. - iElectron.full5x5_e1x5() / iElectron.full5x5_e5x5(), -1.0), 2.0);                                                                    
          else if(*itVarName=="ele_oldr9")                      	TMVAVars[iVar] = std::min(iElectron.full5x5_r9(), float(5.0));

          else if(*itVarName=="ele_scletawidth")                      	TMVAVars[iVar] = iElectron.superCluster()->etaWidth();
          else if(*itVarName=="ele_sclphiwidth")                      	TMVAVars[iVar] = iElectron.superCluster()->phiWidth();
          else if(*itVarName=="ele_he")                      		TMVAVars[iVar] = iElectron.hadronicOverEm();
          else if(*itVarName=="ele_psEoverEraw")                      	TMVAVars[iVar] = iElectron.superCluster()->preshowerEnergy() / iElectron.superCluster()->rawEnergy();


          else if(*itVarName=="ele_kfhits")                     	TMVAVars[iVar] =(validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ;
          else if(*itVarName=="ele_kfchi2")                      	TMVAVars[iVar] = (validKF) ? std::min(myTrackRef->normalizedChi2(), 10.0) : 0;
          else if(*itVarName=="ele_gsfchi2")                      	TMVAVars[iVar] = std::min(iElectron.gsfTrack()->normalizedChi2(), 200.0);

          else if(*itVarName=="ele_fbrem")                      	TMVAVars[iVar] = std::max(iElectron.fbrem(), float(-1.0));
          else if(*itVarName=="ele_gsfhits")                      	TMVAVars[iVar] = iElectron.gsfTrack()->hitPattern().trackerLayersWithMeasurement();
          else if(*itVarName=="ele_expected_inner_hits")                TMVAVars[iVar] = iElectron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

	  
          
          else if(*itVarName=="ele_conversionVertexFitProbability")     TMVAVars[iVar] = vertexFitProbability;
          
          else if(*itVarName=="ele_ep")                      		TMVAVars[iVar] = std::min(iElectron.eSuperClusterOverP(), float(20.0));
          else if(*itVarName=="ele_eelepout")                      	TMVAVars[iVar] = std::min(iElectron.eEleClusterOverPout(), float(20.0));
          else if(*itVarName=="ele_IoEmIop")                      	TMVAVars[iVar] = (1.0/iElectron.ecalEnergy()) - (1.0 / float(iElectron.trackMomentumAtVtx().R()) );
          else if(*itVarName=="ele_deltaetain")                      	TMVAVars[iVar] = std::min(fabs(iElectron.deltaEtaSuperClusterTrackAtVtx()), 0.06);
          else if(*itVarName=="ele_deltaphiin")                      	TMVAVars[iVar] = std::min(fabs(iElectron.deltaPhiSuperClusterTrackAtVtx()), 0.6);
          else if(*itVarName=="ele_deltaetaseed")                      	TMVAVars[iVar] = std::min(fabs(iElectron.deltaEtaSeedClusterTrackAtCalo()), 0.2);
                                                                                                                                         
        }

        if(verbose){
          std::cout << "Electron MVA Variables:" << std::endl;
          for(std::vector<std::string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
            std::cout << "Variable " << itVarName-TMVAVarNames.begin() << ": " << *itVarName << ": " << TMVAVars[itVarName-TMVAVarNames.begin()] << std::endl;
          }
        }

        float TMVAoutput = TMVAReader->EvaluateMVA("ElectronMVA");

        if(verbose) std::cout << "TMVAOutput: " << TMVAoutput << std::endl;

        return TMVAoutput;
}    

