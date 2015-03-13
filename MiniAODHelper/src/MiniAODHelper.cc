#include "../interface/MiniAODHelper.h"

using namespace std;

// Constructor
MiniAODHelper::MiniAODHelper(){

  isSetUp = false;

  vertexIsSet = false;
  rhoIsSet = false;
  jetcorrectorIsSet = false;
  factorizedjetcorrectorIsSet = false;

  // twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging#Preliminary_working_or_operating
  // Preliminary working (or operating) points for CSVv2+IVF
  CSVLwp = 0.423; // 10.1716% DUSG mistag efficiency
  CSVMwp = 0.814; // 1.0623% DUSG mistag efficiency
  CSVTwp = 0.941; // 0.1144% DUSG mistag efficiency

  samplename = "blank";

}

// Destructor
MiniAODHelper::~MiniAODHelper(){

}

// Set up parameters one by one
void MiniAODHelper::SetUp(string iEra, int iSampleNumber, const analysisType::analysisType iAnalysis, bool iIsData){
  // Make sure we don't set up more than once
  if(isSetUp){ ThrowFatalError("Trying to set up 'BEANhelper' for the second time. Check your code."); }

  // Bring in the external values
  era          = iEra;
  sampleNumber = iSampleNumber;
  analysis     = iAnalysis;
  isData       = iIsData;

  // Error checking here
  if((era != "2011") && (era != "2012_52x") && (era != "2012_53x") && (era != "2015_72x") && (era != "2015_73x") && (era != "2015_74x")){ ThrowFatalError("era set to '" + era + "' but it has to be either 2011, 2012_52x, 2012_53x, 2015_72x, 2015_73x, or 2015_74x"); }
  if(sampleNumber==0){ ThrowFatalError("'sampleNumber' cannot be '0'."); }

  // Setup PU reweighing
  //SetUpPUreweighing(iCollisionDS);

  // Setup CSV reweighting
  //SetUpCSVreweighting();

  // Setup jet efficiency scale factors
  //SetUpJetSF();

  // Setup lepton efficiency scale factors
  //SetUpLeptonSF();

  // Awknowledge setup
  isSetUp = true;

}

// Set up parameters one by one
void MiniAODHelper::SetVertex(const reco::Vertex& inputVertex){

  vertex = inputVertex;

  vertexIsSet = true;
}

// Set up parameters one by one
void MiniAODHelper::SetRho(double inputRho){

  useRho = inputRho;

  rhoIsSet = true;
}

// Set up parameters one by one
void MiniAODHelper::SetJetCorrector(const JetCorrector* iCorrector){

  corrector = iCorrector;

  jetcorrectorIsSet = true;
}

// Set up parameters one by one
void MiniAODHelper::SetJetCorrectorUncertainty(){

  std::string inputJECfile = ( isData ) ? string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/Summer13_V5_DATA_Uncertainty_AK5PFchs.txt" : string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/Summer13_V5_MC_Uncertainty_AK5PFchs.txt";

  jecUnc_ = new JetCorrectionUncertainty(inputJECfile);

}

// Set up parameters one by one
void MiniAODHelper::SetFactorizedJetCorrector(){

  // Create the JetCorrectorParameter objects, the order does not matter.
  //JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/"); 
  JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/PLS170_V7AN1_L3Absolute_AK4PFchs.txt");
  JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/PLS170_V7AN1_L2Relative_AK4PFchs.txt");
  JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/PLS170_V7AN1_L1FastJet_AK4PFchs.txt");
  //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
  std::vector<JetCorrectorParameters> vPar;
  vPar.push_back(*L1JetPar);
  vPar.push_back(*L2JetPar);
  vPar.push_back(*L3JetPar);
  //vPar->push_back(ResJetPar);

  useJetCorrector = new FactorizedJetCorrector(vPar);

  std::string inputJECfile = ( isData ) ? string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/Summer13_V5_DATA_Uncertainty_AK5PFchs.txt" : string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/Summer13_V5_MC_Uncertainty_AK5PFchs.txt";

  jecUnc_ = new JetCorrectionUncertainty(inputJECfile);

  factorizedjetcorrectorIsSet = true;
}



std::vector<pat::Muon> 
MiniAODHelper::GetSelectedMuons(const std::vector<pat::Muon>& inputMuons, const float iMinPt, const muonID::muonID iMuonID, const coneSize::coneSize iconeSize, const corrType::corrType icorrType){

  CheckSetUp();

  std::vector<pat::Muon> selectedMuons;

  for( std::vector<pat::Muon>::const_iterator it = inputMuons.begin(), ed = inputMuons.end(); it != ed; ++it ){
    if( isGoodMuon(*it,iMinPt,iMuonID,iconeSize,icorrType) ) selectedMuons.push_back(*it);
  }

  return selectedMuons;
}


std::vector<pat::Electron> 
MiniAODHelper::GetSelectedElectrons(const std::vector<pat::Electron>& inputElectrons, const float iMinPt, const electronID::electronID iElectronID){

  CheckSetUp();

  std::vector<pat::Electron> selectedElectrons;

  for( std::vector<pat::Electron>::const_iterator it = inputElectrons.begin(), ed = inputElectrons.end(); it != ed; ++it ){
    if( isGoodElectron(*it,iMinPt,iElectronID) ) selectedElectrons.push_back(*it);
  }

  return selectedElectrons;
}

std::vector<pat::Tau> 
MiniAODHelper::GetSelectedTaus(const std::vector<pat::Tau>& inputTaus, const float iMinPt, const tau::ID id){

  CheckSetUp();

  std::vector<pat::Tau> selectedTaus;

  for( std::vector<pat::Tau>::const_iterator it = inputTaus.begin(), ed = inputTaus.end(); it != ed; ++it ){
    if( isGoodTau(*it,iMinPt,id) ) selectedTaus.push_back(*it);
  }

  return selectedTaus;
}



std::vector<pat::Jet> 
MiniAODHelper::GetSelectedJets(const std::vector<pat::Jet>& inputJets, const float iMinPt, const float iMaxAbsEta, const jetID::jetID iJetID, const char iCSVwp){

  CheckSetUp();

  std::vector<pat::Jet> selectedJets;

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    if( isGoodJet(*it, iMinPt, iMaxAbsEta, iJetID, iCSVwp) ) selectedJets.push_back(*it);
  }

  return selectedJets;
}


std::vector<pat::Jet> 
MiniAODHelper::GetUncorrectedJets(const std::vector<pat::Jet>& inputJets){

  CheckSetUp();

  std::vector<pat::Jet> outputJets;

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    pat::Jet jet = (*it);
    jet.setP4( it->correctedJet(0).p4() );
    outputJets.push_back(jet);
  }

  return outputJets;
}


std::vector<pat::Jet> 
MiniAODHelper::GetCorrectedJets(const std::vector<pat::Jet>& inputJets, const edm::Event& event, const edm::EventSetup& setup, const sysType::sysType iSysType){

  CheckSetUp();

  std::vector<pat::Jet> outputJets;

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    pat::Jet jet = (*it);
    double scale = 1.;

    if( jetcorrectorIsSet ) scale = corrector->correction(*it, event, setup);
    else std::cout << " !! ERROR !! Trying to use Full Framework GetCorrectedJets without setting jet corrector !" << std::endl;

    jet.scaleEnergy( scale );

    if( iSysType == sysType::JESup || iSysType == sysType::JESdown ){
      jecUnc_->setJetEta(jet.eta());
      jecUnc_->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
      double unc = 1;
      double jes = 1;
      if( iSysType==sysType::JESup ){
	unc = jecUnc_->getUncertainty(true);
	jes = 1 + unc;
      }
      else if( iSysType==sysType::JESdown ){
	unc = jecUnc_->getUncertainty(false);
	jes = 1 - unc;
      }

      jet.scaleEnergy( jes );
    }

    outputJets.push_back(jet);
  }

  return outputJets;
}


std::vector<pat::Jet> 
MiniAODHelper::GetCorrectedJets(const std::vector<pat::Jet>& inputJets, const sysType::sysType iSysType ){

  CheckSetUp();

  std::vector<pat::Jet> outputJets;

  if( !(factorizedjetcorrectorIsSet && rhoIsSet) ){
   std::cout << " !! ERROR !! Trying to use FWLite Framework GetCorrectedJets without setting factorized jet corrector !" << std::endl;

     return inputJets;
  }

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    pat::Jet jet = (*it);
    double scale = 1.;

    useJetCorrector->setJetEta(it->eta());
    useJetCorrector->setJetPt(it->pt());
    useJetCorrector->setJetA(it->jetArea());
    useJetCorrector->setRho(useRho); 

    scale = useJetCorrector->getCorrection();

    jet.scaleEnergy( scale );

    if( iSysType == sysType::JESup || iSysType == sysType::JESdown ){
      jecUnc_->setJetEta(jet.eta());
      jecUnc_->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
      double unc = 1;
      double jes = 1;
      if( iSysType==sysType::JESup ){
	unc = jecUnc_->getUncertainty(true);
	jes = 1 + unc;
      }
      else if( iSysType==sysType::JESdown ){
	unc = jecUnc_->getUncertainty(false);
	jes = 1 - unc;
      }

      jet.scaleEnergy( jes );
    }

    outputJets.push_back(jet);
  }

  return outputJets;
}

bool 
MiniAODHelper::isGoodMuon(const pat::Muon& iMuon, const float iMinPt, const muonID::muonID iMuonID, const coneSize::coneSize iconeSize, const corrType::corrType icorrType){

  CheckVertexSetUp();

  double minMuonPt = iMinPt;

  float maxLooseMuonAbsEta = 2.5;

  float maxTightMuonAbsEta = 2.4;


  // Be skeptical about this muon making it through
  bool passesKinematics	= false;
  bool passesIso        = false;
  bool passesID         = false;
  bool isPFMuon         = false;
  bool passesTrackerID  = false;

  bool passesGlobalTrackID   = false;
  bool passesMuonBestTrackID = false;
  bool passesInnerTrackID    = false;
  bool passesTrackID         = false;


  switch(iMuonID){
  case muonID::muonPreselection:
  case muonID::muonSide:
  case muonID::muonSideLooseMVA:
  case muonID::muonSideTightMVA:
  case muonID::muonPtOnly:
  case muonID::muonPtEtaOnly:
  case muonID::muonPtEtaIsoOnly:
  case muonID::muonPtEtaIsoTrackerOnly:
  case muonID::muonRaw:
  case muonID::muonLooseCutBased:
  case muonID::muonTightCutBased:
  case muonID::muonCutBased:
  case muonID::muonLoose:
    passesKinematics = ((iMuon.pt() >= minMuonPt) && (fabs(iMuon.eta()) <= maxLooseMuonAbsEta));
    passesIso        = (GetMuonRelIso(iMuon,iconeSize,icorrType) < 0.200);
    isPFMuon         = iMuon.isPFMuon();
    passesID         = (( iMuon.isGlobalMuon() || iMuon.isTrackerMuon() ) && isPFMuon);
    break;
  case muonID::muonTight:
    passesKinematics = ((iMuon.pt() >= minMuonPt) && (fabs(iMuon.eta()) <= maxTightMuonAbsEta));
    passesIso        = (GetMuonRelIso(iMuon,iconeSize,icorrType) < 0.12);
    isPFMuon         = iMuon.isPFMuon();

    if( iMuon.globalTrack().isAvailable() ){
      passesGlobalTrackID = ( (iMuon.globalTrack()->normalizedChi2() < 10.) 
			      && (iMuon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0)
			      );
    }
    if( iMuon.muonBestTrack().isAvailable() ){
      passesMuonBestTrackID = ( (fabs(iMuon.muonBestTrack()->dxy(vertex.position())) < 0.2)
				&& (fabs(iMuon.muonBestTrack()->dz(vertex.position())) < 0.5)
				);
    }
    if( iMuon.innerTrack().isAvailable() )
      passesInnerTrackID = (iMuon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0);
    if( iMuon.track().isAvailable() )
      passesTrackID = (iMuon.track()->hitPattern().trackerLayersWithMeasurement() > 5);

    passesTrackerID = ( passesGlobalTrackID && passesMuonBestTrackID && passesInnerTrackID && passesTrackID && (iMuon.numberOfMatchedStations() > 1) );

    passesID        = (iMuon.isGlobalMuon() && isPFMuon && passesTrackerID);
    break;
  }

  return (passesKinematics && passesIso && passesID);
}



bool 
MiniAODHelper::isGoodElectron(const pat::Electron& iElectron, const float iMinPt, const electronID::electronID iElectronID){

  CheckVertexSetUp();

  double minElectronPt = iMinPt;

  float maxLooseElectronAbsEta = 2.4;
  float maxTightElectronAbsEta = 2.4;


  // Be skeptical about this electron making it through
  bool passesKinematics	= false;
  bool passesIso        = false;
  bool passesID         = false;

  double SCeta = (iElectron.superCluster().isAvailable()) ? iElectron.superCluster()->position().eta() : -99;
  double absSCeta = fabs(SCeta);

  bool inCrack = false;
  if( iElectron.superCluster().isAvailable() ) inCrack = ( absSCeta>1.4442 && absSCeta<1.5660 );


  bool myTrigPresel = true;

  //double eleID      = iElectron.electronID("eidTight");
  bool passMVAId53x = true;//( eleID>0.5 );  // For 2012_53x, tighter selection

  bool d02 = false; 
  bool d04 = false;
  bool dZ  = false;
  bool no_exp_inner_trkr_hits = true; //false; // see below
  if( iElectron.gsfTrack().isAvailable() ){
    d02 = ( fabs(iElectron.gsfTrack()->dxy(vertex.position())) < 0.02 );
    d04 = ( fabs(iElectron.gsfTrack()->dxy(vertex.position())) < 0.04 );
    //no_exp_inner_trkr_hits = ( iElectron.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 ); // deprecated in 7_2_0 .. replace with ..?
    dZ = ( fabs(iElectron.gsfTrack()->dz(vertex.position())) < 1. );
  }


  bool notConv = ( iElectron.passConversionVeto() );
  bool id      = ( passMVAId53x && d02 && dZ && notConv );


  switch(iElectronID){
  case electronID::electronPreselection:
  case electronID::electronSide:
  case electronID::electronSideLooseMVA:
  case electronID::electronSideTightMVA:
  case electronID::electronLooseMinusTrigPresel:
  case electronID::electronRaw:
  case electronID::electronLooseCutBased:
  case electronID::electronTightCutBased:
  case electronID::electronCutBased:
  case electronID::electronLoose:
    passesKinematics = ((iElectron.pt() >= minElectronPt) && (fabs(iElectron.eta()) <= maxLooseElectronAbsEta) && !inCrack);
    passesIso        = (GetElectronRelIso(iElectron) < 0.200);
    passesID         = ( passMVAId53x && no_exp_inner_trkr_hits && d04 && notConv && myTrigPresel );
    break;
  case electronID::electronTightMinusTrigPresel:
  case electronID::electronTight:
    passesKinematics = ((iElectron.pt() >= minElectronPt) && (fabs(iElectron.eta()) <= maxTightElectronAbsEta) && !inCrack);
    passesIso        = (GetElectronRelIso(iElectron) < 0.100);
    passesID         = ( id && no_exp_inner_trkr_hits && myTrigPresel );
    break;
  case electronID::electronPhys14L:
  case electronID::electronPhys14M:
  case electronID::electronPhys14T:
    id = PassElectronPhys14Id( iElectron, iElectronID );
    passesIso = id;
    passesID = id;
    passesKinematics = ((iElectron.pt() >= minElectronPt) && (fabs(iElectron.eta()) <= maxTightElectronAbsEta) && !inCrack);
    break;
  }

  

  return (passesKinematics && passesIso && passesID);
}

bool
MiniAODHelper::isGoodTau(const pat::Tau& tau, const float min_pt, const tau::ID id)
{
  CheckVertexSetUp();

  double minTauPt = min_pt;

  bool passesKinematics = false;
  bool passesIsolation = false;
  bool passesID = tau.tauID("decayModeFinding") >= .5;

  switch (id) {
     case tau::nonIso:
        passesID = passesID and \
                   tau.tauID("againstMuonLoose3") >= .5 and \
                   tau.tauID("againstElectronVLooseMVA5") >= .5;
        passesIsolation = true;
        break;
     case tau::loose:
        passesID = passesID and \
                   tau.tauID("againstMuonLoose3") >= .5 and \
                   tau.tauID("againstElectronVLooseMVA5") >= .5;
        passesIsolation = tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") >= .5;
        break;
     case tau::medium:
        passesID = passesID and \
                   tau.tauID("againstMuonLoose3") >= .5 and \
                   tau.tauID("againstElectronLooseMVA5") >= .5;
        passesIsolation = tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") >= .5;
        break;
     case tau::tight:
        passesID = passesID and \
                   tau.tauID("againstMuonTight3") >= .5 and \
                   tau.tauID("againstElectronMediumMVA5") >= .5;
        passesIsolation = tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") >= .5;
        break;
  }

  // systematics are only defined for p_T > 20
  passesKinematics = (tau.pt() >= 20) && (fabs(tau.eta()) <= 2.1) && (tau.pt() > minTauPt);

  return passesKinematics && passesIsolation && passesID;
}

bool 
MiniAODHelper::isGoodJet(const pat::Jet& iJet, const float iMinPt, const float iMaxAbsEta, const jetID::jetID iJetID, const char iCSVworkingPoint){

  CheckVertexSetUp();

  // Transverse momentum requirement
  if( iJet.pt() < iMinPt ) return false;

  // Absolute eta requirement
  if( fabs(iJet.eta()) > iMaxAbsEta ) return false;

  bool loose = (
		iJet.neutralHadronEnergyFraction() < 0.99 &&
		iJet.chargedEmEnergyFraction() < 0.99 &&
		iJet.neutralEmEnergyFraction() < 0.99 &&
		iJet.numberOfDaughters() > 1
		);

  if( fabs(iJet.eta())<2.4 ){
    loose = ( loose &&
	      iJet.chargedHadronEnergyFraction() > 0.0 &&
	      iJet.chargedMultiplicity() > 0
	      );
  }

  // Jet ID
  switch(iJetID){
  case jetID::none:
  case jetID::jetPU:
  case jetID::jetMinimal:
  case jetID::jetLooseAOD:
  case jetID::jetLoose:
  case jetID::jetTight:
    if( !loose ) return false;
    break;
  default:
    break;
  }

  if( !PassesCSV(iJet, iCSVworkingPoint) ) return false;

  return true;
}

float MiniAODHelper::GetMuonRelIso(const pat::Muon& iMuon,const coneSize::coneSize iconeSize, const corrType::corrType icorrType) const
{
  //rho corrections based on phys14
  //details here: https://www.dropbox.com/s/66lzhbro09diksa/effectiveareas-pog-121214.pdf?dl=0
  // !!! NOTE !!! rho used should be: fixedGridRhoFastjetAll
  float result = 9999; 
  
  double correction = 9999.;
  double EffArea = 9999.;
  double Eta = abs(iMuon.eta());
  
  double pfIsoCharged;
  double pfIsoNeutral;
  double pfIsoPUSubtracted;
  
  switch(iconeSize)
    {
    case coneSize::R04:
      if (Eta >= 0. && Eta < 0.8) EffArea = 0.1546;
      else if (Eta >= 0.8 && Eta < 1.3) EffArea = 0.1325;
      else if (Eta >= 1.3 && Eta < 2.0) EffArea = 0.0913;
      else if (Eta >= 2.0 && Eta < 2.2) EffArea = 0.1212;
      else if (Eta >= 2.2 && Eta <= 2.5) EffArea = 0.2085;
      
      pfIsoCharged = iMuon.pfIsolationR04().sumChargedHadronPt;
      pfIsoNeutral = iMuon.pfIsolationR04().sumNeutralHadronEt + iMuon.pfIsolationR04().sumPhotonEt;
      
      switch(icorrType){
      case corrType::rhoEA:
         correction = useRho*EffArea;
         break;
      case corrType::deltaBeta:
         correction =  0.5*iMuon.pfIsolationR04().sumPUPt;
         break;
      }
      
      pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction );
      result = (pfIsoCharged + pfIsoPUSubtracted)/iMuon.pt();
      break;
    case coneSize::R03:
      if (Eta >= 0. && Eta < 0.8) EffArea = 0.0913;
      else if (Eta >= 0.8 && Eta < 1.3) EffArea = 0.0765;
      else if (Eta >= 1.3 && Eta < 2.0) EffArea = 0.0546;
      else if (Eta >= 2.0 && Eta < 2.2) EffArea = 0.0728;
      else if (Eta >= 2.2 && Eta <= 2.5) EffArea = 0.1177;
      
      pfIsoCharged = iMuon.pfIsolationR03().sumChargedHadronPt;
      pfIsoNeutral = iMuon.pfIsolationR03().sumNeutralHadronEt + iMuon.pfIsolationR03().sumPhotonEt;
      
      switch(icorrType){
      case corrType::rhoEA:
         correction = useRho*EffArea;
         break;
      case corrType::deltaBeta:
         correction = 0.5*iMuon.pfIsolationR03().sumPUPt;
         break;
      }
      
      pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction );
      result = (pfIsoCharged + pfIsoPUSubtracted)/iMuon.pt();
      break;
    }
  
  return result;
}

float MiniAODHelper::GetElectronRelIso(const pat::Electron& iElectron,const coneSize::coneSize iconeSize, const corrType::corrType icorrType) const
{
  //rho corrections based on phys14
  //details here: https://www.dropbox.com/s/66lzhbro09diksa/effectiveareas-pog-121214.pdf?dl=0
  // !!! NOTE !!! rho used should be: fixedGridRhoFastjetAll
  float result = 9999; 
  
  double correction = 9999.;
  double EffArea = 9999.;
  double Eta = abs(iElectron.eta());
  
  double pfIsoCharged;
  double pfIsoNeutral;
  double pfIsoPUSubtracted;
  
  switch(iconeSize)
    {
    case coneSize::R04:
    case coneSize::R03:
      if (Eta >= 0. && Eta < 0.8) EffArea = 0.1013;
      else if (Eta >= 0.8 && Eta < 1.3) EffArea = 0.0988;
      else if (Eta >= 1.3 && Eta < 2.0) EffArea = 0.0572;
      else if (Eta >= 2.0 && Eta < 2.2) EffArea = 0.0842;
      else if (Eta >= 2.2 && Eta <= 2.5) EffArea = 0.1530;
      
      pfIsoCharged = iElectron.pfIsolationVariables().sumChargedHadronPt;
      pfIsoNeutral = iElectron.pfIsolationVariables().sumNeutralHadronEt + iElectron.pfIsolationVariables().sumPhotonEt;
      
      switch(icorrType){
      case corrType::rhoEA:
         correction = useRho*EffArea;
         break;
      case corrType::deltaBeta:
         correction = 0.5*iElectron.pfIsolationVariables().sumPUPt;
         break;
      }
      
      pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction );
      result = (pfIsoCharged + pfIsoPUSubtracted)/iElectron.pt();
      break;
    }
  
  return result;
}

bool MiniAODHelper::PassesCSV(const pat::Jet& iJet, const char iCSVworkingPoint){
  CheckSetUp();

  float csvValue = iJet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");

  // CSV b-tagging requirement
  switch(iCSVworkingPoint){
  case 'L':	if(csvValue > CSVLwp){ return true; }	break;
  case 'M':	if(csvValue > CSVMwp){ return true; }	break;
  case 'T':	if(csvValue > CSVTwp){ return true; }	break;
  case '-':	return true;                            break;
  }
  return false;
}


bool MiniAODHelper::PassElectronPhys14Id(const pat::Electron& iElectron, const electronID::electronID iElectronID) const{

  double SCeta = (iElectron.superCluster().isAvailable()) ? iElectron.superCluster()->position().eta() : -99;
  double absSCeta = fabs(SCeta);

  bool isEB = ( absSCeta < 1.479 );

  double pfIsoCharged = iElectron.pfIsolationVariables().sumChargedHadronPt;
  double pfIsoNeutralHadron = iElectron.pfIsolationVariables().sumNeutralHadronEt;
  double pfIsoNeutralPhoton = iElectron.pfIsolationVariables().sumPhotonEt;
  double pfIsoSumPUPt = iElectron.pfIsolationVariables().sumPUPt;

  double relIso = (pfIsoCharged + std::max( pfIsoNeutralHadron + pfIsoNeutralPhoton - 0.5*pfIsoSumPUPt, 0.0 ))/iElectron.pt();


  double full5x5_sigmaIetaIeta = iElectron.full5x5_sigmaIetaIeta();
  double dEtaIn = fabs( iElectron.deltaEtaSuperClusterTrackAtVtx() );
  double dPhiIn = fabs( iElectron.deltaPhiSuperClusterTrackAtVtx() );
  double hOverE = iElectron.hcalOverEcal();

  double ooEmooP = -999;
  if( iElectron.ecalEnergy() == 0 ) ooEmooP = 1e30;
  else if( !std::isfinite(iElectron.ecalEnergy()) ) ooEmooP = 1e30;
  else ooEmooP = fabs(1.0/iElectron.ecalEnergy() - iElectron.eSuperClusterOverP()/iElectron.ecalEnergy() );

  double d0 = -999;
  double dZ = -999;
  double expectedMissingInnerHits = -999;
  if( iElectron.gsfTrack().isAvailable() ){
    d0 = fabs(iElectron.gsfTrack()->dxy(vertex.position()));
    dZ = fabs(iElectron.gsfTrack()->dz(vertex.position()));
    expectedMissingInnerHits = iElectron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
  }

  bool passConversionVeto = ( iElectron.passConversionVeto() );

  bool pass = false;
  switch(iElectronID){
  case electronID::electronPhys14L:
    if( isEB ){
      pass = ( full5x5_sigmaIetaIeta < 0.010557 &&
	       dEtaIn < 0.012442 &&
	       dPhiIn < 0.072624 &&
	       hOverE < 0.121476 &&
	       ooEmooP < 0.221803 &&
	       d0 < 0.022664 &&
	       dZ < 0.173670 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.120026
	       );
    }
    else{
      pass = ( full5x5_sigmaIetaIeta < 0.032602 &&
	       dEtaIn < 0.010654 &&
	       dPhiIn < 0.145129 &&
	       hOverE < 0.131862 &&
	       ooEmooP < 0.142283 &&
	       d0 < 0.097358 &&
	       dZ < 0.198444 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.162914
	       );
    }
    break;
  case electronID::electronPhys14M:
    if( isEB ){
      pass = ( full5x5_sigmaIetaIeta < 0.010399 &&
	       dEtaIn < 0.007641 &&
	       dPhiIn < 0.032643 &&
	       hOverE < 0.060662 &&
	       ooEmooP < 0.153897 &&
	       d0 < 0.011811 &&
	       dZ < 0.070775 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.097213
	       );
    }
    else{
      pass = ( full5x5_sigmaIetaIeta < 0.029524 &&
	       dEtaIn < 0.009285 &&
	       dPhiIn < 0.042447 &&
	       hOverE < 0.104263 &&
	       ooEmooP < 0.137468 &&
	       d0 < 0.051682 &&
	       dZ < 0.180720 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.116708
	       );
    }
    break;
  case electronID::electronPhys14T:
    if( isEB ){
      pass = ( full5x5_sigmaIetaIeta < 0.010181 &&
	       dEtaIn < 0.006574 &&
	       dPhiIn < 0.022868 &&
	       hOverE < 0.037553 &&
	       ooEmooP < 0.131191 &&
	       d0 < 0.009924 &&
	       dZ < 0.015310 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.074355
	       );
    }
    else{
      pass = ( full5x5_sigmaIetaIeta < 0.028766 &&
	       dEtaIn < 0.005681 &&
	       dPhiIn < 0.032046 &&
	       hOverE < 0.081902 &&
	       ooEmooP < 0.106055 &&
	       d0 < 0.027261 &&
	       dZ < 0.147154 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.090185
	       );
    }
    break;
  default:
    break;
  }

  return pass;
}
