#include <iostream>

#include "MiniAOD/MiniAODHelper/interface/CSVHelper.h"

CSVHelper::CSVHelper(std::string hf, std::string lf, int nHFptBins_):nHFptBins(nHFptBins_)
{
    std::string inputFileHF = hf.size() > 0 ? hf : "data/csv_rwt_hf_IT_FlatSF.root";
    std::string inputFileLF = lf.size() > 0 ? lf :"data/csv_rwt_lf_IT_FlatSF.root";

    TFile *f_CSVwgt_HF = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/" + inputFileHF).c_str());
    TFile *f_CSVwgt_LF = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/" + inputFileLF).c_str());

    fillCSVHistos(f_CSVwgt_HF, f_CSVwgt_LF);
}

// fill the histograms (done once)
void
CSVHelper::fillCSVHistos(TFile *fileHF, TFile *fileLF)
{
    for (int iSys = 0; iSys < 9; iSys++) {
        for (int iPt = 0; iPt < 5; iPt++)
            h_csv_wgt_hf[iSys][iPt] = NULL;
        for (int iPt = 0; iPt < 3; iPt++) {
            for (int iEta = 0; iEta < 3; iEta++)
                h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
        }
    }
    for (int iSys = 0; iSys < 5; iSys++) {
        for (int iPt = 0; iPt < 5; iPt++)
            c_csv_wgt_hf[iSys][iPt] = NULL;
    }

    // CSV reweighting /// only care about the nominal ones
    for (int iSys = 0; iSys < 9; iSys++) {
        TString syst_csv_suffix_hf = "final";
        TString syst_csv_suffix_c = "final";
        TString syst_csv_suffix_lf = "final";

        switch (iSys) {
            case 0:
                // this is the nominal case
                break;
            case 1:
                // JESUp
                syst_csv_suffix_hf = "final_JESUp";
                syst_csv_suffix_lf = "final_JESUp";
                syst_csv_suffix_c = "final_cErr1Up";
                break;
            case 2:
                // JESDown
                syst_csv_suffix_hf = "final_JESDown";
                syst_csv_suffix_lf = "final_JESDown";
                syst_csv_suffix_c = "final_cErr1Down";
                break;
            case 3:
                // purity up
                syst_csv_suffix_hf = "final_LFUp";
                syst_csv_suffix_lf = "final_HFUp";
                syst_csv_suffix_c = "final_cErr2Up";
                break;
            case 4:
                // purity down
                syst_csv_suffix_hf = "final_LFDown";
                syst_csv_suffix_lf = "final_HFDown";
                syst_csv_suffix_c = "final_cErr2Down";
                break;
            case 5:
                // stats1 up
                syst_csv_suffix_hf = "final_Stats1Up";
                syst_csv_suffix_lf = "final_Stats1Up";
                break;
            case 6:
                // stats1 down
                syst_csv_suffix_hf = "final_Stats1Down";
                syst_csv_suffix_lf = "final_Stats1Down";
                break;
            case 7:
                // stats2 up
                syst_csv_suffix_hf = "final_Stats2Up";
                syst_csv_suffix_lf = "final_Stats2Up";
                break;
            case 8:
                // stats2 down
                syst_csv_suffix_hf = "final_Stats2Down";
                syst_csv_suffix_lf = "final_Stats2Down";
                break;
        }

        for (int iPt = 0; iPt < nHFptBins; iPt++)
            h_csv_wgt_hf[iSys][iPt] =
                (TH1D *)fileHF->Get(Form("csv_ratio_Pt%i_Eta0_%s", iPt, syst_csv_suffix_hf.Data()));

        if (iSys < 5) {
            for (int iPt = 0; iPt < nHFptBins; iPt++)
                c_csv_wgt_hf[iSys][iPt] =
                    (TH1D *)fileHF->Get(Form("c_csv_ratio_Pt%i_Eta0_%s", iPt, syst_csv_suffix_c.Data()));
        }

        for (int iPt = 0; iPt < 4; iPt++) {
            for (int iEta = 0; iEta < 3; iEta++)
                h_csv_wgt_lf[iSys][iPt][iEta] =
                    (TH1D *)fileLF->Get(Form("csv_ratio_Pt%i_Eta%i_%s", iPt, iEta, syst_csv_suffix_lf.Data()));
        }
    }

    return;
}

double
CSVHelper::getCSVWeight(std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs,
                       std::vector<int> jetFlavors, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF)
{
    int iSysHF = 0;
    switch (iSys) {
        case 7:
            iSysHF = 1;
            break; // JESUp
        case 8:
            iSysHF = 2;
            break; // JESDown
        case 9:
            iSysHF = 3;
            break; // LFUp
        case 10:
            iSysHF = 4;
            break; // LFDown
        case 13:
            iSysHF = 5;
            break; // Stats1Up
        case 14:
            iSysHF = 6;
            break; // Stats1Down
        case 15:
            iSysHF = 7;
            break; // Stats2Up
        case 16:
            iSysHF = 8;
            break; // Stats2Down
        default:
            iSysHF = 0;
            break; // NoSys
    }

    int iSysC = 0;
    switch (iSys) {
        case 21:
            iSysC = 1;
            break;
        case 22:
            iSysC = 2;
            break;
        case 23:
            iSysC = 3;
            break;
        case 24:
            iSysC = 4;
            break;
        default:
            iSysC = 0;
            break;
    }

    int iSysLF = 0;
    switch (iSys) {
        case 7:
            iSysLF = 1;
            break; // JESUp
        case 8:
            iSysLF = 2;
            break; // JESDown
        case 11:
            iSysLF = 3;
            break; // HFUp
        case 12:
            iSysLF = 4;
            break; // HFDown
        case 17:
            iSysLF = 5;
            break; // Stats1Up
        case 18:
            iSysLF = 6;
            break; // Stats1Down
        case 19:
            iSysLF = 7;
            break; // Stats2Up
        case 20:
            iSysLF = 8;
            break; // Stats2Down
        default:
            iSysLF = 0;
            break; // NoSys
    }

    double csvWgthf = 1.;
    double csvWgtC = 1.;
    double csvWgtlf = 1.;

    for (int iJet = 0; iJet < int(jetPts.size()); iJet++) {
        double csv = jetCSVs[iJet];
        double jetPt = jetPts[iJet];
        double jetAbsEta = fabs(jetEtas[iJet]);
        int flavor = jetFlavors[iJet];

        int iPt = -1;
        int iEta = -1;
        if (jetPt >= 19.99 && jetPt < 30)
            iPt = 0;
        else if (jetPt >= 30 && jetPt < 40)
            iPt = 1;
        else if (jetPt >= 40 && jetPt < 60)
            iPt = 2;
        else if (jetPt >= 60 && jetPt < 100)
            iPt = 3;
        else if (jetPt >= 100 && jetPt < 160)
            iPt = 4;
        else if (jetPt >= 160)
            iPt = 5;

        if (jetAbsEta >= 0 && jetAbsEta < 0.8)
            iEta = 0;
        else if (jetAbsEta >= 0.8 && jetAbsEta < 1.6)
            iEta = 1;
        else if (jetAbsEta >= 1.6 && jetAbsEta < 2.41)
            iEta = 2;

        if (iPt < 0 || iEta < 0){
	  //std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt    << ", jetAbsEta = " << jetAbsEta << std::endl;
	}
        else if (abs(flavor) == 5) {
	    // RESET iPt to maximum pt bin (only 5 bins for new SFs)
	    if(iPt>=nHFptBins){
		iPt=nHFptBins-1;
	    }
            int useCSVBin = (csv >= 0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
            double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
            if (iCSVWgtHF != 0)
                csvWgthf *= iCSVWgtHF;

        } else if (abs(flavor) == 4) {
	    // RESET iPt to maximum pt bin (only 5 bins for new SFs)
	    if(iPt>=nHFptBins){
		iPt=nHFptBins-1;
	    }
            int useCSVBin = (csv >= 0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
            double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
            if (iCSVWgtC != 0)
                csvWgtC *= iCSVWgtC;
        } else {
            if (iPt >= 3)
                iPt = 3; /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
            int useCSVBin = (csv >= 0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
            double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
            if (iCSVWgtLF != 0)
                csvWgtlf *= iCSVWgtLF;
        }
    }

    double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

    csvWgtHF = csvWgthf;
    csvWgtLF = csvWgtlf;
    csvWgtCF = csvWgtC;

    return csvWgtTotal;
}
