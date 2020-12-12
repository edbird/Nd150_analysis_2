#ifndef FIT_2E_MAKEHISTOGRAMS_H
#define FIT_2E_MAKEHISTOGRAMS_H


// mode_flag = 0: run function with cuts enabled (normal), write to (normal)
//                location (/rawdata/...)
// mode_flag = 1: run function with cuts disabled, write to location
//                ...
//
//                TODO / NOTE: currently using same NAMES and pointers
//                to hold histogram data, may wish to change behaviour
//                depending on value of the flag such that the NAMES
//                and pointers are always unique (surpress root warnings)
void makeHistograms(
    TString thePath,
    TString sampleName,
    std::ofstream& ofile_cutcount,
    const Int_t mode_flag = 0,
    const uint64_t sample_split_flags = 0)
{

    // enable / disable additional cuts
    // by default the basic cuts are enabled
    // this includes sector and phase cut
    const Int_t mode_flag_2 = 1;

    // swap between raw mode and processed mode
    // raw mode disables all cuts
    TString name_append;
    if(mode_flag == 0)
    {
        //name_append = "";
    }
    else if(mode_flag == 1)
    {
        name_append = "_raw";
    }


    //Declare the histograms.
    TH1D* hRun; // 0
    TH1D* hNElectrons;
    TH1D* hTotalE;
    TH1D* hInternalPullee;
    TH1D* hInternalProbee; 
    TH1D* hExternalPullee; // 5
    TH1D* hExternalProbee;
    TH1D* hCosee;
    TH1D* hCoseeWeighted;
    TH1D* hVertexDZ;
    TH1D* hVertexDR; // 10
    TH1D* hVertexDRPhi;
    TH1D* hNAPromptGgHits;

    TH1D* hEeMax;
    TH1D* hElectronLengthMax;
    TH1D* hVertexZMax; // 15
    TH1D* hVertexSectorMax;
    TH1D* hVertexRMax;
    TH1D* hElectronFirstGgMax;
    TH1D* hElectronLastGgMax;
    TH1D* hVertexMinDistPromptGgMax; // 20
    TH1D* hElectronLDCorrMax;
    TH2D* hVertexZSecMax; // usefull really only for data.

    TH1D* hEeMin;
    TH1D* hElectronLengthMin;
    TH1D* hVertexZMin; // 25
    TH1D* hVertexSectorMin;
    TH1D* hVertexRMin;
    TH1D* hElectronFirstGgMin;
    TH1D* hElectronLastGgMin;
    TH1D* hVertexMinDistPromptGgMin; // 30
    TH1D* hElectronLDCorrMin;
    TH2D* hVertexZSecMin; // usefull really only for data.

    // TODO: obsolete, delete these
    TH1D* hNLowEGammas;
    TH1D* hLowEGammaEnergy;
    TH1D* hSummedLowEGammaE; // 35
    TH1D *hLowEMinDistPromptGg;

    TH2D *hEeMaxVEeMin;

    //Histograms for keeping track of analysis cuts
    TH1D* hNAfterCuts;
    TH1D* hEffAfterCuts;

    // my new histograms
    //
    //
    // missing items?
    TH1D* hNAPromptGgHitsDist2VertexMin; // 40
    TH1D *hTrackSignMin;
    TH1D *hTrackSignMax;
    TH1D *hnGammaClusters;
    TH1D *hnInCluster;
    TH1D *hclusterHitEnergy; // 45
    TH1D *hclusterHitEnergyMin;
    TH1D *hclusterHitEnergyMax;
    TH1D *hnLowEnergyHits;
    TH1D *hclusterEnergy;

    TH1D *hFoilSide;
    TH1D *hTrueVertexR;
    TH2D *hFoilSideTrueVertexR;
    TH1D *hTrueVertexLayer;

    TH1D *hSingleEnergy;
    TH1D *hHighEnergy;
    TH1D *hLowEnergy;
    TH2D *hHighLowEnergy;

    TH1D *hTotalClusterEnergy;
    TH2D *hHotSpotLocation;

    TH1D *hEnergySum;
    TH1D *hEnergyDiff;

    // NOTE: this idea (multuplying histograms) will not work because
    // we lose the information of which true energy values match which
    // reconstructed energy values
    // processing must be done event by event
    //#if TRUTH_ENABLE
    //TH2D *hMCElectronEnergy;
    //TH2D *hMCTrueElectronEnergy;
    //#endif


    TString name_sample_split_additional = "";
    /*
    if(sample_split_flags == 0x01)
    {
        // IN
        std::cout << "sample_split_flags is 1" << std::endl;
        name_sample_split_additional = "_IN";
    }
    else if(sample_split_flags == 0x02)
    {
        // OUT
        std::cout << "sample_split_flags is 2" << std::endl;
        name_sample_split_additional = "_OUT";
    }
    else if(sample_split_flags == 0x04)
    {
        // L0 IN
        std::cout << "sample_split_flags is 4" << std::endl;
        name_sample_split_additional = "_INL0";
    }
    else if(sample_split_flags == 0x08)
    {
        // L0 OUT
        std::cout << "sample_split_flags is 8" << std::endl;
        name_sample_split_additional = "_OUTL0";
    }
    else if(sample_split_flags == 0x10)
    {
        // L>0 IN+OUT
        std::cout << "sample_split_flags is 16" << std::endl;
        name_sample_split_additional = "_Lg0";
    }
    std::cout << "name_sample_split_additional=" << name_sample_split_additional << std::endl;
    */

    hRun                    = new TH1D("hRun_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " runs; run numbers",
                                       8000, 0, 10000); // 1000, 8000

    hNElectrons             = new TH1D("hNElectrons_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " N tracks; N tracks/event",
                                       20, -0.5, 19.5); // limits? all events are 2 for nd150

    hNAPromptGgHits         = new TH1D("hNAPromptGgHits_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " N unassoc. prompt gg hits; N unassoc prompt gg hits",
                                       20, -0.5, 19.5); // limits? blank histogram for nd150?

    hTotalE                 = new TH1D("hTotalE_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " total energy; #SigmaE_{e} (MeV)",
                                       50, 0.0, 5.0); // 0.0, 2.2
                                       // TODO was 4.0

    hInternalPullee         = new TH1D("hInternalPullee_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " internal hypothesis; ee Pull",
                                       50, -40., 40.); // limits? all events within -10, 10 for nd150

    hInternalProbee         = new TH1D("hInternalProbee_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " internal hypothesis; ee Probability",
                                       50, 0.0, 1.); // limits?

    hExternalPullee         = new TH1D("hExternalPullee_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " external hypothesis; ee Pull",
                                       50, -20., 20.); // limits? all events within -30, 10 for nd150

    hExternalProbee         = new TH1D("hExternalProbee_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " external hypothesis; ee Probability",
                                       50, 0.0, 1.); // limits?

    hCosee                  = new TH1D("hCosee_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " angular correlation of electron tracks; cos(#theta)_{ee}",
                                       50, -1., 1.); // limits?

    hCoseeWeighted          = new TH1D("hCoseeWeighted_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + "corrected angular correlation of electron tracks; cos(#theta)_{ee}",
                                       50, -1., 1.); // limits?

    hVertexDZ               = new TH1D("hVertexDZ_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + "distance between recon vertices; #DeltaZ (cm)",
                                       50, -10., 10.); // limits? all events within -10, 10 for nd150

    hVertexDR               = new TH1D("hVertexDR_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + "distance between recon vertices; #DeltaR (cm)",
                                       50, -0.01, 0.01); // limits?

    hVertexDRPhi            = new TH1D("hVertexDRPhi_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + "distance between recon vertices; #DeltaR#phi (cm*rad)",
                                       50, -5., 5.); // limits? all events witnin -10, 10 for nd150

    hEeMax                  = new TH1D("hEeMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy Ee; E_{e} (MeV)"              ,
                                       50, 0.0, 5.0); // limits ok
                                       // TODO was 4.0

    hElectronLengthMax      = new TH1D("hElectronLengthMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy e^{-} track length;  track length (cm)",
                                       50, 0, 600); // limits?

    hVertexZMax             = new TH1D("hVertexZMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy vertex position;  Z (cm)"                  ,
                                       50, -120, 120); // limits?
                                       //50, -150, 150); // limits?

    hVertexSectorMax        = new TH1D("hVertexSectorMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy vertex position;  sector"                 ,
                                       50, 5.7, 5.9); // limits?
                                       //50, 5.5, 6.1); // limits?

    hVertexRMax             = new TH1D("hVertexRMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy vertex position;  R (cm)"                  ,
                                       //50, 154.7, 155.0); // limits? all events within 154.85, 154.9 for nd150
                                       50, 154.8904, 154.8908); // limits? all events within 154.85, 154.9 for nd150

    hElectronFirstGgMax     = new TH1D("hElectronFirstGgMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy e^{-} first gg hit location;first gg layer",
                                       9, -0.5, 8.5); // all in 0 for nd150

    hElectronLastGgMax      = new TH1D("hElectronLastGgMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " higher energy e^{-} last gg hit location;last gg layer",
                                       9, -0.5, 8.5); // all in 7/8 for nd150

    hVertexMinDistPromptGgMax  = new TH1D("hVertexMinDistPromptGgMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " distance from higher energy vertex to closest, prompt unassoc. gg hit;(v_{x}, v_{y}, v_{z}) - (gg_{x},gg_{y},gg_{z}) (cm)",
                                       50, 0, 600); // blank?

    hElectronLDCorrMax      = new TH1D("hElectronLDCorrMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" higher energy e^{-} LD correction;LD corr", 
                                       50, 0.9,1.1); // limits? don't know what this histogram is


    hVertexZSecMax          = new TH2D("hVertexZSecMax_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" higher energy vertex location; sector; Z (cm)",
                                       100, 5.7, 5.9,
                                       //100, 5.5, 6.1,
                                       100, -120, 120); 

    hEeMin                  = new TH1D("hEeMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy Ee; E_{e} (MeV)",
                                       50, 0.0, 5.0);

    hElectronLengthMin      = new TH1D("hElectronLengthMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy e^{-} track length;  track length (cm)",
                                       50, 0, 600);

    hVertexZMin             = new TH1D("hVertexZMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy vertex position;  Z (cm)"                  ,
                                       50, -120, 120);
                                       //50, -150, 150);

    hVertexSectorMin        = new TH1D("hVertexSectorMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy vertex position;  sector"                 ,
                                       50, 5.7, 5.9);
                                       //50, 5.5, 6.1);

    hVertexRMin             = new TH1D("hVertexRMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy vertex position;  R (cm)"                  ,
                                       //50, 154.7, 155.0);
                                       50, 154.8904, 154.8908);

    hElectronFirstGgMin     = new TH1D("hElectronFirstGgMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy e^{-} first gg hit location;first gg layer",
                                       9, -0.5, 8.5);

    hElectronLastGgMin      = new TH1D("hElectronLastGgMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy e^{-} last gg hit location;last gg layer",
                                       9, -0.5, 8.5);

    hVertexMinDistPromptGgMin   = new TH1D("hVertexMinDistPromptGgMin_" + sampleName + name_sample_split_additional + name_append,
                                           "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " distance from lower energy vertex to closest, prompt unassoc. gg hit;(v_{x},v_{y},v_{z}) - (gg_{x},gg_{y},gg_{z}) (cm)",
                                           50, 0, 600);

    hElectronLDCorrMin      = new TH1D("hElectronLDCorrMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy e^{-} LD correction;LD corr"            ,
                                       50, 0.9, 1.1);

    hVertexZSecMin          = new TH2D("hVertexZSecMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" lower energy vertex location; sector; Z (cm)"              ,
                                       150, 5.8904, 5.8907,
                                       100, -120, 120);


    hNLowEGammas            = new TH1D("hNLowEGammas_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" unassoc. scint hits E < 200 keV; N hits"     ,
                                       11, -0.5, 10.5);

    hLowEGammaEnergy        = new TH1D("hLowEGammaEnergy_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" unassoc. scint hits E < 200 keV; Energy (indiv. hits) (MeV)" ,
                                       50, 0, 0.2);

    hSummedLowEGammaE       = new TH1D("hSummedLowEGammaE_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" unassoc. scint hits E < 200 keV; Energy (indiv. hits) (MeV)" ,
                                       50, 0, 0.5);

    hLowEMinDistPromptGg    = new TH1D("hLowEMinDistPromptGg_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase "+Phase+" "+sampleName + name_sample_split_additional + name_append+" distance from unassoc. calo hits to closest, prompt unassoc. gg hit;(v_{x},v_{y},v_{z}) - (gg_{x},gg_{y},gg_{z}) (cm)",
                                       100, 0, 600.);//all low E hits


    hEeMaxVEeMin = new TH2D("hEeMaxVEeMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " electron energies; Ee^{Max} (MeV); Ee^{min} (MeV)",
                                       50, 0, 4, 50, 0, 4);

    hNAPromptGgHitsDist2VertexMin = new TH1D("hNAPromptGgHitsDist2VertexMin_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " min vertex distance to unassociated prompt geiger hit; min vertex distance (cm); y label",
                                       50, 0., 100.);
                                        



    hNAfterCuts      = new TH1D("hNAfterCuts_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Analysis Cut Flow; ; N events",
                                       20, 0, 20); //never know how many we'll make... :)

    hEffAfterCuts      = new TH1D("hEffAfterCuts_" + sampleName + name_sample_split_additional + name_append,
                                       "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Analysis Cut Flow; ; efficiency",
                                       20, 0, 20); //never know how many we'll make... :)

    hTrackSignMax   = new TH1D("hTrackSignMax_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Track Sign (Max)",
                                3, -1, 1);

    hTrackSignMin   = new TH1D("hTrackSignMin_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Track Sign (Min)",
                                3, -1, 1);

    hnGammaClusters = new TH1D("hnGammaClusters_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Num Gamma Clusters",
                                11, 0, 10);

    hnInCluster     = new TH1D("hnInCluster_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Num in Cluster",
                                11, 0, 10);

    hclusterHitEnergy     = new TH1D("hclusterHitEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Cluster Hit Energy",
                                50, 0, 0.5);
    
    hclusterHitEnergyMin     = new TH1D("hclusterHitEnergyMin_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Cluster Hit Energy Min",
                                50, 0, 0.5);

    hclusterHitEnergyMax     = new TH1D("hclusterHitEnergyMax_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Cluster Hit Energy Max",
                                50, 0, 0.5);

    hnLowEnergyHits     = new TH1D("hnLowEnergyHits_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Num in Cluster",
                                11, 0, 10);

    hclusterEnergy     = new TH1D("hclusterEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Cluster Energy",
                                50, 0, 1.0);


    hFoilSide     = new TH1D("hFoilSide_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Foil Side",
                                21, -10, 10);

    hTrueVertexR     = new TH1D("hTrueVertexR_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " True Vertex R",
                                100, -10.0, 400.0);

    hFoilSideTrueVertexR  = new TH2D("hFoilSideTrueVertexR_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Foil Side True Vertex R",
                                21, -10, 10, 100, -10.0, 400.0);

    hTrueVertexLayer     = new TH1D("hTrueVertexLayer_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " True Vertex Layer",
                                20, -5, 14);

    hSingleEnergy     = new TH1D("hSingleEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Single Electron Energy;Energy (MeV)",
                                50, 0.0, 5.0);

    hHighEnergy     = new TH1D("hHighEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " High Energy Electron Energy; Energy (MeV)",
                                50, 0.0, 5.0);

    hLowEnergy     = new TH1D("hLowEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Low Energy Electron Energy;Energy (MeV)",
                                50, 0.0, 5.0);

    hHighLowEnergy     = new TH2D("hHighLowEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Electron Energy;Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                50, 0.0, 5.0, 50, 0.0, 5.0);


    hTotalClusterEnergy     = new TH1D("hTotalClusterEnergy_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Total Cluster Hit Energy",
                                50, 0, 1.0);

    hHotSpotLocation = new TH2D("hHotSpotLocation_" + sampleName + name_sample_split_additional + name_append,
                        "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Hot Spot Locations",
                        50, 5.7371, 5.8706, 50, -120.0, 120.0);

    hEnergySum     = new TH1D("hEnergySum_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Single Electron Energy;Energy (MeV)",
                                50, 0.0, 5.0);
                                
    hEnergyDiff     = new TH1D("hEnergyDiff_" + sampleName + name_sample_split_additional + name_append,
                                "Phase " + Phase + " " + sampleName + name_sample_split_additional + name_append + " Single Electron Energy;Energy (MeV)",
                                50, 0.0, 5.0);




    //#if TRUTH_ENABLE
    // limits used for both are same as hTotalE limits, since these
    // datasets will be used to reweight 150 Nd MC and create new hTotalE
    // histogram

    // 150 Nd MC electron energy, both electrons in 2D histogram
    //hMCElectronEnergy = new TH2D("hMCElectronEnergy_" + sampleName + name_append,
    //                                "Phase " + Phase + " " + sampleName + name_append + " MC Electron Energy",
    //                                50, 0.0, 4.0);

    // 150 Nd MC true electron energy
    //hMCTrueElectronEnergy = new TH2D("hMCTrueElectronEnergy_" + sampleName + name_append,
    //                                "Phase " + Phase + " " + sampleName + name_append + " MC True Electron Energy",
    //                                50, 0.0, 4.0, 50, 0.0, 4.0);
    //#endif



    ///////////////////////////////////////////////////////////////////////////


    // map for data structure for raw data
    std::map<TString, TH1*> histogramPointers_rawdata;
    // map for data structure for processed data
    std::map<TString, TH1*> histogramPointers;


    std::map<TString, TH1*> *hpmap_p = nullptr;
    if(mode_flag == 0)
    {
        hpmap_p = &histogramPointers;
    }
    else if(mode_flag == 1)
    {
        hpmap_p = &histogramPointers_rawdata;
    }
    std::map<TString, TH1*> &hpmap{*hpmap_p};
    hpmap["hRun_"] = hRun;
    hpmap["hNElectrons_"] = hNElectrons;
    hpmap["hTotalE_"] = hTotalE;
    hpmap["hInternalPullee_"] = hInternalPullee;
    hpmap["hInternalProbee_"] = hInternalProbee;
    hpmap["hExternalPullee_"] = hExternalPullee;
    hpmap["hExternalProbee_"] = hExternalProbee;
    hpmap["hCosee_"] = hCosee;
    hpmap["hCoseeWeighted_"] = hCoseeWeighted;
    hpmap["hVertexDZ_"] = hVertexDZ;
    hpmap["hVertexDR_"] = hVertexDR;
    hpmap["hVertexDRPhi_"] = hVertexDRPhi;
    hpmap["hNAPromptGgHits_"] = hNAPromptGgHits;
    hpmap["hEeMax_"] = hEeMax;
    hpmap["hElectronLengthMax_"] = hElectronLengthMax;
    hpmap["hVertexZMax_"] = hVertexZMax;
    hpmap["hVertexSectorMax_"] = hVertexSectorMax;
    hpmap["hVertexRMax_"] = hVertexRMax;
    hpmap["hElectronFirstGgMax_"] = hElectronFirstGgMax;
    hpmap["hElectronLastGgMax_"] = hElectronLastGgMax;
    hpmap["hVertexMinDistPromptGgMax_"] = hVertexMinDistPromptGgMax;
    hpmap["hElectronLDCorrMax_"] = hElectronLDCorrMax;
    hpmap["hVertexZSecMax_"] = hVertexZSecMax;
    hpmap["hEeMin_"] = hEeMin;
    hpmap["hElectronLengthMin_"] = hElectronLengthMin;
    hpmap["hVertexZMin_"] = hVertexZMin;
    hpmap["hVertexSectorMin_"] = hVertexSectorMin;
    hpmap["hVertexRMin_"] = hVertexRMin;
    hpmap["hElectronFirstGgMin_"] = hElectronFirstGgMin;
    hpmap["hElectronLastGgMin_"] = hElectronLastGgMin;
    hpmap["hVertexMinDistPromptGgMin_"] = hVertexMinDistPromptGgMin;
    hpmap["hElectronLDCorrMin_"] = hElectronLDCorrMin;
    hpmap["hVertexZSecMin_"] = hVertexZSecMin;
    hpmap["hNLowEGammas_"] = hNLowEGammas;
    hpmap["hLowEGammaEnergy_"] = hLowEGammaEnergy;
    hpmap["hSummedLowEGammaE_"] = hSummedLowEGammaE;
    hpmap["hLowEMinDistPromptGg_"] = hLowEMinDistPromptGg;
    hpmap["hEeMaxVEeMin_"] = hEeMaxVEeMin;
    hpmap["hNAfterCuts_"] = hNAfterCuts;
    hpmap["hEffAfterCuts_"] = hEffAfterCuts;
    hpmap["hNAPromptGgHitsDist2VertexMin_"] = hNAPromptGgHitsDist2VertexMin;
    hpmap["hTrackSignMax_"] = hTrackSignMax;
    hpmap["hTrackSignMin_"] = hTrackSignMin;

    hpmap["hnGammaClusters_"] = hnGammaClusters;
    hpmap["hnInCluster_"] = hnInCluster;
    hpmap["hclusterHitEnergy_"] = hclusterHitEnergy;
    hpmap["hclusterHitEnergyMin_"] = hclusterHitEnergyMin;
    hpmap["hclusterHitEnergyMax_"] = hclusterHitEnergyMax;
    hpmap["hnLowEnergyHits_"] = hnLowEnergyHits;
    hpmap["hclusterEnergy_"] = hclusterEnergy;

    hpmap["hFoilSide_"] = hFoilSide;
    hpmap["hTrueVertexR_"] = hTrueVertexR;
    hpmap["hFoilSideTrueVertexR_"] = hFoilSideTrueVertexR;
    hpmap["hTrueVertexLayer_"] = hTrueVertexLayer;

    hpmap["hSingleEnergy_"] = hSingleEnergy;
    hpmap["hHighEnergy_"] = hHighEnergy;
    hpmap["hLowEnergy_"] = hLowEnergy;
    hpmap["hHighLowEnergy_"] = hHighLowEnergy;

    hpmap["hTotalClusterEnergy_"] = hTotalClusterEnergy;

    hpmap["hHotSpotLocation_"] = hHotSpotLocation;

    hpmap["hEnergySum_"] = hEnergySum;
    hpmap["hEnergyDiff_"] = hEnergyDiff;


    for(std::map<TString, TH1*>::iterator it{hpmap.begin()}; it != hpmap.end(); ++ it)
    {
        //std::cout << it->second->GetName() << std::endl;
        it->second->Sumw2();
    }

    //#if TRUTH_ENABLE
    //    if((sampleName.CompareTo("nd150_rot_2n2b_m4") == 0) ||
    //       (sampleName.CompareTo("nd150_rot_2b2n_m4") == 0))
    //    {
    //        hpmap["hMCElectronEnergy_"] = hMCElectronEnergy;
    //        hpmap["hMCTrueElectronEnergy_"] = hMCTrueElectronEnergy;
    //    }
    //#endif

    ///////////////////////////////////////////////////////////////////////////


    Long64_t total_event_count = 0;
    Long64_t event_pass_count = 0;


    // this idea didn't work - designed to stop program crashing when file
    // and path is not found
    /*
    TFile *aFile = new TFile;
    for(;;)
    {
        TString filepathname;

        if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
           (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
        {
            filepathname = filePath + thePath + sampleName + "/Nd150_2eNg_output_truth.root";
        }
        else
        {
            filepathname = filePath + thePath + sampleName + "/Nd150_2eNg_output.root";
        }

        //aFile = TFile::Open(filepathname);
        aFile->Open(filepathname);

        if(aFile->IsOpen())
        {
            break;
        }
        else
        {
            std::cout << "Failed to open input file " << filepathname << std::endl;
            std::cin.get();
        }
    }
    */
    TFile *aFile;
    if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
       (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
    {
        #if TRUTH_ENABLE
        //aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output_truth.root");
        aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output_truth_NEW_all_3.root");
        #else
        aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output.root");
        #endif
    }
    else
    {
        aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output.root");
    }
    TTree *theTree = (TTree*)aFile->Get("Nd150_2eNg/Nd150_2eNg");

    #if TRUTH_ENABLE
    TFile *outputFile = nullptr;
    TTree *outputTree = nullptr;
    TFile *outputFile_small = nullptr;
    TTree *outputTree_small = nullptr;
    if(mode_flag == 0)
    {
        outputFile_small = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed_small.root", "recreate");
        outputFile_small->mkdir("Nd150_2eNg");                                                                           
        outputFile_small->cd("Nd150_2eNg");                                                  
        outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
    
        #if 0
        if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
           (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
        {
            // create output tree with information saved
            /*
            outputFile = new TFile("Nd150_2eNg_output_truth_postprocessed.root", "recreate");
            outputFile->mkdir("Nd150_2eNg");                                                                           
            outputFile->cd("Nd150_2eNg");                                                  
            outputTree = new TTree("Nd150_2eNg", "Nd150_2eNg");

            outputFile_small = new TFile("Nd150_2eNg_output_truth_postprocessed_small.root", "recreate");
            outputFile_small->mkdir("Nd150_2eNg");                                                                           
            outputFile_small->cd("Nd150_2eNg");                                                  
            outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
            */
            //#if TRUTH_ENABLE
            /*
            outputFile = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_truth_postprocessed.root", "recreate");
            outputFile->mkdir("Nd150_2eNg");                                                                           
            outputFile->cd("Nd150_2eNg");                                                  
            outputTree = new TTree("Nd150_2eNg", "Nd150_2eNg");
            */

            /*
            outputFile_small = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_truth_postprocessed_small.root", "recreate");
            outputFile_small->mkdir("Nd150_2eNg");                                                                           
            outputFile_small->cd("Nd150_2eNg");                                                  
            outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
            */
            //#else
            /*
            outputFile = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed.root", "recreate");
            outputFile->mkdir("Nd150_2eNg");                                                                           
            outputFile->cd("Nd150_2eNg");                                                  
            outputTree = new TTree("Nd150_2eNg", "Nd150_2eNg");
            */
            // NOTE: changed 2020-08-28: all files uniform naming scheme
            outputFile_small = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed_small.root", "recreate");
            outputFile_small->mkdir("Nd150_2eNg");                                                                           
            outputFile_small->cd("Nd150_2eNg");                                                  
            outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
            //#endif
        }
        else
        {
            /*
            outputFile = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed.root", "recreate");
            outputFile->mkdir("Nd150_2eNg");                                                                           
            outputFile->cd("Nd150_2eNg");                                                  
            outputTree = new TTree("Nd150_2eNg", "Nd150_2eNg");
            */
            outputFile_small = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed_small.root", "recreate");
            outputFile_small->mkdir("Nd150_2eNg");                                                                           
            outputFile_small->cd("Nd150_2eNg");                                                  
            outputTree_small = new TTree("Nd150_2eNg", "Nd150_2eNg");
        }
        #endif
    }
    #endif

    // general
    int         event;
    int         run;
    int         runStatus; // TODO
    int         nElectrons; // TODO

    double      radonWeight;
    double      bi210Weight;
    int         foilSide;
    double      eventTime;
    double      trueVertexR;
    double      trueVertexZ;
    double      trueVertexSector;
    int         trueVertexLayer;

    double      electronEnergy[2];
    double      eTrackLength[2];
    int         electronSide[2];
    double      trackSign[2];
    //double      radiusOfCurvature[2]; // TODO - branch not in Nd150 tree
    double      electronMeasTime[2];
    double      electronDMeasTime[2];
    int         electronBlockType[2];


    double      eeInternalPull;
    double      eeInternalProb;
    double      eeExternalPull;
    double      eeExternalProb;
    double      cosee;
    double      cosee_weight;

    int         electronPMT[2];
    int         electronLDFlag[2];
    double      electronLDCorr[2];
    double      electronLDCorrErr[2];

    double      vertexZ[2];
    double      vertexSec[2];
    double      vertexR[2];

    bool        vertexInHotSpot[2];

    int         firstGgHitLayer[2];
    int         lastGgHitLayer[2];

    int         NAPromptGgHits; // TODO - max value is 20 ? EDIT: APPEARS TO BE 5
    int         NAPromptGgHitsSide[5];
    double      NAPromptGgHitsDist2Vertex[5]; // TODO check number 20 - TODO: NUMBER IS 5 NOT 20
    double      NAPromptGgHitsDist2Calo[5];
    
    
    int         nGammaClusters;
    int         nInCluster[10];                          
    double      clusterEnergy[10];
    double      clusterTimeSpan[10];


    //double clusterLength : clusterLength[nGammaClusters]/D                    *
    //double eMingInternalVertexPull : eMingInternalVertexPull[nGammaClusters]/D*
    //double eMingInternalVertexProb : eMingInternalVertexProb[nGammaClusters]/D*
    //double coseMing  : coseMing[nGammaClusters]/D                             *
    //double eMaxgInternalVertexPull : eMaxgInternalVertexPull[nGammaClusters]/D*
    //double eMaxgInternalVertexProb : eMaxgInternalVertexProb[nGammaClusters]/D*
    //double coseMaxg  : coseMaxg[nGammaClusters]/D                             *
    

    int         nTotalClusterHits; // TODO - check max value is 10 ?
    // TODO: MAX VALUE IS 20 NOT 10!
    double      clusterHitEnergy[20];
    int         clusterHitPMT[20];
    int         clusterHitLDFlag[20];
    double      clusterHitLDCorr[20];
    double      clusterHitLDCorrErr[20];
    double      clusterHitSec[20];
    double      clusterHitZ[20];

    #if TRUTH_ENABLE
    double      trueElectronEnergy[2]; // TODO
    #endif


    // these are the old variable names
    //int         nLowEClusters;
    //double      lowEClusterE[10];
    //double      lowEClusterSector[10];
    //double      lowEClusterZ[10];
    //double      lowEClusterR[10];
    //double      lowEMinDistPromptGg[10];
    //int         lowEClusterPMTs[10];
    // these are the branches they were mapped to
    //theTree->SetBranchAddress("nTotalClusterHits"  , &nLowEClusters);
    //theTree->SetBranchAddress("clusterHitEnergy"   , lowEClusterE);
    //theTree->SetBranchAddress("clusterHitPMT"      , lowEClusterPMTs);
    //theTree->SetBranchAddress("clusterHitSec"      , lowEClusterSector);
    //theTree->SetBranchAddress("clusterHitZ"        , lowEClusterZ);
    // theTree->SetBranchAddress( "lowEClusterR"      , lowEClusterR );
    // theTree->SetBranchAddress( "lowEMinDistPromptGg", lowEMinDistPromptGg );


    theTree->SetBranchAddress("Event"               , &event);
    theTree->SetBranchAddress("Run"                 , &run);
    theTree->SetBranchAddress("runStatus"           , &runStatus);  
    theTree->SetBranchAddress("nElectrons"          , &nElectrons);  
    //theTree->SetBranchAddress( "nGammaClusters"   , &nGammas );  
    //  theTree->SetBranchAddress( "nParticles"     , &nParticles );  

    theTree->SetBranchAddress("radonWeight"         , &radonWeight);
    theTree->SetBranchAddress("bi210Weight"         , &bi210Weight);
    theTree->SetBranchAddress("foilSide"            , &foilSide);
    theTree->SetBranchAddress("eventTime"           , &eventTime);

    theTree->SetBranchAddress("trueVertexR"         , &trueVertexR);
    theTree->SetBranchAddress("trueVertexZ"         , &trueVertexZ);
    theTree->SetBranchAddress("trueVertexSector"    , &trueVertexSector);
    theTree->SetBranchAddress("trueVertexLayer"     , &trueVertexLayer);

    theTree->SetBranchAddress("electronEnergy"      , electronEnergy);
    theTree->SetBranchAddress("eTrackLength"        , eTrackLength);
    theTree->SetBranchAddress("electronSide"        , electronSide);
    theTree->SetBranchAddress("trackSign"           , trackSign);
    //theTree->SetBranchAddress("radiusOfCurvature" , radiusOfCurvature);
    theTree->SetBranchAddress("electronMeasTime"    , electronMeasTime);
    theTree->SetBranchAddress("electronDMeasTime"   , electronDMeasTime);
    theTree->SetBranchAddress("electronBlockType"   , electronBlockType);

    theTree->SetBranchAddress("internalPull"        , &eeInternalPull);
    theTree->SetBranchAddress("internalProb"        , &eeInternalProb);
    theTree->SetBranchAddress("externalPull"        , &eeExternalPull);
    theTree->SetBranchAddress("externalProb"        , &eeExternalProb);
    theTree->SetBranchAddress("cosee"               , &cosee);
    theTree->SetBranchAddress("cosee_weight"        , &cosee_weight);

    theTree->SetBranchAddress("electronPMT"         , electronPMT);
    theTree->SetBranchAddress("electronLDFlag"      , electronLDFlag);
    theTree->SetBranchAddress("electronLDCorr"      , electronLDCorr);
    theTree->SetBranchAddress("electronLDCorrErr"   , electronLDCorrErr);

    theTree->SetBranchAddress("vertexZ"             , vertexZ);
    theTree->SetBranchAddress("vertexSec"           , vertexSec);
    theTree->SetBranchAddress("vertexR"             , vertexR);

    theTree->SetBranchAddress("vertexInHotSpot"     , vertexInHotSpot);

    theTree->SetBranchAddress("firstGgHitLayer"     , firstGgHitLayer);
    theTree->SetBranchAddress("lastGgHitLayer"      , lastGgHitLayer);

    // number of unassociated prompt geiger hits
    theTree->SetBranchAddress("NAPromptGgHits"              , &NAPromptGgHits);
    theTree->SetBranchAddress("NAPromptGgHitsSide"          , NAPromptGgHitsSide);
    theTree->SetBranchAddress("NAPromptGgHitsDist2Vertex"   , NAPromptGgHitsDist2Vertex);
    theTree->SetBranchAddress("NAPromptGgHitsDist2Calo"     , NAPromptGgHitsDist2Calo);

    // not sure if these branches are in the tree?
    //theTree->SetBranchAddress("electronCaloZ"     ,     electronCaloZ);
    //theTree->SetBranchAddress("electronCaloR"     ,     electronCaloR);
    //theTree->SetBranchAddress("electronCaloSector",     electronCaloSector);
    //  theTree->SetBranchAddress( "electronPMT"       , electronPMT );

    theTree->SetBranchAddress("nGammaClusters"              , &nGammaClusters);
    theTree->SetBranchAddress("nInCluster"                 , &nInCluster);
    theTree->SetBranchAddress("clusterEnergy"               , clusterEnergy);
    theTree->SetBranchAddress("clusterTimeSpan"             , clusterTimeSpan);
    //double clusterLength : clusterLength[nGammaClusters]/D                    *
    //double eMingInternalVertexPull : eMingInternalVertexPull[nGammaClusters]/D*
    //double eMingInternalVertexProb : eMingInternalVertexProb[nGammaClusters]/D*
    //double coseMing  : coseMing[nGammaClusters]/D                             *
    //double eMaxgInternalVertexPull : eMaxgInternalVertexPull[nGammaClusters]/D*
    //double eMaxgInternalVertexProb : eMaxgInternalVertexProb[nGammaClusters]/D*
    //double coseMaxg  : coseMaxg[nGammaClusters]/D                             *
    //

    theTree->SetBranchAddress("nTotalClusterHits"           , &nTotalClusterHits);
    theTree->SetBranchAddress("clusterHitEnergy"            , clusterHitEnergy);
    theTree->SetBranchAddress("clusterHitPMT"               , clusterHitPMT);
    theTree->SetBranchAddress("clusterHitLDFlag"            , clusterHitLDFlag);
    theTree->SetBranchAddress("clusterHitLDCorr"            , clusterHitLDCorr);
    theTree->SetBranchAddress("clusterHitLDCorrErr"         , clusterHitLDCorrErr);
    theTree->SetBranchAddress("clusterHitSec"               , clusterHitSec);
    theTree->SetBranchAddress("clusterHitZ"                 , clusterHitZ);
    // NOTE: these were already removed by summer
    // theTree->SetBranchAddress( "lowEClusterR"      , lowEClusterR );
    // theTree->SetBranchAddress( "lowEMinDistPromptGg", lowEMinDistPromptGg );
   
    #if TRUTH_ENABLE
    if(mode_flag == 0)
    {
        if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
           (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
        {
            theTree->SetBranchAddress("trueElectronEnergy" , trueElectronEnergy);
        }
        else
        {
            trueElectronEnergy[0] = 0.;
            trueElectronEnergy[1] = 0.;
        }

        ///////////////////////////////////////////////////////////////////////
        // do this for all sample types
        ///////////////////////////////////////////////////////////////////////

        // big tree, all information

        /*
        outputTree->Branch("Event"                      , &event                        , "event/I");
        outputTree->Branch("Run"                        , &run                          , "run/I");
        outputTree->Branch("runStatus"                  , &runStatus                    , "runStatus/I");
        outputTree->Branch("nElectrons"                 , &nElectrons                   , "nElectrons/I");

        outputTree->Branch("radonWeight"                , &radonWeight                  , "radonWeight/D");
        outputTree->Branch("bi210Weight"                , &bi210Weight                  , "bi210Weight/D");
        outputTree->Branch("foilSide"                   , &foilSide                     , "foilSide/I");
        outputTree->Branch("eventTime"                  , &eventTime                    , "eventTime/D");

        outputTree->Branch("trueVertexR"                , &trueVertexR                  , "trueVertexR/D");
        outputTree->Branch("trueVertexZ"                , &trueVertexZ                  , "trueVertexZ/D");
        outputTree->Branch("trueVertexSector"           , &trueVertexSector             , "trueVertexSector/D");
        outputTree->Branch("trueVertexLayer"            , &trueVertexLayer              , "trueVertexLayer/I");

        outputTree->Branch("electronEnergy"             , electronEnergy                , "electronEnergy[2]/D");
        outputTree->Branch("eTrackLength"               , eTrackLength                  , "eTrackLength[2]/D");
        outputTree->Branch("electronSide"               , electronSide                  , "electronSide[2]/I");
        outputTree->Branch("trackSign"                  , trackSign                     , "trackSign[2]/D");
        outputTree->Branch("electronMeasTime"           , electronMeasTime              , "electronMeasTime[2]/D");
        outputTree->Branch("electronDMeasTime"          , electronDMeasTime             , "electronDMeasTime[2]/D");
        outputTree->Branch("electronBlockType"          , electronBlockType             , "electronBlockType[2]/I");

        outputTree->Branch("internalPull"               , &eeInternalPull               , "eeInternalPull/D");
        outputTree->Branch("internalProb"               , &eeInternalProb               , "eeInternalProb/D");
        outputTree->Branch("externalPull"               , &eeExternalPull               , "eeExternalPull/D");
        outputTree->Branch("externalProb"               , &eeExternalProb               , "eeExternalProb/D");
        outputTree->Branch("cosee"                      , &cosee                        , "cosee/D");
        outputTree->Branch("cosee_weight"               , &cosee_weight                 , "cosee_weight/D");

        outputTree->Branch("electronPMT"                , electronPMT                   , "electronPMT[2]/I");
        outputTree->Branch("electronLDFlag"             , electronLDFlag                , "electronLDFlag[2]/I");
        outputTree->Branch("electronLDCorr"             , electronLDCorr                , "electronLDCorr[2]/D");
        outputTree->Branch("electronLDCorrErr"          , electronLDCorrErr             , "electronLDCorrErr[2]/D");

        outputTree->Branch("vertexZ"                    , vertexZ                       , "vertexZ[2]/D");
        outputTree->Branch("vertexSec"                  , vertexSec                     , "vertexSec[2]/D");
        outputTree->Branch("vertexR"                    , vertexR                       , "vertexR[2]/D");

        outputTree->Branch("vertexInHotSpot"            , vertexInHotSpot               , "vertexInHotSpot[2]/O");

        outputTree->Branch("firstGgHitLayer"            , firstGgHitLayer               , "firstGgHitLayer[2]/I");
        outputTree->Branch("lastGgHitLayer"             , lastGgHitLayer                , "lastGgHitLayer[2]/I");

        outputTree->Branch("NAPromptGgHits"             , &NAPromptGgHits               , "NAPromptGgHits/I");
        outputTree->Branch("NAPromptGgHitsSide"         , NAPromptGgHitsSide            , "NAPromptGgHitsSide/I"); // TODO: should be array?
        outputTree->Branch("NAPromptGgHitsDist2Vertex"  , NAPromptGgHitsDist2Vertex     , "NAPromptGgHitsDist2Vertex[NAPromptGgHits]/D");
        outputTree->Branch("NAPromptGgHitsDist2Calo"    , NAPromptGgHitsDist2Calo       , "NAPromptGgHitsDist2Calo[NAPromptGgHits]/D");

        outputTree->Branch("nGammaClusters"             , &nGammaClusters               , "nGammaClusters/I");
        outputTree->Branch("nInCluster"                 , &nInCluster                   , "nInCluster[nGammaClusters]/I");
        outputTree->Branch("clusterEnergy"              , clusterEnergy                 , "clusterEnergy[nGammaClusters]/D");
        outputTree->Branch("clusterTimeSpan"            , clusterTimeSpan               , "clusterTimeSpan[nGammaClusters]/D");

        outputTree->Branch("nTotalClusterHits"          , &nTotalClusterHits            , "nTotalClusterHits/I");
        outputTree->Branch("clusterHitEnergy"           , clusterHitEnergy              , "clusterHitEnergy[nTotalClusterHits]/D");
        outputTree->Branch("clusterHitPMT"              , clusterHitPMT                 , "clusterHitPMT[nTotalClusterHits]/I");
        outputTree->Branch("clusterHitLDFlag"           , clusterHitLDFlag              , "clusterHitLDFlag[nTotalClusterHits]/I");
        outputTree->Branch("clusterHitLDCorr"           , clusterHitLDCorr              , "clusterHitLDCorr[nTotalClusterHits]/D");
        outputTree->Branch("clusterHitLDCorrErr"        , clusterHitLDCorrErr           , "clusterHitLDCorrErr[nTotalClusterHits]/D");
        outputTree->Branch("clusterHitSec"              , clusterHitSec                 , "clusterHitSec[nTotalClusterHits]/D");
        outputTree->Branch("clusterHitZ"                , clusterHitZ                   , "clusterHitZ[nTotalClusterHits]/D");

        outputTree->Branch("trueElectronEnergy"         , trueElectronEnergy            , "trueElectronEnergy[2]/D"); 
        */

        // small tree, essential information

        outputTree_small->Branch("Run"                          , &run                      , "run/I");
        outputTree_small->Branch("eventTime"                    , &eventTime                , "eventTime/D");
        //outputTree_small->Branch("nElectrons"                 , &nElectrons               , "nElectrons/I"); // could remove
        outputTree_small->Branch("radonWeight"                  , &radonWeight              , "radonWeight/D");
        outputTree_small->Branch("bi210Weight"                  , &bi210Weight              , "bi210Weight/D");
        outputTree_small->Branch("electronEnergy"               , electronEnergy            , "electronEnergy[2]/D");
        // truth information, save in all files to maintain homogenaity
        outputTree_small->Branch("trueElectronEnergy"           , trueElectronEnergy        , "trueElectronEnergy[2]/D");
    }
    #else
    if(mode_flag == 0)
    {
        // TODO:
        // don't really care about filling this in because never run code
        // without TRUTH_ENABLE now
    }
    #endif

    if(mode_flag == 0)
    {
        std::cout << "Processing: " << sampleName << std::endl;
    }

    Int_t electron_cut_A = 0;
    //Int_t electron_cut_B = 0;

    int nCuts = 17; // 26;
    Int_t cut_counter[nCuts];
    //Int_t cut_counter_index = 0;
    for(int i = 0; i < nCuts; ++ i)
    {
        cut_counter[i] = 0;
    }

    Long_t events = (Long_t)theTree->GetEntries();
    for(Long_t event_i = 0; event_i < events; event_i++)
    {

        bool cut = false;
        //if(event_i % 1000 == 0)
        //{
        //    std::cout << "\r Processing complete : " << 100 * event_i / events << "%" << std::flush;
        //}
        theTree->GetEvent(event_i);


        if((electronEnergy[0] < 0.3) || (electronEnergy[1] < 0.3))
        {
            ++ electron_cut_A;
        }
        
        ++ total_event_count;







        double weight = 1.;
        if(sampleName.Contains("bi214_swire"))
        {
            weight = radonWeight;
        }
        // TODO: due to secular equlibrium, the following MAY be required
        // (unsure at present time)
        else if(sampleName.Contains("pb214_swire"))
        {
            weight = radonWeight;
        }
        else if(sampleName.Contains("bi207"))
        {
            // halflife
            const double T12 = 31557600. * 32.97818366312;
            // https://periodictable.com/Isotopes/083.207/index2.full.dm.html
            
            // lambda
            const double lambda = std::log(2.0) / T12;
            double arg = -lambda * eventTime;
            // I assume event time is in SI units
            weight = std::exp(arg);
        }
        else if(sampleName.Contains("bi210"))
        {
            // halflife
            // https://periodictable.com/Isotopes/083.210/index.p.full.dm.html
            //const double T21 = 86400. * 5.011574074074;

            //https://periodictable.com/Isotopes/082.210/index.p.full.dm.html
            const double T21 = 31557600. * 22.22856418062;
            
            // lambda
            const double lambda = std::log(2.0) / T21;
            double arg = -lambda * eventTime;
            // I assume event time is in SI units
            weight = bi210Weight * std::exp(arg);
        }
        else if(sampleName.Contains("co60"))
        {
            // halflife
            // https://periodictable.com/Isotopes/027.60/index.p.html
            const double T21 = 31557600. * 5.274923896499;
            // lambda
            const double lambda = std::log(2.0) / T21;
            double arg = -lambda * eventTime;
            // I assume event time is in SI units
            weight = std::exp(arg);
        }
        else if(sampleName.Contains("eu152"))
        {
            // halflife
            // https://periodictable.com/Isotopes/063.152/index.p.html
            const double T21 = 31557600. * 13.537;
            // lambda
            const double lambda = std::log(2.0) / T21;
            double arg = -lambda * eventTime;
            // I assume event time is in SI units
            weight = std::exp(arg);
        }
        else if(sampleName.Contains("eu154"))
        {
            // halflife
            // https://periodictable.com/Isotopes/063.154/index.html
            const double T21 = 31557600. * 8.593;
            // lambda
            const double lambda = std::log(2.0) / T21;
            double arg = -lambda * eventTime;
            // I assume event time is in SI units
            weight = std::exp(arg);
        }
        // TODO: can fill 208Tl with 208Tl weight here instead of changing weight arbitrarily later
        // TODO: I get completely different numbers
        /*    else if ( histName.CompareTo("bi210_swire") == 0 ) {weight = bi210Weight * pow(0.5,eventTime/433123.);}
        else if ( (histName.CompareTo("bi210_swire") == 0) || (histName.CompareTo("bi210_sscin") == 0) || (histName.CompareTo("bi210_sfoil") == 0)) {weight = bi210Weight * pow(0.5,eventTime/433123.);}
        else if ((histName.CompareTo("co60_cuTower")==0) || (histName.CompareTo("co60_steel")==0) || (histName.CompareTo("co60_muMetal")==0) || (histName.CompareTo("co60_cuPetals")==0) )weight = pow(0.5,(eventTime/166344192.));
        else if ( (histName.CompareTo("eu152_sscin")==0) || (histName.CompareTo("eu152_int")==0)) {weight = pow(0.5,(eventTime/4.265e8));}
        else if ((histName.CompareTo("eu154_int")==0)) {weight = pow(0.5,(eventTime/2.71169e8));}
        */
        
        // TODO: re-implement this
        // http://nile.hep.utexas.edu/DocDB/ut-nemo/docs/0037/003716/005/position_paper.pdf
        // Page 20. Section on short-lived isotopes.
        // Top of page 31 has Bi events being weighted by half-life of Pb
        // (Doing a ctrl-f on weight works quite well, but those are the main two.)
        // Table 9.1
        // Some isotopes contributing to the background rates have relatively short half-lives
        // on the time scale of the NEMO-3 experiment. The events from these isotopes are
        // therefore weighted according to their known half-lives and the time of the event.
        // I assume that in order to do this, one must have a measurement of the activity of
        // these isotopes at a particular time, and then the weight is given by
        // A * exp( -k (t - t0) ), where A is the initial weight measured at time t0 (A=1), k is a decay rate
        // On Page 31 I find this "Events from 210Bi are weighted by the half life of 210Pb as described in Section 3.4.1.
        // Therefore, the activities reported here correspond to the activity as of 01/02/03"
        // I guess 01/02/03 is 1st Feb 2003, and is t0 for isotopes 210Bi and 210Pb
        // Similar info from Table 4.8 "Events from 207Bi and 152,154Eu have been weighted by their respective
        // half-lives as described in Section 3.4.1, and therefore the activities reported here correspond to their activities as of 01/02/03."
        //


        // count number of events before cuts
        // count the initial number of events
        int cc = 0;
        ++ cut_counter[cc]; // cut 0 - input count
        ++ cc;






        // additional sample splitting
        // side = 1: OUT
        // side = 0: IN
        // side = -1: FAIL
        // layer: 0 - 8 inclusive = valid
        // -1 = fail
        #if 0
        if(mode_flag_2 == 1)
        {
            if(sample_split_flags == 0x01)
            {
                // IN
                if(foilSide != 0)
                {
                    continue;
                }
            }
            else if(sample_split_flags == 0x02)
            {
                // OUT
                if(foilSide != 1)
                {
                    continue;
                }
            }
            else if(sample_split_flags == 0x04)
            {
                // L0 IN
                name_sample_split_additional = "_INL0";
                if(foilSide != 0)
                {
                    continue;
                }
                if(trueVertexLayer != 0)
                {
                    continue;
                }
            }
            else if(sample_split_flags == 0x08)
            {
                // L0 OUT
                name_sample_split_additional = "_OUTL0";
                if(foilSide != 1)
                {
                    continue;
                }
                if(trueVertexLayer != 0)
                {
                    continue;
                }
            }
            else if(sample_split_flags == 0x10)
            {
                // L>0 IN+OUT
                name_sample_split_additional = "_Lg0";
                if(trueVertexLayer < 1 || trueVertexLayer > 8)
                {
                    continue;
                }
            }
        }
        #endif


        ++ cut_counter[cc]; // cut 0 - input count
        ++ cc;






        // 1: Accepted Run
        // (1) Accepted Run

        // is that the same as Phase cut?
        

        ///////////////////////////////////////////////////////////////////////
        // TODO: probably do NOT want to fill this histo with "weight" variable
        // do not expect to see these used in the fit
        // TODO: this is filled before all cuts, including phase cut!
        // in Summers code it was filled at the end!
        //hRun->Fill(run, weight);
        ///////////////////////////////////////////////////////////////////////


        // phase cut moved to end


        // maybe 5.87
        //if(vertexSec[0] < 5.74 || vertexSec[0] > 5.88 || vertexSec[1] < 5.74 || vertexSec[1] > 5.88)
        // The following values taken from Position Paper Section 2,
        // seems to match data from 2e directory
        // No cut in 2eNg, do not know why
        if(mode_flag_2 == 1)
        {
            //if(vertexSec[0] < 5.7371 + 0.001 || vertexSec[0] > 5.8706 - 0.001 ||
            //   vertexSec[1] < 5.7571 + 0.001 || vertexSec[1] > 5.8706 - 0.001)
            //if(vertexSec[0] < 5.7371 + 0.000 || vertexSec[0] > 5.8706 - 0.001 ||
            //   vertexSec[1] < 5.7571 + 0.000 || vertexSec[1] > 5.8706 - 0.001)
            if(vertexSec[0] < 5.7371 || vertexSec[0] > 5.8706 || vertexSec[1] < 5.7371 || vertexSec[1] > 5.8706)
            // when working with 2e, this should not cut any events
            //if(vertexSec[0] < 5.7 || vertexSec[0] > 5.9 || vertexSec[1] < 5.7 || vertexSec[1] > 5.9)
            {
                continue;
            }
        }
        ++ cut_counter[cc]; // cut 2 - vertex sector
        ++ cc;


        ///////////////////////////////////////////////////////////////////////
        // TODO: probably do NOT want to fill this histo with "weight" variable
        // do not expect to see these used in the fit
        //hNElectrons->Fill(nElectrons, weight);
        //moved to fill after cuts
        ///////////////////////////////////////////////////////////////////////

    /*

        // 2: Two electrons from 150Nd foil
        // (not 1, 2) - this is number not 1 -> 2 in position paper
        // (2) Two Electrons from 150Nd foil
        // 1. Exactly two electrons with reconstructed vertices inside of the 150Nd source foil boundaries.
        if(nElectrons != 2) continue;
       
        ++ cut_counter[cc];
        ++ cc;

    */

        // Set the high and low energy index
        int highE_index = -1;
        int lowE_index = -1;
        if(electronEnergy[0] > electronEnergy[1])
        {
            highE_index = 0;
            lowE_index = 1;
        }
        else
        {
            highE_index = 1;
            lowE_index = 0;
        }



        if(vertexInHotSpot[0])
        {
            hHotSpotLocation->Fill(vertexSec[0], vertexZ[0], 1.0);
        }
        if(vertexInHotSpot[1])
        {
            hHotSpotLocation->Fill(vertexSec[1], vertexZ[1], 1.0);
        }


    /*
     *
     * TODO: some of these were making cuts that removed events, but they should not have!

        ///////////////////////////////////////////////////////
        // here is a summary of all the cuts from references //
        ///////////////////////////////////////////////////////
        
        // *******************************
        // *** Position Paper, Page 61 ***
        // *******************************
        //
        // 1. Exactly two electrons with reconstructed vertices inside of the 150Nd source foil boundaries.
        // 2. The separation between each individually reconstructed vertex must be <= 4 cm radially
        // and <= 8 cm longitudinally to ensure that the electrons originate from a common vertex.
        // 3. Electron tracks with vertices located in hot spot regions are not accepted in order to
        // improve the signal to background ratio in this channel.
        // 4. Each of the electron tracks must be at least 30 cm long, have a hit in the first Gg cell layer,
        // and be associated to unique scintillator blocks.
        // 5. Tracks extrapolated to the petal blocks in the calorimeter layer closest to the source foil
        // are not accepted as their energy calibration is less precise than the other block types.
        // 6. The impact of the tracks must be on the front face of the scintillator block for proper light
        // collection.
        // 7. The energy of each electron is required to be at least 0.3 MeV, as single electron energies
        // below this threshold are not well described by the simulation.
        // 8. No alpha candidates are allowed within 15 cm of the electron vertices in order to reduce the
        // event rate coming from 214Bi decays.
        // 9. To further improve this rejection of 214Bi, no more than one prompt Geiger hit which is
        // not associated to a track is allowed within 15 cm of the track vertices. If the tracks on the
        // same side of the foil, no unassociated prompt Geiger hits are allowed on the other side of
        // the foil.
        // This allows for some level of Geiger cell noise which is not simulated, while rejecting a large
        // number of 214Bi decays.
        // 10. No unassociated Geiger hits are allowed within 15 cm of the calorimeter hits as this is
        // indicative of a reflected electron, which therefore has not deposited its full energy in the
        // calorimeter.
        // 11. Additional, isolated calorimeter hits which do not have a track associated to them are
        // allowed if the single hit energies are less than 0.2 MeV (indicative of noise which is not
        // simulated).
        // 12. The internal probability between the two electrons must be at least 1% and the external
        // probability must be no more than 1%

        // *****************************************
        // *** Position Paper, Table 6.1 Page 63 ***
        // *****************************************
        //
        // (1) Accepted Run
        // (2) Two electrons from 150Nd foil
        // (3) No other tracks in event
        // (4) Tracks intersect foil
        // (5) No alpha candidates near vertex
        // (6) Ee >= 0.2 MeV
        // (7) <= 20 total free scintillators and <= 10 gamma clusters (data management)
        // (8) All PMTs have good status
        // (9) Isolated electron calorimeter hits
        // (10) Electron hits in different scintillator blocks
        // (11) Impacts on block faces
        // (12) Electron LD flags <= 2
        // (13) Electrons have LD corrections
        // (14) Electron HS flags = 0
        // (15) Gamma hits have LD flags <= 2
        // (16) Gamma hits have LD corrections
        // (17) No gamma hits with HS = 1
        // (18) Tracks start in Gg layer 0 or 1
        // (19) Tracks miss <= 1 gg layer
        // (20) Track lengths each > 20 cm
        // (21) <= 5 NA prompt Gg hits
        // (22) Electron block types != L2/L3
        // (23) No gamma with E_gamma > 200 keV
        // (24) Tracks have negative curvature
        // (25) Track lengths > 30 cm
        // (26) Tracks not from hot spots
        // (27) No NA prompt Gg hits within 15 cm of scint hits
        // (28) <= 1 NA prompt Gg hit within 15 cm of vertex
        // (29) No NA prompt Gg hits on opposite side of foil if tracks are on same side
        // (30) Both tracks have hit in Gg L0
        // (31) Pint > 0.01 and Pext < 0.01
        // (32) Delta R <= 4 cm and delta Z <= 8 cm
        // (33) Each E_e >= 300 keV

        // *********************
        // *** Summer Thesis ***
        // *********************
        //
        // 1: Accepted Run
        // 2: Two electrons from 150Nd foil
        // 3: No other tracks in event
        // 4: No alpha candidates near vertex
        // 5: E_e >= 0.2 MeV
        // 6: <= 20 total free scintillators and <= 10 gamma clusters
        // 7: All PMTs have good status
        // 8: Isolated electron calorimeter hits
        // 9: Electron hits in different scintillator blocks
        // 10: Impacts on block faces
        // 11: All PMTs have good laser survey and counting rates identified by LD and HS flags
        // 12: Tracks start in Gg layer 0 or 1
        // 13: Tracks miss <= 1 Gg layer before impact on calorimeter
        // 14: Each L_e > 20 cm
        // 15: <= 5 unassociated prompt Gg hits
        // 16: Electron block types != L2/L3
        // 17: No gamma rays in event
        // 18: Tracks have negative curvature
        // 19: Each L_e ≥ 30 cm
        // 20: Tracks not from hot spots
        // 21: <= 1 unassociated, prompt Gg hit within 15 cm in XY coordinate of event vertex
        // 22: No unassociated, prompt Gg hits on opposite side of foil if tracks are on same side
        // 23: Both tracks have hit in Gg L0
        // 24: P_int >= 0.01 and P_ext <= 0.01
        // 25: Delta R <= 4 cm and delta Z <= 8 cm
        // 26: Each E_e >= 300 keV

        // Here is a map of how all the above points connect
        //
        // 1: (1) Accepted Run
        // 2: (2) Two electrons
        // 2: (2) from Nd150 foil
        //
        // 25: (32) 2. Delta R <= 4cm and Delta Z <= 8 cm
        //
        //


        // Numbers followed by colon : are from Summers Thesis
        // Numbers in (parenthesys) are from Table 6.1 of Page 63 from Position Paper
        // http://nile.hep.utexas.edu/DocDB/ut-nemo/docs/0037/003716/005/position_paper.pdf
        // Numbers followed by a period xx. are from List on Page 61 of Position paper

        // 3: No other tracks in event
        // (3) No other tracks in event
        // (4) Tracks intersect foil
        ++ cut_counter[cc];
        ++ cc;
        
        // 4: No alpha candidates near vertex
        // (5) No alpha candidates near vertex
        // TODO: what does near mean?
        ++ cut_counter[cc];
        ++ cc
        
        // 5: E_e >= 0.2 MeV
        // (6) E_e >= 0.2 MeV
        // (23) No gamma with E_gamma > 200 keV
        if( (electronEnergy[0] < 0.2) || (electronEnergy[1] < 0.2) ) continue;
        ++ cut_counter[cc];
        ++ cc;

        // 6: <= 20 total free scintillators and <= 10 gamma clusters
        // (7) <= 20 total free scintillators and <= 10 gamma clusters (data management)
        // NOTE: nothing to be done here, this cut is done in data management
        // by file /home/blotsd/NEMO3/release_1.0.0/SDBanalysis/src/Nd150_2e.cpp
        // by file /home/blotsd/NEMO3/release_1.0.0/SDBanalysis/src/Nd150_2eNg.cpp
        if(nGammaClusters > 10) continue; 
        ++ cut_counter[cc];
        ++ cc;

        // 7: All PMTs have good status
        // (8) All PMTs have good status
        // TODO : not a clue how to implement this
        ++ cut_counter[cc];
        ++ cc;

        // 8: Isolated electron calorimeter hits
        // (9) Isolated electron calorimeter hits
        // TODO: again not a clue
        ++ cut_counter[cc];
        ++ cc;

        // 9: Electron hits in different scintillator blocks
        // (10) Electron hits in different scintillator blocks
        // [3]
        // NOTE: This is the same as cut number 27 and [3] 
        // 27: Each electron track must be associated to a unique scinitllator block
        // [3] this index appears multiple times
        if(electronPMT[0] == electronPMT[1]) continue;
        ++ cut_counter[cc];
        ++ cc

        // 10: Impacts on block faces
        // (11) Impacts on block faces
        // TODO: again not a clue
        ++ cut_counter[cc];
        ++ cc;

        // 11: All PMT have good laser survey and counting rates identified by LD and HS flags
        // (12) Electron LD flags <= 2
        // (13) Electrons have LD corrections
        // (14) Electron HS flags = 0
        // (15) Gamma hits have LD flags <= 2
        // (16) Gamma hits have LD corrections
        // (17) No gamma hits with HS = 1
        // TODO: not sure about this, this might only be 1 part of it
        if((electronLDFlag[0] > 0) || (electronLDFlag[1] > 0)) continue;
        ++ cut_counter[cc];
        ++ cc
        
        // 12: Tracks start in Gg layer 0 or 1
        // (18) Tracks start in Gg later 0 or 1
        if((firstGgHitLayer[0] > 1) || (firstGgHitLayer[1] > 1)) continue;
        ++ cut_counter[cc];
        ++ cc;
        
        // 13: Tracks miss <= 1 Gg layer before impact on calorimeter
        // (19) Tracks miss <= 1 gg layer
        // TODO: not a clue
        ++ cut_counter[cc];
        ++ cc;
        
        // 14: Each L_e > 20 cm
        // (20) Track lengths each > 20 cm
        if((eTrackLength[0] <= 20. ) || (eTrackLength[1] <= 20.)) continue;
        ++ cut_counter[cc];
        ++ cc;
        
        // 15: <= 5 unassociated prompt Gg hits
        // (21) <= 5 NA prompt Gg hits
        if(NAPromptGgHits > 5) continue;
        ++ cut_counter[cc];
        ++ cc;

    */  


        // Fill histograms with relevant variables which we want to plot pre-cuts

        ///////////////////////////////////////////////////////////////////////
        /*
        Double_t NAPromptGgHitsDist2VertexMin = -1.;
        for(int i{0}; i < NAPromptGgHits; ++ i)
        {
            if(NAPromptGgHitsDist2VertexMin == -1.)
            {
                NAPromptGgHitsDist2VertexMin = NAPromptGgHitsDist2Vertex[i];
            }
            else if(NAPromptGgHitsDist2Vertex[i] < NAPromptGgHitsDist2VertexMin)
            {
                NAPromptGgHitsDist2VertexMin = NAPromptGgHitsDist2Vertex[i];
            }

            // TODO: changed to fill all values
            hNAPromptGgHitsDist2VertexMin->Fill(NAPromptGgHitsDist2Vertex[i], weight);
        }
        if(NAPromptGgHitsDist2VertexMin != -1.)
        {
            //hNAPromptGgHitsDist2VertexMin->Fill(NAPromptGgHitsDist2VertexMin, weight);
        }
        */
        // TODO: I disabled this because I couldn't be bothered to move it
        // to the end of the code to fill after all cuts with the rest of the
        // variables. Check Summers code to see how this was filled. May need
        // to move the Fill() function call with the for loop (most likely)
        // or just the Fill() by itself? (Unlikely)
        ///////////////////////////////////////////////////////////////////////



    /* WORKING HERE */

        // file: /unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/CutFlow.txt
        // number of events: 690220
        // last line in file, line 24
        // matches number 15 in table from summers thesis
        // this is the first cut which should have an effect...TODO check

        // 16: Elecron block types != L2 / L3
        // (22) Electron block types != L2 / L3
        // 5. Tracks extrapolated to the petal blocks in the calorimeter layer closest to the source foil
        // are not accepted as their energy calibration is less precise than the other block types.
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                
                //if( (electronBlockType[0] == 4) || (electronBlockType[1] == 4) ) continue; // not in petals
                //if( (electronBlockType[0] == 5) || (electronBlockType[1] == 5) ) continue; // not in petals
                //if( (electronBlockType[0] < 0) || (electronBlockType[1] < 0) ) continue; // does not do anything (makes no difference to cut count)
                //if( (electronBlockType[0] > 5) || (electronBlockType[1] > 5) ) continue; // cuts way to many events
                                                                                         // disabling cut on 4 doesn't give enough events back (cuts too many events in total)
                                                                                         // disabling cut on 5 gives too many events back (doesn't cut enough events in total)

                if( (electronBlockType[0] == 4) || (electronBlockType[1] == 4) ) continue; // not in petals
                if( (electronBlockType[0] == 5) || (electronBlockType[1] == 5) ) continue; // not in petals

                //if( (electronBlockType[0] > 2) || (electronBlockType[1] > 2) ) continue;
                
                //if( (electronBlockType[0] == 3) || (electronBlockType[1] == 3) ) continue;
                //if( (electronBlockType[0] == 4) || (electronBlockType[1] == 4) ) continue; // nope.

            }
        }
        ++ cut_counter[cc]; // cut 3 - block type
        ++ cc;
        
        // 17: No gamma rays in event
        //if(nGammaClusters > 0) continue;
        // Think this is the correct way to cut any gamma rays
        // Summer removed this for Position paper
        // replacing it with no gamma rays with energy > 200 keV
        // check this TODO
        //if(mode_flag == 0)
        //{
        //    if(mode_flag_2 == 1)
        //    {
        //        int cluster_hit_index = 0;
        //        int cesum = 0;
        //        for(int i{0}; i < nGammaClusters; ++ i)
        //        {
        //            /*
        //            if(clusterEnergy[i] > 0.2)
        //            {
        //                ++ cesum;
        //            }
        //            */
        //            for(int j{0}; j < nInCluster[i]; ++ j)
        //            {
        //                if(clusterHitEnergy[cluster_hit_index] > 0.2)
        //                {
        //                    ++ cesum;
        //                }
        //                ++ cluster_hit_index;
        //            }
        //        }
        //        if(cesum > 0) continue;
        //    }
        //}
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                double cluster_total_energy = 0.;
                int cluster_hit_index = 0;
                int cesum = 0;

                //if(nGammaClusters > 0) continue;

                for(int i{0}; i < nGammaClusters; ++ i)
                {
                    /*
                    if(clusterEnergy[i] > 0.2)
                    {
                        ++ cesum;
                    }
                    */
                    /*
                    for(int j{0}; j < nInCluster[i]; ++ j)
                    {
                        if(clusterHitEnergy[cluster_hit_index] > 0.2)
                        {
                            ++ cesum;
                            cluster_total_energy += clusterHitEnergy[cluster_hit_index];
                        }
                        ++ cluster_hit_index;
                    }
                    */
                    //if(clusterEnergy[i] > 0.2) continue;
                    cluster_total_energy += clusterEnergy[i];
                }
                //if(cesum > 0) continue;
                //if(cluster_total_energy > 0.2) continue;
                if(cluster_total_energy > 0.0) continue;
                
            }
        }
        ++ cut_counter[cc]; // cut 4 - gamma energy
        ++ cc;
        // Check that this does not cut any events - it shouldn't
        // TODO
        // Justin says there may be gamma rays in addition to 2e
        // but then what does 2eNg contain?
        // In her position paper analysis, she doesn’t require the no-γ cut, so I think we’re OK to not use it.
        // Instead, see table 6.1 on page 63 of her position paper, she just gets rid of gamma with E > 200 keV.

        // 18: Tracks have negative curvature
        // (24) Track have negative curvature
        // This appears to have been removed from Summers Position Paper
        // TODO CHECK THIS
        // how many events does this cut?
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1 || mode_flag_2 == 0)
            {
                if((trackSign[0] >= 0) || (trackSign[1] >= 0)) continue;
                // 2020-05-03 changed TODO change back?
                // result: makes NO DIFFERENCE ?!
                //if((trackSign[0] > 0) || (trackSign[1] > 0)) continue;
            }
        }
        ++ cut_counter[cc]; // cut 5 - track negative curvature
        ++ cc;
        
        // 19: Each L_e >= 30 cm
        // (25) Track lengths > 30 cm
        // 4. - position paper number/index, this index appears multiple times
        // 4. Each of the electron tracks must be at least 30 cm long, [have a hit in the first Gg cell layer,
        // and be associated to unique scintillator blocks.]
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                if((eTrackLength[0] < 30. ) || (eTrackLength[1] < 30.)) continue;
            }
        }
        ++ cut_counter[cc]; // cut 6 - track length
        ++ cc;

        ///////////////////////////////////////////////////////////////////////
        /*
        hElectronFirstGgMax->Fill(firstGgHitLayer[highE_index], weight);
        hElectronLastGgMax->Fill(lastGgHitLayer[highE_index], weight);
        hElectronFirstGgMin->Fill(firstGgHitLayer[lowE_index], weight);
        hElectronLastGgMin->Fill(lastGgHitLayer[lowE_index], weight);
        */
        // moved to fill with others after cuts
        ///////////////////////////////////////////////////////////////////////

        // 4. [Each of the electron tracks must be at least 30 cm long,] have a hit in the first Gg cell layer,
        // and be associated to unique scintillator blocks.
        //if((firstGgHitLayer[0] > 0) || (firstGgHitLayer[1] > 0)) continue;
        //++ cut_counter[cc]; // cut 7
        //++ cc;
        // TODO unique scintillator blocks
        // from below
        // 23: Both tracks have hit in Gg L0
        // (not 4, 18) - position paper number/index, this index appears multiple times
        // (30) Both tracks have hit in Gg L0
        //if(electronPMT[0] == electronPMT[1]) continue;
        //++ cut_counter[cc]; // cut 8
        //++ cc;
        // 27: Each electron track must be associated to a unique scinitllator block
        // (3) - position paper number/index, this index appears multiple times

        
        // 20: Tracks not from hot spots
        // (26) Tracks not from hot spots
        // 3. Electron tracks with vertices located in hot spot regions are not accepted in order to
        // improve the signal to background ratio in this channel.
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                //if(vertexInHotSpot[0] || vertexInHotSpot[1]) continue;

                double vertexZ_1 = vertexZ[0];
                double vertexZ_2 = vertexZ[1];
                double vertexSec_1 = vertexSec[0];
                double vertexSec_2 = vertexSec[1];

                bool flag_continue = false;

                for(int elix = 0; elix < N_HOTSPOT_ELIPSE; ++ elix)
                {

                    double x1 = vertexSec_1;
                    double y1 = vertexZ_1;
                    double h1 = hotspot_elipse[elix][0];
                    double k1 = hotspot_elipse[elix][1];
                    double a1 = hotspot_elipse[elix][2];
                    double b1 = hotspot_elipse[elix][3];

                    double p1 = (std::pow(x1 - h1, 2.0) / pow(a1, 2.0)) +
                               (std::pow(y1 - k1, 2.0) / pow(b1, 2.0));

                    if(p1 <= 1.0)
                    {
                        flag_continue = true;
                        break;
                    }

                    double x2 = vertexSec_2;
                    double y2 = vertexZ_2;
                    double h2 = hotspot_elipse[elix][0];
                    double k2 = hotspot_elipse[elix][1];
                    double a2 = hotspot_elipse[elix][2];
                    double b2 = hotspot_elipse[elix][3];

                    double p2 = (std::pow(x2 - h2, 2.0) / pow(a2, 2.0)) +
                               (std::pow(y2 - k2, 2.0) / pow(b2, 2.0));

                    if(p2 <= 1.0)
                    {
                        flag_continue = true;
                        break;
                    }
                               
                }

                if(flag_continue == true) continue;

            }
        }
        ++ cut_counter[cc]; // cut 7 (was cut 9) - track hot spot
        ++ cc;
        
        // 10. No unassociated Geiger hits are allowed within 15 cm of the calorimeter hits as this is
        // indicative of a reflected electron, which therefore has not deposited its full energy in the
        // calorimeter.
        // (27) No NA prompt Gg hits within 15 cm of scintillator hits
        // This does not appear to be included in the Table from Summers thesis
        // TODO
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                int ssum = 0;
                for(int i{0}; i < NAPromptGgHits; ++ i)
                {
                    if(NAPromptGgHitsDist2Calo[i] < 15.)
                    {
                        ++ ssum;
                    }
                }
                if(ssum > 0) continue;
            }
        }
        ++ cut_counter[cc]; // cut 8 - Gg hit scint dist
        ++ cc;

        // 21: <= 1 unassociated, prompt Gg hit within 15 cm in XY coordinate of event vertex
        // (28) <= 1 NA prompt Gg hit within 15 cm of vertex
        // 9. [To further improve this rejection of 214Bi, no more than one prompt Geiger hit which is
        // not associated to a track is allowed within 15 cm of the track vertices.] If the tracks on the
        // same side of the foil, no unassociated prompt Geiger hits are allowed on the other side of
        // the foil.
        // This allows for some level of Geiger cell noise which is not simulated, while rejecting a large
        // number of 214Bi decays.
        // NAPromptGgHitsDist2Vertex is distance between the Geiger hits and the average
        // electron vertex in the X,Y plane
        // The average electron vertex is obtained by adding the 2 electron vertices and
        // dividing by 2
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                int ssum = 0;
                for(int i{0}; i < NAPromptGgHits; ++ i)
                {
                    // assuming units of cm TODO
                    if(NAPromptGgHitsDist2Vertex[i] < 15.)
                    {
                        ++ ssum;
                    }
                }
                if(ssum > 1) continue;
            }
        }
        ++ cut_counter[cc];
        ++ cc; // cut 9 (was cut 10) - Gg hit vertex dist
        
        // 22: No unassociated, prompt Gg hits on opposite side of foil if tracks are on same side
        // 9. To further improve this rejection of 214Bi, no more than one prompt Geiger hit which is
        // not associated to a track is allowed within 15 cm of the track vertices. [If the tracks on the
        // same side of the foil, no unassociated prompt Geiger hits are allowed on the other side of
        // the foil.]
        // This allows for some level of Geiger cell noise which is not simulated, while rejecting a large
        // number of 214Bi decays.
        // (29) No NA prompt Gg hits on opposite side of foil if tracks are on same side
        // TODO: not sure if this is right
        // this is clearly not right
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                int ssum = 0;
                if(electronSide[0] == electronSide[1])
                {
                    // tracks on same side
                    //if((electronFirstGG[0] > 1) && (electronFirstGG[1] > 1)) cut = true;
                    for(int i{0}; i < NAPromptGgHits; ++ i)
                    {
                        if(electronSide[0] != NAPromptGgHitsSide[i])
                        {
                            ++ ssum;
                        }
                    }
                }
                else if(electronSide[0] != electronSide[1])
                {
                    // tracks on opposite side
                    // do nothing in this case
                }
                if(ssum != 0) continue;
            }
        }
        ++ cut_counter[cc]; // cut 10 (was cut 11) - Gg hit opposite side
        ++ cc;


        
        
        // 23: Both tracks have hit in Gg L0
        // TODO: copied from above, but this looks like something which should
        // go here - might need to check against layer 1 and 11 ?
        //if((electronFirstGG[0] > 1) || (electronFirstGG[1] > 1)) cut = true; TODO: changed from this?
        // (not 4, 18) - position paper number/index, this index appears multiple times
        // (30) Both tracks have hit in Gg L0
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                if((firstGgHitLayer[0] > 0) || (firstGgHitLayer[1] > 0)) continue;
            }
        }
        ++ cut_counter[cc]; // cut 11 - Gg L0
        ++ cc;
        // moved above

        // more position paper stuff TODO
        // 6. The impact of the tracks must be on the front face of the scintillator block for proper light
        // collection.
        // 8. No alpha candidates are allowed within 15 cm of the electron vertices in order to reduce the
        // event rate coming from 214Bi decays.
        // 11. Additional, isolated calorimeter hits which do not have a track associated to them are
        // allowed if the single hit energies are less than 0.2 MeV (indicative of noise which is not
        // simulated).


        ///////////////////////////////////////////////////////////////////////
        // TODO: why are these being filled with different weight
        // is it because these are not used in the fit?
        /*
        hInternalPullee->Fill(eeInternalPull, 1.0);
        hInternalProbee->Fill(eeInternalProb, 1.0);
        hExternalPullee->Fill(eeExternalPull, 1.0);
        hExternalProbee->Fill(eeExternalProb, 1.0);
        */
        // moved to fill after all cuts
        ///////////////////////////////////////////////////////////////////////

        // 24: P_{int} >= 0.01 and P_{ext} <= 0.01
        // (31) P_int > 0.01 and P_ext < 0.01
        // 12. The internal probability between the two electrons must be at least 1% and the external
        // probability must be no more than 1%
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1 || mode_flag_2 == 0)
            {
                if(eeInternalProb < 0.01) continue;
                if(eeExternalProb > 0.01) continue;
            }
        }
        ++ cut_counter[cc]; // cut 12 (was cut 13) - prob
        ++ cc;
        

        double rPhi0 = (vertexR[0] * vertexSec[0] * TMath::TwoPi() / 20.);
        double rPhi1 = (vertexR[1] * vertexSec[1] * TMath::TwoPi() / 20.);
        ///////////////////////////////////////////////////////////////////////
        /*
        hVertexDZ->Fill(vertexZ[highE_index] - vertexZ[lowE_index], weight);
        hVertexDR->Fill(vertexR[highE_index] - vertexR[lowE_index], weight);
        hVertexDRPhi->Fill(rPhi0 - rPhi1, weight);
        */
        // moved to fill after all cuts
        //TODO DSec
        ///////////////////////////////////////////////////////////////////////

        // 25: Delta R <= 4cm and Delta Z <= 8 cm
        // (32) Delta R <= 4 cm and Delta Z <= 8 cm
        // 2. The separation between each individually reconstructed vertex must be ≤ 4 cm radially
        // and ≤ 8 cm longitudinally to ensure that the electrons originate from a common vertex.
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                if(std::abs(vertexR[0] - vertexR[1]) > 4.) continue;
                if(std::abs(vertexZ[0] - vertexZ[1]) > 8.) continue;

                //if(std::abs(vertexR[0] - vertexR[1]) > 4.0e-2) continue;
                //if(std::abs(vertexZ[0] - vertexZ[1]) > 8.0e-2) continue;
            }
        }
        ++ cut_counter[cc]; // cut 13 (was cut 14) - delta
        ++ cc;
        
        // save data to file if electronEnergy >= 0.2 MeV
        // and phase is P1 or P2
        if(mode_flag == 0)
        {
            if((electronEnergy[0] < 0.2) || (electronEnergy[1] < 0.2))
            {
                // do not save to file
            }
            else
            {
                if(thePhase == 0)
                {
                    //if((1869 <= run) && (run <= 3395)) cut = false;
                    
                    #if TRUTH_ENABLE
                    //outputTree->Fill();
                    outputTree_small->Fill();
                    #endif
                }
                else if(thePhase == 1)
                {
                    //if((3396 <= run) && (run <= 9186)) cut = false;

                    #if TRUTH_ENABLE
                    //outputTree->Fill();
                    outputTree_small->Fill();
                    #endif
                }
            }
        }



        // 26: Each E_e >= 300 keV
        // (33) Each E_e >= 300 keV
        // 7. The energy of each electron is required to be at least 0.3 MeV, as single electron energies
        // below this threshold are not well described by the simulation.
        if(mode_flag == 0)
        {
            if(mode_flag_2 == 1)
            {
                if((electronEnergy[0] < 0.3) || (electronEnergy[1] < 0.3)) continue;
            }
        }
        ++ cut_counter[cc]; // cut 14 (was cut 15) - energy
        ++ cc;

        // 27: Each electron track must be associated to a unique scinitllator block
        // (3) - position paper number/index, this index appears multiple times
        //if(electronPMT[0] == electronPMT[1]) continue;
        //++ cut_counter[cc];
        //++ cc;
        // moved above
        

        // Phase cut moved to be final cut
        // 1869 - 3395 is P1
        // 3396 - 9186 is P2
        if(mode_flag == 0)
        {
            cut = true;
            if(thePhase == 0)
            {
                if((1869 <= run) && (run <= 3395)) cut = false;
            }
            else if(thePhase == 1)
            {
                if((3396 <= run) && (run <= 9186)) cut = false;
            }
            else if(nocutonphase == 1)
            {
                cut = false;
            }
            else
            {
                cut = true;
            }

            // accept P1 and P2 for testing purposes
            //cut = false;
            if(cut == true) continue;
        }
        ++ cut_counter[cc]; // cut 1 - phase
        ++ cc;


        //std::cout << "Total number of cuts recorded in cut_counter: " << cc << std::endl;
/**/

        //for(int i = 0; i < nUnAssocPromptGg; i++)
        /*
        if(thePath.CompareTo("nd150/") == 0)
        {
            if(NAPromptGgHits < 1)
            {
                std::cout << "THERE ARE NO PROMPT GG HITS FOR THIS SAMPLE" << std::endl;
            }
        }
        */
        // for Nd150, there do not appear to be any unassociated prompt geiger
        // hits
        for(int i = 0; i < NAPromptGgHits; i++)
        {
/*
            if(highE_index == 0)
            {
                hVertexMinDistPromptGgMax->Fill(dist2VertexPromptGg0[i]);
                hVertexMinDistPromptGgMin->Fill(dist2VertexPromptGg1[i]);
            }
            else
            {
                hVertexMinDistPromptGgMax->Fill(dist2VertexPromptGg1[i]);
                hVertexMinDistPromptGgMin->Fill(dist2VertexPromptGg0[i]);
            }
*/

            /*
            if(thePath.CompareTo("nd150/") == 0)
            {
                std::cout << "fill, value = " << NAPromptGgHitsDist2Vertex[i] << std::endl;
            }
            hVertexMinDistPromptGgMax->Fill(NAPromptGgHitsDist2Vertex[i]);
            */

        }
/*


*/
        //hNAPromptGgHits->Fill(NAPromptGgHits);
        // moved below
/*
        double totalLowE = 0.;
        for(int i = 0; i < nLowEClusters; i++)
        {
            totalLowE += lowEClusterE[i];
            hLowEGammaEnergy->Fill(lowEClusterE[i], weight);
            hLowEMinDistPromptGg->Fill(lowEMinDistPromptGg[i], weight);
        }
        hNLowEGammas->Fill(nLowEClusters, weight);
        hSummedLowEGammaE->Fill(totalLowE, weight);

        //if ( nLowEClusters > 0 )continue;
*/

        ++ event_pass_count;

        // moved back here from above
        hRun->Fill(run, weight);
        hNElectrons->Fill(nElectrons, weight);

        //#if TRUTH_ENABLE
        //if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
        //   (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
        //{
        //    hMCElectronEnergy->Fill(electronEnergy[0], electronEnergy[1]);
        //    hMCTrueElectronEnergy->Fill(trueElectronEnergy[0], trueElectronEnergy[1]);
        //}
        //#endif

        // moved from above
        hElectronFirstGgMax->Fill(firstGgHitLayer[highE_index], weight);
        hElectronLastGgMax->Fill(lastGgHitLayer[highE_index], weight);
        hElectronFirstGgMin->Fill(firstGgHitLayer[lowE_index], weight);
        hElectronLastGgMin->Fill(lastGgHitLayer[lowE_index], weight);

        hInternalPullee->Fill(eeInternalPull, 1.0);
        hInternalProbee->Fill(eeInternalProb, 1.0);
        hExternalPullee->Fill(eeExternalPull, 1.0);
        hExternalProbee->Fill(eeExternalProb, 1.0);

        hVertexDZ->Fill(vertexZ[highE_index] - vertexZ[lowE_index], weight);
        hVertexDR->Fill(vertexR[highE_index] - vertexR[lowE_index], weight);
        hVertexDRPhi->Fill(rPhi0 - rPhi1, weight);

        hNAPromptGgHits->Fill(NAPromptGgHits);
        // these were moved from above


        hTotalE->Fill(electronEnergy[0] + electronEnergy[1], weight);
        hCosee->Fill(cosee, weight);
        hCoseeWeighted->Fill(cosee, weight * cosee_weight);

        hEeMax->Fill(electronEnergy[highE_index], weight);
        hElectronLengthMax->Fill(eTrackLength[highE_index], weight);
        hVertexSectorMax->Fill(vertexSec[highE_index], weight);
        hVertexZMax->Fill(vertexZ[highE_index], weight);
        hVertexRMax->Fill(vertexR[highE_index], weight);
        hVertexZSecMax->Fill(vertexSec[highE_index], vertexZ[highE_index], weight);
        hElectronLDCorrMax->Fill(electronLDCorr[highE_index], weight);

        hEeMin->Fill(electronEnergy[lowE_index], weight);
        hElectronLengthMin->Fill(eTrackLength[lowE_index], weight);
        hVertexSectorMin->Fill(vertexSec[lowE_index], weight);
        hVertexZMin->Fill(vertexZ[lowE_index], weight);
        hVertexRMin->Fill(vertexR[lowE_index], weight);
        hVertexZSecMin->Fill(vertexSec[lowE_index], vertexZ[lowE_index], weight);
        hElectronLDCorrMin->Fill(electronLDCorr[lowE_index], weight);

        hEeMaxVEeMin->Fill(electronEnergy[highE_index], electronEnergy[lowE_index], weight);

        hTrackSignMax->Fill(trackSign[highE_index], weight);
        hTrackSignMin->Fill(trackSign[lowE_index], weight);


        hnGammaClusters->Fill(nGammaClusters, weight);

        hFoilSide->Fill(foilSide, weight);
        hTrueVertexR->Fill(trueVertexR, weight);
        hFoilSideTrueVertexR->Fill(foilSide, trueVertexR, weight);
        hTrueVertexLayer->Fill(trueVertexLayer, weight);
        
        hSingleEnergy->Fill(electronEnergy[lowE_index], weight);
        hSingleEnergy->Fill(electronEnergy[highE_index], weight);
        hHighEnergy->Fill(electronEnergy[highE_index], weight);
        hLowEnergy->Fill(electronEnergy[lowE_index], weight);
        hHighLowEnergy->Fill(electronEnergy[lowE_index], electronEnergy[highE_index], weight);
        hEnergySum->Fill(electronEnergy[highE_index] + electronEnergy[lowE_index], weight);
        hEnergyDiff->Fill(electronEnergy[highE_index] - electronEnergy[lowE_index], weight);


        int hit_counter = 0;
        int nLowEnergyHits = 0;
        double clusterHitEnergyMin;
        double clusterHitEnergyMax;
        double total_cluster_energy;
        for(int i = 0; i < nGammaClusters; ++ i)
        {
            hnInCluster->Fill(nInCluster[i], weight);
            hclusterEnergy->Fill(clusterEnergy[i], weight);

            clusterHitEnergyMin = 0.;
            clusterHitEnergyMax = 0.;
            total_cluster_energy = 0.0;
            for(int j = 0; j < nInCluster[i]; ++ j)
            {
                total_cluster_energy += clusterHitEnergy[hit_counter];

                if(j == 0)
                {
                    clusterHitEnergyMin = clusterHitEnergy[hit_counter];
                    clusterHitEnergyMax = clusterHitEnergy[hit_counter];
                }

                hclusterHitEnergy->Fill(clusterHitEnergy[hit_counter], weight);
                ++ hit_counter;

                if(clusterHitEnergy[hit_counter] < .2)
                {
                    ++ nLowEnergyHits;
                }
            }
            hTotalClusterEnergy->Fill(total_cluster_energy, weight);
            if(nInCluster[i] >= 2)
            {
                // don't store case where no data was set (0)
                // or case where Min = Max (1)
                hclusterHitEnergyMin->Fill(clusterHitEnergyMin, weight);
                hclusterHitEnergyMax->Fill(clusterHitEnergyMax, weight);
            }
        }
        hnLowEnergyHits->Fill(nLowEnergyHits, weight);
        if(hit_counter != nTotalClusterHits)
        {
            std::cout << "hit_counter=" << hit_counter << ", nTotalClusterHits=" << nTotalClusterHits << std::endl;
            std::cin.get();
        }



    } //~event_i
     

    if(mode_flag == 0)
    {
        //std::cout << hSingleEnergy->GetName() << " -> " << hSingleEnergy->GetBinContent(10) << std::endl;
    }


    if(mode_flag == 0)
    {
        ofile_cutcount  << sampleName << "    " << "" << total_event_count << "    " << "" << event_pass_count << std::endl;
        std::cout       << sampleName << "    " << "total_event_count=" << total_event_count << "    " << "event_pass_count=" << event_pass_count << std::endl;
    }
    //std::cout << "electron_cut_A=" << electron_cut_A << std::endl;
    std::cout << std::endl;

    #if TRUTH_ENABLE
    if(mode_flag == 0)
    {
        //if((sampleName.CompareTo("nd150_rot_2b2n_m4") == 0)   ||
        //   (sampleName.CompareTo("nd150_rot_2n2b_m4") == 0))
        //{
        //outputTree->Write();
        //outputFile->Close();

        outputTree_small->Write();
        outputFile_small->Close();
        //}
    }
    #endif

    aFile->Close();
    
    std::string cut_description[nCuts];
    int cc = 0;
    cut_description[cc] = "Input count";
    ++ cc;
    cut_description[cc] = "IN/OUT/Layer Sample Split";
    ++ cc;
    cut_description[cc] = "Vertex sector";
    ++ cc;
    cut_description[cc] = "Block type";
    ++ cc;
    cut_description[cc] = "No gamma with E > 0.2 MeV";
    ++ cc;
    cut_description[cc] = "Tracks negative curvature";
    ++ cc;
    cut_description[cc] = "Track length > 30 cm";
    ++ cc;
    cut_description[cc] = "Tracks not from hot spots";
    ++ cc;
    cut_description[cc] = "0 Prompt Gg hit within 15 cm scint";
    ++ cc;
    cut_description[cc] = "1 Prompt Gg hit within 15 cm vertex";
    ++ cc;
    cut_description[cc] = "Unassoc prompt Gg on opposite side of foil";
    ++ cc;
    cut_description[cc] = "Tracks hits in Gg L0";
    ++ cc;
    cut_description[cc] = "Pint / Pext";
    ++ cc;
    cut_description[cc] = "Delta R / Delta Z";
    ++ cc;
    cut_description[cc] = "Energy 300 keV";
    ++ cc;
    cut_description[cc] = "Phase";
    ++ cc;
    cut_description[cc] = "none";

    //std::cout << "Here are the cut counts" << std::endl;
    std::string of_cutcount_filename;
    if(mode_flag == 0)
    {
        of_cutcount_filename = "of_cutcount_processeddata_" + std::string(sampleName) + ".txt";
    }
    else if(mode_flag == 1)
    {
        of_cutcount_filename = "of_cutcount_rawdata_" + std::string(sampleName) + ".txt";
    }
    std::ofstream of_cutcount(of_cutcount_filename, std::ios::out);
    of_cutcount << sampleName << std::endl;
    for(int i = 0; i < nCuts; i++)
    {
        //std::cout << "Passed cut " << i << ": " << cut_counter[i] << " events. Description=\"" << cut_description[i] << "\"" << std::endl;
        of_cutcount << "Passed cut " << i << ": " << cut_counter[i] << " events. Description=\"" << cut_description[i] << "\"" << std::endl; 
    }
    of_cutcount << std::endl;
    

    for(int i = 0; i < nCuts; i++ )
    {
        // NOTE: bug, should be i + 1, bin 0 is underflow, bin -1 does not exist
        //hNAfterCuts->SetBinContent(i - 1, cut_counter[i]);
        hNAfterCuts->SetBinContent(i + 1, cut_counter[i]);
    }


    TFile *myFile = new TFile("Nd150_2e_P" + Phase + ".root", "UPDATE");
    TDirectory *dir = gDirectory;
    //dir->cd("singleHistos/unscaled");
    //dir->cd("rawdata");
    //dir->mkdir(sampleName);
    //dir->cd(sampleName);
    // TODO: don't want dirs organized by sample name, want them organized
    // by histogram name

    //dir->cd("rawdata/" + sampleName);

    //dir->pwd();
    //dir->print();

    //std::cout << "writing histograms into \"singleHistos/unscaled\'" << std::endl;
    //std::cout << "writing: singleHistos/unscaled/" << hTotalE->GetName() << std::endl;
    //std::cout << "writing histograms into \"rawdata" << "\" or \"processeddata\"" << std::endl;
    //std::cout << "writing: rawdata or processeddata ... " << hTotalE->GetName() << std::endl;

    //for(int i = 0; i < numHistograms; ++ i)
    for(std::map<TString, TH1*>::iterator it{hpmap.begin()}; it != hpmap.end(); ++ it)
    {
        //std::cout << "i=" << i << std::endl;

        //std::string histogram_name = std::string(histogramNames[i].Data());
        //TString histogram_name = histogramNames[i];
        TString histogram_name = it->first;
        //std::cout << "histogram_name=" << histogram_name << std::endl;
        
        // only create the directories once
        // this functions called with DataFile first
        // bit of a dodgy hack
        if(sampleName.CompareTo(DataFile) == 0)
        {
            if(mode_flag == 0)
            {
                TDirectory *dir_histogram_name = dir->mkdir("processeddata/" + histogram_name);
                TDirectory *dir_histogram_name_2 = dir->mkdir("scaled/" + histogram_name);
//                TDirectory *dir_histogram_name_3 = dir->mkdir("scaled2/" + histogram_name);
            }
            else if(mode_flag == 1)
            {
//                TDirectory *dir_histogram_name = dir->mkdir("rawdata/" + histogram_name);
            }
        }

        if(mode_flag == 0)
        {
            myFile->cd("processeddata/" + histogram_name);
            it->second->Write();
        }
        else if(mode_flag == 1)
        {
            myFile->cd("rawdata/" + histogram_name);
// 2020-12-11 disabled            it->second->Write();
        }
        //myFile->pwd();
        
        //hpmap.at(histogram_name)->Write();
        // same as
        
    }


    myFile->Close();

}

#endif // FIT_2E_MAKEHISTOGRAMS_H
