#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>


#include "../include/InputFileDef.h"


TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
//TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/";




void splitsamples()
{

    // map of flags, used for selecting events which are split by
    // L0, L>0, IN, OUT
    // 0 = default
    // 1 = select IN events
    // 2 = select OUT events
    // 4 = select L0 IN events
    // 8 = select L0 OUT events
    // 16 = select L>0 (IN+OUT) events
    std::map<std::string, uint64_t> sample_split_flags;
    sample_split_flags["bi214_sfoil"] = (uint64_t)1 | (uint64_t)2;
    sample_split_flags["pb214_sfoil"] = (uint64_t)1 | (uint64_t)2;
    sample_split_flags["bi214_swire"] = (uint64_t)4 | (uint64_t)8 | (uint64_t)16;
    sample_split_flags["pb214_swire"] = (uint64_t)4 | (uint64_t)8 | (uint64_t)16;

    
    std::cout << ">>>>> Reading in Rn222 backgrounds to split into output files" << std::endl;
    for(int i = 0; i < nRn222BkgsInput; i++)
    {
        TString sample_name = Rn222BkgFilesInput[i];
        std::cout << "sample_name=" << sample_name << std::endl;
        std::string std_string_sample_name = std::string(sample_name.Data());
        if(sample_split_flags.count(std_string_sample_name) > 0)
        {
            //std::cout << "catch: flags: " << std_string_sample_name << std::endl;
            //std::cin.get();

            uint64_t flags = sample_split_flags.at(std_string_sample_name);
            //std::cout << "flags=" << flags << std::endl;
            int counter = 0;
            //               |  ||  ||  ||  |
            while(flags & (0xFFFFFFFFFFFFFFFF << counter))
            {
                //std::cout << "while" << std::endl;

                //uint64_t flag_bit = flags & (uint64_t)0x01;
                //                             |  ||  ||  ||  |
                uint64_t flag_bit = flags & (0x0000000000000001 << counter);
                //std::cout << "flag_bit=" << flag_bit << std::endl;
                if(flag_bit)
                {
                    //uint64_t f = flag_bit >> counter;
                    //std::cout << "calling makeHistograms with flag_bit=" << flag_bit << std::endl;
                    //std::cin.get();
                    //makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, 0, flag_bit);
                    //makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, 1, flag_bit);
                    TString thePath = "externals/";
                    TString sampleName = Rn222BkgFilesInput[i];
                    const Int_t mode_flag = 0;
                    const uint64_t sample_split_flags = flag_bit;


                    TString name_sample_split_additional = "";
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

                    
                    TFile *aFile;
                    aFile = TFile::Open(filePath + thePath + sampleName + "/Nd150_2eNg_output.root");
                    TTree *theTree = (TTree*)aFile->Get("Nd150_2eNg/Nd150_2eNg");

                    TString sampleName_output = sampleName + name_sample_split_additional;
                    std::cout << "sampleName_output=" << sampleName_output << std::endl;

                    std::cout << "writing to file: " << filePath + thePath + sampleName_output + "/Nd150_2eNg_output.root" << std::endl;
                    //std::cin.get();

                    TFile *outputFile = nullptr;
                    TTree *outputTree = nullptr;
                    outputFile = new TFile(filePath + thePath + sampleName_output + "/Nd150_2eNg_output.root", "recreate");
                    outputFile->mkdir("Nd150_2eNg");                                                                           
                    outputFile->cd("Nd150_2eNg");                                                  
                    outputTree = new TTree("Nd150_2eNg", "Nd150_2eNg");


                    ///////////////////////////////////////////////////////////
                    //
                    ///////////////////////////////////////////////////////////

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


                    int         nTotalClusterHits; // TODO - check max value is 10 ?
                    // TODO: MAX VALUE IS 20 NOT 10!
                    double      clusterHitEnergy[20];
                    int         clusterHitPMT[20];
                    int         clusterHitLDFlag[20];
                    double      clusterHitLDCorr[20];
                    double      clusterHitLDCorrErr[20];
                    double      clusterHitSec[20];
                    double      clusterHitZ[20];



                    ///////////////////////////////////////////////////////////
                    //
                    ///////////////////////////////////////////////////////////

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
                    theTree->SetBranchAddress("nInCluster"                  , nInCluster);
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


                    ///////////////////////////////////////////////////////////
                    //
                    ///////////////////////////////////////////////////////////

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
                    outputTree->Branch("nInCluster"                 , nInCluster                    , "nInCluster[nGammaClusters]/I");
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

                    Long_t fill_count = 0;

                    Long_t events = (Long_t)theTree->GetEntries();
                    std::cout << "events=" << events << std::endl;
                    for(Long_t event_i = 0; event_i < events; event_i++)
                    {

                        theTree->GetEvent(event_i);

                        /*
                        if(sample_split_flags > 2)
                        {
                            std::cout << "foilSide=" << foilSide << std::endl;
                            std::cout << "trueVertexLayer=" << trueVertexLayer << std::endl;
                        }
                        */

                        // additional sample splitting
                        // side = 1: OUT
                        // side = 0: IN
                        // side = -1: FAIL
                        // layer: 0 - 8 inclusive = valid
                        // -1 = fail
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
                            if((trueVertexLayer < 1) || (trueVertexLayer > 8))
                            {
                                continue;
                            }
                        }


                        /*
                        if(sample_split_flags > 2)
                        {
                            std::cin.get();
                        }
                        */

                        ++ fill_count;
                        outputTree->Fill();
                    }

                    std::cout << "fill_count=" << fill_count << " / " << theTree->GetEntries() << std::endl;

                    outputTree->Write();
                    outputFile->Close();

                    aFile->Close();

                } // if flagbit

                //flags_shift << 1;
                ++ counter;

                //std::cout << std::hex;
                //std::cout << "0xFFFFFFFFFFFFFFFF << counter = " << (0xFFFFFFFFFFFFFFFF << counter) << " flags=" << flags << std::endl;
                //std::cout << std::dec;

            }

        }
        else
        {
            // do nothing
        }

        std::cout << std::endl;

    }


    


}
