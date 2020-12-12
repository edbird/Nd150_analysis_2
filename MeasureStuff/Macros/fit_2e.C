#include "TH1D.h"
#include "TFractionFitter.h"
#include "TArray.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"

#include "THStack.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include <string>
#include "TLegend.h"
#include <vector>
#include <exception>

#include "../include/fit_2e.h"

#include "fit_2e_scale.h"
#include "fit_2e_makeHistograms.h"

// TODO:
//
// - change to double precision floating point

//-----------------------------------------
//    Functions
//----------------------------------------
void loadFiles();

//void makeHistograms(TString thePath, TString sampleName, std::ofstream &ofile_cutcount, const Int_t mode_flag, const uint64_t sample_split_flags = 0);
//void  makeHistograms(TTree *tree, TString sampleName, Int_t isotope);
void fitHistograms();
//void drawPlots(TH1D *data,TH1D *mc, THStack *hs);
double getActErr(double Npass, double Ngen, double Ndata, double sf_err );
double getChi2(TH1D *data, TH1D* mc, Int_t &ndof);

void loadFiles()
{

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    //First, create a root file to hold all of the histograms
    TFile *myFile = new TFile("Nd150_2e_P" + Phase + ".root", "RECREATE");
    TDirectory *dir = gDirectory;
    //TDirectory *dir_sh = dir->mkdir("singleHistos"); 
    //dir_sh->mkdir("unscaled");
    
    // raw data, no cuts
//TDirectory *dir_rawdata = dir->mkdir("rawdata");
    // processed data with cuts
    TDirectory *dir_processeddata = dir->mkdir("processeddata");
    // processed data with cuts and scaling
    TDirectory *dir_scaled = dir->mkdir("scaled");
//    TDirectory *dir_scaled2 = dir->mkdir("scaled2");

    
    myFile->Close();


    // map of flags, used for selecting events which are split by
    // L0, L>0, IN, OUT
    // 0 = default
    // 1 = select IN events
    // 2 = select OUT events
    // 4 = select L0 IN events
    // 8 = select L0 OUT events
    // 16 = select L>0 (IN+OUT) events
    /*
    std::map<std::string, uint64_t> sample_split_flags;
    sample_split_flags["bi214_sfoil"] = (uint64_t)1 | (uint64_t)2;
    sample_split_flags["pb214_sfoil"] = (uint64_t)1 | (uint64_t)2;
    sample_split_flags["bi214_swire"] = (uint64_t)4 | (uint64_t)8 | (uint64_t)16;
    sample_split_flags["pb214_swire"] = (uint64_t)4 | (uint64_t)8 | (uint64_t)16;
    */


    // open and reset file for storing the number of events that pass cuts
    std::ofstream ofile_cutcount("cutcount.txt", std::ofstream::out | std::ofstream::trunc);
    // TODO: need to count the number of events that pass cuts and store
    // into variable ExternalBkgPassCut[thePhase][i]
    // cannot pass that variable here because this would require
    // loadFiles() to be called before fitHistograms()
    // which is slow
    // instead, store results into a file

    // Read in all files and make the cuts you want.
    //First we will read in the data
    std::cout << ">>>>> Reading in data file...";
    makeHistograms("betabeta/", DataFile, ofile_cutcount, 0);
//    makeHistograms("betabeta/", DataFile, ofile_cutcount, 1);
    std::cout << "... DONE." << std::endl;

    // Read in all the backgrounds: externals, Rn222, Rn220, internals and neighbouring foils
    std::cout << std::endl;
    std::cout << ">>>>> Reading in external backgrounds..." << std::endl;
    for(int i = 0; i < nExternalBkgs; i++)
    {
        makeHistograms("externals/", ExternalBkgFiles[i], ofile_cutcount, 0);
//        makeHistograms("externals/", ExternalBkgFiles[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;


    std::cout << std::endl;
    std::cout << ">>>>> Reading in Rn220 backgrounds..." << std::endl;
    for(int i = 0; i < nRn220Bkgs; i++)
    {
        makeHistograms("externals/", Rn220BkgFiles[i], ofile_cutcount, 0);
//        makeHistograms("externals/", Rn220BkgFiles[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;

    #if 0
    std::cout << std::endl;
    std::cout << ">>>>> Reading in Rn222 backgrounds..." << std::endl;
    for(int i = 0; i < nRn222BkgsInput; i++)
    {
        TString sample_name = Rn222BkgFilesInput[i];
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
                    makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, 0, flag_bit);
                    makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, 1, flag_bit);
                }

                //flags_shift << 1;
                ++ counter;

                //std::cout << std::hex;
                //std::cout << "0xFFFFFFFFFFFFFFFF << counter = " << (0xFFFFFFFFFFFFFFFF << counter) << " flags=" << flags << std::endl;
                //std::cout << std::dec;
            }
        }
        else
        {
            makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, 0);
            makeHistograms("externals/", Rn222BkgFilesInput[i], ofile_cutcount, 1);
        }

    /*
        if(sample_name == "bi214_sfoil" ||
           sample_name == "pb214_sfoil")
        {
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 0, sample_split_flags[sample_name]);
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 1, sample_split_flags[sample_name]);
        }
        else if(sample_name == "bi214_swire" ||
                sample_name == "pb214_swire")
        {
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 0, sample_split_flags[sample_name]);
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 1, sample_split_flags[sample_name]);
        }
        else
        {
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 0);
            makeHistograms("externals/", Rn222BkgFiles[i], ofile_cutcount, 1);
        }
    */
    }
    std::cout << "... DONE." << std::endl;
    #endif




    std::cout << std::endl;
    std::cout << ">>>>> Reading in Rn222 backgrounds..." << std::endl;
    for(int i = 0; i < nRn222Bkgs; i++)
    {
        makeHistograms("externals/", Rn222BkgFilesNew[i], ofile_cutcount, 0);
//        makeHistograms("externals/", Rn222BkgFilesNew[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;




    std::cout << std::endl;
    std::cout << ">>>>> Reading in internal backgrounds..." << std::endl;
    for(int i = 0; i < nInternalBkgs; i++)
    {
        makeHistograms("internals/", InternalBkgFiles[i], ofile_cutcount, 0);
//        makeHistograms("internals/", InternalBkgFiles[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;


    std::cout << std::endl;
    std::cout << ">>>>> Reading in neighbouring foil backgrounds..." << std::endl;
    for(int i = 0; i < nNeighbours; i++)
    {
        makeHistograms("neighbours/", NeighbourFiles[i], ofile_cutcount, 0);
//        makeHistograms("neighbours/", NeighbourFiles[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;


    std::cout << std::endl;
    std::cout << ">>>>> Reading in Nd150 backgrounds..." << std::endl;
    for(int i = 0; i < nNd150Samples; i++)
    {
        makeHistograms("nd150/", Nd150Files[i], ofile_cutcount, 0);
//        makeHistograms("nd150/", Nd150Files[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;

    std::cout << std::endl;
    std::cout << ">>>>> Reading in Nd150 systematics..." << std::endl;
    for(int i = 0; i < nNd150Systematics; i++)
    {
        makeHistograms("nd150/", Nd150Systematics[i], ofile_cutcount, 0);
        //makeHistograms("nd150/", Nd150Systematics[i], ofile_cutcount, 1);
    }
    std::cout << "... DONE." << std::endl;
    


    ofile_cutcount.close();

}





// allSamples are individual histograms, index order: histogram type, background type
// index order: array index ([]), followed by ->At() function call index for TObjArray
// hAllSamples is a single histogram, index order: histogram type
// index order: array index ([])
// hAllMc is a single histogram, index order: histogram type
// index order: array index ([])
// same as hAllSamples, but contains all the different types of background/mc/data in a single
// histogram rather than splitting between several different objects (arrays of pointers) one
// for each type of background/mc/data
void stackfunction(int i,                           // INPUT: histogram index, should possibly be moved inside function (can be done by moving calls to new for hMinorStacks inside, and calling when allMC_first == true?)
                   int nSamples,                    // INPUT: number of samples in list (number of objects in lists below)
                   THStack **hMinorStacks,          // OUTPUT: ?
                   TObjArray **allSamples,          // INPUT: all samples
                   THStack **hMajorStacks,          // OUTPUT: ?
                   const TString& sampleHumanName,  // INPUT:
                   const TString& sampleName,       // INPUT: "Nd150", "Externals", "Internals", "Radon", used to when calling Clone() to change histogram name
                   TH1D **hAllSamples,              // OUTPUT: single histogram stack composed of samples components added together by function call to TH1::Add()
                   const TString *const sampleNames,
                   const TString *const sampleFiles,
                   TH1D **hAllMC,
                   bool allMC_first,
                   bool hAllSamples_addtoexisting,  // INPUT: flag, if true objects in allSamples are added to hAllSamples without call to Clone to change name
                   Color_t majorStackColor,
                   bool hMajorStacks_add)
{

    // for printing statements
    TString h_name;
    try
    {
        h_name = ((TH1D*)allSamples[0]->At(i))->GetName();
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "h_name is not set" << std::endl;
    }

    // i = histogram index (type, eg; electron energy min/max, total energy, track length...)
    // j = background/data/mc index

    TH1D *tmpHist; //we'll be using this later.

    // this block relating to hAllSamples
    // which is one of hAllExternals, hAllInternals, etc...
    // add all samples of a type into a histogram for the group hAllSamples
    for(int j = 0; j < nSamples; j++)
    {
        if(h_name.Contains("hTotalE"))
        {
            std::cout << "stackfunction: j=" << j << std::endl;
        }

        if(hAllSamples_addtoexisting == false)
        {
            if(j == 0)
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllSamples[i] = (TH1D*)allSamples[j]->At(i)->Clone(sampleName);
                if(h_name.Contains("hTotalE"))
                {
                    std::cout << "adding to hAllSamples[" << i << "]: " << ((TH1D*)allSamples[j]->At(i))->GetName() << " (clone)" << std::endl;
                }
            }
            else
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllSamples[i]->Add((TH1D*)allSamples[j]->At(i));
                if(h_name.Contains("hTotalE"))
                {
                    std::cout << "adding to hAllSamples[" << i << "]: " << ((TH1D*)allSamples[j]->At(i))->GetName() << std::endl;
                }
            }
        }
        else
        {
            // NOTE: allSamples is either allExternals, allInternals, etc...
            hAllSamples[i]->Add((TH1D*)allSamples[j]->At(i));
            if(h_name.Contains("hTotalE"))
            {
                std::cout << "adding to hAllSamples[" << i << "]: " << ((TH1D*)allSamples[j]->At(i))->GetName() << std::endl;
            }
        }

        // 2020-03-31: moved from below block to this location
        tmpHist = (TH1D*)allSamples[j]->At(i);
        if(tmpHist->Integral() > 0)
        {
            //std::cout << ExternalBkgFiles[j] << " " << tmpHist->Integral() << std::endl;
            tmpHist->SetTitle(sampleNames[j]); // TODO: want this above?
            hMinorStacks[i]->Add((TH1D*)allSamples[j]->At(i)->Clone(sampleFiles[j]));
        }
        else
        {
            #if 0
            std::cout << "MC: " << sampleNames[j] << ", HIST: " << tmpHist->GetName() << ", Integral() == 0" << std::endl;
            #endif
        }


        // below block all related to hAllMC

        // allMC goes here
        if(allMC_first == true)
        {
            if(j == 0)
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllMC[i] = (TH1D*)allSamples[j]->At(i)->Clone("Total MC");
            }
            else
            {
                // NOTE: allSamples is either allExternals, allInternals, etc...
                hAllMC[i]->Add((TH1D*)allSamples[j]->At(i));
            }
        }
        else
        {
            // NOTE: allSamples is either allExternals, allInternals, etc...
            hAllMC[i]->Add((TH1D*)allSamples[j]->At(i));
        }

   
        //(TH1D*)allSamples[j]->At(i)->Write();
        // 2020-03-30: removed
        

    } // generic


    TString events;
    //std::cout << "150Nd: Add()" << std::endl;
    //events.Form("%i",(int)hAllNd150[i]->Integral());
    //hAllNd150[i]->SetTitle("^{150}Nd 2#nu2#beta ("+events+")");
    //hMajorStacks[i]->Add((TH1D*)hAllNd150[i]);
    
    if(hMajorStacks_add == true)
    {
        if(h_name.Contains("hTotalE"))
        {
            std::cout << sampleName << ": Add()" << std::endl;
        }
        events.Form("%i", (int)hAllSamples[i]->Integral());
        TString title_string = sampleHumanName + " (" + events + ")";
        //tmpHist->SetFillColor(dataColors[i]);
        //tmpHist->SetLineColor(dataColors[i]);
        hAllSamples[i]->SetFillColor(majorStackColor);
        hAllSamples[i]->SetLineColor(majorStackColor);
        hAllSamples[i]->SetLineWidth(0);
        hAllSamples[i]->SetTitle(title_string);
        hMajorStacks[i]->Add((TH1D*)hAllSamples[i]);
        
        //std::cout << "Externals: Add()" << std::endl;
        //events.Form("%i",(int)hAllExternals[i]->Integral());
        //hAllExternals[i]->SetTitle("Externals ("+events+")");
        //hMajorStacks[i]->Add((TH1D*)hAllExternals[i]);
    }
    else
    {
        // do nothing
    }
    

}


// TODO re-enable RAWENABLE ?
// need to save into file in different location, or not save at all
void fitHistograms()
{
    // TODO: draw the 150nd 2D histograms for hMCElectronEnergy and
    // hMCTrueElectronEnergy


    std::cout << "fitHistograms()" << std::endl;

    TFile *myFile = TFile::Open("Nd150_2e_P"+Phase+".root", "UPDATE");
    TDirectory* dir = gDirectory; 
    //std::cout << "dir->ls():" << std::endl;
    //dir->ls();
    TObjArray *allDataHistograms = new TObjArray();    
    TObjArray *allDataHistograms_rawdata = new TObjArray();    
    
    //Get the data histograms
    for(int i = 0; i < numHistograms; i++)
    {

        //TString name_append;
        //if(mode_flag == 0)
        //{
        //    //name_append = "";
        //}
        //else if(mode_flag == 1)
        //{
        //    name_append = "_raw";
        //}

        TString histogram_name = histogramNames[i];

        // Note: DataFile = "data_2e"
        // processed
        {
            //std::string name("singleHistos/unscaled/" + histogramNames[i] + DataFile);
            //std::string name("rawdata/" + histogram_name  + "/" + histogramNames[i] + DataFile);
            //std::string name("processeddata/" + histogram_name  + "/" + histogramNames[i] + DataFile + name_append);
            std::string name("processeddata/" + histogram_name  + "/" + histogramNames[i] + DataFile);
            //std::cout << "name=" << name << std::endl;
            if(histogram_name.CompareTo("hTotalE_") == 0)
            {
                std::cout << "Get() : " << name << " from file, Clone() : " << histogramNames[i] + "data" << std::endl;
            }
            allDataHistograms->Add((TH1D*)myFile->Get(name.c_str())->Clone(histogramNames[i] + "data"));
            // 2020-04-02: changed "_data" to "data"
        }

#if RAWENABLE
        // raw
        {
            //std::string name_rawdata("rawdata/" + histogram_name  + "/" + histogramNames[i] + DataFile + name_append);
            std::string name_rawdata("rawdata/" + histogram_name  + "/" + histogramNames[i] + DataFile + "_raw");
            if(histogram_name.CompareTo("hTotalE_") == 0)
            {
                std::cout << "Get() : " << name_rawdata << " from file, Clone() : " << histogramNames[i] + "data_raw" << std::endl;
            }
            allDataHistograms_rawdata->Add((TH1D*)myFile->Get(name_rawdata.c_str())->Clone(histogramNames[i] + "data_raw"));
        }
#endif
    }





#if 0

            std::cout << "DONT DO THIS" << std::endl;
            std::cin.get();

            TString histogram_name = histogramNames[j];
            TString new_histogram_name = histogram_name + sampleFiles[i] + "_fit";
            //TString name = "rawdata/" + histogram_name + "/" + histogram_name + sampleFiles[i];

            if(mode_flag == 0)
            {
                TString name = "processeddata/" + histogram_name + "/" + histogram_name + sampleFiles[i]; // + name_append
                std::cout << "Get() : " << name << " from file, Clone() : " << new_histogram_name << std::endl;
                tmpHist = (TH1D*)myFile->Get(name)->Clone(new_histogram_name);
            }
            else if(mode_flag == 1)
            {
                TString name = "rawdata/" + histogram_name + "/" + histogram_name + sampleFiles[i] + "_raw"; // + name_append;
                std::cout << "Get() : " << name << " from file, Clone() : " << new_histogram_name << std::endl;
                tmpHist = (TH1D*)myFile->Get(name)->Clone(new_histogram_name);
            }
            // NOTE: will have 2 new histograms here
            // 1: scaled, processeddata histogram (save this one to file)
            // 2: scaled, rawdata histogram (do not save this one to file)

            //tmpHist->Sumw2();
            tmpHist->SetFillColor(sampleColors[i]);
            tmpHist->SetLineColor(sampleColors[i]);
            tmpHist->SetLineWidth(0);
            tmpHist->SetTitle(sampleNames[i]);

            tmpHist->Scale(TotalTime / sampleNGenMC[i]);
            if(mode_flag == 0)
            {
                tmpHist->Write();
            }
            tmpHist->Scale(sampleActivity[thePhase * nSamples + i]);

            // NOTE: done AFTER call to Write()
            // so sample is scaled with activity value if used again from allSamples
            allSamples[i]->Add(tmpHist);


#endif















    // TODO: TotalTime is 33859178. + 133770114.
    // is this correct? perhaps this is wrong?
    // considering dead_time_scale_factor
    // unless dead_time_scale_factor is actually a re-scaling for the entire
    // detector P1 + P2 time
    AcceptedTime[0] = 33859178.;
    AcceptedTime[1] = 133770114.;
    AcceptedTime[2] = TotalTime;
                              
    std::ofstream oftemp("oftemp1.txt", std::ofstream::out | std::ofstream::trunc);
    oftemp.close();

/*
void scale(TFile* myFile,                       // INPUT: unscaled histograms are loaded from this root file
           const std::string& activityfile,     // INPUT: activity values read from this file
           const std::string& typedir,          // INPUT: forms part of path string, eg: "nd150/", "internals/", "externals/", "neighbours/", all radon files are in "externals/"
           const int nData,                     // INPUT: number of objects in dataFiles/dataNames
           const TString *const dataFiles,      // INPUT: list of samples, eg: bi214_sfoil_rot
           const TString *const dataNames,      // INPUT: corresponding human readable names, Latex format, eg: ^{214}Bi Foil Surface
           const Color_t *const dataColor,      // INPUT: (unused)
           const Color_t *const dataColors,     // INPUT: data color set using values in this array
           const double *const AcceptedTime,    // INPUT: accepted time for each phase is passed in here 
           double *const dataNGenMC,            // OUTPUT: number of generated MC events are written to this location
           double *const dataActivity,          // OUTPUT: activity values written to this location
           double *const dataEfficiency,        // OUTPUT: efficiency values written to this location
           TObjArray** allSamples)              // OUTPUT: each histogram has colors set and is scaled (using scaling with activity), before being stored as an output here
*/
    // documenation for scale() function:
    //
    // myFile is a TFile containing all raw histograms (unscaled)
    // "../include/activities.txt" -> activityfile is a path and filename for activity file which contains all sample activties (pre-fit)
    // "internals/" -> typedir forms part of path string to MC files
    // nInternalBkgs -> nData contains number of MC histograms to process
    // InternalBkgFiles -> dataFiles contains list of histogram names, eg: bi214_int_rot
    // InternalBkgNames -> dataNames contains corresponding human readable names, eg: ^{214}Bi Internal
    // InternalBkgColor -> dataColor is unused
    // InternalBkgColors -> dataColors is used to set the line and fill color for each MC sample
    // AcceptedTune Phase 1 and Phase 2 accepted time
    // InternalBkgNGenMC -> dataNGenMC is an output containing number of generated MC
    // InternalBkgActivity -> dataActivity is an output containing the activity for Phase 1 and Phase 2
    // InternalBkgEfficiency -> dataEfficiency is an output containing the efficiency for Phase 1 and Phase 2
    // allInternals -> allSamples is an output, format TObjArray, contains each MC histogram (as a clone) in a TObjArray

    // after scale function, these contain
    // TH1D objects
    // the named object below is an array of TObjArray pointers
    // the TObjArray is allocated with new in scale
    // there is one for each background type
    // a single TH1D* object is added to the TObjArray
    // so it is an array of length 1
    TObjArray *allExternals[nExternalBkgs];
    TObjArray *allRn220Bkgs[nRn220Bkgs];
    TObjArray *allRn222Bkgs[nRn222Bkgs];
    TObjArray *allInternals[nInternalBkgs];
    TObjArray *allNeighbours[nNeighbours];
    TObjArray *allNd150Samples[nNd150Samples];
    // rawdata versions
    TObjArray *allExternals_rawdata[nExternalBkgs];
    TObjArray *allRn220Bkgs_rawdata[nRn220Bkgs];
    TObjArray *allRn222Bkgs_rawdata[nRn222Bkgs];
    TObjArray *allInternals_rawdata[nInternalBkgs];
    TObjArray *allNeighbours_rawdata[nNeighbours];
    TObjArray *allNd150Samples_rawdata[nNd150Samples];

    // clear this output file
    std::ofstream of_numberofevents("of_numberofevents.txt", std::ofstream::out | std::ofstream::app);

#if INTERNALS_ON
    std::cout << "scale: internals" << std::endl;

    scale(myFile,
          "../include/activities.txt",
          "internals/",
          nInternalBkgs,
          InternalBkgFiles,
          InternalBkgNames,
          InternalBkgColors,
          AcceptedTime,
          InternalBkgNGenMC,
          &InternalBkgActivity[0][0],
          &InternalBkgEfficiency[0][0],
          allInternals);
    
    #if RAWENABLE

    std::cout << "RAWENABLE" << std::endl;
    std::cin.get();

    scale(myFile,
          "../include/activities.txt",
          "internals/",
          nInternalBkgs,
          InternalBkgFiles,
          InternalBkgNames,
          InternalBkgColors,
          AcceptedTime,
          InternalBkgNGenMC_rawdata,
          &InternalBkgActivity_rawdata[0][0],
          &InternalBkgEfficiency_rawdata[0][0],
          allInternals_rawdata, 1);
    #endif
#endif

#if NEIGHBOURS_ON
    std::cout << "scale: neighbours" << std::endl;

    scale(myFile,
          "../include/activities.txt",
          "neighbours/",
          nNeighbours,
          NeighbourFiles,
          NeighbourNames,
          NeighbourColors,
          AcceptedTime,
          NeighbourNGenMC,
          &NeighbourActivity[0][0],
          &NeighbourEfficiency[0][0],
          allNeighbours);
    
    #if RAWENABLE
    scale(myFile,
          "../include/activities.txt",
          "neighbours/",
          nNeighbours,
          NeighbourFiles,
          NeighbourNames,
          NeighbourColors,
          AcceptedTime,
          NeighbourNGenMC_rawdata,
          &NeighbourActivity_rawdata[0][0],
          &NeighbourEfficiency_rawdata[0][0],
          allNeighbours_rawdata, 1);
    #endif
#endif

#if ND_ON
    std::cout << "scale: 150Nd" << std::endl;

    scale(myFile,
          "../include/activities.txt",
          "nd150/",
          nNd150Samples,
          Nd150Files,
          Nd150Names,
          Nd150Colors,
          AcceptedTime,
          Nd150NGenMC,
          &Nd150Activity[0][0],
          &Nd150Efficiency[0][0],
          allNd150Samples);

    #if RAWENABLE
    scale(myFile,
          "../include/activities.txt",
          "nd150/",
          nNd150Samples,
          Nd150Files,
          Nd150Names,
          Nd150Colors,
          AcceptedTime,
          Nd150NGenMC_rawdata,
          &Nd150Activity_rawdata[0][0],
          &Nd150Efficiency_rawdata[0][0],
          allNd150Samples_rawdata, 1);
    #endif
#endif

#if EXTERNALS_ON
    std::cout << "scale: externals" << std::endl;

    scale(myFile,
          "../include/activities.txt",
          "externals/",
          nExternalBkgs,
          ExternalBkgFiles,
          ExternalBkgNames,
          ExternalBkgColors,
          AcceptedTime,
          ExternalBkgNGenMC,
          &ExternalBkgActivity[0][0],
          &ExternalBkgEfficiency[0][0],
          allExternals);

    #if RAWENABLE
    scale(myFile,
          "../include/activities.txt",
          "externals/",
          nExternalBkgs,
          ExternalBkgFiles,
          ExternalBkgNames,
          ExternalBkgColors,
          AcceptedTime,
          ExternalBkgNGenMC_rawdata,
          &ExternalBkgActivity_rawdata[0][0],
          &ExternalBkgEfficiency_rawdata[0][0],
          allExternals_rawdata, 1);
    #endif
#endif

#if RADON_ON
    std::cout << "scale: radon" << std::endl;

    #if 0
        scale(myFile,
              "../include/activities.txt",
              "externals/",
              nRn222Bkgs,
              Rn222BkgFiles,
              Rn222BkgNames,
              Rn222BkgColors,
              AcceptedTime,
              Rn222BkgNGenMC,
              &Rn222BkgActivity[0][0],
              &Rn222BkgEfficiency[0][0],
              allRn222Bkgs);

        #if RAWENABLE
        scale(myFile,
              "../include/activities.txt",
              "externals/",
              nRn222Bkgs,
              Rn222BkgFiles,
              Rn222BkgNames,
              Rn222BkgColors,
              AcceptedTime,
              Rn222BkgNGenMC_rawdata,
              &Rn222BkgActivity_rawdata[0][0],
              &Rn222BkgEfficiency_rawdata[0][0],
              allRn222Bkgs_rawdata, 1);
        #endif
    #else
        scale(myFile,
              "../include/activities.txt",
              "externals/",
              nRn222Bkgs,
              Rn222BkgFilesNew,
              //Rn222BkgNamesNew,
              Rn222BkgNames,
              //Rn222BkgColorsNew,
              Rn222BkgColors,
              AcceptedTime,
              Rn222BkgNGenMC,
              &Rn222BkgActivity[0][0],
              &Rn222BkgEfficiency[0][0],
              allRn222Bkgs);

        #if RAWENABLE
        scale(myFile,
              "../include/activities.txt",
              "externals/",
              nRn222Bkgs,
              Rn222BkgFilesNew,
              //Rn222BkgNamesNew,
              Rn222BkgNames,
              //Rn222BkgColorsNew,
              Rn222BkgColors,
              AcceptedTime,
              Rn222BkgNGenMC_rawdata,
              &Rn222BkgActivity_rawdata[0][0],
              &Rn222BkgEfficiency_rawdata[0][0],
              allRn222Bkgs_rawdata, 1);
        #endif
    #endif

    scale(myFile,
          "../include/activities.txt",
          "externals/",
          nRn220Bkgs,
          Rn220BkgFiles,
          Rn220BkgNames,
          Rn220BkgColors,
          AcceptedTime,
          Rn220BkgNGenMC,
          &Rn220BkgActivity[0][0],
          &Rn220BkgEfficiency[0][0],
          allRn220Bkgs);

    #if RAWENABLE
    scale(myFile,
          "../include/activities.txt",
          "externals/",
          nRn220Bkgs,
          Rn220BkgFiles,
          Rn220BkgNames,
          Rn220BkgColors,
          AcceptedTime,
          Rn220BkgNGenMC_rawdata,
          &Rn220BkgActivity_rawdata[0][0],
          &Rn220BkgEfficiency_rawdata[0][0],
          allRn220Bkgs_rawdata, 1);
    #endif
#endif

    // myFile->Close();

    std::cout << "Stack" << std::endl;
    // hAllExternals, hAllRadon, hAllInternals, hAllNeighbours, hAllNd150 are
    // passed into stackfunction() as argument "hAllSamples"
    // stackfunction calls Add() on object hAllSamples to add each histogram
    // component from each "category" (externals/radon/internals/neighbours
    // /Nd150) to hAllSamples
    // each histogram is pulled from TObjArray allSamples, which is an argument
    // passed in from one of allExternals, allRadon, allInternals,
    // AllNeighbours, allNd150
    // 
    
    // single histogram containing sum of all sub-samples, created using
    // TH1::Add()
    TH1D *hAllExternals[numHistograms];
    TH1D *hAllRadon[numHistograms];
    TH1D *hAllInternals[numHistograms];
    TH1D *hAllNeighbours[numHistograms];
    TH1D *hAllNd150[numHistograms];
    // rawdata versions
    TH1D *hAllExternals_rawdata[numHistograms];
    TH1D *hAllRadon_rawdata[numHistograms];
    TH1D *hAllInternals_rawdata[numHistograms];
    TH1D *hAllNeighbours_rawdata[numHistograms];
    TH1D *hAllNd150_rawdata[numHistograms];

    TH1D *hAllMC[numHistograms];
    TH1D *hAllMC_rawdata[numHistograms];
    // for all MC samples to calculate chi2
    // also used to draw black line around stack of all MC samples

    // stacks for drawing: hMajorStacks: samples are grouped by type,
    // eg: radon, external, internal, neighbour, 150Nd
    // hMinorStacks: samples are not grouped, each individual sample
    // is displayed seperatly
    THStack *hMajorStacks[numHistograms];
    THStack *hMinorStacks[numHistograms];
    THStack *hMajorStacks_rawdata[numHistograms];
    THStack *hMinorStacks_rawdata[numHistograms];

    // I don't know why this is here
    //dir->cd("singleHistos/");
    
    // -1 because we don't want the 2D hist
    // TODO: think this -1 is wrong, because there are several 2d hists
    //for(int i = 0; i < numHistograms-1; i++)
    for(int i = 0; i < numHistograms; i++)
    {
        TString i_str;
        i_str.Form("%i", i);

        hMajorStacks[i] = new THStack("majorstack" + i_str, histogramNames[i]);
        hMinorStacks[i] = new THStack("minorstack" + i_str, histogramNames[i]);
    #if RAWENABLE
        hMajorStacks_rawdata[i] = new THStack("majorstack_rawdata" + i_str, histogramNames[i]);
        hMinorStacks_rawdata[i] = new THStack("minorstack_rawdata" + i_str, histogramNames[i]);
    #endif



#if EXTERNALS_ON
        if(i == 0)
        {
            std::cout << "stack: externals" << std::endl;
        }

        stackfunction(i, nExternalBkgs,
                      hMinorStacks, allExternals, // TObjArray *allExternals[nExternalBkgs]
                      hMajorStacks, "External Backgrounds",
                      "Externals", hAllExternals, // TH1D *hAllExt[numHistograms]
                      ExternalBkgNames, ExternalBkgFiles, hAllMC,
                      true, false,
                      ExternalBkgColor,
                      true);

    #if RAWENABLE
        stackfunction(i, nExternalBkgs,
                      hMinorStacks_rawdata, allExternals_rawdata, // TObjArray *allExternals[nExternalBkgs]
                      hMajorStacks_rawdata, "External Backgrounds",
                      "Externals", hAllExternals_rawdata, // TH1D *hAllExt[numHistograms]
                      ExternalBkgNames, ExternalBkgFiles, hAllMC_rawdata,
                      true, false,
                      ExternalBkgColor,
                      true);
    #endif
#endif

#if RADON_ON
        if(i == 0)
        {
            std::cout << "stack: radon" << std::endl;
        }

        // first
        // add to existing
        // add to majors
        #if 0
        stackfunction(i, nRn222Bkgs,
                      hMinorStacks, allRn222Bkgs,
                      hMajorStacks, "Radon Backgrounds",
                      "Radon", hAllRadon,
                      Rn222BkgNames, Rn222BkgFiles, hAllMC,
                      false, false,
                      RadonBkgColor,
                      false);

            #if RAWENABLE
        stackfunction(i, nRn222Bkgs,
                      hMinorStacks_rawdata, allRn222Bkgs_rawdata,
                      hMajorStacks_rawdata, "Radon Backgrounds",
                      "Radon", hAllRadon_rawdata,
                      Rn222BkgNames, Rn222BkgFiles, hAllMC_rawdata,
                      false, false,
                      RadonBkgColor,
                      false);
            #endif
        #else
        if(i == 0)
        {
            std::cout << "stackfunction Rn222" << std::endl;
        }
        stackfunction(i, nRn222Bkgs,
                      hMinorStacks, allRn222Bkgs,
                      hMajorStacks, "Radon Backgrounds",
                      "Radon", hAllRadon,
                      Rn222BkgNames, Rn222BkgFilesNew, hAllMC,
                      #if EXTERNALS_ON
                      false, false,
                      #else
                      true, false,
                      #endif
                      RadonBkgColor,
                      false);

            #if RAWENABLE
        if(i == 0)
        {   
            std::cout << "stackfunction Rn222 (raw)" << std::endl;
        }
        stackfunction(i, nRn222Bkgs,
                      hMinorStacks_rawdata, allRn222Bkgs_rawdata,
                      hMajorStacks_rawdata, "Radon Backgrounds",
                      "Radon", hAllRadon_rawdata,
                      Rn222BkgNames, Rn222BkgFilesNew, hAllMC_rawdata,
                      #if EXTERNALS_ON
                      false, false,
                      #else
                      true, false,
                      #endif
                      RadonBkgColor,
                      false);
            #endif
        #endif


        if(i == 0)
        {
            std::cout << "stackfunction Rn220" << std::endl;
        }
        stackfunction(i, nRn220Bkgs,
                      hMinorStacks, allRn220Bkgs,
                      hMajorStacks, "Radon Backgrounds",
                      "Radon", hAllRadon,
                      Rn220BkgNames, Rn220BkgFiles, hAllMC,
                      false, true,
                      RadonBkgColor,
                      true);

            #if RAWENABLE
        if(i == 0)
        {
            std::cout << "stackfunction Rn220 (raw)" << std::endl;
        }
        stackfunction(i, nRn220Bkgs,
                      hMinorStacks_rawdata, allRn220Bkgs_rawdata,
                      hMajorStacks_rawdata, "Radon Backgrounds",
                      "Radon", hAllRadon_rawdata,
                      Rn220BkgNames, Rn220BkgFiles, hAllMC_rawdata,
                      false, true,
                      RadonBkgColor,
                      true);
            #endif
#endif

#if INTERNALS_ON
        if(i == 0)
        {
            std::cout << "stack: internals" << std::endl;
        }

        stackfunction(i, nInternalBkgs,
                      hMinorStacks, allInternals,
                      hMajorStacks, "Internal Backgrounds",
                      "Internals", hAllInternals,
                      InternalBkgNames, InternalBkgFiles, hAllMC,
                      false, false,
                      InternalBkgColor,
                      true);

    #if RAWENABLE
        stackfunction(i, nInternalBkgs,
                      hMinorStacks_rawdata, allInternals_rawdata,
                      hMajorStacks_rawdata, "Internal Backgrounds",
                      "Internals", hAllInternals_rawdata,
                      InternalBkgNames, InternalBkgFiles, hAllMC_rawdata,
                      false, false,
                      InternalBkgColor,
                      true);
    #endif
#endif

#if NEIGHBOURS_ON
        if(i == 0)
        {
            std::cout << "stack: neighbours" << std::endl;
        }

        stackfunction(i, nNeighbours,
                      hMinorStacks, allNeighbours,
                      hMajorStacks, "Neighbouring Foils",
                      "Neighbours", hAllNeighbours,
                      NeighbourNames, NeighbourFiles, hAllMC,
                      false, false,
                      NeighbourColor,
                      true);

    #if RAWENABLE
        stackfunction(i, nNeighbours,
                      hMinorStacks_rawdata, allNeighbours_rawdata,
                      hMajorStacks_rawdata, "Neighbouring Foils",
                      "Neighbours", hAllNeighbours_rawdata,
                      NeighbourNames, NeighbourFiles, hAllMC_rawdata,
                      false, false,
                      NeighbourColor,
                      true);
    #endif
#endif

#if ND_ON
        if(i == 0)
        {
            std::cout << "stack: 150Nd" << std::endl;
        }

        stackfunction(i, nNd150Samples,
                      hMinorStacks, allNd150Samples,
                      hMajorStacks, "^{150}Nd 2#nu2#beta",
                      "Nd150", hAllNd150,
                      Nd150Names, Nd150Files, hAllMC,
                      false, false,
                      Nd150Color,
                      true);

    #if RAWENABLE
        stackfunction(i, nNd150Samples,
                      hMinorStacks_rawdata, allNd150Samples_rawdata,
                      hMajorStacks_rawdata, "^{150}Nd 2#nu2#beta",
                      "Nd150", hAllNd150_rawdata,
                      Nd150Names, Nd150Files, hAllMC_rawdata,
                      false, false,
                      Nd150Color,
                      true);
    #endif
#endif

        
    }

    std::cout << "Completed the numHistograms loop" << std::endl;

    // stack functions for 2D histograms
/*
    TH2D *hAllMC2D[numHistograms]; //for all MC samples to calculate chi2
    TH2D *tmp2D;
    for(int i = numHistograms-1; i < numHistograms; i++)
    {
        //stack doesnt support 3D hist
        TString i_str;
        i_str.Form("%i",i);

        #if EXTERNALS_ON
        for( int j = 0; j < nExternalBkgs; j++)
        {
            (TH2D*)allExternals[j]->At(i)->Write();
            if(j == 0)
            {
                // TODO: this one different because of clone option
                hAllMC2D[i] = (TH2D*)allExternals[j]->At(i)->Clone("allMC"+histogramNames[i]);
            }
            else
            {
                hAllMC2D[i]->Add((TH2D*)allExternals[j]->At(i)->Clone());
            }
        }//externals
        #endif
    
        #if RADON_ON
        for(int j = 0; j < nRn222Bkgs; j++)
        {
            (TH2D*)allRn222Bkgs[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2D*)allRn222Bkgs[j]->At(i)->Clone());
        }//radon222

        for(int j = 0; j < nRn220Bkgs; j++)
        {
            (TH2D*)allRn220Bkgs[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2D*)allRn220Bkgs[j]->At(i)->Clone());
        }//radon220
        #endif

        #if INTERNALS_ON
        for(int j = 0; j < nInternalBkgs; j++)
        {
            (TH2D*)allInternals[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2D*)allInternals[j]->At(i)->Clone());
        }//internals
        #endif

        #if NEIGHBOURS_ON
        for(int j = 0; j < nNeighbours; j++)
        {
            (TH2D*)allNeighbours[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2D*)allNeighbours[j]->At(i)->Clone());
        }//neighbours
        #endif

        #if ND_ON
        for(int j = 0; j < nNd150Samples; j++)
        {
            (TH2D*)allNd150Samples[j]->At(i)->Write();
            hAllMC2D[i]->Add((TH2D*)allNd150Samples[j]->At(i)->Clone());
        }//150Nd
        #endif

        hAllMC2D[i]->Write();

        tmp2D = (TH2D*)allDataHistograms->At(i);
        tmp2D->SetName(histogramNames[i]+DataName);
        tmp2D->Write();
    }
    std::cout << "Finished the cloning loop" << std::endl;
*/



    // loop to write out histograms before running canvas loop
    // also sets hist structure for use in following loop
    TH1D* hist[numHistograms];
    for(int i = 0; i < numHistograms; i++)
    {

        hist[i] = (TH1D*)allDataHistograms->At(i);
        //hist[i]->SetName(histogramNames[i] + DataName);

        // inserted from below loop (2)
        TString Ndata_str;
        Ndata_str.Form("%i", (int)hist[i]->Integral());
        hist[i]->SetTitle("Data (" + Ndata_str + ")");
        hist[i]->SetLineWidth(2);
        hist[i]->SetMarkerStyle(20);

        // 2020-04-02: removed
        // 2020-04-02: re-enabled, this is required to write the data
        // histogram to file in the root directory
        // 2020-04-02: disabled again, will read data histogram from
        // /singleHistos/unscaled directory in fitting code
        // 2020-04-02: requires fiddling around with GetListOfKeys() and
        // cd() to different TDirectory, can't be bothered with that
        // so will re-enable
        std::cout << "writing histogram with name " << hist[i]->GetName() << " at end of code which was removed" << std::endl;
// disabled 2020-12-11       hist[i]->Write();
        hist[i]->Write();
        // TODO: don't like that other functions write out the equivalent
        // MC histograms to file in a different place from this Write() call
        // expecially since this call to Write() is inside a block of code
        // which creates a TCanvas - move this and SetName() elsewhere

    }
    #if RAWENABLE
    // same for rawdata versions
    TH1D* hist_rawdata[numHistograms];
    for(int i = 0; i < numHistograms; i++)
    {

        hist_rawdata[i] = (TH1D*)allDataHistograms_rawdata->At(i);
        hist_rawdata[i]->SetName(histogramNames[i] + DataName);
        // TODO: bug here: same name

        // inserted from below loop (2)
        TString Ndata_str;
        Ndata_str.Form("%i", (int)hist_rawdata[i]->Integral());
        hist_rawdata[i]->SetTitle("Data (" + Ndata_str + ")");
        hist_rawdata[i]->SetLineWidth(2);
        hist_rawdata[i]->SetMarkerStyle(20);

        // 2020-04-02: removed
        // 2020-04-02: re-enabled, this is required to write the data
        // histogram to file in the root directory
        // 2020-04-02: disabled again, will read data histogram from
        // /singleHistos/unscaled directory in fitting code
        // 2020-04-02: requires fiddling around with GetListOfKeys() and
        // cd() to different TDirectory, can't be bothered with that
        // so will re-enable
        std::cout << "writing histogram with name " << hist_rawdata[i]->GetName() << " at end of code which was removed" << std::endl;
        // TODO: fix bug by not writing this histogram to file
        // TODO: should it go in the rawdata folder?
        // TODO: it already is
        //hist_rawdata[i]->Write();
        // TODO: don't like that other functions write out the equivalent
        // MC histograms to file in a different place from this Write() call
        // expecially since this call to Write() is inside a block of code
        // which creates a TCanvas - move this and SetName() elsewhere

    }
    #endif


/*
    for(std::map<TString, TH1*>::iterator it{histogramPointers.begin();
        it != histogramPointers.end(); ++ it)
    {

        TH1D *histtmp = (TH1D*)it->second;

        std::ofstream of_hTotalE_numberofevents("hTotalE_numberofevents.txt")
        if(TString(histtmp->GetName()).Contains("hTotalE"))
        {
            std::cout << histtmp->GetName() << " number of events " << histtmp->Integral() << std::endl;
            of_hTotalE_numberofevents << histtmp->GetName() << ", " << histtmp->Integral() << std::endl;
        }
    }
*/


    // TODO: draw raw data histograms 


    TCanvas *c_processeddata[numHistograms];
    TLegend *leg_processeddata[numHistograms];

    #if RAWENABLE
    TCanvas *c_rawdata[numHistograms];
    TLegend *leg_rawdata[numHistograms];
    #endif
    
    std::cout << "About to start the canvas loop" << std::endl;
    for(int j = 0; j < numHistograms; j++)
    {

        ////////////////////
        // processed data //
        ////////////////////
        
        {
            TCanvas *&c = c_processeddata[j];
            TLegend *&leg = leg_processeddata[j];

            if(histogramDrawFlag[j] == 0) continue;

            std::cout << "At start of canvas loop" << endl;

            TString i_str;
            TString Nmc_str;
            
            i_str.Form("%i", j);
            c = new TCanvas("c_processeddata" + i_str + "_");
            c->SetFillColor(kWhite);
            // c[i]->SetGrid(1,1);
            //c[i]->SetTicks(1,1);
            //c[i]->Divide(2,2);

            //c[i]->cd(1);
            //c0 = (TPad*)c[i]->GetPad(1);
            // c1_1->SetBottomMargin(0.1);
            //c0->SetGrid(1,1);
            //c0->SetTicks(1,1);
            c->SetGrid(1, 1);
            c->SetTicks(1, 1);

            // moved to above (1)

            Nmc_str.Form("%i", (int)hAllMC[j]->Integral());
            // moved to above (2)
            // draw twice to set max
            if(j == 58)
            {
            //    hist[j]->Draw("colz");
            }
            else
            {
                hist[j]->Draw("PE");
            }

            //std::cout << " data: " << hist[i]->Integral() << std::endl;
            //hSimpleStacks[j]->Draw("hist same");
            if(j == 58)
            {
            }
            else
            {
                //hMinorStacks[j]->Draw("hist");
                hMinorStacks[j]->Draw("hist same");
                //hMajorStacks[j]->Draw("hist same");
            }


            //hAllMC[j]->SetMaximum(1.0e+04);
            //hAllMC[j]->GetYaxis()->SetMaximum(1.0e+04);
            //hAllMC[j]->GetYaxis()->SetRangeUser(0.0, 1.0e+04);

            hAllMC[j]->SetLineWidth(2);
            hAllMC[j]->SetLineColor(kBlack);
            hAllMC[j]->SetFillColor(kBlack);
            hAllMC[j]->SetFillStyle(0);
            //hAllMC[j]->Sumw2();
            if(j == 58)
            {
            }
            else
            {
                hAllMC[j]->Draw("hist same");
            }

            // moved from above
            if(j == 58)
            {
                hist[j]->Draw("colz");
            }
            else
            {
                hist[j]->Draw("PE same");
            }

            double chi2;
            int ndf, igood;
            TString chi2_str, ndf_str;
            double prob = hist[j]->Chi2TestX(hAllMC[j], chi2, ndf, igood, "UW");
            chi2_str.Form("%4.3f", chi2);
            ndf_str.Form("%i", ndf);
            hAllMC[j]->SetTitle("total MC (" + Nmc_str + ")");
            // TODO: moved Ndata_str out to above loop, consider moving out Nmc_str

            //hAllMC[i]->Write();
            // 2020-03-30: removed

            // TODO: forget difference between hAllMC and hMinorStacks
            // why is only 1 written to file? does this write all the
            // histograms to file? cannot be the case, as some
            // histos would be skipped by the continue statement above

            //leg[i] = c0->BuildLegend();
            leg = c->BuildLegend();
            leg->SetName("leg" + i_str + "_");
            // TODO: build legend using parameter_names.list rather than
            // individual MC samples
            leg->SetFillColor(kWhite);
            leg->AddEntry((TObject*)0, "#chi^{2}/ndf = " + chi2_str + "/" + ndf_str, "");
            //c->cd(2);
            TLegend *tmpLeg = (TLegend*)leg->Clone("leg_copy");
            tmpLeg->Draw();
            //c0->cd();
            leg->Delete();
            hist[j]->SetTitle(histogramNames[j]); //switch to titles once you've written them
            /*
            if(j == 58)
            {
            }
            else
            {
                hist[j]->Draw("PEsame");
            }
            */

            //c[i]->cd(3);
            //c2 = (TPad*)c[i]->GetPad(3);
            // c1_1->SetBottomMargin(0.1);
            //c2->SetGrid(1,1);
            //c2->SetTicks(1,1);

            //make ratio plot
            //TH1D *ratio = (TH1D*)hist[i]->Clone("ratio_" + histogramNames[i]);
            //ratio->Sumw2();
            //ratio->Divide((TH1D*)allMC[i]);
            //ratio->SetLineWidth(2);
            //ratio->SetMarkerStyle(7);
            //ratio->SetTitle("");
            //ratio->GetYaxis()->SetTitle("data / MC");
            //ratio->GetYaxis()->CenterTitle();

            // ratio->GetYaxis()->SetTitleSize(.1);
            // ratio->GetYaxis()->SetLabelSize(.11);
            // ratio->GetYaxis()->SetNdivisions(9+(100*5));
            // ratio->GetXaxis()->SetTitleSize(.15);

            //  ratio->GetXaxis()->SetTitleOffset(.8);
            // ratio->GetXaxis()->SetLabelOffset(.01);
            //    ratio->GetXaxis()->SetLabelSize(.11);
            //ratio->Draw("PE");
            std::cout << "Iteration of canvasy loop finished" << std::endl;
        }


        //////////////
        // raw data //
        //////////////
        
    #if RAWENABLE
        if(histogramDrawFlag_rawdata == 1)
        {
            TCanvas *&c = c_rawdata[j];
            TLegend *&leg = leg_rawdata[j];

            if(histogramDrawFlag[j] == 0) continue;

            std::cout << "At start of canvas loop" << endl;

            TString i_str;
            TString Nmc_str;
            
            i_str.Form("%i", j);
            c = new TCanvas("c_rawdata" + i_str + "_");
            c->SetFillColor(kWhite);
            // c[i]->SetGrid(1,1);
            //c[i]->SetTicks(1,1);
            //c[i]->Divide(2,2);

            //c[i]->cd(1);
            //c0 = (TPad*)c[i]->GetPad(1);
            // c1_1->SetBottomMargin(0.1);
            //c0->SetGrid(1,1);
            //c0->SetTicks(1,1);
            c->SetGrid(1, 1);
            c->SetTicks(1, 1);

            // moved to above (1)

            Nmc_str.Form("%i", (int)hAllMC_rawdata[j]->Integral());
            // moved to above (2) 
            if(j == 58)
            {
                //hist_rawdata[j]->Draw("colz");
            }
            else
            {
                hist_rawdata[j]->Draw("PE");
            }

            //std::cout << " data: " << hist[i]->Integral() << std::endl;
            //hSimpleStacks[i]->Draw("hist same");
            if(j == 58)
            {
            }
            else
            {
                //hMinorStacks_rawdata[j]->Draw("hist");
                //hMajorStacks_rawdata[j]->Draw("hist same");
                hMinorStacks_rawdata[j]->Draw("hist same");
            }


            //hAllMC_rawdata[j]->SetMaximum(1.0e+04);
            //hAllMC_rawdata[j]->GetYaxis()->SetRangeUser(0.0, 1.0e+04);

            hAllMC_rawdata[j]->SetLineWidth(2);
            hAllMC_rawdata[j]->SetLineColor(kBlack);
            hAllMC_rawdata[j]->SetFillColor(kBlack);
            hAllMC_rawdata[j]->SetFillStyle(0);
            //hAllMC_rawdata[j]->Sumw2();
            if(j == 58)
            {
            }
            else
            {
                hAllMC_rawdata[j]->Draw("hist same");
                //hAllMC_rawdata[j]->Draw("hist");
            }

            // moved to here
            if(j == 58)
            {
                hist_rawdata[j]->Draw("colz");
            }
            else
            {
                hist_rawdata[j]->Draw("PE same");
            }

            double chi2;
            int ndf, igood;
            TString chi2_str, ndf_str;
            double prob = hist_rawdata[j]->Chi2TestX(hAllMC_rawdata[j], chi2, ndf, igood, "UW");
            chi2_str.Form("%4.3f", chi2);
            ndf_str.Form("%i", ndf);
            hAllMC_rawdata[j]->SetTitle("total MC (" + Nmc_str + ")");
            // TODO: moved Ndata_str out to above loop, consider moving out Nmc_str

            //hAllMC[i]->Write();
            // 2020-03-30: removed

            // TODO: forget difference between hAllMC and hMinorStacks
            // why is only 1 written to file? does this write all the
            // histograms to file? cannot be the case, as some
            // histos would be skipped by the continue statement above

            //leg[i] = c0->BuildLegend();
            leg = c->BuildLegend();
            leg->SetName("leg_rawdata" + i_str + "_");
            leg->SetFillColor(kWhite);
            leg->AddEntry((TObject*)0, "#chi^{2}/ndf = " + chi2_str + "/" + ndf_str, "");
            //c[i]->cd(2);
            TLegend *tmpLeg = (TLegend*)leg->Clone("leg_rawdata_copy");
            tmpLeg->Draw();
            //c0->cd();
            leg->Delete();
            hist_rawdata[j]->SetTitle(histogramNames[j]); //switch to titles once you've written them
            /*
            if()
            {
                hist_rawdata[j]->Draw("colz");
            }
            else
            {
                hist_rawdata[j]->Draw("PEsame");
            }
            */

            //make ratio plot
            //TH1D *ratio = (TH1D*)hist[i]->Clone("ratio_" + histogramNames[i]);
            //ratio->Sumw2();
            //ratio->Divide((TH1D*)allMC[i]);
            //ratio->SetLineWidth(2);
            //ratio->SetMarkerStyle(7);
            //ratio->SetTitle("");
            //ratio->GetYaxis()->SetTitle("data / MC");
            //ratio->GetYaxis()->CenterTitle();

            // ratio->GetYaxis()->SetTitleSize(.1);
            // ratio->GetYaxis()->SetLabelSize(.11);
            // ratio->GetYaxis()->SetNdivisions(9+(100*5));
            // ratio->GetXaxis()->SetTitleSize(.15);

            //  ratio->GetXaxis()->SetTitleOffset(.8);
            // ratio->GetXaxis()->SetLabelOffset(.01);
            //    ratio->GetXaxis()->SetLabelSize(.11);
            //ratio->Draw("PE");
            std::cout << "Iteration of canvasy loop finished" << std::endl;
        }
    #endif
   

/*
    TCanvas *c[numHistograms];
    TH1D* hist[numHistograms];
    TLegend *leg[numHistograms];
    TPad *c0;
    TPad *c1;
    TPad *c2;
    std::cout << "About to start the canvasy loop" << std::endl;
    // for ( int i = 2; i < 3; i++) {
    for( int i = 2; i < 4; i++)
    {
        std::cout << "At start of canvasy loop" << endl;

        TString i_str, Ndata_str, Nmc_str;
        i_str.Form("%i",i);
        c[i] = new TCanvas("c"+i_str+"_");
        c[i]->SetFillColor(kWhite);
        // c[i]->SetGrid(1,1);
        //c[i]->SetTicks(1,1);
        c[i]->Divide(2,2);

        c[i]->cd(1);
        c0 = (TPad*)c[i]->GetPad(1);
        // c1_1->SetBottomMargin(0.1);
        c0->SetGrid(1,1);
        c0->SetTicks(1,1);

        hist[i] = (TH1D*)allDataHistograms->At(i);
        hist[i]->SetName(histogramNames[i]+DataName);
        hist[i]->Write();

        Ndata_str.Form("%i",(int)hist[i]->Integral());
        Nmc_str.Form("%i",(int)allMC[i]->Integral());
        hist[i]->SetTitle("Data ("+Ndata_str+")");
        hist[i]->SetLineWidth(2);
        hist[i]->SetMarkerStyle(20);
        hist[i]->Draw("PE");

        //Justin printing out integrals:
        std::cout << "Justin printing out integrals:" << std::endl;
        std::cout << "Data integral: " << hist[i]->Integral() << endl;
        std::cout << "MC integral: " << allMC[i]->Integral() << endl;
        
        std::cout << " data: " << hist[i]->Integral() << std::endl;
        //hSimpleStacks[i]->Draw("hist same");
        hMinorStacks[i]->Draw("hist same");

        allMC[i]->SetLineWidth(2);
        allMC[i]->SetLineColor(kBlack);
        allMC[i]->SetFillColor(kWhite);
        allMC[i]->SetFillStyle(0);
        allMC[i]->Sumw2();
        allMC[i]->Draw("hist sames");

        double chi2;
        int ndf, igood;
        TString chi2_str, ndf_str;
        double prob = hist[i]->Chi2TestX(allMC[i], chi2, ndf, igood, "UW");
        chi2_str.Form("%4.3f", chi2);
        ndf_str.Form("%i", ndf);
        allMC[i]->SetTitle("total MC (" + Nmc_str + ")");

        allMC[i]->Write();

        leg[i] = c0->BuildLegend();
        leg[i]->SetName("leg" + i_str + "_");
        leg[i]->SetFillColor(kWhite);
        leg[i]->AddEntry((TObject*)0, "#chi^{2}/ndf = " + chi2_str + "/" + ndf_str, "");
        c[i]->cd(2);
        TLegend *tmpLeg = (TLegend*)leg[i]->Clone("leg_copy");
        tmpLeg->Draw();
        c0->cd();
        leg[i]->Delete();
        hist[i]->SetTitle(histogramNames[i]); //switch to titles once you've written them
        hist[i]->Draw("PEsame");

        c[i]->cd(3);
        c2 = (TPad*)c[i]->GetPad(3);
        // c1_1->SetBottomMargin(0.1);
        c2->SetGrid(1,1);
        c2->SetTicks(1,1);

        //make ratio plot
        TH1D *ratio = (TH1D*)hist[i]->Clone("ratio_" + histogramNames[i]);
        ratio->Sumw2();
        ratio->Divide((TH1D*)allMC[i]);
        ratio->SetLineWidth(2);
        ratio->SetMarkerStyle(7);
        ratio->SetTitle("");
        ratio->GetYaxis()->SetTitle("data / MC");
        ratio->GetYaxis()->CenterTitle();

        // ratio->GetYaxis()->SetTitleSize(.1);
        // ratio->GetYaxis()->SetLabelSize(.11);
        // ratio->GetYaxis()->SetNdivisions(9+(100*5));
        // ratio->GetXaxis()->SetTitleSize(.15);

        //  ratio->GetXaxis()->SetTitleOffset(.8);
        // ratio->GetXaxis()->SetLabelOffset(.01);
        //    ratio->GetXaxis()->SetLabelSize(.11);
        ratio->Draw("PE");
        std::cout << "Iteration of canvasy loop finished" << std::endl;

    */
    }
    std::cout << "Finished with canvas" << std::endl;
    // drawPlots(hEe_data, hEe_totalMC, hs_hEe);


    // myFile->Close();
}


#if 0
void drawPlots(TH1D *data,TH1D *mc, THStack *hs)
{
    TString whichHist = mc->GetName();
    gStyle->SetTitleFillColor(kWhite);
    gROOT->SetStyle("Plain");

    TFile *myFile = TFile::Open("final"+whichHist+"_2e_P"+Phase+".root", "RECREATE");  

    std::cout << whichHist << std::endl;
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1_"+whichHist,"c1_"+whichHist,700,700);

    c1->SetFillColor(kWhite);
    c1->SetGrid(1,1);
    c1->Divide(2,2);

    c1->cd(1);
    TPad *c1_1;
    c1_1 = (TPad*)c1->GetPad(1);
    // c1_1->SetBottomMargin(0.1);
    c1_1->SetGrid(1,1);
    c1_1->SetTicks(1,1);

    data->SetMarkerStyle(20);
    data->SetLineWidth(2);
    data->SetTitle("");
    data->GetXaxis()->SetLabelOffset(2.);
    data->GetYaxis()->SetLabelSize(0.05);
    data->Draw("PE");

    hs->Draw("hist same");
    mc->SetLineWidth(2);
    mc->Draw("same");

    c1_1->Update();
    // TLegend *leg;
    // leg = c1_1->BuildLegend();
    // leg->SetName("leg_"+whichHist);
    // leg->SetFillColor(kWhite);
    data->Draw("PEsame");

    c1->cd(2);
    TPad *c1_2;
    c1_2 = (TPad*)c1->GetPad(2);
    TLegend *leg;
    leg = c1_1->BuildLegend();
    leg->SetName("leg_"+whichHist);
    leg->SetFillColor(kWhite);
    leg->Draw();

    c1->cd(3);
    TPad *c1_3;
    c1_3 = (TPad*)c1->GetPad(3);
    //  c1_2->SetTopMargin(0.001);
    // c1_2->SetBottomMargin(2.);
    c1_3->SetGrid(1,1);
    c1_3->SetTicks(1,1);
    double ymin;// = c1_1->GetY1();
    double ymax;// = c1_1->GetY2();
    double xmin;// = c1_1->GetX1();
    double xmax;// = c1_1->GetX2();

    c1_3->GetPadPar(xmin, ymin, xmax, ymax);
    c1_3->SetPad(xmin, ymin + ((ymax - ymin) / 4.), xmax, ymax - ((ymax - ymin) / 4.));
    // c1_2->SetBottomMargin(3.);

    TH1D *ratio = (TH1D*)data->Clone("ratio_" + whichHist); 
    ratio->Sumw2();
    ratio->Divide(mc);
    ratio->SetLineWidth(1);
    ratio->SetMarkerStyle(7);
    ratio->SetTitle("");
    ratio->GetYaxis()->SetTitle("data / MC");
    ratio->GetYaxis()->CenterTitle();
    ratio->GetYaxis()->SetTitleOffset(.2);
    ratio->GetYaxis()->SetTitleSize(.1);
    ratio->GetYaxis()->SetLabelSize(.11);
    ratio->GetYaxis()->SetNdivisions(9 + (100 * 5));
    ratio->GetXaxis()->SetTitleSize(.15);
    ratio->GetXaxis()->SetTitleOffset(.8);
    ratio->GetXaxis()->SetLabelOffset(.01);
    ratio->GetXaxis()->SetLabelSize(.11);
    ratio->Draw("PE");

    c1_3->Update();

    //  TDirectory* dir = gDirectory;
    //  dir->cd();
    //TODO this is wrong
    data->SetTitle("Phase 1OR2"); // TODO ?
    data->SetName("data_" + whichHist);
    data->Write();
    //  ratio->Write();
    hs->SetName("hs_" + whichHist);
    hs->Write();
    leg->SetName("leg_" + whichHist);
    leg->Write();
    mc->SetName("mcSum_" + whichHist);
    mc->Write();

    nHistograms++;

}
#endif


double getActErr(double Npass, double Ngen, double Ndata, double sf_err)
{
    double err = TMath::Sqrt( (1.0 / Npass) + (1.0 / Ngen) + (1.0 / Ndata) + (sf_err * sf_err) );
    return err;
}

double getChi2(TH1D *data, TH1D* mc, Int_t &ndof)
{

    size_t N_bins =  data->GetXaxis()->GetNbins();
    double chi2 = 0.;
    ndof = 0;

    for(size_t i=0; i<N_bins ; i++)
    {

        double nData = data->GetBinContent(i + 1) ;

        double nMC = mc->GetBinContent(i + 1);

        if(nData + nMC == 0) continue;

        ndof ++;

        double diff = (nMC - nData)/sqrt(nMC + nData); // gaussian prob

        chi2 += diff*diff;
    }

    //  cout << " chi2 " << chi2 << " ndof " << ndof << endl;
    //fflush(stdout);
    std::cout.flush();

    double P = TMath::Prob(chi2,ndof);

    //  cout << " Prob " << P  << endl;
    //fflush(stdout);
    std::cout.flush();

    return chi2;
}

