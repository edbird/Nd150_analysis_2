#ifndef NEWLOGLIKFITTER_REBUILD_BACKGROUND_MC_H
#define NEWLOGLIKFITTER_REBUILD_BACKGROUND_MC_H


void reload_background_MC_helper(TString sample_name, TString sample_dir)
{

    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        
        TString histname = TString(channel_histname_1D[channel]);

        TString search_object_P1 = histname + sample_name + "_P1" + "_fit";
        TString search_object_P2 = histname + sample_name + "_P2" + "_fit";

        if(rebuild_background_MC_first_run == false)
        {
            std::cout << "search_object_P1=" << search_object_P1 << std::endl;
            std::cout << "search_object_P2=" << search_object_P2 << std::endl;
        }

        TObject *tmpobj_P1 = allMCSamples1D[channel]->FindObject(search_object_P1);
        TObject *tmpobj_P2 = allMCSamples1D[channel]->FindObject(search_object_P2);

        if(tmpobj_P1 == nullptr)
        {
            if(rebuild_background_MC_first_run == false)
            {
                std::cout << "could not find object: " << search_object_P1 << std::endl;
            }
        }
        else
        {

        }

        if(tmpobj_P2 == nullptr)
        {
            if(rebuild_background_MC_first_run == false)
            {
                std::cout << "could not find object: " << search_object_P2 << std::endl;
            }
        }
        else
        {

        }

        if(allMCSamples1D[channel]->Remove(tmpobj_P1) == nullptr)
        {
            if(rebuild_background_MC_first_run == false)
            {
                std::cout << "FAILED TO REMOVE" << std::endl;
            }
        }
        allMCSamples1D[channel]->Compress();

        if(allMCSamples1D[channel]->Remove(tmpobj_P2) == nullptr)
        {
            if(rebuild_background_MC_first_run == false)
            {
                std::cout << "FAILED TO REMOVE" << std::endl;
            }
        }
        allMCSamples1D[channel]->Compress();
    }

    for(int channel = 0; channel < number2DHists; ++ channel)
    {
        
        TString histname = TString(channel_histname_2D[channel]);

        TString search_object_P1 = histname + sample_name + "_P1" + "_fit";
        TString search_object_P2 = histname + sample_name + "_P2" + "_fit";

        if(rebuild_background_MC_first_run == false)
        {
            std::cout << "search_object_P1=" << search_object_P1 << std::endl;
            std::cout << "search_object_P2=" << search_object_P2 << std::endl;
        }

        TObject *tmpobj_P1 = allMCSamples2D[channel]->FindObject(search_object_P1);
        TObject *tmpobj_P2 = allMCSamples2D[channel]->FindObject(search_object_P2);

        if(tmpobj_P1 == nullptr)
        {
            if(rebuild_background_MC_first_run == false)
            {
                std::cout << "could not find object: " << search_object_P1 << std::endl;
            }
        }
        else
        {

        }

        if(tmpobj_P2 == nullptr)
        {
            if(rebuild_background_MC_first_run == false)
            {
                std::cout << "could not find object: " << search_object_P2 << std::endl;
            }
        }
        else
        {

        }

        if(allMCSamples2D[channel]->Remove(tmpobj_P1) == nullptr)
        {
            if(rebuild_background_MC_first_run == false)
            {
                std::cout << "FAILED TO REMOVE" << std::endl;
            }
        }
        allMCSamples2D[channel]->Compress();

        if(allMCSamples2D[channel]->Remove(tmpobj_P2) == nullptr)
        {
            if(rebuild_background_MC_first_run == false)
            {
                std::cout << "FAILED TO REMOVE" << std::endl;
            }
        }
        allMCSamples2D[channel]->Compress();
    }


    ///////////////////////////////////////////////////////////////////////
    // start load from file
    ///////////////////////////////////////////////////////////////////////

    TString tinput_filename = "Nd150_2eNg_output_postprocessed_small.root";

    const TString sampleName = sample_name;
    const TString name_append = "";

    // P1
    deleter_helper("hTotalE_" + sampleName + name_append + "_P1" + "_fit");
    deleter_helper("hSingleEnergy_" + sampleName + name_append + "_P1" + "_fit");
    deleter_helper("hHighEnergy_" + sampleName + name_append + "_P1" + "_fit");
    deleter_helper("hLowEnergy_" + sampleName + name_append + "_P1" + "_fit");
    deleter_helper("hEnergySum_" + sampleName + name_append + "_P1" + "_fit");
    deleter_helper("hEnergyDiff_" + sampleName + name_append + "_P1" + "_fit");
    deleter_helper("hHighLowEnergy_" + sampleName + name_append + "_P1" + "_fit");

    // P2
    deleter_helper("hTotalE_" + sampleName + name_append + "_P2" + "_fit");
    deleter_helper("hSingleEnergy_" + sampleName + name_append + "_P2" + "_fit");
    deleter_helper("hHighEnergy_" + sampleName + name_append + "_P2" + "_fit");
    deleter_helper("hLowEnergy_" + sampleName + name_append + "_P2" + "_fit");
    deleter_helper("hEnergySum_" + sampleName + name_append + "_P2" + "_fit");
    deleter_helper("hEnergyDiff_" + sampleName + name_append + "_P2" + "_fit");
    deleter_helper("hHighLowEnergy_" + sampleName + name_append + "_P2" + "_fit");


    ///////////////////////////////////////////////////////////////////////
    // create new histograms
    ///////////////////////////////////////////////////////////////////////

    // P1
    TH1D *hTotalE_P1;
    TH1D *hSingleEnergy_P1;
    TH1D *hHighEnergy_P1;
    TH1D *hLowEnergy_P1;
    TH1D *hEnergySum_P1;
    TH1D *hEnergyDiff_P1;
    TH2D *hHighLowEnergy_P1;
    // P2
    TH1D *hTotalE_P2;
    TH1D *hSingleEnergy_P2;
    TH1D *hHighEnergy_P2;
    TH1D *hLowEnergy_P2;
    TH1D *hEnergySum_P2;
    TH1D *hEnergyDiff_P2;
    TH2D *hHighLowEnergy_P2;

    ///////////////////////////////////////////////////////////////////
    // phase 1 histograms (output)
    ///////////////////////////////////////////////////////////////////

    hTotalE_P1          = new TH1D("hTotalE_" + sampleName + name_append + "_P1" + "_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                                50, 0.0, 5.0);

    hSingleEnergy_P1    = new TH1D("hSingleEnergy_" + sampleName + name_append + "_P1" + "_fit",
                                TString("Phase ") + "P1" + " " + sampleName  + name_append + " Single Energy",
                                50, 0.0, 5.0);
    
    //std::cout << "High new" << std::endl;
    hHighEnergy_P1     = new TH1D("hHighEnergy_" + sampleName + name_append + "_P1" + "_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " High Energy; Energy (MeV)",
                                50, 0.0, 5.0);

    //std::cout << "Low new" << std::endl;
    hLowEnergy_P1     = new TH1D("hLowEnergy_" + sampleName + name_append + "_P1" + "_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Sum new" << std::endl;
    hEnergySum_P1     = new TH1D("hEnergySum_" + sampleName + name_append + "_P1" + "_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Diff new" << std::endl;
    hEnergyDiff_P1     = new TH1D("hEnergyDiff_" + sampleName + name_append + "_P1" + "_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "HighLow new" << std::endl;
    hHighLowEnergy_P1     = new TH2D("hHighLowEnergy_" + sampleName + name_append + "_P1" + "_fit",
                                TString("Phase ") + "P1" + " " + sampleName + name_append + ";Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                12, 0.0, 3.6, 12, 0.0, 3.6);


    ///////////////////////////////////////////////////////////////////
    // phase 2 histograms (output)
    ///////////////////////////////////////////////////////////////////

    hTotalE_P2           = new TH1D("hTotalE_" + sampleName + name_append + "_P2" + "_fit",
                                //"Phase " + Phase + " " + sampleName + name_append + " total energy; #SigmaE_{e} (MeV)",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " total energy; Total Energy #SigmaE_{e} (MeV)",
                                50, 0.0, 5.0);

    //std::cout << "Single new" << std::endl;
    hSingleEnergy_P2    = new TH1D("hSingleEnergy_" + sampleName + name_append + "_P2" + "_fit",
                                TString("Phase ") + "P2" + " " + sampleName  + name_append + " Single Energy",
                                50, 0.0, 5.0);
    
    //std::cout << "High new" << std::endl;
    hHighEnergy_P2     = new TH1D("hHighEnergy_" + sampleName + name_append + "_P2" + "_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " High Energy; Energy (MeV)",
                                50, 0.0, 5.0);

    //std::cout << "Low new" << std::endl;
    hLowEnergy_P2     = new TH1D("hLowEnergy_" + sampleName + name_append + "_P2" + "_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Sum new" << std::endl;
    hEnergySum_P2     = new TH1D("hEnergySum_" + sampleName + name_append + "_P2" + "_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "Diff new" << std::endl;
    hEnergyDiff_P2     = new TH1D("hEnergyDiff_" + sampleName + name_append + "_P2" + "_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + " Low Energy",
                                50, 0.0, 5.0);

    //std::cout << "HighLow new" << std::endl;
    hHighLowEnergy_P2     = new TH2D("hHighLowEnergy_" + sampleName + name_append + "_P2" + "_fit",
                                TString("Phase ") + "P2" + " " + sampleName + name_append + ";Low Energy Electron Energy (MeV);High Energy Electron Energy (MeV)",
                                12, 0.0, 3.6, 12, 0.0, 3.6);


    
    hTotalE_P1->Sumw2();
    hSingleEnergy_P1->Sumw2();
    hHighEnergy_P1->Sumw2();
    hLowEnergy_P1->Sumw2();
    hHighLowEnergy_P1->Sumw2();
    hEnergySum_P1->Sumw2();
    hEnergyDiff_P1->Sumw2();


    hTotalE_P2->Sumw2();
    hSingleEnergy_P2->Sumw2();
    hHighEnergy_P2->Sumw2();
    hLowEnergy_P2->Sumw2();
    hHighLowEnergy_P2->Sumw2();
    hEnergySum_P2->Sumw2();
    hEnergyDiff_P2->Sumw2();
    


    ///////////////////////////////////////////////////////////////////
    // load data from file
    ///////////////////////////////////////////////////////////////////

    TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    TString thePath = sample_dir; //MCNameToPathMap.at(sampleName);

    TFile *inputFile = new TFile(filePath + thePath + sampleName + "/Nd150_2eNg_output_postprocessed_small.root");
    if(!inputFile->IsOpen())
    {
        std::cout << "failed to open" << std::endl;
    }
    TTree *inputTree = (TTree*)inputFile->Get("Nd150_2eNg/Nd150_2eNg");

    int run;
    double eventTime;
    double radonWeight;
    double bi210Weight;
    double electronEnergy[2];
    double trueElectronEnergy[2];
    double gauss[2];

    inputTree->SetBranchAddress("Run", &run);
    inputTree->SetBranchAddress("eventTime", &eventTime);
    inputTree->SetBranchAddress("radonWeight", &radonWeight);
    inputTree->SetBranchAddress("bi210Weight", &bi210Weight);
    inputTree->SetBranchAddress("electronEnergy", electronEnergy);
    inputTree->SetBranchAddress("trueElectronEnergy", trueElectronEnergy);
    inputTree->SetBranchAddress("gauss", gauss);

    double el_energy_0;
    double el_energy_1;

    double weight_total_150Nd = 0.0;
    double event_pass_count_150Nd = 0.0;

    Long_t events = (Long_t)inputTree->GetEntries();
    for(Long_t event_i = 0; event_i < events; ++ event_i)
    {
    
        inputTree->GetEvent(event_i);

        
        ///////////////////////////////////////////////////////////////
        // see file: fit_2e.C for documentation
        ///////////////////////////////////////////////////////////////
        double weight = 1.0;
        if(sampleName.Contains("bi214_swire"))
        {
            weight = radonWeight;
        }
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



        if(sampleName.Contains("nd150"))
        {
            const Double_t trueT0{trueElectronEnergy[highE_index]};
            const Double_t trueT1{trueElectronEnergy[lowE_index]};

//                Double_t _weight{ReWeight3(trueT0, trueT1, xi_31_baseline, xi_31, h_nEqNull, h_nEqTwo, psiN0, psiN2, "true")}; // TODO remove true?

            /*
            if(_weight != 1.0)
            {
                std::cout << "_weight=" << _weight << std::endl;
            }
            */

            const double _weight = 1.0; // HSD

            weight *= _weight;

            //weight_total_150Nd += _weight;
        }


        // set phase
        // 1869 - 3395 is P1
        // 3396 - 9186 is P2
        bool cut = true;
        Int_t thePhase = -1;
        if((1869 <= run) && (run <= 3395))
        {
            cut = false;
            thePhase = 0;
        }
        else if((3396 <= run) && (run <= 9186))
        {
            cut = false;
            thePhase = 1;
        }
        if(cut == true) continue;


        ///////////////////////////////////////////////////////////////
        // apply systematics
        ///////////////////////////////////////////////////////////////

        Double_t el_energy_0{electronEnergy[highE_index]};
        Double_t el_energy_1{electronEnergy[lowE_index]};

        // optical correction
        if(OPTICAL_CORRECTION_ENABLE)
        {
            double cf0 = gSystematics.g_systematic_optical_correction_n->Eval(el_energy_0);
            double cf1 = gSystematics.g_systematic_optical_correction_n->Eval(el_energy_1);
            el_energy_0 *= (cf0);
            el_energy_1 *= (cf1); 
        }


        ///////////////////////////////////////////////////////////////
        // standard cuts not implemented in preprocessing
        ///////////////////////////////////////////////////////////////

        if((el_energy_0 < 0.3) || (el_energy_1 < 0.3))
        {
            continue;
        }

        ///////////////////////////////////////////////////////////////
        // fill histograms
        ///////////////////////////////////////////////////////////////

        if(thePhase == 0)
        {
            hTotalE_P1->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hSingleEnergy_P1->Fill(el_energy_1, 1.0 * weight);
            hSingleEnergy_P1->Fill(el_energy_0, 1.0 * weight);
            hHighEnergy_P1->Fill(el_energy_0, 1.0 * weight);
            hLowEnergy_P1->Fill(el_energy_1, 1.0 * weight);
            hHighLowEnergy_P1->Fill(el_energy_1, el_energy_0, 1.0 * weight);
            hEnergySum_P1->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hEnergyDiff_P1->Fill(el_energy_0 - el_energy_1, 1.0 * weight);
        }
        else if(thePhase == 1)
        {
            hTotalE_P2->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hSingleEnergy_P2->Fill(el_energy_0, 1.0 * weight);
            hSingleEnergy_P2->Fill(el_energy_1, 1.0 * weight);
            hHighEnergy_P2->Fill(el_energy_0, 1.0 * weight);
            hLowEnergy_P2->Fill(el_energy_1, 1.0 * weight);
            hHighLowEnergy_P2->Fill(el_energy_1, el_energy_0, 1.0 * weight);
            hEnergySum_P2->Fill(el_energy_0 + el_energy_1, 1.0 * weight);
            hEnergyDiff_P2->Fill(el_energy_0 - el_energy_1, 1.0 * weight);
        }
        else
        {
            std::cout << "Error: thePhase=" << thePhase << " invalid value" << std::endl;
        }


    }


    inputFile->Close();


    
    // scale
    std::ifstream inFile;
    const TString typedir = thePath;
    const TString sample_file_name = sampleName;
    inFile.open(filePath + typedir + sample_file_name + "/JobSummary.txt");
    //std::cout << "file: " << filePath + typedir + sample_file_name + "/JobSummary.txt" << std::endl;
    if(!inFile.is_open())
    {
        std::cout << "file not open" << std::endl;
        std::cin.get();
    }
    std::string dummy;
    double sampleNGenMC;
    inFile >> dummy >> sampleNGenMC;
    inFile.close();

    const double TotalTime = 167629292.;

    // P1
    hTotalE_P1->Scale(TotalTime / sampleNGenMC);
    hSingleEnergy_P1->Scale(TotalTime / sampleNGenMC);
    hHighEnergy_P1->Scale(TotalTime / sampleNGenMC);
    hLowEnergy_P1->Scale(TotalTime / sampleNGenMC);
    hHighLowEnergy_P1->Scale(TotalTime / sampleNGenMC);
    hEnergySum_P1->Scale(TotalTime / sampleNGenMC);
    hEnergyDiff_P1->Scale(TotalTime / sampleNGenMC);

    // P2
    hTotalE_P2->Scale(TotalTime / sampleNGenMC);
    hSingleEnergy_P2->Scale(TotalTime / sampleNGenMC);
    hHighEnergy_P2->Scale(TotalTime / sampleNGenMC);
    hLowEnergy_P2->Scale(TotalTime / sampleNGenMC);
    hHighLowEnergy_P2->Scale(TotalTime / sampleNGenMC);
    hEnergySum_P2->Scale(TotalTime / sampleNGenMC);
    hEnergyDiff_P2->Scale(TotalTime / sampleNGenMC);


    
    if(sampleName.Contains("tl208"))
    {
        if(!sampleName.Contains("air"))
        {
            const double scale_factor = 0.36;

            // P1
            hTotalE_P1->Scale(scale_factor);
            hSingleEnergy_P1->Scale(scale_factor);
            hHighEnergy_P1->Scale(scale_factor);
            hLowEnergy_P1->Scale(scale_factor);
            hHighLowEnergy_P1->Scale(scale_factor);
            hEnergySum_P1->Scale(scale_factor);
            hEnergyDiff_P1->Scale(scale_factor);

            // P2
            hTotalE_P2->Scale(scale_factor);
            hSingleEnergy_P2->Scale(scale_factor);
            hHighEnergy_P2->Scale(scale_factor);
            hLowEnergy_P2->Scale(scale_factor);
            hHighLowEnergy_P2->Scale(scale_factor);
            hEnergySum_P2->Scale(scale_factor);
            hEnergyDiff_P2->Scale(scale_factor);
        }
    }

    if(sampleName.Contains("bi207"))
    {
        const double scale_factor = 1.846;
        
        // P1
        hTotalE_P1->Scale(scale_factor);
        hSingleEnergy_P1->Scale(scale_factor);
        hHighEnergy_P1->Scale(scale_factor);
        hLowEnergy_P1->Scale(scale_factor);
        hHighLowEnergy_P1->Scale(scale_factor);
        hEnergySum_P1->Scale(scale_factor);
        hEnergyDiff_P1->Scale(scale_factor);

        // P2
        hTotalE_P2->Scale(scale_factor);
        hSingleEnergy_P2->Scale(scale_factor);
        hHighEnergy_P2->Scale(scale_factor);
        hLowEnergy_P2->Scale(scale_factor);
        hHighLowEnergy_P2->Scale(scale_factor);
        hEnergySum_P2->Scale(scale_factor);
        hEnergyDiff_P2->Scale(scale_factor);
    }


    double scale_factor_P1 = 0.0;
    double scale_factor_P2 = 0.0;
    int param_number_P1 = 0;
    int param_number_P2 = 0;

    std::string mc_name = std::string(sample_name);

    bool success = g_pg.convert_MC_name_to_param_number(mc_name, param_number_P1, param_number_P2);
    if(success == true)
    {
        if(param_number_P1 != -1)
        {
            scale_factor_P1 = g_pg.file_params.at(param_number_P1).paramInitValue;
        }
        else
        { 
            std::cout << "WARNING: param_number_P1=" << param_number_P1 << std::endl;
            std::cout << "mc_name=" << mc_name << std::endl;
        }

        if(param_number_P2 != -1)
        {
            scale_factor_P2 = g_pg.file_params.at(param_number_P2).paramInitValue;
        }
        else
        { 
            std::cout << "WARNING: param_number_P2=" << param_number_P2 << std::endl;
            std::cout << "mc_name=" << mc_name << std::endl;
        }
    }

    if((param_number_P1 == 1) || (param_number_P2 == 1))
    {
        std::cout << "ERROR" << std::endl;
        throw "error";
    }


    // P1

    // scale to get default number of events
    hTotalE_P1->Scale(scale_factor_P1);
    hSingleEnergy_P1->Scale(scale_factor_P1);
    hHighEnergy_P1->Scale(scale_factor_P1);
    hLowEnergy_P1->Scale(scale_factor_P1);
    hHighLowEnergy_P1->Scale(scale_factor_P1);
    hEnergySum_P1->Scale(scale_factor_P1);
    hEnergyDiff_P1->Scale(scale_factor_P1);
    
    // P2
    
    // scale to get default number of events
    hTotalE_P2->Scale(scale_factor_P2);
    hSingleEnergy_P2->Scale(scale_factor_P2);
    hHighEnergy_P2->Scale(scale_factor_P2);
    hLowEnergy_P2->Scale(scale_factor_P2);
    hHighLowEnergy_P2->Scale(scale_factor_P2);
    hEnergySum_P2->Scale(scale_factor_P2);
    hEnergyDiff_P2->Scale(scale_factor_P2);

    // P1
    allMCSamples1D[0]->Add(hTotalE_P1);
    allMCSamples1D[1]->Add(hSingleEnergy_P1);
    allMCSamples1D[2]->Add(hHighEnergy_P1);
    allMCSamples1D[3]->Add(hLowEnergy_P1);
    allMCSamples1D[4]->Add(hEnergySum_P1);
    allMCSamples1D[5]->Add(hEnergyDiff_P1);

    allMCSamples2D[0]->Add(hHighLowEnergy_P1);

    // P2
    allMCSamples1D[0]->Add(hTotalE_P2);
    allMCSamples1D[1]->Add(hSingleEnergy_P2);
    allMCSamples1D[2]->Add(hHighEnergy_P2);
    allMCSamples1D[3]->Add(hLowEnergy_P2);
    allMCSamples1D[4]->Add(hEnergySum_P2);
    allMCSamples1D[5]->Add(hEnergyDiff_P2);

    allMCSamples2D[0]->Add(hHighLowEnergy_P2);


    


}


void reload_background_MC()
{

    ///////////////////////////////////////////////////////////////////////////
    // init map
    // convert mc name to path
    ///////////////////////////////////////////////////////////////////////////

    std::map<TString, TString> MCNameToPathMap;
    for(Int_t index = 0; index < nRn222Bkgs; ++ index)
    {
        TString mcname = Rn222BkgFilesNew[index];
        TString mcpath = "externals/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }
    for(Int_t index = 0; index < nRn220Bkgs; ++ index)
    {
        TString mcname = Rn220BkgFiles[index];
        TString mcpath = "externals/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }
    for(Int_t index = 0; index < nExternalBkgs; ++ index)
    {
        TString mcname = ExternalBkgFiles[index];
        TString mcpath = "externals/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }
    for(Int_t index = 0; index < nInternalBkgs; ++ index)
    {
        TString mcname = InternalBkgFiles[index];
        TString mcpath = "internals/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }
    for(Int_t index = 0; index < nNeighbours; ++ index)
    {
        TString mcname = NeighbourFiles[index];
        TString mcpath = "neighbours/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }
    for(Int_t index = 0; index < nNd150Samples; ++ index)
    {
        TString mcname = Nd150Files[index];
        TString mcpath = "nd150/";
        MCNameToPathMap.insert(std::make_pair(mcname, mcpath));
    }


    std::vector<TString> bkg_sample_name;
    std::vector<TString> bkg_sample_dir;


    for(int i = 0; i < nExternalBkgs; ++ i)
    {
        TString sample_name = ExternalBkgFiles[i];
        bkg_sample_name.push_back(sample_name);
        bkg_sample_dir.push_back(MCNameToPathMap.at(sample_name));
    }


    for(int i = 0; i < nInternalBkgs; ++ i)
    {
        TString sample_name = InternalBkgFiles[i];
        bkg_sample_name.push_back(sample_name);
        bkg_sample_dir.push_back(MCNameToPathMap.at(sample_name));
    }


    for(int i = 0; i < nRn222Bkgs; ++ i)
    {
        TString sample_name = Rn222BkgFilesNew[i];
        bkg_sample_name.push_back(sample_name);
        bkg_sample_dir.push_back(MCNameToPathMap.at(sample_name));
    }


    for(int i = 0; i < nRn220Bkgs; ++ i)
    {
        TString sample_name = Rn220BkgFiles[i];
        bkg_sample_name.push_back(sample_name);
        bkg_sample_dir.push_back(MCNameToPathMap.at(sample_name));
    }


    for(int i = 0; i < nNeighbours; ++ i)
    {
        TString sample_name = NeighbourFiles[i];
        bkg_sample_name.push_back(sample_name);
        bkg_sample_dir.push_back(MCNameToPathMap.at(sample_name));
    }


    // ! should not call more than once !
    for(int i = 0; i < nNd150Samples; ++ i)
    {
        TString sample_name = Nd150Files[i];
        bkg_sample_name.push_back(sample_name);
        bkg_sample_dir.push_back(MCNameToPathMap.at(sample_name));
    }

  
    for(int i = 0; i < bkg_sample_name.size(); ++ i)
    {
        TString sample_name = bkg_sample_name.at(i);
        TString sample_dir = bkg_sample_dir.at(i);

        reload_background_MC_helper(sample_name, sample_dir);
    }



    TString thePhase_arg_P1 = "P1";
    TString thePhase_arg_P2 = "P2";
    TString theChannel = "2e_";
    TFile *aFileP1 = new TFile("Nd150_" + theChannel + thePhase_arg_P1 + ".root");
    TFile *aFileP2 = new TFile("Nd150_" + theChannel + thePhase_arg_P2 + ".root");


    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        TString theHistogram = channel_histname_1D[channel];

        std::string directory("processeddata/" + theHistogram + "/");
        std::string name(theHistogram + DataFile);
        std::string fullname = directory + name;

        std::string new_name_P1(theHistogram + DataFile + "_P1");
        std::string new_name_P2(theHistogram + DataFile + "_P2");
        TH1D *tmpHist_P1 = (TH1D*)aFileP1->Get(fullname.c_str())->Clone(new_name_P1.c_str());
        TH1D *tmpHist_P2 = (TH1D*)aFileP2->Get(fullname.c_str())->Clone(new_name_P2.c_str());
        if(tmpHist_P1 != nullptr)
        {
            allDataSamples1D->Add((TH1D*)tmpHist_P1);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << fullname << std::endl;
            std::cout << __func__ << " could not find histogram in file: " << fullname << " - disabling parameter number TODO fix this error" << std::endl;
        }

        if(tmpHist_P2 != nullptr)
        {
            allDataSamples1D->Add((TH1D*)tmpHist_P2);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << fullname << std::endl;
            std::cout << __func__ << " could not find histogram in file: " << fullname << " - disabling parameter number TODO fix this error" << std::endl;
        }
    }

    for(int channel = 0; channel < number2DHists; ++ channel)
    {
        TString theHistogram = channel_histname_2D[channel];

        std::string directory("processeddata/" + theHistogram + "/");
        std::string name(theHistogram + DataFile);
        std::string fullname = directory + name;

        std::string new_name_P1(theHistogram + DataFile + "_P1");
        std::string new_name_P2(theHistogram + DataFile + "_P2");
        TH1D *tmpHist_P1 = (TH1D*)aFileP1->Get(fullname.c_str())->Clone(new_name_P1.c_str());
        TH1D *tmpHist_P2 = (TH1D*)aFileP2->Get(fullname.c_str())->Clone(new_name_P2.c_str());
        if(tmpHist_P1 != nullptr)
        {
            allDataSamples2D->Add((TH1D*)tmpHist_P1);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << fullname << std::endl;
            std::cout << __func__ << " could not find histogram in file: " << fullname << " - disabling parameter number TODO fix this error" << std::endl;
        }

        if(tmpHist_P2 != nullptr)
        {
            allDataSamples2D->Add((TH1D*)tmpHist_P2);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << fullname << std::endl;
            std::cout << __func__ << " could not find histogram in file: " << fullname << " - disabling parameter number TODO fix this error" << std::endl;
        }
    }

    aFileP1->Close();
    delete aFileP1;
    aFileP2->Close();
    delete aFileP2;

    
    if(rebuild_background_MC_first_run == true)
    {
        rebuild_background_MC_first_run = false;
    }


}


#endif // NEWLOGLIKFITTER_REBUILD_BACKGROUND_MC_H
