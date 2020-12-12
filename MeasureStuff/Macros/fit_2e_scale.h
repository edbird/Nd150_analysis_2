#ifndef FIT_2E_SCALE_H
#define FIT_2E_SCALE_H


void scale(TFile* myFile,                       // INPUT: unscaled histograms are loaded from this root file
           const std::string& activityfile,     // INPUT: activity values read from this file
           const std::string& typedir,          // INPUT: forms part of path string,
           // eg: "nd150/", "internals/", "externals/", "neighbours/",
           // all radon files are in "externals/"
           const int nSamples,                     // INPUT: number of objects in dataFiles/dataNames
           const TString *const sampleFiles,      // INPUT: list of samples, eg: bi214_sfoil_rot
           const TString *const sampleNames,      // INPUT: corresponding human readable names, Latex format, eg: ^{214}Bi Foil Surface
           const Color_t *const sampleColors,     // INPUT: sample color set using values in this array
           const double *const AcceptedTime,    // INPUT: accepted time for each phase is passed in here 
           double *const sampleNGenMC,            // OUTPUT: number of generated MC events are written to this location
           double *const sampleActivity,          // OUTPUT: activity values written to this location
           double *const sampleEfficiency,        // OUTPUT: efficiency values written to this location
           TObjArray** allSamples,              // OUTPUT: each histogram has colors set and is scaled (using scaling with activity), before being stored as an output here
           const Int_t mode_flag = 0)
{
    std::cout << "scale" << std::endl;

    std::ofstream oftemp("oftemp1.txt", std::ofstream::out | std::ofstream::app);




    TH1D *tmpHist; //we'll be using this later.

    //Let's get the internals/externals/data and add them to an array
    //TObjArray *allSamples[nSamples];
    for(int i = 0; i < nSamples; i++)
    {

        // create new object array
        //std::cout << sampleFiles[i] << std::endl;
        allSamples[i] = new TObjArray();


        // read efficiency from file
        std::ifstream ifile_cutcount("cutcount.txt");
        if(!ifile_cutcount.is_open())
        {
            std::cout << "failed to open file " << "cutcount.txt" << std::endl;
        }

        //double sample_total_event_count;
        double sample_pass_count;
        double sample_efficiency;
        while(!ifile_cutcount.eof())
        {
            std::string sample_name;
            Long64_t total_events;
            Long64_t pass_count;

            ifile_cutcount >> sample_name >> total_events >> pass_count;

            // TODO: removed std::string conversion?
            if(sample_name == std::string(sampleFiles[i]))
            {
                //sample_total_event_count = (double)total_events;
                sample_pass_count = (double)pass_count;
                //sample_efficiency = (double)pass_count / (double)total_events;
            }
        }


        // read in the information needed for calculating efficiencies
        // read in number of generated MC events
        std::ifstream inFile;
        std::string sample_file_name = sampleFiles[i].Data();
        static const std::vector<std::string> search_objects = {"_Lg0", "_INL0", "_OUTL0", "_IN", "_OUT"};
        for(int i = 0; i < search_objects.size(); ++ i)
        {
            std::string search_object = search_objects.at(i);
            if(sample_file_name.find(search_object) != std::string::npos)
            {
                //std::cout << "sample_file_name=" << sample_file_name << std::endl;
                sample_file_name = sample_file_name.substr(0, sample_file_name.size() - search_object.size());
                //std::cout << "sample_file_name=" << sample_file_name << std::endl;
                break;
            }
        }
        inFile.open(filePath + typedir + sample_file_name + "/JobSummary.txt");
        if(!inFile.is_open())
        {
            std::cout << "failed to open " << filePath + typedir + sample_file_name + "/JobSummary.txt" << std::endl;
        }
        // cout << "JobSummary stuff is here: "
        // 	 << filePath << "internals/" << InternalBkgFiles[i] << "/JobSummary.txt"
        // 	 << endl;
        std::string dummy;
        inFile >> dummy >> sampleNGenMC[i];
        // cout << "Read in as dummy: " << dummy << endl;
        // cout << "Read in as InternalBkgNGenMC[i]: " << InternalBkgNGenMC[i] << endl;
        inFile.close();

        // if we want to use a pre-measured background activity (which we do for the time being) then we
        // have to scale the histogram appropriately.  So now we need to fetch the activity from the txt file
        std::ifstream inFile2;
        inFile2.open(activityfile.c_str());
        if(!inFile2.is_open())
        {
            std::cout << "failed to open " << activityfile << std::endl;
        }

        bool found = false;
        while((!found) && (!inFile2.eof()))
        {
            // phase 1 and phase 2 activity
            double activityPhase1, activityPhase2; 
            // name of MC sample
            TString mcname;
            
            std::stringstream ss;
            std::string s;
            std::getline(inFile2, s);

            if(s.size() > 0)
            {
                ss << s;
                ss >> mcname >> activityPhase1 >> activityPhase2;
                //std::cout << mcname << ", " << activityPhase1 << ", " << activityPhase2 << std::endl;
                
                if(mcname.CompareTo(sampleFiles[i]) == 0)
                {
                    sampleActivity[0 * nSamples + i] = activityPhase1;
                    sampleActivity[1 * nSamples + i] = activityPhase2;
                    found = true;
                }

            }
        }
        if(!found)
        {
            std::cout << "Fault: sampleFiles[" << i << "]=" << sampleFiles[i] << std::endl;
            std::cin.get();
        }

        inFile2.close();

        for(int j = 0; j < numHistograms; j++)
        {

            /*
            TString name_append;
            if(mode_flag == 0)
            {
                //name_append = "";
            }
            else if(mode_flag == 1)
            {
                name_append = "_raw";
            }
            */
            // don't need

            TString histogram_name = histogramNames[j];
            TString new_histogram_name = histogram_name + sampleFiles[i] + "_fit";
            //TString name = "rawdata/" + histogram_name + "/" + histogram_name + sampleFiles[i];

            TString name;
            if(mode_flag == 0)
            {
                /*TString*/ name = "processeddata/" + histogram_name + "/" + histogram_name + sampleFiles[i]; // + name_append
                if(histogram_name.CompareTo("hTotalE_") == 0)
                {
                    std::cout << "Get() : " << name << " from file, Clone() : " << new_histogram_name << std::endl;
                }
                //tmpHist = (TH1D*)myFile->Get(name)->Clone(new_histogram_name);
            }
            else if(mode_flag == 1)
            {
                /*TString*/ name = "rawdata/" + histogram_name + "/" + histogram_name + sampleFiles[i] + "_raw"; // + name_append;
                if(histogram_name.CompareTo("hTotalE_") == 0)
                {
                    std::cout << "Get() : " << name << " from file, Clone() : " << new_histogram_name << std::endl;
                }
                //tmpHist = (TH1D*)myFile->Get(name)->Clone(new_histogram_name);
            }
            try
            {
                tmpHist = (TH1D*)myFile->Get(name)->Clone(new_histogram_name);
            }
            catch(std::exception &e)
            {
                std::cout << "name=" << name << std::endl;
                std::cout << e.what() << std::endl;
                throw e;
            }
            // NOTE: will have 2 new histograms here
            // 1: scaled, processeddata histogram (save this one to file)
            // 2: scaled, rawdata histogram (do not save this one to file)

            // this code was replaced with above: 2020-05-03
            //std::string name("singleHistos/unscaled/" + histogramNames[j] + sampleFiles[i]);
            //std::cout << "Get() : " << name << " from file, Clone() : " << histogramNames[j] + sampleFiles[i] + "_fit" << std::endl;
            //tmpHist = (TH1D*)myFile->Get(name.c_str())->Clone(histogramNames[j] + sampleFiles[i] + "_fit");
            
            //calculate efficiency
            // TODO: this is calculated using the 0th histogram, should perhaps
            // change to the hTotalE ?
            /*
            if(j == 0)
            {
                double NPass = tmpHist->Integral();
                // NPass is phase specific
                // sampleNGenMC[i] is for both phases
                std::cout << "NPass: " << NPass << std::endl;
                std::cout << "sampleNGenMC[i=" << i << "]=" << sampleNGenMC[i] << std::endl;

                sample_efficiency = sample_pass_count / sampleNGenMC[i];
                
                // TODO: don't fully understand this
                // efficiency = (number of events which pass cuts * total time) / 
                // (number of generated events * accepted time)
                // accepted time is like total time - dead time?
                //double eff = NPass / (double)(sampleNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime) );
                // this efficiency is phase specific
                double eff = (double)NPass / (double)sampleNGenMC[i];
                //sampleEfficiency[thePhase][i] = eff;
            }
            */
            
            //tmpHist->Sumw2();
            tmpHist->SetFillColor(sampleColors[i]);
            tmpHist->SetLineColor(sampleColors[i]);
            tmpHist->SetLineWidth(0);
            tmpHist->SetTitle(sampleNames[i]);

            //tmpHist->Scale(sampleActivity[thePhase][i] * AcceptedTime[thePhase] * 0.001 / (sampleNGenMC[i] * (AcceptedTime[thePhase]/(double)TotalTime)) );
            //tmpHist->Scale(sampleActivity[thePhase * nSamples + i] * 0.001 * (double)TotalTime / sampleNGenMC[i] );
            //std::cout << "scaling, index=" << thePhase + nSamples + i << std::endl;
            
            // In Position Paper
            // Section 3.4.1, Equation 3.7
            // A = N / (t * eps)
            // A: activity, N: number of events, t: time, eps: efficiency
            // Can calculate number of events backwards from activity, efficiency and time

            
            

                // TODO: working here, the activity can change between phases
                // not sure what the best way is to work out the number
                // of expected events on each phase and compare this to
                // the histogram integrals
                // also the number of events appears to be different depending
                // on whether hTotalE is used or whether hRun is used
                // does hRun contain events for both phases?
                // what are the limits on hRun? how is it filled?
                // copied from filling code:
                //hRun->Fill(run, weight);
                //hTotalE->Fill(electronEnergy[0] + electronEnergy[1], weight);

            
            // moved from stackFunction()
            // 2020-03-30: 
            //(TH1D*)allSamples[j]->At(i)->Write();
            // saving the histogram without scaling it by the activity will
            // cause the fit output (and input) to be in units of activty
            // is this actually the case, what about the factor of the detector
            // livetime which I have left in above?
            /*
            Double_t scale_factor_without_activity = (double)TotalTime / sampleNGenMC[i]; //scale_factor / the_activity;
            tmpHist->Scale(scale_factor_without_activity);
            std::cout << "writing: " << tmpHist->GetName() << std::endl;
            tmpHist->Write();

            tmpHist->Scale(scale_factor / scale_factor_without_activity);
            */



            if(TString(tmpHist->GetName()).Contains("hTotalE"))
            {
                sampleEfficiency[thePhase * nSamples + i] = tmpHist->Integral() / sampleNGenMC[i];
                std::ofstream of_sampleEfficiency("of_sampleEfficiency.txt", std::ios::app);
                of_sampleEfficiency << tmpHist->GetName() << " " << tmpHist->Integral() << " " << sampleNGenMC[i] << " " << sampleEfficiency[0 * nSamples + i] << " " << sampleEfficiency[1 * nSamples + i] << std::endl;
            }

            if(TString(tmpHist->GetName()).Contains("hTotalE"))
            {
                std::cout << "TotalTime=" << TotalTime
                          << ", sampleNGenMC[" << i << "]=" << sampleNGenMC[i]
                          << ", sampleActivity[thePhase="
                              << thePhase
                              << " * nSamples="
                              << nSamples
                              << " + i="
                              << i << "]="
                          << sampleActivity[thePhase * nSamples + i] << std::endl;
                //std::cin.get();
                std::cout << std::endl;
            }



            tmpHist->Scale(TotalTime / sampleNGenMC[i]);

            if(mode_flag == 0)
            {
                if(TString(tmpHist->GetName()).Contains("hSingleEnergy_"))
                {
                    std::cout << tmpHist->GetName() << " -> " << tmpHist->GetBinContent(10) << std::endl;
                    //std::cin.get();
                }
                //std::cout << tmpHist->GetName() << " -> scale=" << TotalTime / sampleNGenMC[i] << std::endl;
            }
            
            if(mode_flag == 0)
            {
                // write into root directory to be read back by newLogLikFitter.C
                myFile->cd("/");
// disabled 2020-12-11                tmpHist->Write();

                // write into scaled / histogram_name subdirectory
                myFile->cd("scaled/" + histogram_name);
                TString clone_name = TString(tmpHist->GetName()) + TString("_scaled");
                TH1D *tmpHistClone = (TH1D*)tmpHist->Clone(clone_name);
                tmpHistClone->Write();
            }
            //if(sampleActivity[thePhase * nSamples + i] <= 0.0)
            //{
            //    std::cout << "fault" << std::endl;
            //    std::cout << "the name is: " << tmpHist->GetName() << std::endl;
            //    std::cout << sampleFiles[i] << std::endl;
            //    std::cin.get();
            //    //throw "fault";
            //}
            tmpHist->Scale(sampleActivity[thePhase * nSamples + i]);
            if(mode_flag == 0)
            {
                if(TString(tmpHist->GetName()).Contains("hSingleEnergy_"))
                {
                    std::cout << "after scaling by activity A=" << sampleActivity[thePhase * nSamples + i] << " the 10th bin contains:" << std::endl;
                    std::cout << tmpHist->GetName() << " -> " << tmpHist->GetBinContent(10) << std::endl;
                    //std::cin.get();
                }
                //std::cout << tmpHist->GetName() << " -> scale=" << TotalTime / sampleNGenMC[i] << std::endl;
            }
            if(mode_flag == 0)
            {
                // write into scaled2 / histogram_name subdirectory
//                myFile->cd("scaled2/" + histogram_name);
//                TString clone_name = TString(tmpHist->GetName()) + TString("_scaled2");
//                TH1D *tmpHistClone = (TH1D*)tmpHist->Clone(clone_name);
// disabled 2020-12-11                tmpHistClone->Write();
            }

            std::ofstream of_numberofevents("of_numberofevents.txt", std::ofstream::out | std::ofstream::app);
            if(mode_flag == 0 && TString(tmpHist->GetName()).Contains("hTotalE"))
            {
                std::cout << tmpHist->GetName() << " number of events " << tmpHist->Integral() << std::endl;
                of_numberofevents << tmpHist->GetName() << " number of events " << tmpHist->Integral() << std::endl;
            }
            of_numberofevents.close();

            // NOTE: done AFTER call to Write()
            // so sample is scaled with activity value if used again from allSamples
            allSamples[i]->Add(tmpHist);

        } // all histos
    }   //all of sample type (internal/external/...)

    oftemp.close();
}

#endif // FIT_2E_SCALE_H
