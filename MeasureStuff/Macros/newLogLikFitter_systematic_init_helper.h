#ifndef NEWLOGLIKFITTER_SYSTEMATIC_INIT_HELPER_H
#define NEWLOGLIKFITTER_SYSTEMATIC_INIT_HELPER_H


void systematic_n_init_helper
(
    std::vector<double> *systematic_n_highlow_1D_P1[][number1DHists],
    std::vector<double> *systematic_n_highlow_1D_P2[][number1DHists],
    const int systematic_n_index
)
{

    // loop over all channels
    for(int channel = 0; channel < number1DHists; ++ channel)
    {

        std::string histname = std::string(channel_histname_1D[channel]);
        std::string search_object_P1;
        std::string search_object_P2;

        search_object_P1 = histname + std::string("fakedata") + "_P1";
        search_object_P2 = histname + std::string("fakedata") + "_P2";

        TH1D *tmpDataHist1D_P1 = nullptr;
        TH1D *tmpDataHist1D_P2 = nullptr;
    
        tmpDataHist1D_P1 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P1.c_str());
        tmpDataHist1D_P2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());

        // phase 1
        if(tmpDataHist1D_P1 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }

        // phase 2
        if(tmpDataHist1D_P2 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
            throw "problem";
        }

        // phase 1
        for(Int_t bin_ix{1}; bin_ix <= tmpDataHist1D_P1->GetNbinsX(); ++ bin_ix)
        {
            Double_t content = tmpDataHist1D_P1->GetBinContent(bin_ix);
            //systematic_offset_high_1D_P1[channel]->push_back(content);
            systematic_n_highlow_1D_P1[systematic_n_index][channel]->push_back(content);
        }

        // phase 2
        for(Int_t bin_ix{1}; bin_ix <= tmpDataHist1D_P2->GetNbinsX(); ++ bin_ix)
        {
            Double_t content = tmpDataHist1D_P2->GetBinContent(bin_ix);
            //systematic_offset_high_1D_P2[channel]->push_back(content);
            systematic_n_highlow_1D_P2[systematic_n_index][channel]->push_back(content);
        }

    }
}


void systematic_init_helper
(
    // pointer to first object in series in matrix
    std::vector<double> **systematic_array_1D_P1, // nominal / high / low
    std::vector<double> **systematic_array_1D_P2, // nominal / high / low
    std::vector<double> **systematic_array_2D_P1, // nominal / high / low
    std::vector<double> **systematic_array_2D_P2  // nominal / high / low
)
{
    
    // loop over all channels
    for(int channel = 0; channel < number1DHists; ++ channel)
    {

        std::string histname = std::string(channel_histname_1D[channel]);
        std::string search_object_P1;
        std::string search_object_P2;

        search_object_P1 = histname + std::string("fakedata") + "_P1";
        search_object_P2 = histname + std::string("fakedata") + "_P2";

        TH1D *tmpDataHist1D_P1 = nullptr;
        TH1D *tmpDataHist1D_P2 = nullptr;
    
        tmpDataHist1D_P1 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P1.c_str());
        tmpDataHist1D_P2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());

        // phase 1
        if(tmpDataHist1D_P1 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }

        // phase 2
        if(tmpDataHist1D_P2 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
            throw "problem";
        }

        // phase 1
        for(Int_t bin_ix{1}; bin_ix <= tmpDataHist1D_P1->GetNbinsX(); ++ bin_ix)
        {
            Double_t content = tmpDataHist1D_P1->GetBinContent(bin_ix);
            //systematic_array_1D_P1[channel]->push_back(content);
            (*(systematic_array_1D_P1 + channel))->push_back(content);
        }

        // phase 2
        for(Int_t bin_ix{1}; bin_ix <= tmpDataHist1D_P2->GetNbinsX(); ++ bin_ix)
        {
            Double_t content = tmpDataHist1D_P2->GetBinContent(bin_ix);
            (*(systematic_array_1D_P2 + channel))->push_back(content);
        }

    }
    


    // loop over all channels
    for(int channel = 0; channel < number2DHists; ++ channel)
    {

        std::string histname = std::string(channel_histname_2D[channel]);
        std::string search_object_P1;
        std::string search_object_P2;

        search_object_P1 = histname + std::string("fakedata") + "_P1";
        search_object_P2 = histname + std::string("fakedata") + "_P2";

        TH2D *tmpDataHist2D_P1 = nullptr;
        TH2D *tmpDataHist2D_P2 = nullptr;
    
        tmpDataHist2D_P1 = (TH2D*)allFakeDataSamples2D->FindObject(search_object_P1.c_str());
        tmpDataHist2D_P2 = (TH2D*)allFakeDataSamples2D->FindObject(search_object_P2.c_str());

        // phase 1
        if(tmpDataHist2D_P1 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }

        // phase 2
        if(tmpDataHist2D_P2 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
            throw "problem";
        }

        // phase 1
        for(Int_t bin_iy{1}; bin_iy <= tmpDataHist2D_P1->GetNbinsY(); ++ bin_iy)
        {
            for(Int_t bin_ix{1}; bin_ix <= tmpDataHist2D_P1->GetNbinsX(); ++ bin_ix)
            {
                Double_t content = tmpDataHist2D_P1->GetBinContent(bin_ix, bin_iy);
                //systematic_array_2D_P1[channel]->push_back(content);
                (*(systematic_array_2D_P1 + channel))->push_back(content);
            }
        }

        // phase 2
        for(Int_t bin_iy{1}; bin_iy <= tmpDataHist2D_P2->GetNbinsY(); ++ bin_iy)
        {
            for(Int_t bin_ix{1}; bin_ix <= tmpDataHist2D_P2->GetNbinsX(); ++ bin_ix)
            {
                Double_t content = tmpDataHist2D_P2->GetBinContent(bin_ix, bin_iy);
                (*(systematic_array_2D_P2 + channel))->push_back(content);
            }
        }

    }

}


#endif // NEWLOGLIKFITTER_SYSTEMATIC_INIT_HELPER_H
