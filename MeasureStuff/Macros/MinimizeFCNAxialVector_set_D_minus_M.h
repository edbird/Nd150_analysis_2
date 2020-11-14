#ifndef MINIMIZEFCNAXIALVECTOR_SET_D_MINUS_M
#define MINIMIZEFCNAXIALVECTOR_SET_D_MINUS_M


void
MinimizeFCNAxialVector::set_D_minus_M() const
{
    // M is set, D is set

    
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(channel_enable_1D[channel] == 1)
        {
            // set D_minus_M
            //for(Int_t binx{1}; binx <= D_1D_P1[channel]->GetNbinsX(); ++ binx)
            for(Int_t binx{0}; binx < D_1D_P1_data[channel]->size(); ++ binx)
            {
                // P1
                {
                    //Double_t content_D = D_1D_P1[channel]->GetBinContent(binx, 1);
                    //Double_t content_M = M_1D_P1[channel]->GetBinContent(binx, 1);
                    #if VECTOR_RANGE_CHECK
                    Double_t content_D = D_1D_P1_data[channel]->at(binx);
                    Double_t content_M = M_1D_P1_data[channel]->at(binx);
                    #else
                    Double_t content_D = D_1D_P1_data[channel]->operator[](binx);
                    Double_t content_M = M_1D_P1_data[channel]->operator[](binx);
                    #endif
                    /*
                    if(content_D != 0.0)
                    {
                        std::cout << "binx=" << binx << " ~> " << content_D << " " << content_M << std::endl;
                    }
                    */
                    Double_t content_D_minus_M = content_D - content_M;
                    //D_minus_M_1D_P1[channel]->SetBinContent(binx, 1, content_D_minus_M);
                    #if VECTOR_RANGE_CHECK
                    D_minus_M_1D_P1_data[channel]->at(binx) = content_D_minus_M;
                    #else
                    D_minus_M_1D_P1_data[channel]->operator[](binx) = content_D_minus_M;
                    #endif
                }

                // P2
                {
                    //Double_t content_D = D_1D_P2[channel]->GetBinContent(binx, 1);
                    //Double_t content_M = M_1D_P2[channel]->GetBinContent(binx, 1);
                    #if VECTOR_RANGE_CHECK
                    Double_t content_D = D_1D_P2_data[channel]->at(binx);
                    Double_t content_M = M_1D_P2_data[channel]->at(binx);
                    #else
                    Double_t content_D = D_1D_P2_data[channel]->operator[](binx);
                    Double_t content_M = M_1D_P2_data[channel]->operator[](binx);
                    #endif
                    /*
                    if(content_D != 0.0)
                    {
                        std::cout << "binx=" << binx << " ~> " << content_D << " " << content_M << std::endl;
                    }
                    */
                    Double_t content_D_minus_M = content_D - content_M;
                    //D_minus_M_1D_P2[channel]->SetBinContent(binx, 1, content_D_minus_M);
                    #if VECTOR_RANGE_CHECK
                    D_minus_M_1D_P2_data[channel]->at(binx) = content_D_minus_M;
                    #else
                    D_minus_M_1D_P2_data[channel]->operator[](binx) = content_D_minus_M;
                    #endif
                }
            }
        }
    }
}


#if 0
void MinimizeFCNAxialVector::set_D_minus_M(const std::vector<double> &param) const
{

    // loop over all channels
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(debuglevel >= 4)
        {
            std::cout << "channel=" << channel << std::endl;
        }

        // check channel enabled
        if(channel_enable_1D[channel] == 0)
        {
            if(debuglevel >= 5)
            {
                std::cout << "1D: channel " << channel << " disabled, skip" << std::endl;
            }
            continue;
        }

        // TODO: this no longer works because there are multiple phases
        // in the tmpData1D arrays
        // copy code from below to find correct object
        // TODO: phase 1 and phase 2 data objects
        //Double_t *tmpData1D_P1 = nullptr; //(TH1D*)allDataSamples1D->At(channel);
        //Double_t *tmpData1D_P2 = nullptr;
        //Double_t *tmpFakeData1D_P1 = nullptr; //(TH1D*)allFakeDataSamples1D->At(channel);
        //Double_t *tmpFakeData1D_P2 = nullptr;

        std::string histname = std::string(channel_histname_1D[channel]);
        std::string search_object_P1;
        std::string search_object_P2;
        if(g_mode_fake_data == false)
        {
            search_object_P1 = histname + std::string(DataFile) + "_P1";
            search_object_P2 = histname + std::string(DataFile) + "_P2";
        }
        else if(g_mode_fake_data == true)
        {
            search_object_P1 = histname + std::string("fakedata") + "_P1";
            search_object_P2 = histname + std::string("fakedata") + "_P2";
        }
        TH1D *tmpDataHist1D_P1 = nullptr;
        TH1D *tmpDataHist1D_P2 = nullptr;
        
        if(debuglevel >= 6)
        {
            std::cout << "search_object_P1=" << search_object_P1
                      << " search_object_P2=" << search_object_P2 << std::endl;
        }

        if(g_mode_fake_data == false)
        {
            tmpDataHist1D_P1 = (TH1D*)allDataSamples1D->FindObject(search_object_P1.c_str());
            tmpDataHist1D_P2 = (TH1D*)allDataSamples1D->FindObject(search_object_P2.c_str());
        }
        else if(g_mode_fake_data == true)
        {
            tmpDataHist1D_P1 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P1.c_str());
            tmpDataHist1D_P2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());
        }

        if(tmpDataHist1D_P1 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }
        if(tmpDataHist1D_P2 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }

        ///////////////////////////////////////////////////////////////////////
        // init P1

        // initialize D_minus_M with D value
        for(Int_t bin_x{0}; bin_x < tmpDataHist1D_P1->GetNbinsX(); ++ bin_x)
        {
            Double_t content = tmpDataHist1D_P1->GetBinContent(bin_x + 1);
            #if VECTOR_RANGE_CHECK
            D_minus_M_1D_P1_data[channel]->at(bin_x) = content;
            #else
            D_minus_M_1D_P1_data[channel]->operator[](bin_x) = content;
            #endif
        }

        ///////////////////////////////////////////////////////////////////////
        // init P2
        
        //std::cout << "LIST OF P2 DATA channel=" << channel << std::endl;
        //std::cin.get();
        for(Int_t bin_x{0}; bin_x < tmpDataHist1D_P2->GetNbinsX(); ++ bin_x)
        {
            Double_t content = tmpDataHist1D_P2->GetBinContent(bin_x + 1);
            #if VECTOR_RANGE_CHECK
            D_minus_M_1D_P2_data[channel]->at(bin_x) = content;
            #else
            D_minus_M_1D_P2_data[channel]->operator[](bin_x) = content;
            #endif
            //std::cout << "bin_x=" << bin_x + 1 << " content=" << content_output << std::endl;
        }

        ///////////////////////////////////////////////////////////////////////
        // subtract MC

        // loop over all the parameters
        std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
        for(; it != g_pg.file_params.end(); ++ it)
        {
            int paramNumberInt = -1;

            int paramNumber = it->second.paramNumber;
            //std::cout << "paramNumber=" << paramNumber << std::endl;
            bool paramEnabled = it->second.paramEnabled;
            bool paramEnabledP1 = it->second.paramEnabledP1;
            bool paramEnabledP2 = it->second.paramEnabledP2;
            double paramInitValue = it->second.paramInitValue;
            double paramInitError = it->second.paramInitError;
            int paramConstraintMode = it->second.paramConstraintMode;

            if(debuglevel >= 5)
            {
                std::cout << "paramNumber=" << paramNumber << std::endl;
            }

            // NOTE: notes on the V MATRIX method and the enable/disable
            // of parameters
            //
            // 1: how can parameters be enabled/disabled
            //
            // a:
            // There are two sources, the first is the
            // parameter_names.lst file, which has options to
            // enable/disable the parameter totally, or enable/disable
            // it for P1 and P2
            // 
            // b:
            // The other source is the channel enable/disable switch
            // in the header file
            // This affects the V_CHEN matrix because it switches
            // on/off groups of BINS in the chisquare calculation
            //
            // By contrast, (a) switches on and off parameters
            // which does not affect the number of bins only the bin
            // content

            if(paramEnabled == false)
            {
                if(debuglevel > 0)
                {
                    std::cout << __func__ << " channel [channel TODO] is disabled (overall) skip" << std::endl;
                    std::cin.get();
                }
                continue;
            }
            if((paramEnabledP1 || paramEnabledP2) == false)
            {
                // if both are false
                if(debuglevel > 0)
                {
                    std::cout << __func__ << " channel [channel TODO] is disabled for P1 and P2, skip" << std::endl;
                    std::cin.get();
                }
                continue;
            }
            // TODO: this correctly ignores any parameter which is disabled such that
            // paramEnabled == false
            // however, it also ignores parameters which are set as disabled for P1
            // and P2, and when these are irrelevent due to the value of gEnablePhaseX
            // so... the internal and external index will not match
            // need to add some code to fix this when the parameters are read from
            // file, (probably)
            // unless I just ignore that here... perhaps paramEnabled dictates
            // whether parameter is drawn and the phase1/phase2 enable flag
            // is to decide whether minuit does the fit or not (in which case
            // the param may still contribute to chisquare but may not be minimized
            // by minuit)


            // looping over parameters...
            // get a paramNumber
            // corresponds to a paramNumberInternal (minuit param number)
            // which is the index of the vector param where the current
            // activity can be found
            // also get a paramName, which is something like ac228_bi212_tl208_int
            // which is not very useful
            // but get a list of mc names, which for this example will be
            // ac228_int_rot, bi212_int_rot, tl208_int_rot
            // the corresponding MC samples names will be
            // theHistogram + BkgFiles[i] + "_Px_fit"
            // aka:
            // hTotalE_ + bi212_int_rot + _P1_fit
            // hTotalE_ + bi212_int_rot + _P2_fit
            // now, the paramNumber can be used to look up the parameter in
            // g_pg, from which we can obtain paramEnabled, paramEnabledP1
            // and paramEnabledP2, as well as the constraint mode, constraint
            // value and error
            // so, have to strip off the "_fit" ending, and strip off the
            // histogram (channel) name from the front, aka remove "hTotalE_"
            // then left with something like
            // bi212_int_rot_P1 and a channel number which comes from the loop
            // over allDataSamples1D->GetEntries()
            // EXCEPT this will now be wrong, because we have data for P1 and P2
            // so what we want to loop over is an index channel which runs from
            // 0 to number1DHists
            // so we have a channel index which is used to select the channel
            // which selects the index of allMCSamples1D[channel]
            // we then have to somehow "get" from this TObject
            // the histogram with the name bi212_int_rot_P1
            // and bi212_int_rot_P2
            // once we have done that we switch on paramEnabled, and paramEnabledP1
            // and paramEnabledP2 to decide whether or not we skip (in the case
            // where paramEnabled == false), apply a scaling with param[index]
            // for P1 and P2 if they are enabled, or simply add a value
            // where the scaling factor is given by paramInitValue for P1 and P2
            // if they are disabled
            // TODO: does this last part make sense? no it doesn't. skip if
            // P1/P2 not enabled

            std::vector<std::string>::iterator mc_name_it{it->second.MCNameList.begin()};
            for(; mc_name_it != it->second.MCNameList.end(); ++ mc_name_it)
            {
                std::string mc_name = *mc_name_it;
                std::string histname = std::string(channel_histname_1D[channel]);
                std::string search_object_P1 = histname + mc_name + "_P1_fit";
                std::string search_object_P2 = histname + mc_name + "_P2_fit";
                TH1D *tmpHist1D_P1 = nullptr;
                TH1D *tmpHist1D_P2 = nullptr;

                if(debuglevel >= 6)
                {
                    std::cout << "search_object_P1=" << search_object_P1
                              << " search_object_P2=" << search_object_P2 << std::endl;
                }

                paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);
                if(debuglevel >= 5)
                {
                    std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
                }

                // phase 1
                if(paramEnabledP1 == true)
                {
                    tmpHist1D_P1 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P1.c_str());

                    if(tmpHist1D_P1 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                        throw "problem";
                    }

                    Double_t scale_factor_P1 = param.at(paramNumberInt);
                    if(debuglevel >= 6)
                    {
                        std::cout << "enabled P1: scale factor P1: " << scale_factor_P1 << std::endl;
                    }

                    for(Int_t bin_x{0}; bin_x < tmpHist1D_P1->GetNbinsX(); ++ bin_x)
                    {
                        #if VECTOR_RANGE_CHECK
                        Double_t content_input = D_minus_M_1D_P1_data[channel]->at(bin_x);
                        #else
                        Double_t content_input = D_minus_M_1D_P1_data[channel]->operator[](bin_x);
                        #endif

                        Double_t content_add = scale_factor_P1 * tmpHist1D_P1->GetBinContent(bin_x + 1);
                        
                        Double_t content_output = content_input + content_add;

                        #if VECTOR_RANGE_CHECK
                        D_minus_M_1D_P1_data[channel]->at(bin_x) = content_output;
                        #else
                        D_minus_M_1D_P1_data[channel]->operator[](bin_x) = content_output;
                        #endif
                        //std::cout << "content_input=" << content_input << " content_output=" << content_output << " content_add=" << content_add << std::endl;
                    }
                }
                else
                {
                    if(debuglevel >= 4)
                    {
                        std::cout << "disabled P1" << std::endl;
                    }
                }

                // phase 2
                if(paramEnabledP2 == true)
                {
                    tmpHist1D_P2 = (TH1D*)allMCSamples1D[channel]->FindObject(search_object_P2.c_str());
                
                    if(tmpHist1D_P2 == nullptr)
                    {
                        std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                        throw "problem";
                    }

                    Double_t scale_factor_P2 = param.at(paramNumberInt);
                    if(debuglevel >= 4)
                    {
                        std::cout << "enabled P2: scale factor P2: " << scale_factor_P2 << std::endl;
                    }

                    //std::cout << "super_index (start) : " << channel * 2 * 50 + 50 + 0 << std::endl;
                    for(Int_t bin_x{0}; bin_x < tmpHist1D_P2->GetNbinsX(); ++ bin_x)
                    {
                        #if VECTOR_RANGE_CHECK
                        Double_t content_input = D_minus_M_1D_P2_data[channel]->at(bin_x);
                        #else
                        Double_t content_input = D_minus_M_1D_P2_data[channel]->operator[](bin_x);
                        #endif
                        
                        Double_t content_add = scale_factor_P2 * tmpHist1D_P2->GetBinContent(bin_x + 1);
                        
                        Double_t content_output = content_input + content_add;

                        #if VECTOR_RANGE_CHECK
                        D_minus_M_1D_P2_data[channel]->at(bin_x) = content_output;
                        #else
                        D_minus_M_1D_P2_data[channel]->operator[](bin_x) = content_output;
                        #endif
                        //std::cout << "debug: " << "super_index=" << super_index << " content_input=" << content_input << " content_add=" << content_add << " content_output=" << content_output << " M:" << M->GetBinContent(super_index + 1, 1) << std::endl;

                    }
                }
                else
                {
                    if(debuglevel >= 4)
                    {
                        std::cout << "disabled P2" << std::endl;
                    }
                }




                // TODO: what is parameter is NOT enabled

            } // mc sample name iterator

        } // file_param iterator
    }
}
#endif


#endif // MINIMIZEFCNAXIALVECTOR_SET_D_MINUS_M
