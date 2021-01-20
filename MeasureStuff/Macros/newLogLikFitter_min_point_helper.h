#ifndef NEWLOGLIKFITTER_MIN_POINT_HELPER_H
#define NEWLOGLIKFITTER_MIN_POINT_HELPER_H


int
min_point_helper(
    const std::string &fit_channel_mode,
    const bool enable_systematics_in_fit,
    const bool mode_fake_data,
    const int numberFreeParams,
    const bool fix_xi_31_param,
    const double xi_31_fixed_value,
    MinPoint & min_point,
    const std::string& min_point_fname,
    const std::string& draw_output_fname,
    const std::string& draw_output_fdir,
    const double xi_31_systematics_reweight_value
)
{

    int restore_channel_enable_1D[number1DHists];
    for(int i = 0; i < number1DHists; ++ i)
    {
        restore_channel_enable_1D[i] = channel_enable_1D[i];
    }

    // switch channels
    if(fit_channel_mode == "CH0")
    {
        channel_enable_1D[0] = 1;
        channel_enable_1D[1] = 0;
        channel_enable_1D[2] = 0;
        channel_enable_1D[3] = 0;
    }
    else if(fit_channel_mode == "CH1")
    {
        channel_enable_1D[0] = 0;
        channel_enable_1D[1] = 1;
        channel_enable_1D[2] = 0;
        channel_enable_1D[3] = 0;
    }
    else if(fit_channel_mode == "CH2CH3")
    {
        channel_enable_1D[2] = 1;
        channel_enable_1D[3] = 1;
        channel_enable_1D[0] = 0;
        channel_enable_1D[1] = 0;
    }
    else
    {
        return -1;
    }
    
    // switch systematics
    if(enable_systematics_in_fit == false)
    {
        V_ENABLE_SYS_stack_push();

        V_ENABLE_SYSALL = false;
        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            V_ENABLE_SYSn[i] = false;
        }
    }
    else
    {
        V_ENABLE_SYS_stack_push();

        V_ENABLE_SYSALL = true;
    }

    // switch g_mode_fake_data
    const bool restore_g_mode_fake_data = g_mode_fake_data;
    g_mode_fake_data = mode_fake_data;



    int ret = min_point_load(min_point_fname, min_point);
    if(ret == 0)
    {
        // min point is loaded, do not calculate
    }
    else
    {
        // min point is not loaded, calculate and save to file

        if(mode_fake_data == true)
        {
            rebuild_fake_data_systematics(xi_31_systematics_reweight_value, xi_31_baseline);
        }

        recalculate_V_PHYS_SYS = true;

        int restore_gNumberFreeParams = gNumberFreeParams;
        gNumberFreeParams = numberFreeParams;

            // create minimizer
            ROOT::Minuit2::MnUserParameterState theParameterStateBefore;
            ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
            MinimizeFCNAxialVector theFCN;

            // initialize fit
            const int xi_31_param_number = g_pg.get_xi_31_ext_param_number();
            const double xi_31_value = g_pg.file_params.at(xi_31_param_number).paramInitValue;
            const double xi_31_error = g_pg.file_params.at(xi_31_param_number).paramInitError;
            std::cout << "xi_31_param_number=" << xi_31_param_number
                      << " xi_31=" << xi_31_value << " +- " << xi_31_error << std::endl;
            fitBackgrounds_init(theParameterStateBefore, theMinimizer, xi_31_value, xi_31_error);

            TString minuit_param_name;
            if(fix_xi_31_param == true)
            {
                // fix xi_31 parameter
                TString i_str;
                i_str.Form("%i", 1);
                TString minuit_param_number_str;
                minuit_param_number_str.Form("%i", 1);
                minuit_param_name = "_" + i_str + "_" + minuit_param_number_str + "_";
                theParameterStateBefore.Fix(std::string(minuit_param_name));
                theParameterStateBefore.SetValue(std::string(minuit_param_name), xi_31_fixed_value);
            }

            // get parameters and chi2 value before fit
            std::vector<double> params_before = theParameterStateBefore.Params();
            std::vector<double> param_errs_before = theParameterStateBefore.Errors();
            double fval_before = theFCN.operator()(params_before);
            int nch = theFCN.nch;
            int nfp = gNumberFreeParams;
            int ndf = nch - nfp;

            // draw before fit
            draw_input_data drawinputdata;
            drawinputdata.chi2 = fval_before;
            drawinputdata.nch = nch;
            drawinputdata.nfp = nfp;
            drawinputdata.serial_dir = draw_output_fdir;
            drawinputdata.saveas_filename = std::string(draw_output_fname) + std::string("_before");
            drawinputdata.saveas_png = true;
           
            draw(drawinputdata,
                 params_before,
                 param_errs_before);
            
            // exec fit
            // this will fit backgrounds and the 150Nd amplitude parameter
            // but xi_31 is fixed
            ROOT::Minuit2::FunctionMinimum FCN_min =
                fitBackgrounds_exec(
                    theParameterStateBefore,
                    theMinimizer,
                    theFCN);

            // get result
            ROOT::Minuit2::MnUserParameterState theParameterStateAfter = FCN_min.UserParameters();

            // get parameters and chi2 value after fit
            std::vector<double> params_after = theParameterStateAfter.Params();
            std::vector<double> param_errs_after = theParameterStateAfter.Errors();
            double fval_after = theFCN.operator()(params_after);
            nch = theFCN.nch;
            nfp = gNumberFreeParams;
            ndf = nch - nfp;

            // draw after fit
            drawinputdata.chi2 = fval_after;
            drawinputdata.nch = nch;
            drawinputdata.nfp = nfp;
            drawinputdata.saveas_filename = std::string(draw_output_fname) + std::string("_after");
           
            draw(drawinputdata,
                 params_after,
                 param_errs_after);

            if(fix_xi_31_param == true)
            {
                // release fixed param
                theParameterStateBefore.Release(std::string(minuit_param_name));
            }

        // save min point to memory
        min_point.xi_31 = params_after.at(1);
        min_point.A = params_after.at(0);
        min_point.xi_31_err = param_errs_after.at(1);
        min_point.A_err = param_errs_after.at(0);
        min_point.fval = fval_after;
        /* not implemented
        min_point.eff = params_after.at(g_pg.get_efficiency_int_param_number());
        min_point.eff_err = param_errs_after.at(g_pg.get_efficiency_int_param_number());
        */
        min_point.eff = 1.0;
        min_point.eff_err = -1.0;
        //std::cout << "efficiency parameter number is " << g_pg.get_efficiency_int_param_number() << std::endl;
        //std::cout << "CHECK,STOP" << std::endl;
        //std::cin.get();

        std::cout << "min_point: " << min_point.xi_31 << " " << min_point.A << std::endl;

        // save min point memory to file
        min_point_save(min_point_fname, min_point);

        std::cout << "??? Fit: NEMO3 ??? (CH?)" << " min_point_fname=" << min_point_fname << std::endl;
        gSystematics_print();

        std::cout << "Result: " << std::endl;
        std::cout << "fval_before=" << fval_before << std::endl;
        std::cout << "fval_after=" << fval_after
                  << " for params_after[0]=" << params_after[0]
                  << " +- " << param_errs_after[0]
                  << " params_after[1]=" << params_after[1]
                  << " +- " << param_errs_after[1]
                  << std::endl;
        std::cout << std::endl;

        std::cout << "---- NUMBER OF DECAYS ----" << std::endl;
        std::string mc_name = "nd150_rot_2n2b_m4";
        std::string histname = std::string(channel_histname_1D[1]);
        std::string search_object_P1 = histname + mc_name + "_P1_fit";
        std::string search_object_P2 = histname + mc_name + "_P2_fit";
        TH1D *tmpHist1D_P1 = nullptr;
        TH1D *tmpHist1D_P2 = nullptr;

        std::cout << "search_object_P1=" << search_object_P1
                  << " search_object_P2=" << search_object_P2 << std::endl;

        double number_of_decays = 0.0;
        double number_of_decays_err = 0.0;

        int param_number_P1 = -2;
        int param_number_P2 = -2;
        double scale_factor_P1 = 0.0;
        double scale_factor_P2 = 0.0;
        double scale_factor_P1_err = 0.0;
        double scale_factor_P2_err = 0.0;
        bool success = g_pg.convert_MC_name_to_param_number(mc_name, param_number_P1, param_number_P2);
        if(success == true)
        {
            if(param_number_P1 != -1)
            {
                scale_factor_P1 = params_after.at(g_pg.ExtToIntParamNumberMap.at(param_number_P1));
                scale_factor_P1_err = param_errs_after.at(g_pg.ExtToIntParamNumberMap.at(param_number_P1));

                tmpHist1D_P1 = (TH1D*)allMCSamples1D[1]->FindObject(search_object_P1.c_str());
                if(tmpHist1D_P1 == nullptr)
                {
                    std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
                    throw "problem";
                }
                number_of_decays += scale_factor_P1 * tmpHist1D_P1->Integral();
                std::cout << "scale_factor_P1=" << scale_factor_P1 << std::endl;
                std::cout << "Integral P1: " << tmpHist1D_P1->Integral() << std::endl;
                std::cout << "P1 number of decays: " << scale_factor_P1 * tmpHist1D_P1->Integral() << std::endl;
                double error = scale_factor_P1_err * tmpHist1D_P1->Integral();
                number_of_decays_err = std::sqrt(std::pow(number_of_decays_err, 2.0) + std::pow(error, 2.0));
            }
            else
            { 
                std::cout << "WARNING: param_number_P1=" << param_number_P1 << std::endl;
                std::cout << "mc_name=" << mc_name << std::endl;
            }

            if(param_number_P2 != -1)
            {
                scale_factor_P2 = params_after.at(g_pg.ExtToIntParamNumberMap.at(param_number_P2));
                scale_factor_P2_err = param_errs_after.at(g_pg.ExtToIntParamNumberMap.at(param_number_P2));

                tmpHist1D_P2 = (TH1D*)allMCSamples1D[1]->FindObject(search_object_P2.c_str());
                if(tmpHist1D_P2 == nullptr)
                {
                    std::cout << "ERROR: Could not find object " << search_object_P2 << std::endl;
                    throw "problem";
                }
                number_of_decays += scale_factor_P2 * tmpHist1D_P2->Integral();
                std::cout << "scale_factor_P2=" << scale_factor_P2 << std::endl;
                std::cout << "Integral P2: " << tmpHist1D_P2->Integral() << std::endl;
                std::cout << "P2 number of decays: " << scale_factor_P2 * tmpHist1D_P2->Integral() << std::endl;
                double error = scale_factor_P2_err * tmpHist1D_P2->Integral();
                number_of_decays_err = std::sqrt(std::pow(number_of_decays_err, 2.0) + std::pow(error, 2.0));
            }
            else
            { 
                std::cout << "WARNING: param_number_P2=" << param_number_P2 << std::endl;
                std::cout << "mc_name=" << mc_name << std::endl;
            }
        }

        std::cout << "number_of_decays=" << number_of_decays << " +- " << number_of_decays_err << std::endl;


        int gNumberFreeParams = restore_gNumberFreeParams;
    }


    // unswitch g_mode_fake_data
    g_mode_fake_data = restore_g_mode_fake_data;

    // unswitch systematics
    if(enable_systematics_in_fit == false)
    {
        V_ENABLE_SYS_stack_pop();
    }
    else
    {
        V_ENABLE_SYS_stack_pop();
    }

    // unswitch channels
    for(int i = 0; i < number1DHists; ++ i)
    {
        channel_enable_1D[i] = restore_channel_enable_1D[i];
    }

    return 0;
}


#endif // NEWLOGLIKFITTER_MIN_POINT_HELPER_H
