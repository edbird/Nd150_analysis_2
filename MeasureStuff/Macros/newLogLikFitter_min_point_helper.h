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
    double *const min_point, double *const min_point_err, double &min_point_fval,
    const std::string& min_point_fname,
    const std::string& draw_output_fname,
    const std::string& draw_output_fdir,
    const double xi_31_systematics_reweight_value
)
{

    // switch channels
    if(fit_channel_mode == "CH0")
    {
        channel_enable_1D[0] = 1;
        channel_enable_1D[1] = 0;
    }
    else if(fit_channel_mode == "CH1")
    {
        channel_enable_1D[0] = 0;
        channel_enable_1D[1] = 1;
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



    int ret = min_point_load(min_point_fname, min_point, min_point_err, min_point_fval);
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
        min_point[0] = params_after.at(1);
        min_point[1] = params_after.at(0);
        min_point_err[0] = param_errs_after.at(1);
        min_point_err[1] = param_errs_after.at(0);
        min_point_fval = fval_after;

        std::cout << "min_point: " << min_point[0] << " " << min_point[1] << std::endl;

        // save min point memory to file
        min_point_save(min_point_fname, min_point, min_point_err, min_point_fval);

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
    channel_enable_1D[0] = 0;
    channel_enable_1D[1] = 1;

    return 0;
}


#endif // NEWLOGLIKFITTER_MIN_POINT_HELPER_H
