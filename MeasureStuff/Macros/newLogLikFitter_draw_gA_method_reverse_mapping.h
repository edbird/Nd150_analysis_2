#ifndef NEWLOGLIKFITTER_DRAW_GA_METHOD_REVERSE_MAPPING_H
#define NEWLOGLIKFITTER_DRAW_GA_METHOD_REVERSE_MAPPING_H

void draw_gA_method_reverse_mapping()
{


    std::string output_name = "mps_resultsmatrix";
    
    ///////////////////////////////////////////////////////////////////////////
    // OPEN FILE
    ///////////////////////////////////////////////////////////////////////////

    const bool mode_fake_data_flag = g_mode_fake_data;

    std::string after_string = "after";

    mpsdrawdata mps_draw_data_after_sysall;
    mps_draw_data_after_sysall.read(output_name, after_string, true, mode_fake_data_flag);


    ///////////////////////////////////////////////////////////////////////////
    // get phase space object
    ///////////////////////////////////////////////////////////////////////////

    TH2D *h_mps = mps_draw_data_after_sysall.h_mps;
    const double psiN0 = G0_ps_integral_yrinv;
    const double psiN2 = G2_ps_integral_yrinv;

    const double A0_150Nd = 3.45e-04;
    const double NA = 6.022e+23;
    const double m = 36.6;
    const double sec_to_yr = 1.0 / 31557600.0;


    ///////////////////////////////////////////////////////////////////////////
    // convert to a new phase space but loop over the new phase space
    // first, and then use interpolation in the original phase space
    ///////////////////////////////////////////////////////////////////////////

    TH2D *h_ps_loop_primary = new TH2D("h_ps_loop_primary", "h_ps_loop_primary",
                                       100, 0.0, 0.02, 100, 0.0, 1.2);
    for(Int_t j = 1; j <= h_ps_loop_primary->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ps_loop_primary->GetNbinsX(); ++ i)
        {
            const double gA = h_ps_loop_primary->GetYaxis()->GetBinCenter(j);
            const double MGT3 = h_ps_loop_primary->GetXaxis()->GetBinCenter(i);
            const double tmp_value = std::pow(gA, 4.0) * std::pow(MGT3, 2.0);

            // tmp_value = (xi_31 ** 2) / (T12 * (G0 + xi_31 * G2))

            // my origin phase space contains xi_31 in the range x_min:x_max
            // and T12 in the range function(y_min:y_max)

            // first do an experiment by fixing T12 using A0 = 3.45e-4
            //double T12 = (std::log(2.0) * NA * m) / (150.0 * A0_150Nd * A_150Nd) * sec_to_yr;
            //double A_150Nd = 1.0;

            for(Int_t ii = 1; ii < h_mps->GetNbinsX(); ++ ii)
            {
                double xi_31 = h_mps->GetXaxis()->GetBinCenter(ii);

                // convert tmp_value into T12
                const double G0 = G0_ps_integral_yrinv;
                const double G2 = G2_ps_integral_yrinv;
                double T12 = std::pow(xi_31, 2.0) / (tmp_value * (G0 + xi_31 * G2));

                // convert T21 into activity
                const double A0_150Nd = 3.45e-04;
                const double NA = 6.022e+23;
                const double m = 36.6;
                //const double sec_to_yr = 1.0 / 31557600.0;
                const double yr_to_sec = 31557600.0;
                //double A_150Nd = (T12 * yr_to_sec * 150.0 * A0_150Nd) / (std::log(2.0) * NA * m);
                double A_150Nd = (std::log(2.0) * NA * m) / (150.0 * T12 * yr_to_sec * A0_150Nd);

                Int_t jj = h_mps->GetYaxis()->FindBin(A_150Nd);
                if((1 <= jj) && (jj <= 50))
                {
                    //double chi2 = h_mps->GetBinContent(ii, jj);
                    double chi2 = h_mps->Interpolate(xi_31, A_150Nd);

                    double content = h_ps_loop_primary->GetBinContent(i, j);
                    if((content == 0.0) || (content == -1.0))
                    {
                        h_ps_loop_primary->SetBinContent(i, j, chi2);
                    }
                    else if(chi2 < content)
                    {
                        h_ps_loop_primary->SetBinContent(i, j, chi2);
                        //std::cout << "bin i=" << i << " j=" << j << " chi2: " << content << " -> " << chi2 << std::endl;
                    }
                }
                else
                {
                    //std::cout << "i=" << i << " j=" << j << " does not map to any value inside original phase space" << std::endl;

                    double content = h_ps_loop_primary->GetBinContent(i, j);
                    if(content == 0.0)
                    {
                        h_ps_loop_primary->SetBinContent(i, j, -1.0);
                    }
                    else
                    {
                        // do nothing
                    }
                    //else if(chi2 < content)
                    //{
                    //    h_ps_loop_primary->SetBinContent(i, j, chi2);
                    //    std::cout << "bin i=" << i << " j=" << j << " chi2: " << content << " -> " << chi2 << std::endl;
                    //}
                }
            }
            for(Int_t jj = 1; jj <= h_ps_loop_primary->GetNbinsY(); ++ jj)
            {
                for(Int_t ii = 1; ii <= h_ps_loop_primary->GetNbinsX(); ++ ii)
                {
                    if(h_ps_loop_primary->GetBinContent(ii, jj) == -1.0)
                    {
                        //h_ps_loop_primary->SetBinContent(ii, jj, 40.01);
                        h_ps_loop_primary->SetBinContent(ii, jj, 0.01);
                    }
                }
            }
        }
    }

    TCanvas *c_ps_loop_primary = new TCanvas("c_ps_loop_primary", "c_ps_loop_primary");
    h_ps_loop_primary->SetContour(1000);
//    h_ps_loop_primary->GetZaxis()->SetRangeUser(40.0, 60.0);
    /*
            for(Int_t jj = 1; jj <= h_ps_loop_primary->GetNbinsY(); ++ jj)
            {
                for(Int_t ii = 1; ii <= h_ps_loop_primary->GetNbinsX(); ++ ii)
                {
                        double x = h_ps_loop_primary->GetXaxis()->GetBinCenter(ii);
                        double y = h_ps_loop_primary->GetYaxis()->GetBinCenter(jj);
                        h_ps_loop_primary->SetBinContent(ii, jj, std::pow(x, 2.0) * std::pow(y, 4.0));
                }
            }
    */
    h_ps_loop_primary->Draw("colz");


    /*
    TGraph *gtmp = new TGraph(100);
    for(int i = 0; i < 100; ++ i)
    {
        const double xi_31 = xi_31_max[0];
        const double G0 = G0_ps_integral_yrinv;
        const double G2 = G2_ps_integral_yrinv;
        const double x = 0.2 * ((double)(i + 0.5) / (double)(100));
        const double MGT3 = x;
        double value = std::pow(MGT3, 2.0) / std::pow(xi_31, 2.0) * (G0 + xi_31 * G2);
        double T12 = 7.5e+18;
        double value2 = 1.0 / (T12 * value);
        double gA = std::pow(value2, 0.25);
        gtmp->SetPoint(i, MGT3, gA);
    }

    gtmp->SetLineWidth(2.0);
    gtmp->Draw("Lsame");
    */


}

#endif // NEWLOGLIKFITTER_DRAW_GA_METHOD_REVERSE_MAPPING_H
