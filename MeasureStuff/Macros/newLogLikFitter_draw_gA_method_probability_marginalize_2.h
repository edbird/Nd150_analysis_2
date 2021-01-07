#ifndef NEWLOGLIKFITTER_DRAW_GA_METHOD_PROBABILITY_MARGINALIZE_2_H
#define NEWLOGLIKFITTER_DRAW_GA_METHOD_PROBABILITY_MARGINALIZE_2_H


// uses exterior integral of xi_31 to find CL limits
//


void draw_gA_method_probability_marginalize_2()
{

    //std::string output_name = "mps_resultsmatrix";
    std::string output_name = "noparallel";
    
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

    
    ///////////////////////////////////////////////////////////////////////////
    // check minimum
    ///////////////////////////////////////////////////////////////////////////

    double min_chi2 = 0.0;
    double min_x = 0.0;
    double min_y = 0.0;
    Int_t min_bin_x = 0;
    Int_t min_bin_y = 0;
    for(Int_t j = 1; j <= h_mps->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_mps->GetNbinsX(); ++ i)
        {
            double content = h_mps->GetBinContent(i, j);
            if(((i == 1) && (j == 1)) || content < min_chi2)
            {
                min_chi2 = content;
                min_x = h_mps->GetXaxis()->GetBinCenter(i);
                min_y = h_mps->GetYaxis()->GetBinCenter(j);
                min_bin_x = i;
                min_bin_y = j;
            }
        }
    }
    std::cout << std::endl;
    std::cout << "****************" << std::endl;
    std::cout << "minimum chi2 point:" << std::endl;
    std::cout << "chi2=" << min_chi2 << std::endl;
    std::cout << "i, j = " << min_bin_x << ", " << min_bin_y << std::endl;
    std::cout << "x, y = " << min_x << ", " << min_y << std::endl;
    std::cout << "bin width x, y = " << h_mps->GetXaxis()->GetBinWidth(1) << ", " << h_mps->GetYaxis()->GetBinWidth(1) << std::endl;
    std::cout << "****************" << std::endl;


    ///////////////////////////////////////////////////////////////////////////
    // convert to probability space
    ///////////////////////////////////////////////////////////////////////////

    TH2D *h_prob = new TH2D("h_prob", "h_prob",
                            h_mps->GetNbinsX(),
                            h_mps->GetXaxis()->GetBinLowEdge(1),
                            h_mps->GetXaxis()->GetBinUpEdge(h_mps->GetNbinsX()),
                            h_mps->GetNbinsY(),
                            h_mps->GetYaxis()->GetBinLowEdge(1),
                            h_mps->GetYaxis()->GetBinUpEdge(h_mps->GetNbinsY()));

    std::cout << "h_prob: x axis: limits: "
              << h_mps->GetXaxis()->GetBinLowEdge(1) << " "
              << h_mps->GetXaxis()->GetBinUpEdge(h_mps->GetNbinsX())
              << std::endl;

    for(Int_t j = 1; j <= h_mps->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_mps->GetNbinsX(); ++ i)
        {
            const double content = h_mps->GetBinContent(i, j);
            const double value = std::exp(-0.5 * content);
            h_prob->SetBinContent(i, j, value);
        }
    }

    // normalize
    double sum = 0.0;
    for(Int_t j = 1; j <= h_prob->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_prob->GetNbinsX(); ++ i)
        {
            sum += h_prob->GetBinContent(i, j);
        }
    }

    const double PROB_X_MIN = h_prob->GetXaxis()->GetBinLowEdge(1);
    const double PROB_X_MAX = h_prob->GetXaxis()->GetBinUpEdge(h_prob->GetNbinsX());
    const double PROB_X_RANGE = PROB_X_MAX - PROB_X_MIN;
    const double PROB_Y_MIN = h_prob->GetYaxis()->GetBinLowEdge(1);
    const double PROB_Y_MAX = h_prob->GetYaxis()->GetBinUpEdge(h_prob->GetNbinsY());
    const double PROB_Y_RANGE = PROB_Y_MAX - PROB_Y_MIN;

    std::cout << "Scale Factor: " << 1.0 / sum << std::endl;
    const double bin_area = (PROB_Y_RANGE / 301.0) * (PROB_X_RANGE / 301.0);
    std::cout << "bin_area=" << bin_area << std::endl;
    h_prob->Scale(1.0 / (sum * bin_area));
    sum = 0.0;
    for(Int_t j = 1; j <= h_prob->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_prob->GetNbinsX(); ++ i)
        {
            sum += h_prob->GetBinContent(i, j);
        }
    }
    std::cout << "Integral of 2D P-space after scale: " << sum * bin_area << std::endl;


    ///////////////////////////////////////////////////////////////////////////
    // check maximum
    ///////////////////////////////////////////////////////////////////////////

    double max_prob = 0.0;
    double max_x = 0.0;
    double max_y = 0.0;
    Int_t max_bin_x = 0;
    Int_t max_bin_y = 0;
    for(Int_t j = 1; j <= h_prob->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_prob->GetNbinsX(); ++ i)
        {
            double content = h_prob->GetBinContent(i, j);
            if(((i == 1) && (j == 1)) || content > max_prob)
            {
                max_prob = content;
                max_x = h_prob->GetXaxis()->GetBinCenter(i);
                max_y = h_prob->GetYaxis()->GetBinCenter(j);
                max_bin_x = i;
                max_bin_y = j;
            }
        }
    }
    std::cout << std::endl;
    std::cout << "****************" << std::endl;
    std::cout << "maximum prob point:" << std::endl;
    std::cout << "prob=" << max_prob << std::endl;
    std::cout << "i, j = " << max_bin_x << ", " << max_bin_y << std::endl;
    std::cout << "x, y = " << max_x << ", " << max_y << std::endl;
    std::cout << "bin width x, y = " << h_prob->GetXaxis()->GetBinWidth(1) << ", " << h_prob->GetYaxis()->GetBinWidth(1) << std::endl;
    std::cout << "****************" << std::endl;

    
    ///////////////////////////////////////////////////////////////////////////
    // draw p-space
    ///////////////////////////////////////////////////////////////////////////

    TCanvas *c_prob = new TCanvas("c_prob", "c_prob");
    c_prob->SetTicks(2, 2);
    c_prob->SetRightMargin(0.15);
    c_prob->SetBottomMargin(0.15);
    h_prob->SetTitle("");
    h_prob->SetStats(0);
    h_prob->GetZaxis()->SetLabelOffset(0.005);
    h_prob->GetXaxis()->SetLabelSize(17.0);
    h_prob->GetXaxis()->SetLabelFont(43);
    h_prob->GetYaxis()->SetLabelSize(17.0);
    h_prob->GetYaxis()->SetLabelFont(43);
    h_prob->GetZaxis()->SetLabelSize(17.0);
    h_prob->GetZaxis()->SetLabelFont(43);
    h_prob->GetXaxis()->SetTitleSize(18.0);
    h_prob->GetXaxis()->SetTitleFont(43);
    h_prob->GetYaxis()->SetTitleSize(18.0);
    h_prob->GetYaxis()->SetTitleFont(43);
    h_prob->GetXaxis()->SetTitleOffset(1.5);
    h_prob->GetYaxis()->SetTitleOffset(1.2);
    h_prob->GetXaxis()->SetLabelOffset(0.01);
    h_prob->GetYaxis()->SetLabelOffset(0.01);
    h_prob->SetXTitle("#xi^{2#nu#beta#beta}_{31}");
    h_prob->SetYTitle("^{150}Nd Amplitude Scale Factor (#times 3.45#times10^{-4} Bq)");
    h_prob->Draw("colz");
    c_prob->Update();

    TPaletteAxis *palette = (TPaletteAxis*)h_prob->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.89 + 0.03);
    palette->SetX2NDC(0.93 + 0.03);
    palette->SetY1NDC(0.15);
    palette->SetY2NDC(0.9);
    palette->Draw();
    gPad->Modified();
    gPad->Update();
    c_prob->Modified();

    c_prob->SaveAs("c_prob.eps");
    c_prob->SaveAs("c_prob.pdf");
    c_prob->SaveAs("c_prob.png");


    ///////////////////////////////////////////////////////////////////////////
    // marginal_xi_31
    ///////////////////////////////////////////////////////////////////////////

    TH1D *h_marginal_xi_31 = new TH1D("h_marginal_xi_31", "h_marginal_xi_31",
                                h_prob->GetNbinsX(),
                                h_prob->GetXaxis()->GetBinLowEdge(1),
                                h_prob->GetXaxis()->GetBinUpEdge(h_prob->GetNbinsX()));

        const double bin_width_x = PROB_X_RANGE / 301.0;
        const double bin_width_y = PROB_Y_RANGE / 301.0;

        double peak_y = 0.0;
        double peak_x = 0.0;
        Int_t peak_i = 0;
        for(Int_t i = 1; i <= h_prob->GetNbinsX(); ++ i)
        {
            double sum = 0.0;
            for(Int_t j = 1; j <= h_prob->GetNbinsY(); ++ j)
            {
                sum += h_prob->GetBinContent(i, j) * bin_width_y;
            }
            h_marginal_xi_31->SetBinContent(i, sum);

            if((i == 1) || (sum > peak_y))
            {
                peak_y = sum;
                peak_x = h_prob->GetXaxis()->GetBinCenter(i);
                peak_i = i;
            }
        }
        //std::cout << "sum=" << sum << " ->Scale(" << 1.0 / (sum * bin_width) << ")" << std::endl;
        //h_marginal_xi_31->Scale(1.0 / (sum * bin_width));

        std::cout << "peak: " << peak_x << " " << peak_y << std::endl;
        std::cout << "peak_i=" << peak_i << std::endl;

        sum = 0.0;
        for(Int_t i = 1; i <= h_marginal_xi_31->GetNbinsX(); ++ i)
        {
            sum += h_marginal_xi_31->GetBinContent(i);
        }
        std::cout << "integral in 1D P-space: " << sum * bin_width_x << std::endl;


        ///////////////////////////////////////////////////////////////////////////
        // exteriour (canvas)
        ///////////////////////////////////////////////////////////////////////////

    {



        ///////////////////////////////////////////////////////////////////////
        // cumulative exteriour
        ///////////////////////////////////////////////////////////////////////


        TH1D *h_cumulative_xi_31_low_ext = new TH1D("h_cumulative_xi_31_low_ext", "h_cumulative_xi_31_low_ext",
                                      h_marginal_xi_31->GetNbinsX(),
                                      h_marginal_xi_31->GetXaxis()->GetBinLowEdge(1),
                                      h_marginal_xi_31->GetXaxis()->GetBinUpEdge(h_marginal_xi_31->GetNbinsX()));
        h_cumulative_xi_31_low_ext->SetTitle(0);
        h_cumulative_xi_31_low_ext->SetStats(0);

        sum = 0.0;
        for(Int_t i = 1; i <= h_marginal_xi_31->GetNbinsX(); ++ i)
        {
            sum += h_marginal_xi_31->GetBinContent(i) * bin_width_x;
            h_cumulative_xi_31_low_ext->SetBinContent(i, sum);
        }
/*
        TCanvas *c_cumulative_xi_31_low_ext = new TCanvas("c_cumulative_xi_31_low_ext", "c_cumulative_xi_31_low_ext");
        h_cumulative_xi_31_low_ext->Draw("hist");
*/

        TH1D *h_cumulative_xi_31_high_ext = new TH1D("h_cumulative_xi_31_high_ext", "h_cumulative_xi_31_high_ext",
                                      h_marginal_xi_31->GetNbinsX(),
                                      h_marginal_xi_31->GetXaxis()->GetBinLowEdge(1),
                                      h_marginal_xi_31->GetXaxis()->GetBinUpEdge(h_marginal_xi_31->GetNbinsX()));
        h_cumulative_xi_31_high_ext->SetTitle(0);
        h_cumulative_xi_31_high_ext->SetStats(0);

        sum = 0.0;
        for(Int_t i = h_marginal_xi_31->GetNbinsX(); i >= 1; -- i)
        {
            sum += h_marginal_xi_31->GetBinContent(i) * bin_width_x;
            h_cumulative_xi_31_high_ext->SetBinContent(i, sum);
        }
/*
        TCanvas *c_cumulative_xi_31_high_ext = new TCanvas("c_cumulative_xi_31_high_ext", "c_cumulative_xi_31_high_ext");
        h_cumulative_xi_31_high_ext->Draw("hist");
*/



        TGraph *g_cumulative_xi_31_low_ext_reverse = new TGraph(h_cumulative_xi_31_low_ext->GetNbinsX());
        g_cumulative_xi_31_low_ext_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_xi_31_low_ext_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_xi_31_low_ext->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_xi_31_low_ext->GetBinContent(i + 1);
            g_cumulative_xi_31_low_ext_reverse->SetPoint(i, y, x);
        }
/*
        TCanvas *c_cumulative_xi_31_low_ext_reverse = new TCanvas("c_cumulative_xi_31_low_ext_reverse", "c_cumulative_xi_31_low_ext_reverse");
        g_cumulative_xi_31_low_ext_reverse->Draw("al");
*/


        TGraph *g_cumulative_xi_31_high_ext_reverse = new TGraph(h_cumulative_xi_31_high_ext->GetNbinsX());
        g_cumulative_xi_31_high_ext_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_xi_31_high_ext_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_xi_31_high_ext->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_xi_31_high_ext->GetBinContent(i + 1);
            g_cumulative_xi_31_high_ext_reverse->SetPoint(i, y, x);
        }
/*
        TCanvas *c_cumulative_xi_31_high_ext_reverse = new TCanvas("c_cumulative_xi_31_high_ext_reverse", "c_cumulative_xi_31_high_ext_reverse");
        g_cumulative_xi_31_high_ext_reverse->Draw("al");
*/

        const double xi_31_1sigma_low = g_cumulative_xi_31_low_ext_reverse->Eval((1.0 - 0.683) / 2.0);
        const double xi_31_1sigma_high = g_cumulative_xi_31_high_ext_reverse->Eval((1.0 - 0.683) / 2.0);
        const double xi_31_90CL_low = g_cumulative_xi_31_low_ext_reverse->Eval((1.0 - 0.90) / 2.0);
        const double xi_31_90CL_high = g_cumulative_xi_31_high_ext_reverse->Eval((1.0 - 0.90) / 2.0);

        std::cout << "1 sigma (low ext): " << xi_31_1sigma_low << std::endl;
        std::cout << "90 % (low ext): " << xi_31_90CL_low << std::endl;
        std::cout << "1 sigma (high ext): " << xi_31_1sigma_high << std::endl;
        std::cout << "90 % (high ext): " << xi_31_90CL_high << std::endl;

        const int xi_31_1sigma_low_index = h_marginal_xi_31->GetXaxis()->FindBin(xi_31_1sigma_low);
        const int xi_31_1sigma_high_index = h_marginal_xi_31->GetXaxis()->FindBin(xi_31_1sigma_high);
        const int xi_31_90CL_low_index = h_marginal_xi_31->GetXaxis()->FindBin(xi_31_90CL_low);
        const int xi_31_90CL_high_index = h_marginal_xi_31->GetXaxis()->FindBin(xi_31_90CL_high);

        std::cout << "*** bin index ***" << std::endl;
        std::cout << xi_31_1sigma_low_index << std::endl;
        std::cout << xi_31_90CL_low_index << std::endl;
        std::cout << xi_31_1sigma_high_index << std::endl;
        std::cout << xi_31_90CL_high_index << std::endl;

//        std::cout << "1 sigma: " << g_cumulative_xi_31_reverse->Eval(0.683) << std::endl;
//        std::cout << "90 %: " << g_cumulative_xi_31_reverse->Eval(0.90) << std::endl;
//        std::cout << "95 %: " << g_cumulative_xi_31_reverse->Eval(0.95) << std::endl;
//        std::cout << "99 %: " << g_cumulative_xi_31_reverse->Eval(0.99) << std::endl;

//        std::cout << "*** bin index ***" << std::endl;
//        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_reverse->Eval(0.683)) << std::endl;
//        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_reverse->Eval(0.90)) << std::endl;
//        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_reverse->Eval(0.95)) << std::endl;
//        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_reverse->Eval(0.99)) << std::endl;

/*
        double sum0 = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sum3 = 0.0;
        double sum4 = 0.0;
        for(Int_t i = 1; i <= h_marginal_xi_31->GetNbinsX(); ++ i)
        {
            const double xi_31 = h_marginal_xi_31->GetXaxis()->GetBinCenter(i);
            const double content = h_marginal_xi_31->GetBinContent(i);

            sum0 += content;
            if(xi_31 < g_cumulative_xi_31_reverse->Eval(0.683))
            {
                sum1 += content;
            }
            if(xi_31 < g_cumulative_xi_31_reverse->Eval(0.90))
            {
                sum2 += content;
            }
            if(xi_31 < g_cumulative_xi_31_reverse->Eval(0.95))
            {
                sum3 += content;
            }
            if(xi_31 < g_cumulative_xi_31_reverse->Eval(0.99))
            {
                sum4 += content;
            }
        }
        std::cout << "sum0=" << sum0 << " sum1=" << sum1 << " sum2=" << sum2 << " sum3=" << sum3 << " sum4=" << sum4 << std::endl;
*/

        TCanvas *c_marginal_xi_31_ext = new TCanvas("c_marginal_xi_31_ext", "c_marginal_xi_31_ext");
        c_marginal_xi_31_ext->SetTicks(1, 2);
        c_marginal_xi_31_ext->SetBottomMargin(0.15);
        c_marginal_xi_31_ext->SetRightMargin(0.1);
        c_marginal_xi_31_ext->SetLeftMargin(0.12);
        h_marginal_xi_31->SetTitle("");
        h_marginal_xi_31->SetStats(0);
        h_marginal_xi_31->GetZaxis()->SetLabelOffset(0.005);
        h_marginal_xi_31->GetXaxis()->SetLabelSize(17.0);
        h_marginal_xi_31->GetXaxis()->SetLabelFont(43);
        h_marginal_xi_31->GetYaxis()->SetLabelSize(17.0);
        h_marginal_xi_31->GetYaxis()->SetLabelFont(43);
        h_marginal_xi_31->GetZaxis()->SetLabelSize(17.0);
        h_marginal_xi_31->GetZaxis()->SetLabelFont(43);
        h_marginal_xi_31->GetXaxis()->SetTitleSize(18.0);
        h_marginal_xi_31->GetXaxis()->SetTitleFont(43);
        h_marginal_xi_31->GetYaxis()->SetTitleSize(18.0);
        h_marginal_xi_31->GetYaxis()->SetTitleFont(43);
        h_marginal_xi_31->GetXaxis()->SetTitleOffset(1.5);
        h_marginal_xi_31->GetYaxis()->SetTitleOffset(1.2);
        h_marginal_xi_31->GetXaxis()->SetLabelOffset(0.01);
        h_marginal_xi_31->GetYaxis()->SetLabelOffset(0.01);
        h_marginal_xi_31->SetXTitle("#xi^{2#nu#beta#beta}_{31}");
        h_marginal_xi_31->SetYTitle("Probability Density Function p(#xi^{2#nu#beta#beta}_{31})");
        std::cout << "bin width bin(1): " << h_marginal_xi_31->GetBinWidth(1) << std::endl;
        std::cout << "INT=" << h_marginal_xi_31->Integral() * h_marginal_xi_31->GetBinWidth(1) << std::endl;
        TH1D *h_marginal_xi_31_clone = (TH1D*)h_marginal_xi_31->Clone();
    //    h_marginal_xi_31_clone->Scale(1.0 / (2.0 / 301.0));
    //    for(Int_t i = 0; i <= h_marginal_xi_31->GetNbinsX() + 1; ++ i)
        for(Int_t i = 1; i <= h_marginal_xi_31->GetNbinsX(); ++ i)
        {
            const double content = h_marginal_xi_31_clone->GetBinContent(i);
    //        h_marginal_xi_31_clone->SetBinContent(i, 301.0 / 2.0 * content);
    //        std::cout << "i=" << i << " " << h_marginal_xi_31->GetBinContent(i) << std::endl;
        }
//        h_marginal_xi_31_clone->SetMaximum(1.7);
//        h_marginal_xi_31_clone->SetMinimum(0.0);
        h_marginal_xi_31_clone->SetMaximum(1.6);
        h_marginal_xi_31_clone->Draw("axis");
//        h_marginal_xi_31_clone->Draw("hist");

        TH1D *h_marginal_xi_31_clone1 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone1->SetFillColor(kRed);
        h_marginal_xi_31_clone1->SetFillStyle(3002);
        h_marginal_xi_31_clone1->SetLineColor(kRed);
        h_marginal_xi_31_clone1->SetLineWidth(2);
        h_marginal_xi_31_clone1->GetXaxis()->SetRange(xi_31_90CL_low_index, xi_31_1sigma_low_index);
        h_marginal_xi_31_clone1->Draw("][histsame");

        TH1D *h_marginal_xi_31_clone2 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone2->SetFillColor(kOrange - 3);
        h_marginal_xi_31_clone2->SetFillStyle(3002);
        h_marginal_xi_31_clone2->SetLineColor(kOrange - 3);
        h_marginal_xi_31_clone2->SetLineWidth(2);
        h_marginal_xi_31_clone2->GetXaxis()->SetRange(1, xi_31_90CL_low_index);
        h_marginal_xi_31_clone2->Draw("][histsame");

        TH1D *h_marginal_xi_31_clone3 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone3->SetFillColor(kBlue - 4);
        h_marginal_xi_31_clone3->SetFillStyle(3002);
        h_marginal_xi_31_clone3->SetLineColor(kBlue - 4);
        h_marginal_xi_31_clone3->SetLineWidth(2);
        h_marginal_xi_31_clone3->GetXaxis()->SetRange(xi_31_1sigma_high_index, xi_31_90CL_high_index);
        h_marginal_xi_31_clone3->Draw("][histsame");

        TH1D *h_marginal_xi_31_clone4 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone4->SetFillColor(kViolet - 3);
        h_marginal_xi_31_clone4->SetFillStyle(3002);
        h_marginal_xi_31_clone4->SetLineColor(kViolet - 3);
        h_marginal_xi_31_clone4->SetLineWidth(2);
        h_marginal_xi_31_clone4->GetXaxis()->SetRange(xi_31_90CL_high_index, 301);
        h_marginal_xi_31_clone4->Draw("][histsame");

        TH1D *h_marginal_xi_31_clone5 = (TH1D*)h_marginal_xi_31->Clone();
//        h_marginal_xi_31_clone5->SetFillColor(kBlack);
//        h_marginal_xi_31_clone5->SetFillStyle(3002);
        h_marginal_xi_31_clone5->SetLineColor(kBlack);
        h_marginal_xi_31_clone5->SetLineWidth(2);
        h_marginal_xi_31_clone5->GetXaxis()->SetRange(xi_31_1sigma_low_index, xi_31_1sigma_high_index);
        h_marginal_xi_31_clone5->Draw("][histsame");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(63);
        latex.SetTextSize(15);

        TString sext = "Exterior Integral";
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.17, 0.83, sext);

        TString lstring0;
        lstring0.Form("#int_{-0.1}^{1.9} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 1");
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.17, 0.72, lstring0);
        
        TString lstring1;
        lstring1.Form("#int_{-0.1}^{%.3f} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 68.3 %%", xi_31_1sigma_low);
        latex.SetTextColor(kRed);
        latex.DrawLatex(0.55, 0.78, lstring1);
        
        TString lstring2;
        lstring2.Form("#int_{-0.1}^{%.3f} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 90 %%", xi_31_90CL_low);
        latex.SetTextColor(kOrange - 3);
        latex.DrawLatex(0.59, 0.65, lstring2);
        
        TString lstring3;
        lstring3.Form("#int_{%.3f}^{1.9} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 68.3 %%", xi_31_1sigma_high);
        latex.SetTextColor(kBlue - 4);
        latex.DrawLatex(0.61, 0.50, lstring3);
        
        TString lstring4;
        lstring4.Form("#int_{%.3f}^{1.9} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 90 %%", xi_31_90CL_high);
        latex.SetTextColor(kViolet - 3);
        latex.DrawLatex(0.65, 0.35, lstring4);
    
        c_marginal_xi_31_ext->SaveAs("c_marginal_xi_31_ext.eps");
        c_marginal_xi_31_ext->SaveAs("c_marginal_xi_31_ext.pdf");
        c_marginal_xi_31_ext->SaveAs("c_marginal_xi_31_ext.png");

    }



    ///////////////////////////////////////////////////////////////////////////
    // interior (canvas)
    ///////////////////////////////////////////////////////////////////////////

    // note this doesn't work at all


#if 0
    {

        TCanvas *c_marginal_xi_31_int = new TCanvas("c_marginal_xi_31_int", "c_marginal_xi_31_int");
        c_marginal_xi_31_int->SetTicks(1, 2);
        c_marginal_xi_31_int->SetBottomMargin(0.15);
        c_marginal_xi_31_int->SetRightMargin(0.1);
        c_marginal_xi_31_int->SetLeftMargin(0.12);

        TH1D *h_marginal_xi_31_clone2 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone2->Draw("axis");

        TH1D *h_marginal_xi_31_clone1 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone1->SetFillColor(kRed);
        h_marginal_xi_31_clone1->SetFillStyle(3002);
        h_marginal_xi_31_clone1->SetLineColor(kRed);
        h_marginal_xi_31_clone1->SetLineWidth(2);
    //    h_marginal_xi_31_clone1->GetXaxis()->SetRangeUser(-0.1, 0.864839);
        h_marginal_xi_31_clone1->GetXaxis()->SetRange(70, 90);
        h_marginal_xi_31_clone1->Draw("][histsame");

        TH1D *h_marginal_xi_31_clone2 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone2->SetFillColor(kOrange - 3);
        h_marginal_xi_31_clone2->SetFillStyle(3002);
        h_marginal_xi_31_clone2->SetLineColor(kOrange - 3);
        h_marginal_xi_31_clone2->SetLineWidth(2);
    //    h_marginal_xi_31_clone2->GetXaxis()->SetRangeUser(0.864839, 1.10726);
        h_marginal_xi_31_clone2->GetXaxis()->SetRange(1, 70);
        h_marginal_xi_31_clone2->Draw("][histsame");

        TH1D *h_marginal_xi_31_clone3 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone3->SetFillColor(kBlue - 4);
        h_marginal_xi_31_clone3->SetFillStyle(3002);
        h_marginal_xi_31_clone3->SetLineColor(kBlue - 4);
        h_marginal_xi_31_clone3->SetLineWidth(2);
    //    h_marginal_xi_31_clone3->GetXaxis()->SetRangeUser(1.10726, 1.2279);
        h_marginal_xi_31_clone3->GetXaxis()->SetRange(170, 201);
        h_marginal_xi_31_clone3->Draw("][histsame");

        TH1D *h_marginal_xi_31_clone4 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone4->SetFillColor(kViolet - 3);
        h_marginal_xi_31_clone4->SetFillStyle(3002);
        h_marginal_xi_31_clone4->SetLineColor(kViolet - 3);
        h_marginal_xi_31_clone4->SetLineWidth(2);
    //    h_marginal_xi_31_clone4->GetXaxis()->SetRangeUser(1.2279, 1.47319);
        h_marginal_xi_31_clone4->GetXaxis()->SetRange(201, 301);
        h_marginal_xi_31_clone4->Draw("][histsame");

        TH1D *h_marginal_xi_31_clone5 = (TH1D*)h_marginal_xi_31->Clone();
//        h_marginal_xi_31_clone5->SetFillColor(kBlack);
//        h_marginal_xi_31_clone5->SetFillStyle(3002);
        h_marginal_xi_31_clone5->SetLineColor(kBlack);
        h_marginal_xi_31_clone5->SetLineWidth(2);
    //    h_marginal_xi_31_clone5->GetXaxis()->SetRangeUser(1.47319, 1.9);
        h_marginal_xi_31_clone5->GetXaxis()->SetRange(90, 170);
        h_marginal_xi_31_clone5->Draw("][histsame");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(63);
        latex.SetTextSize(15);

        TString sext = "Exterior Integral";
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.17, 0.83, sext);

        TString lstring0;
        lstring0.Form("#int_{-0.1}^{1.9} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 1");
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.17, 0.72, lstring0);
        
        TString lstring1;
        lstring1.Form("#int_{-0.1}^{0.50} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 68.3 %%");
        latex.SetTextColor(kRed);
        latex.DrawLatex(0.55, 0.78, lstring1);
        
        TString lstring2;
        lstring2.Form("#int_{-0.1}^{0.36} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 90 %%");
        latex.SetTextColor(kOrange - 3);
        latex.DrawLatex(0.59, 0.65, lstring2);
        
        TString lstring3;
        lstring3.Form("#int_{1.03}^{1.9} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 68.3 %%");
        latex.SetTextColor(kBlue - 4);
        latex.DrawLatex(0.61, 0.50, lstring3);
        
        TString lstring4;
        lstring4.Form("#int_{1.23}^{1.9} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 90 %%");
        latex.SetTextColor(kViolet - 3);
        latex.DrawLatex(0.65, 0.35, lstring4);

        c_marginal_xi_31_int->SaveAs("c_marginal_xi_31_int.eps");
        c_marginal_xi_31_int->SaveAs("c_marginal_xi_31_int.pdf");
        c_marginal_xi_31_int->SaveAs("c_marginal_xi_31_int.png");






        ///////////////////////////////////////////////////////////////////////
        // 2 sided limit
        ///////////////////////////////////////////////////////////////////////

        TH1D *h_cumulative_xi_31_low_int = nullptr;
        TH1D *h_cumulative_xi_31_high_int = nullptr;

        double x_low, x_mid, x_high;
        double x_range;
        Int_t n_bins;

        if(is_even(h_marginal_xi_31->GetNbinsX()))
        {
            std::cout << "EVEN" << std::endl;

            x_low = h_marginal_xi_31->GetXaxis()->GetBinLowEdge(1);
            x_mid = h_marginal_xi_31->GetXaxis()->GetBinUpEdge(1 + (h_marginal_xi_31->GetNbinsX() / 2));
            x_high = h_marginal_xi_31->GetXaxis()->GetBinUpEdge(h_marginal_xi_31->GetNbinsX());
            n_bins = h_marginal_xi_31->GetNbinsX() / 2;

            std::cout << "n_bins=" << n_bins << " x_low=" << x_low << " x_mid=" << x_mid << " x_high=" << x_high << std::endl;

            h_cumulative_xi_31_low_int = 
                new TH1D("h_cumulative_xi_31_low_int", "h_cumulative_xi_31_low_int",
                         n_bins, x_low, x_mid);

            h_cumulative_xi_31_high_int = 
                new TH1D("h_cumulative_xi_31_high_int", "h_cumulative_xi_31_high_int",
                         n_bins, x_mid, x_high);
        }
        else
        {
            std::cout << "ODD" << std::endl;

            x_low = h_marginal_xi_31->GetXaxis()->GetBinLowEdge(1);
            Int_t index_mid = peak_i;
            //Int_t index_mid = max_bin_y;
            std::cout << ">>>>>>>> " << "peak_i=" << peak_i << " max_bin_x=" << max_bin_x << std::endl;
            x_mid = h_marginal_xi_31->GetXaxis()->GetBinCenter(index_mid);
            x_high = h_marginal_xi_31->GetXaxis()->GetBinUpEdge(h_marginal_xi_31->GetNbinsX());
            x_range = x_high - x_low;

            std::cout << "x_low=" << x_low << " x_mid=" << x_mid << " x_high=" << x_high << std::endl;

            h_cumulative_xi_31_low_int = 
                new TH1D("h_cumulative_xi_31_low_int", "h_cumulative_xi_31_low_int",
                         h_marginal_xi_31->GetNbinsX(), x_mid - x_range, x_mid);
            h_cumulative_xi_31_low_int->SetTitle(0);
            h_cumulative_xi_31_low_int->SetStats(0);

            h_cumulative_xi_31_high_int = 
                new TH1D("h_cumulative_xi_31_high_int", "h_cumulative_xi_31_high_int",
                         h_marginal_xi_31->GetNbinsX(), x_mid, x_mid + x_range);
            h_cumulative_xi_31_high_int->SetTitle(0);
            h_cumulative_xi_31_high_int->SetStats(0);


            double sum_low = 0.0;
            double sum_high = 0.0;
            //Int_t index_mid = (h_marginal_xi_31->GetNbinsX() + 1) / 2;
            std::cout << "index_mid=" << index_mid << std::endl;
            sum_low += 0.5 * h_marginal_xi_31->GetBinContent(index_mid) * bin_width_x;
            sum_high += 0.5 * h_marginal_xi_31->GetBinContent(index_mid) * bin_width_x;
            for(Int_t i = 1; i <= h_marginal_xi_31->GetNbinsX(); ++ i)
            {
                Int_t index_low = index_mid - i;
                Int_t index_high = index_mid + i;

                //std::cout << "index: " << index_low << " " << index_high << std::endl;
                if(index_low >= 1)
                {
                    sum_low += h_marginal_xi_31->GetBinContent(index_low) * bin_width_x;
                }
                if(i <= h_cumulative_xi_31_low_int->GetNbinsX())
                {   
                    h_cumulative_xi_31_low_int->SetBinContent(h_cumulative_xi_31_low_int->GetNbinsX() - (i - 1), sum_low);
                }

                if(index_high <= h_marginal_xi_31->GetNbinsX())
                {
                    sum_high += h_marginal_xi_31->GetBinContent(index_high) * bin_width_x;
                }
                if(i <= h_cumulative_xi_31_high_int->GetNbinsX())
                {
                    h_cumulative_xi_31_high_int->SetBinContent(i, sum_high);
                }
            }
        }

        TCanvas *c_cumulative_xi_31_low_int = new TCanvas("c_cumulative_xi_31_low_int", "c_cumulative_xi_31_low_int");
        h_cumulative_xi_31_low_int->Draw("hist");

        TCanvas *c_cumulative_xi_31_high_int = new TCanvas("c_cumulative_xi_31_high_int", "c_cumulative_xi_31_high_int");
        h_cumulative_xi_31_high_int->Draw("hist");

        // cumulative reverse low

        TGraph *g_cumulative_xi_31_low_int_reverse = new TGraph(h_cumulative_xi_31_low_int->GetNbinsX());
        g_cumulative_xi_31_low_int_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_xi_31_low_int_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_xi_31_low_int->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_xi_31_low_int->GetBinContent(i + 1);
            g_cumulative_xi_31_low_int_reverse->SetPoint(i, y, x);
        }

        TCanvas *c_cumulative_xi_31_low_int_reverse = new TCanvas("c_cumulative_xi_31_low_int_reverse", "c_cumulative_xi_31_low_int_reverse");
        //TCanvas *c_cumulative_xi_31_high_intlow_reverse = new TCanvas("c_cumulative_xi_31_high_intlow_reverse", "c_cumulative_xi_31_high_intlow_reverse");
        g_cumulative_xi_31_low_int_reverse->Draw("al");

        std::cout << "1 sigma: " << g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.683) << std::endl;
        std::cout << "90 %: " << g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.90) << std::endl;
        //std::cout << "95 %: " << g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.95) << std::endl;
        //std::cout << "99 %: " << g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.99) << std::endl;

        std::cout << "*** bin index ***" << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.683)) << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.90)) << std::endl;
        //std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.95)) << std::endl;
        //std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.99)) << std::endl;


        // cumulative reverse high

        TGraph *g_cumulative_xi_31_high_int_reverse = new TGraph(h_cumulative_xi_31_high_int->GetNbinsX());
        g_cumulative_xi_31_high_int_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_xi_31_high_int_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_xi_31_high_int->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_xi_31_high_int->GetBinContent(i + 1);
            g_cumulative_xi_31_high_int_reverse->SetPoint(i, y, x);
        }

        TCanvas *c_cumulative_xi_31_high_int_reverse = new TCanvas("c_cumulative_xi_31_high_int_reverse", "c_cumulative_xi_31_high_int_reverse");
        g_cumulative_xi_31_high_int_reverse->Draw("al");
        //g_cumulative_xi_31_high_int_reverse->Draw("lsame");

        std::cout << "1 sigma: " << g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.683) << std::endl;
        std::cout << "90 %: " << g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.90) << std::endl;
        //std::cout << "95 %: " << g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.95) << std::endl;
        //std::cout << "99 %: " << g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.99) << std::endl;

        std::cout << "*** bin index ***" << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.683)) << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.90)) << std::endl;
        //std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.95)) << std::endl;
        //std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.99)) << std::endl;

        std::cout << "1 sigma: xi_31 = " << peak_x << " " 
                  << "+ " << (g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.683) - peak_x) << " "
                  << "- " << (peak_x - g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.683)) << std::endl;
        std::cout << "90 %: xi_31 = " << peak_x << " " 
                  << "+ " << (g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.90) - peak_x) << " "
                  << "- " << (peak_x - g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.90)) << std::endl;
        //std::cout << "95 %: xi_31 = " << peak_x << " " 
        //          << "+ " << (g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.95) - peak_x) << " "
        //          << "- " << (peak_x - g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.95)) << std::endl;
        //std::cout << "99 %: xi_31 = " << peak_x << " " 
        //          << "+ " << (g_cumulative_xi_31_high_int_reverse->Eval(0.5 * 0.99) - peak_x) << " "
        //          << "- " << (peak_x - g_cumulative_xi_31_low_int_reverse->Eval(0.5 * 0.99)) << std::endl;


        ///////////////////////////////////////////////////////////////////////
        // cumulative interiour
        ///////////////////////////////////////////////////////////////////////

#if 0

        TH1D *h_cumulative_xi_31_low_int = new TH1D("h_cumulative_xi_31_low_int", "h_cumulative_xi_31_low_int",
                                      h_marginal_xi_31->GetNbinsX(),
                                      h_marginal_xi_31->GetXaxis()->GetBinLowEdge(1),
                                      h_marginal_xi_31->GetXaxis()->GetBinUpEdge(h_marginal_xi_31->GetNbinsX()));
        h_cumulative_xi_31_low_int->SetTitle(0);
        h_cumulative_xi_31_low_int->SetStats(0);

        sum = 0.0;
        for(Int_t i = 1; i <= h_marginal_xi_31->GetNbinsX(); ++ i)
        {
            sum += h_marginal_xi_31->GetBinContent(i) * bin_width_x;
            h_cumulative_xi_31_low_int->SetBinContent(i, sum);
        }
/*
        TCanvas *c_cumulative_xi_31_low_int = new TCanvas("c_cumulative_xi_31_low_int", "c_cumulative_xi_31_low_int");
        h_cumulative_xi_31_low_int->Draw("hist");
*/

        TH1D *h_cumulative_xi_31_high_int = new TH1D("h_cumulative_xi_31_high_int", "h_cumulative_xi_31_high_int",
                                      h_marginal_xi_31->GetNbinsX(),
                                      h_marginal_xi_31->GetXaxis()->GetBinLowEdge(1),
                                      h_marginal_xi_31->GetXaxis()->GetBinUpEdge(h_marginal_xi_31->GetNbinsX()));
        h_cumulative_xi_31_high_int->SetTitle(0);
        h_cumulative_xi_31_high_int->SetStats(0);

        sum = 0.0;
        for(Int_t i = h_marginal_xi_31->GetNbinsX(); i >= 1; -- i)
        {
            sum += h_marginal_xi_31->GetBinContent(i) * bin_width_x;
            h_cumulative_xi_31_high_int->SetBinContent(i, sum);
        }
/*
        TCanvas *c_cumulative_xi_31_high_int = new TCanvas("c_cumulative_xi_31_high_int", "c_cumulative_xi_31_high_int");
        h_cumulative_xi_31_high_int->Draw("hist");
*/



        TGraph *g_cumulative_xi_31_low_int_reverse = new TGraph(h_cumulative_xi_31_low_int->GetNbinsX());
        g_cumulative_xi_31_low_int_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_xi_31_low_int_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_xi_31_low_int->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_xi_31_low_int->GetBinContent(i + 1);
            g_cumulative_xi_31_low_int_reverse->SetPoint(i, y, x);
        }
/*
        TCanvas *c_cumulative_xi_31_low_int_reverse = new TCanvas("c_cumulative_xi_31_low_int_reverse", "c_cumulative_xi_31_low_int_reverse");
        g_cumulative_xi_31_low_int_reverse->Draw("al");
*/


        TGraph *g_cumulative_xi_31_high_int_reverse = new TGraph(h_cumulative_xi_31_high_int->GetNbinsX());
        g_cumulative_xi_31_high_int_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_xi_31_high_int_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_xi_31_high_int->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_xi_31_high_int->GetBinContent(i + 1);
            g_cumulative_xi_31_high_int_reverse->SetPoint(i, y, x);
        }
/*
        TCanvas *c_cumulative_xi_31_high_int_reverse = new TCanvas("c_cumulative_xi_31_high_int_reverse", "c_cumulative_xi_31_high_int_reverse");
        g_cumulative_xi_31_high_int_reverse->Draw("al");
*/

#endif

        std::cout << "1 sigma (low int): " << g_cumulative_xi_31_low_int_reverse->Eval((1.0 - 0.683) / 2.0) << std::endl;
        std::cout << "90 % (low int): " << g_cumulative_xi_31_low_int_reverse->Eval((1.0 - 0.90) / 2.0) << std::endl;
        std::cout << "1 sigma (high int): " << g_cumulative_xi_31_high_int_reverse->Eval((1.0 - 0.683) / 2.0) << std::endl;
        std::cout << "90 % (high int): " << g_cumulative_xi_31_high_int_reverse->Eval((1.0 - 0.90) / 2.0) << std::endl;

        std::cout << "*** bin index ***" << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_low_int_reverse->Eval((1.0 - 0.683) / 2.0)) << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_low_int_reverse->Eval((1.0 - 0.90) / 2.0)) << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_high_int_reverse->Eval((1.0 - 0.683) / 2.0)) << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_high_int_reverse->Eval((1.0 - 0.90) / 2.0)) << std::endl;







    }
#endif



    ///////////////////////////////////////////////////////////////////////////
    // marginal T12
    ///////////////////////////////////////////////////////////////////////////

    TH1D *h_marginal_T12 = new TH1D("h_marginal_T12", "h_marginal_T12",
                                h_prob->GetNbinsY(),
                                h_prob->GetYaxis()->GetBinLowEdge(1),
                                h_prob->GetYaxis()->GetBinUpEdge(h_prob->GetNbinsY()));

    {
        const double bin_width_x = PROB_X_RANGE / 301.0;
        const double bin_width_y = PROB_Y_RANGE / 301.0;

        double peak_y = 0.0;
        double peak_x = 0.0;
        Int_t peak_j = 0;
        for(Int_t j = 1; j <= h_prob->GetNbinsY(); ++ j)
        {
            double sum = 0.0;
            for(Int_t i = 1; i <= h_prob->GetNbinsX(); ++ i)
            {
                sum += h_prob->GetBinContent(i, j) * bin_width_x;
            }
            h_marginal_T12->SetBinContent(j, sum);

            if((j == 1) || (sum > peak_y))
            {
                peak_y = sum;
                peak_x = h_prob->GetYaxis()->GetBinCenter(j);
                peak_j = j;
            }
        }
        //std::cout << "sum=" << sum << " ->Scale(" << 1.0 / (sum * bin_width) << ")" << std::endl;
        //h_marginal_T12->Scale(1.0 / (sum * bin_width));

        std::cout << "peak: " << peak_x << " " << peak_y << std::endl;
        std::cout << "peak_j=" << peak_j << std::endl;

        sum = 0.0;
        for(Int_t i = 1; i <= h_marginal_T12->GetNbinsX(); ++ i)
        {
            sum += h_marginal_T12->GetBinContent(i);
        }
        std::cout << "integral in 1D P-space (T21): " << sum * bin_width_y << std::endl;



        ///////////////////////////////////////////////////////////////////////
        // compute mean and variance of marginal T12
        ///////////////////////////////////////////////////////////////////////

        const double marginal_T12_mean = peak_x;
        std::cout << "Mean: " << marginal_T12_mean << std::endl;
        const double marginal_T12_var = variance_x_fixed_mean(h_marginal_T12, marginal_T12_mean);
        std::cout << "Var: " << marginal_T12_var << std::endl;
        const double marginal_T12_se = std::sqrt(marginal_T12_var);
        std::cout << "Standard Error: " << marginal_T12_se << std::endl;
        std::cout << "Alternative mean: " << mean_x(h_marginal_T12) << std::endl;



        TCanvas *c_marginal_T12 = new TCanvas("c_marginal_T12", "c_marginal_T12");
        c_marginal_T12->SetTicks(1, 2);
        c_marginal_T12->SetBottomMargin(0.15);
        c_marginal_T12->SetRightMargin(0.1);
        c_marginal_T12->SetLeftMargin(0.12);
        h_marginal_T12->SetTitle("");
        h_marginal_T12->SetStats(0);
        h_marginal_T12->GetZaxis()->SetLabelOffset(0.005);
        h_marginal_T12->GetXaxis()->SetLabelSize(17.0);
        h_marginal_T12->GetXaxis()->SetLabelFont(43);
        h_marginal_T12->GetYaxis()->SetLabelSize(17.0);
        h_marginal_T12->GetYaxis()->SetLabelFont(43);
        h_marginal_T12->GetZaxis()->SetLabelSize(17.0);
        h_marginal_T12->GetZaxis()->SetLabelFont(43);
        h_marginal_T12->GetXaxis()->SetTitleSize(18.0);
        h_marginal_T12->GetXaxis()->SetTitleFont(43);
        h_marginal_T12->GetYaxis()->SetTitleSize(18.0);
        h_marginal_T12->GetYaxis()->SetTitleFont(43);
        h_marginal_T12->GetXaxis()->SetTitleOffset(1.5);
        h_marginal_T12->GetYaxis()->SetTitleOffset(1.2);
        h_marginal_T12->GetXaxis()->SetLabelOffset(0.01);
        h_marginal_T12->GetYaxis()->SetLabelOffset(0.01);
        h_marginal_T12->SetXTitle("^{150}Nd Amplitude Scale Factor (#times 3.45#times10^{-4} Bq)");
        h_marginal_T12->SetYTitle("Probability Density Function p(A^{^{150}Nd})");
        std::cout << h_marginal_T12->GetBinWidth(1) << std::endl;
        std::cout << "INT (T12) = " << h_marginal_T12->Integral() * h_marginal_T12->GetBinWidth(1) << std::endl;
        TH1D *h_marginal_T12_clone = (TH1D*)h_marginal_T12->Clone();
    //    h_marginal_T12_clone->Scale(1.0 / (2.0 / 301.0));
    //    for(Int_t i = 0; i <= h_marginal_T12->GetNbinsX() + 1; ++ i)
        for(Int_t i = 1; i <= h_marginal_T12->GetNbinsX(); ++ i)
        {
            const double content = h_marginal_T12_clone->GetBinContent(i);
    //        h_marginal_T12_clone->SetBinContent(i, 301.0 / 2.0 * content);
    //        std::cout << "i=" << i << " " << h_marginal_T12->GetBinContent(i) << std::endl;
        }
        h_marginal_T12_clone->SetLineWidth(1);
        h_marginal_T12_clone->SetLineColor(kBlack);
        h_marginal_T12_clone->SetMaximum(6.0);
        h_marginal_T12_clone->Draw("axis");


//        const double sigma_estimated = std::sqrt(1.0 / (2.0 * 4.0 * std::atan(1.0) * std::pow(peak_y, 2.0)));
        //const double sigma_estimated = 0.0707099; // measured at end
        const double sigma_measured = marginal_T12_se;
        {
            double sum = 0.0;
            const double width = h_marginal_T12_clone->GetXaxis()->GetBinWidth(1);
            std::cout << "The estimate of sigma from the height is: "
                      << 1.0 / (std::sqrt(2.0 * 4.0 * std::atan(1.0)) * peak_y) << std::endl;
            //int low = h_marginal_T12_clone->GetXaxis()->FindBin(peak_x - 0.0694279);
            //int high = h_marginal_T12_clone->GetXaxis()->FindBin(peak_x + 0.0694279);
            int low = h_marginal_T12_clone->GetXaxis()->FindBin(peak_x - sigma_measured);
            int high = h_marginal_T12_clone->GetXaxis()->FindBin(peak_x + sigma_measured);

            //std::cout << "low=" << low << std::endl;
            //std::cout << "bin limits: " << h_marginal_T12_clone->GetXaxis()->GetBinLowEdge(low) 
            //          << " " << h_marginal_T12_clone->GetXaxis()->GetBinUpEdge(low)
            //          << " value=" << peak_x - sigma_estimated << std::endl;

            for(Int_t i = low; i <= high; ++ i)
            {
                if((i == low) || (i == high))
                {
                    sum += 0.5 * h_marginal_T12_clone->GetBinContent(i) * width;
                }
                else
                {
                    sum += h_marginal_T12_clone->GetBinContent(i) * width;
                }
            }
            //std::cout << "Integration between " << peak_x << " +- 0.0694279: " << sum << std::endl;
            std::cout << "Integration between " << peak_x << " +- " << sigma_measured << ": " << sum << std::endl;

            double sum2 = 0.0;
            int low2 = h_marginal_T12_clone->GetXaxis()->FindBin(peak_x - 1.64 * sigma_measured);
            int high2 = h_marginal_T12_clone->GetXaxis()->FindBin(peak_x + 1.64 * sigma_measured);
            for(Int_t i = low2; i <= high2; ++ i)
            {
                if((i == low2) || (i == high2))
                {
                    sum2 += 0.5 * h_marginal_T12_clone->GetBinContent(i) * width;
                }
                else
                {
                    sum2 += h_marginal_T12_clone->GetBinContent(i) * width;
                }
            }
            std::cout << "Integration between " << peak_x << " +- " << 1.64 * sigma_measured << ": " << sum2 << std::endl;


            double sum3 = 0.0;
            int low3 = h_marginal_T12_clone->GetXaxis()->FindBin(peak_x - 2.58 * sigma_measured);
            int high3 = h_marginal_T12_clone->GetXaxis()->FindBin(peak_x + 2.58 * sigma_measured);
            for(Int_t i = low3; i <= high3; ++ i)
            {
                if((i == low3) || (i == high3))
                {
                    sum3 += 0.5 * h_marginal_T12_clone->GetBinContent(i) * width;
                }
                else
                {
                    sum3 += h_marginal_T12_clone->GetBinContent(i) * width;
                }
            }
            std::cout << "Integration between " << peak_x << " +- " << 2.58 * sigma_measured << ": " << sum3 << std::endl;
            
        }


        /*
        TH1D *h_marginal_T12_estimated = new TH1D("h_marginal_T12_estimated", "h_marginal_T12_estimated",
                                                  h_marginal_T12->GetNbinsX(),
                                                  h_marginal_T12->GetXaxis()->GetBinLowEdge(1),
                                                  h_marginal_T12->GetXaxis()->GetBinUpEdge(h_marginal_T12->GetNbinsX()));
        */
        TGraph *g_marginal_T12_estimated = new TGraph(h_marginal_T12->GetNbinsX());

        std::vector<std::pair<double, double>> p_marginal_T12_estimated_A;
        std::vector<std::pair<double, double>> p_marginal_T12_estimated_B;
        std::vector<std::pair<double, double>> p_marginal_T12_estimated_C;
        std::vector<std::pair<double, double>> p_marginal_T12_estimated_D;
        std::vector<std::pair<double, double>> p_marginal_T12_estimated_E;
        std::vector<std::pair<double, double>> p_marginal_T12_estimated_F;

        //std::cout << "peak_y=" << peak_y << " predicted peak: " << 1.0 / std::sqrt(2.0 * 4.0 * std::atan(1.0) * std::pow(0.0694279, 2.0)) << std::endl;
        std::cout << "peak_y=" << peak_y << " predicted peak: " << 1.0 / std::sqrt(2.0 * 4.0 * std::atan(1.0) * std::pow(sigma_measured, 2.0)) << std::endl;
        for(Int_t i = 0; i < g_marginal_T12_estimated->GetN(); ++ i)
        {
            double x = h_marginal_T12_clone->GetXaxis()->GetBinCenter(i + 1);
//            double content = peak_y * std::exp(-std::pow((x - 1.19003) / (0.240864 * 0.4), 2.0));
//            double content = peak_y * std::exp(-std::pow((x - 1.19003) / 0.0580154, 2.0));
//            double content = peak_y * std::exp(-std::pow((x - 1.19003) / 0.0694279, 2.0));
            //double content = peak_y * std::exp(-std::pow((x - 1.19003) / (0.112576 * 0.9), 2.0));
            //double content = peak_y * std::exp(-std::pow((x - 1.19333) / (0.112576 * 0.88), 2.0));
//            double content = peak_y * std::exp(-0.5 * std::pow((x - 1.19333) / (0.0694279), 2.0));
            //double content = (1.0 / (std::sqrt(2.0 * 4.0 * std::atan(1.0) * std::pow(0.0694279, 2.0)))) * std::exp(-0.5 * std::pow((x - 1.19333) / (0.0694279), 2.0));
            double content = (1.0 / (std::sqrt(2.0 * 4.0 * std::atan(1.0)) * sigma_measured))
                * std::exp(-0.5 * std::pow((x - peak_x) / (sigma_measured), 2.0));


            //h_marginal_T12_estimated->SetBinContent(i, content);
            g_marginal_T12_estimated->SetPoint(i, x, content);

            if((peak_x - sigma_measured < x) && (x < peak_x + sigma_measured))
            {
                p_marginal_T12_estimated_A.push_back(std::make_pair(x, content));
            }

            if((peak_x - 1.64 * sigma_measured  < x) && (x <= peak_x - sigma_measured))
            {
                p_marginal_T12_estimated_C.push_back(std::make_pair(x, content));
            }

            if((peak_x + sigma_measured <= x) && (x < peak_x + 1.64 * sigma_measured))
            {
                p_marginal_T12_estimated_D.push_back(std::make_pair(x, content));
            }

            if((peak_x - 2.58 * sigma_measured  < x) && (x <= peak_x - 1.64 * sigma_measured))
            {
                p_marginal_T12_estimated_E.push_back(std::make_pair(x, content));
            }

            if((peak_x + 1.64 * sigma_measured <= x) && (x < peak_x + 2.58 * sigma_measured))
            {
                p_marginal_T12_estimated_F.push_back(std::make_pair(x, content));
            }


        }

        TGraph *g_marginal_T12_estimated_A = new TGraph(p_marginal_T12_estimated_A.size() + 2);
        TGraph *g_marginal_T12_estimated_B = new TGraph(p_marginal_T12_estimated_A.size() + 2);
        TGraph *g_marginal_T12_estimated_C = new TGraph(p_marginal_T12_estimated_C.size() + 2);
        TGraph *g_marginal_T12_estimated_D = new TGraph(p_marginal_T12_estimated_D.size() + 2);
        TGraph *g_marginal_T12_estimated_E = new TGraph(p_marginal_T12_estimated_E.size() + 2);
        TGraph *g_marginal_T12_estimated_F = new TGraph(p_marginal_T12_estimated_F.size() + 2);

        g_marginal_T12_estimated_A->SetPoint(0, p_marginal_T12_estimated_A.front().first, 0.0);
        g_marginal_T12_estimated_C->SetPoint(0, p_marginal_T12_estimated_C.front().first, 0.0);
        g_marginal_T12_estimated_D->SetPoint(0, p_marginal_T12_estimated_D.front().first, 0.0);
        g_marginal_T12_estimated_E->SetPoint(0, p_marginal_T12_estimated_E.front().first, 0.0);
        g_marginal_T12_estimated_F->SetPoint(0, p_marginal_T12_estimated_F.front().first, 0.0);

        for(int i = 0; i < p_marginal_T12_estimated_A.size(); ++ i)
        {
            g_marginal_T12_estimated_A->SetPoint(i + 1,
                p_marginal_T12_estimated_A.at(i).first,
                p_marginal_T12_estimated_A.at(i).second);
        }

        for(int i = 0; i < p_marginal_T12_estimated_C.size(); ++ i)
        {
            g_marginal_T12_estimated_C->SetPoint(i + 1,
                p_marginal_T12_estimated_C.at(i).first,
                p_marginal_T12_estimated_C.at(i).second);
        }

        for(int i = 0; i < p_marginal_T12_estimated_D.size(); ++ i)
        {
            g_marginal_T12_estimated_D->SetPoint(i + 1,
                p_marginal_T12_estimated_D.at(i).first,
                p_marginal_T12_estimated_D.at(i).second);
        }

        for(int i = 0; i < p_marginal_T12_estimated_E.size(); ++ i)
        {
            g_marginal_T12_estimated_E->SetPoint(i + 1,
                p_marginal_T12_estimated_E.at(i).first,
                p_marginal_T12_estimated_E.at(i).second);
        }

        for(int i = 0; i < p_marginal_T12_estimated_F.size(); ++ i)
        {
            g_marginal_T12_estimated_F->SetPoint(i + 1,
                p_marginal_T12_estimated_F.at(i).first,
                p_marginal_T12_estimated_F.at(i).second);
        }
 
        g_marginal_T12_estimated_A->SetPoint(g_marginal_T12_estimated_A->GetN() - 1
            , p_marginal_T12_estimated_A.back().first, 0.0);

        g_marginal_T12_estimated_C->SetPoint(g_marginal_T12_estimated_C->GetN() - 1
            , p_marginal_T12_estimated_C.back().first, 0.0);
        
        g_marginal_T12_estimated_D->SetPoint(g_marginal_T12_estimated_D->GetN() - 1
            , p_marginal_T12_estimated_D.back().first, 0.0);

        g_marginal_T12_estimated_E->SetPoint(g_marginal_T12_estimated_E->GetN() - 1
            , p_marginal_T12_estimated_E.back().first, 0.0);
        
        g_marginal_T12_estimated_F->SetPoint(g_marginal_T12_estimated_F->GetN() - 1
            , p_marginal_T12_estimated_F.back().first, 0.0);


        {
            double x1, x2, y1, y2, x, y;

            g_marginal_T12_estimated_C->GetPoint(g_marginal_T12_estimated_C->GetN() - 2, x1, y1);
            g_marginal_T12_estimated_A->GetPoint(1, x2, y2);
            //std::cout << "x1=" << x1 << " y1=" << y1 << std::endl;
            //std::cout << "x2=" << x2 << " y2=" << y2 << std::endl;
            x = 0.5 * (x1 + x2);
            y = 0.5 * (y1 + y2);
            g_marginal_T12_estimated_C->SetPoint(g_marginal_T12_estimated_C->GetN() - 2, x, y);
            g_marginal_T12_estimated_A->SetPoint(1, x, y);


            g_marginal_T12_estimated_C->GetPoint(g_marginal_T12_estimated_C->GetN() - 1, x1, y1);
            g_marginal_T12_estimated_A->GetPoint(0, x2, y2);
            //std::cout << "x1=" << x1 << " y1=" << y1 << std::endl;
            //std::cout << "x2=" << x2 << " y2=" << y2 << std::endl;
            x = 0.5 * (x1 + x2);
            y = 0.5 * (y1 + y2);
            g_marginal_T12_estimated_C->SetPoint(g_marginal_T12_estimated_C->GetN() - 1, x, y);
            g_marginal_T12_estimated_A->SetPoint(0, x, y);


            g_marginal_T12_estimated_A->GetPoint(g_marginal_T12_estimated_A->GetN() - 2, x1, y1);
            g_marginal_T12_estimated_D->GetPoint(1, x2, y2);
            //std::cout << "x1=" << x1 << " y1=" << y1 << std::endl;
            //std::cout << "x2=" << x2 << " y2=" << y2 << std::endl;
            x = 0.5 * (x1 + x2);
            y = 0.5 * (y1 + y2);
            g_marginal_T12_estimated_A->SetPoint(g_marginal_T12_estimated_A->GetN() - 2, x, y);
            g_marginal_T12_estimated_D->SetPoint(1, x, y);


            g_marginal_T12_estimated_A->GetPoint(g_marginal_T12_estimated_A->GetN() - 1, x1, y1);
            g_marginal_T12_estimated_D->GetPoint(0, x2, y2);
            //std::cout << "x1=" << x1 << " y1=" << y1 << std::endl;
            //std::cout << "x2=" << x2 << " y2=" << y2 << std::endl;
            x = 0.5 * (x1 + x2);
            y = 0.5 * (y1 + y2);
            g_marginal_T12_estimated_A->SetPoint(g_marginal_T12_estimated_A->GetN() - 1, x, y);
            g_marginal_T12_estimated_D->SetPoint(0, x, y);


            g_marginal_T12_estimated_E->GetPoint(g_marginal_T12_estimated_E->GetN() - 2, x1, y1);
            g_marginal_T12_estimated_C->GetPoint(1, x2, y2);
            //std::cout << "x1=" << x1 << " y1=" << y1 << std::endl;
            //std::cout << "x2=" << x2 << " y2=" << y2 << std::endl;
            x = 0.5 * (x1 + x2);
            y = 0.5 * (y1 + y2);
            g_marginal_T12_estimated_E->SetPoint(g_marginal_T12_estimated_E->GetN() - 2, x, y);
            g_marginal_T12_estimated_C->SetPoint(1, x, y);


            g_marginal_T12_estimated_E->GetPoint(g_marginal_T12_estimated_E->GetN() - 1, x1, y1);
            g_marginal_T12_estimated_C->GetPoint(0, x2, y2);
            //std::cout << "x1=" << x1 << " y1=" << y1 << std::endl;
            //std::cout << "x2=" << x2 << " y2=" << y2 << std::endl;
            x = 0.5 * (x1 + x2);
            y = 0.5 * (y1 + y2);
            g_marginal_T12_estimated_E->SetPoint(g_marginal_T12_estimated_E->GetN() - 1, x, y);
            g_marginal_T12_estimated_C->SetPoint(0, x, y);


            g_marginal_T12_estimated_D->GetPoint(g_marginal_T12_estimated_D->GetN() - 2, x1, y1);
            g_marginal_T12_estimated_F->GetPoint(1, x2, y2);
            //std::cout << "x1=" << x1 << " y1=" << y1 << std::endl;
            //std::cout << "x2=" << x2 << " y2=" << y2 << std::endl;
            x = 0.5 * (x1 + x2);
            y = 0.5 * (y1 + y2);
            g_marginal_T12_estimated_D->SetPoint(g_marginal_T12_estimated_D->GetN() - 2, x, y);
            g_marginal_T12_estimated_F->SetPoint(1, x, y);


            g_marginal_T12_estimated_D->GetPoint(g_marginal_T12_estimated_D->GetN() - 1, x1, y1);
            g_marginal_T12_estimated_F->GetPoint(0, x2, y2);
            //std::cout << "x1=" << x1 << " y1=" << y1 << std::endl;
            //std::cout << "x2=" << x2 << " y2=" << y2 << std::endl;
            x = 0.5 * (x1 + x2);
            y = 0.5 * (y1 + y2);
            g_marginal_T12_estimated_D->SetPoint(g_marginal_T12_estimated_D->GetN() - 1, x, y);
            g_marginal_T12_estimated_F->SetPoint(0, x, y);
        }

        g_marginal_T12_estimated_A->SetFillColor(kRed);
        g_marginal_T12_estimated_A->SetFillStyle(3002);
        g_marginal_T12_estimated_A->SetLineColor(kRed);
        g_marginal_T12_estimated_A->SetLineWidth(0);
        g_marginal_T12_estimated_A->Draw("Lfsame");

        g_marginal_T12_estimated_C->SetFillColor(kViolet);
        g_marginal_T12_estimated_C->SetFillStyle(3002);
        g_marginal_T12_estimated_C->SetLineColor(kViolet);
        g_marginal_T12_estimated_C->SetLineWidth(0);
        g_marginal_T12_estimated_C->Draw("Lfsame");

        g_marginal_T12_estimated_D->SetFillColor(kViolet);
        g_marginal_T12_estimated_D->SetFillStyle(3002);
        g_marginal_T12_estimated_D->SetLineColor(kViolet);
        g_marginal_T12_estimated_D->SetLineWidth(0);
        g_marginal_T12_estimated_D->Draw("Lfsame");

        g_marginal_T12_estimated_E->SetFillColor(kBlue - 4);
        g_marginal_T12_estimated_E->SetFillStyle(3002);
        g_marginal_T12_estimated_E->SetLineColor(kBlue - 4);
        g_marginal_T12_estimated_E->SetLineWidth(0);
        g_marginal_T12_estimated_E->Draw("Lfsame");

        g_marginal_T12_estimated_F->SetFillColor(kBlue - 4);
        g_marginal_T12_estimated_F->SetFillStyle(3002);
        g_marginal_T12_estimated_F->SetLineColor(kBlue - 4);
        g_marginal_T12_estimated_F->SetLineWidth(0);
        g_marginal_T12_estimated_F->Draw("Lfsame");

        //h_marginal_T12_estimated->SetLineColor(kRed);
        //h_marginal_T12_estimated->Draw("histsame");
        g_marginal_T12_estimated->SetLineColor(kRed);
        g_marginal_T12_estimated->SetLineWidth(2);
        g_marginal_T12_estimated->Draw("lsame");

        h_marginal_T12_clone->Draw("histsame");

        TH2D *h_marginal_T12_clone_clone = (TH2D*)h_marginal_T12_clone->Clone();
        TGraph *g_marginal_T12_estimated_clone = (TGraph*)g_marginal_T12_estimated->Clone();
        h_marginal_T12_clone_clone->SetLineWidth(5);
        g_marginal_T12_estimated_clone->SetLineWidth(5);

        TLegend *leg = new TLegend(0.6, 0.75, 0.86, 0.86);
        leg->AddEntry(h_marginal_T12_clone_clone, "Marginal Distribution", "L");
        leg->AddEntry(g_marginal_T12_estimated_clone, "Parameterized Model", "L");
        leg->SetBorderSize(5);
        leg->SetShadowColor(kGray + 2);
        leg->SetTextFont(63);
        leg->SetTextSize(15);
        leg->SetMargin(0.2);
        leg->Draw("BR");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextColor(kRed);
        latex.SetTextFont(43);
        latex.SetTextSize(20);

        TString s2;
        s2.Form("#mu = %.4f", peak_x);
        TString s1;
        s1.Form("#sigma = %.4f", sigma_measured);

        latex.DrawLatex(0.72, 0.55, s2);
        latex.DrawLatex(0.72, 0.48, s1);

        TString s3;
        s3.Form("#int_{0.95}^{1.45} P(A^{^{150}Nd}) dA^{^{150}Nd} = 1");
        latex.SetTextColor(kBlack);
        latex.SetTextSize(15);
        latex.DrawLatex(0.65, 0.65, s3);

        TString s4;
        s4.Form("#int_{%.3f}^{%.3f} P(A^{^{150}Nd}) dA^{^{150}Nd} = 1#sigma",
            peak_x - sigma_measured,
            peak_x + sigma_measured);
        latex.SetTextColor(kRed - 4);
        latex.DrawLatex(0.15, 0.55, s4);

        TString s5;
        s5.Form("#int_{%.3f}^{%.3f} P(A^{^{150}Nd}) dA^{^{150}Nd} = 90 %%",
            peak_x - 1.64 * sigma_measured,
            peak_x + 1.64 * sigma_measured);
        latex.SetTextColor(kViolet);
        latex.DrawLatex(0.15, 0.67, s5);

        TString s6;
        s6.Form("#int_{%.3f}^{%.3f} P(A^{^{150}Nd}) dA^{^{150}Nd} = 99 %%",
            peak_x - 2.58 * sigma_measured,
            peak_x + 2.58 * sigma_measured);
        latex.SetTextColor(kBlue - 4);
        latex.DrawLatex(0.15, 0.79, s6);

        c_marginal_T12->SaveAs("c_marginal_T12_2.eps");
        c_marginal_T12->SaveAs("c_marginal_T12_2.pdf");
        c_marginal_T12->SaveAs("c_marginal_T12_2.png");


        std::cout << "table of results for latex" << std::endl;
        std::cout << peak_x << std::endl;
        std::cout << peak_x - sigma_measured << " & " << peak_x + sigma_measured << " & "
                  << sigma_measured << std::endl;
        std::cout << peak_x - 1.66 * sigma_measured << " & " << peak_x + 1.66 * sigma_measured << " & "
                  << 1.66 * sigma_measured << std::endl;
        std::cout << peak_x - 1.96 * sigma_measured << " & " << peak_x + 1.96 * sigma_measured << " & "
                  <<  1.96 *sigma_measured << std::endl;
        std::cout << peak_x - 2.58 * sigma_measured << " & " << peak_x + 2.58 * sigma_measured << " & "
                  <<  2.58 *sigma_measured << std::endl;


        /*
        TH1D *h_marginal_T12_clone1 = (TH1D*)h_marginal_T12->Clone();
        h_marginal_T12_clone1->SetFillColor(kRed);
        h_marginal_T12_clone1->SetFillStyle(3002);
        h_marginal_T12_clone1->SetLineColor(kRed);
        h_marginal_T12_clone1->SetLineWidth(2);
    //    h_marginal_T12_clone1->GetXaxis()->SetRangeUser(-0.1, 0.864839);
        h_marginal_T12_clone1->GetXaxis()->SetRange(1, 146);
        h_marginal_T12_clone1->Draw("][histsame");
        TH1D *h_marginal_T12_clone2 = (TH1D*)h_marginal_T12->Clone();
        h_marginal_T12_clone2->SetFillColor(kOrange + 8);
        h_marginal_T12_clone2->SetFillStyle(3002);
        h_marginal_T12_clone2->SetLineColor(kOrange + 8);
        h_marginal_T12_clone2->SetLineWidth(2);
    //    h_marginal_T12_clone2->GetXaxis()->SetRangeUser(0.864839, 1.10726);
        h_marginal_T12_clone2->GetXaxis()->SetRange(146, 182);
        h_marginal_T12_clone2->Draw("][histsame");
        TH1D *h_marginal_T12_clone3 = (TH1D*)h_marginal_T12->Clone();
        h_marginal_T12_clone3->SetFillColor(kBlue - 4);
        h_marginal_T12_clone3->SetFillStyle(3002);
        h_marginal_T12_clone3->SetLineColor(kBlue - 4);
        h_marginal_T12_clone3->SetLineWidth(2);
    //    h_marginal_T12_clone3->GetXaxis()->SetRangeUser(1.10726, 1.2279);
        h_marginal_T12_clone3->GetXaxis()->SetRange(182, 200);
        h_marginal_T12_clone3->Draw("][histsame");
        TH1D *h_marginal_T12_clone4 = (TH1D*)h_marginal_T12->Clone();
        h_marginal_T12_clone4->SetFillColor(kViolet - 3);
        h_marginal_T12_clone4->SetFillStyle(3002);
        h_marginal_T12_clone4->SetLineColor(kViolet - 3);
        h_marginal_T12_clone4->SetLineWidth(2);
    //    h_marginal_T12_clone4->GetXaxis()->SetRangeUser(1.2279, 1.47319);
        h_marginal_T12_clone4->GetXaxis()->SetRange(200, 237);
        h_marginal_T12_clone4->Draw("][histsame");
        TH1D *h_marginal_T12_clone5 = (TH1D*)h_marginal_T12->Clone();
        h_marginal_T12_clone5->SetFillColor(kBlack);
        h_marginal_T12_clone5->SetFillStyle(3002);
        h_marginal_T12_clone5->SetLineColor(kBlack);
        h_marginal_T12_clone5->SetLineWidth(2);
    //    h_marginal_T12_clone5->GetXaxis()->SetRangeUser(1.47319, 1.9);
        h_marginal_T12_clone5->GetXaxis()->SetRange(237, 301);
        h_marginal_T12_clone5->Draw("][histsame");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(63);
        latex.SetTextSize(15);
        std::cout << kRed + 1 << std::endl;
        std::cout << kRed << std::endl;
        std::cout << kRed - 4 << std::endl;
        std::cout << kRed - 7 << std::endl;

        TString lstring0;
        lstring0.Form("#int_{-0.1}^{1.9} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 1");
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.17, 0.78, lstring0);
        
        TString lstring1;
        lstring1.Form("#int_{-0.1}^{0.86} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 1 #sigma (68.3 %%)");
        latex.SetTextColor(kRed);
        latex.DrawLatex(0.55, 0.78, lstring1);
        
        TString lstring2;
        lstring2.Form("#int_{-0.1}^{1.11} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 90 %%");
        latex.SetTextColor(kOrange + 8);
        latex.DrawLatex(0.61, 0.65, lstring2);
        
        TString lstring3;
        lstring3.Form("#int_{-0.1}^{1.23} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 95 %%");
        latex.SetTextColor(kBlue - 4);
        latex.DrawLatex(0.63, 0.50, lstring3);
        
        TString lstring4;
        lstring4.Form("#int_{-0.1}^{1.47} P(#xi^{2#nu#beta#beta}_{31}) d#xi^{2#nu#beta#beta}_{31} = 99 %%");
        latex.SetTextColor(kViolet - 3);
        latex.DrawLatex(0.65, 0.35, lstring4);

        c_marginal_T12->SaveAs("c_marginal_T12.eps");
        c_marginal_T12->SaveAs("c_marginal_T12.pdf");
        c_marginal_T12->SaveAs("c_marginal_T12.png");
        */

        ///////////////////////////////////////////////////////////////////////
        // 1 sided limit
        ///////////////////////////////////////////////////////////////////////

        TH1D *h_cumulative_T12 = new TH1D("h_cumulative_T12", "h_cumulative_T12",
                                      h_marginal_T12->GetNbinsX(),
                                      h_marginal_T12->GetXaxis()->GetBinLowEdge(1),
                                      h_marginal_T12->GetXaxis()->GetBinUpEdge(h_marginal_T12->GetNbinsX()));
        h_cumulative_T12->SetTitle(0);
        h_cumulative_T12->SetStats(0);

        sum = 0.0;
        for(Int_t i = 1; i <= h_marginal_T12->GetNbinsX(); ++ i)
        {
            sum += h_marginal_T12->GetBinContent(i) * bin_width_y;
            h_cumulative_T12->SetBinContent(i, sum);
        }

        /*
        TCanvas *c_cumulative_T12 = new TCanvas("c_cumulative_T12", "c_cumulative_T12");
        h_cumulative_T12->Draw("hist");
        */

        TGraph *g_cumulative_T12_reverse = new TGraph(h_cumulative_T12->GetNbinsX());
        g_cumulative_T12_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_T12_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_T12->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_T12->GetBinContent(i + 1);
            g_cumulative_T12_reverse->SetPoint(i, y, x);
        }

        /*
        TCanvas *c_cumulative_T12_reverse = new TCanvas("c_cumulative_T12_reverse", "c_cumulative_T12_reverse");
        g_cumulative_T12_reverse->Draw("al");
        */

        std::cout << "1 sigma: " << g_cumulative_T12_reverse->Eval(0.683) << std::endl;
        std::cout << "90 %: " << g_cumulative_T12_reverse->Eval(0.90) << std::endl;
        std::cout << "95 %: " << g_cumulative_T12_reverse->Eval(0.95) << std::endl;
        std::cout << "99 %: " << g_cumulative_T12_reverse->Eval(0.99) << std::endl;

        std::cout << "*** bin index ***" << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_reverse->Eval(0.683)) << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_reverse->Eval(0.90)) << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_reverse->Eval(0.95)) << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_reverse->Eval(0.99)) << std::endl;

        double sum0 = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sum3 = 0.0;
        double sum4 = 0.0;
        for(Int_t i = 1; i <= h_marginal_T12->GetNbinsX(); ++ i)
        {
            const double xi_31 = h_marginal_T12->GetXaxis()->GetBinCenter(i);
            const double content = h_marginal_T12->GetBinContent(i);

            sum0 += content;
            if(xi_31 < g_cumulative_T12_reverse->Eval(0.683))
            {
                sum1 += content;
            }
            if(xi_31 < g_cumulative_T12_reverse->Eval(0.90))
            {
                sum2 += content;
            }
            if(xi_31 < g_cumulative_T12_reverse->Eval(0.95))
            {
                sum3 += content;
            }
            if(xi_31 < g_cumulative_T12_reverse->Eval(0.99))
            {
                sum4 += content;
            }
        }
        std::cout << "sum0=" << sum0 << " sum1=" << sum1 << " sum2=" << sum2 << " sum3=" << sum3 << " sum4=" << sum4 << std::endl;




    }


}

#endif // NEWLOGLIKFITTER_DRAW_GA_METHOD_PROBABILITY_MARGINALIZE_2_H
