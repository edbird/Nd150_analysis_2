#ifndef NEWLOGLIKFITTER_DRAW_GA_H
#define NEWLOGLIKFITTER_DRAW_GA_H


bool is_even(const int var)
{
    if(var % 2 == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


double mean_x(const TH1D* const hist)
{
    double sum = 0.0;
    double psum = 0.0;
    for(Int_t i = 1; i <= hist->GetNbinsX(); ++ i)
    {
        // bin width = dx
        double width = hist->GetXaxis()->GetBinWidth(1);

        // probability p(x)
        double p = hist->GetBinContent(i);
        psum += p * width;

        // x
        double x = hist->GetXaxis()->GetBinCenter(i);

        sum += x * p * width;
        //sum += (x - mean) * p * width;
    }
//    std::cout << "psum=" << psum << std::endl;
    return sum;
}

double variance_x_fixed_mean(const TH1D* const hist, const double mean)
{
    double sum = 0.0;
    double psum = 0.0;
    for(Int_t i = 1; i <= hist->GetNbinsX(); ++ i)
    {
        // bin width = dx
        double width = hist->GetXaxis()->GetBinWidth(1);

        // probability p(x)
        double p = hist->GetBinContent(i);
        psum += p * width;

        // x
        double x = hist->GetXaxis()->GetBinCenter(i);

        // method 1
        sum += std::pow(x - mean, 2.0) * p * width;

        // method 2
        //sum += x * x * p * width;

        //sum += (x - mean) * p * width;
    }
//    std::cout << "psum=" << psum << std::endl;

    // method 2
    //sum -= mean * mean;

    return sum;
}


void draw_gA_method_probability_marginalize()
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
    std::cout << "Scale Factor: " << 1.0 / sum << std::endl;
    const double bin_area = (2.0 / 301.0) * ((1.45 - 0.95) / 301.0);
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
    h_prob->SetYTitle("^{150}Nd Amplitude Scale Factor (A_{0} = 3.45#times10^{-4} Bq)");
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

    {
        const double bin_width_x = (1.9 + 0.1) / 301.0;
        const double bin_width_y = (1.45 - 0.95) / 301.0;

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

        TCanvas *c_marginal_xi_31 = new TCanvas("c_marginal_xi_31", "c_marginal_xi_31");
        c_marginal_xi_31->SetTicks(1, 2);
        c_marginal_xi_31->SetBottomMargin(0.15);
        c_marginal_xi_31->SetRightMargin(0.1);
        c_marginal_xi_31->SetLeftMargin(0.12);
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
//        h_marginal_xi_31_clone->Draw("axis");
        h_marginal_xi_31_clone->Draw("hist");

    /*
        TH1D *h_marginal_xi_31_clone1 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone1->SetFillColor(kRed);
        h_marginal_xi_31_clone1->SetFillStyle(3002);
        h_marginal_xi_31_clone1->SetLineColor(kRed);
        h_marginal_xi_31_clone1->SetLineWidth(2);
    //    h_marginal_xi_31_clone1->GetXaxis()->SetRangeUser(-0.1, 0.864839);
        h_marginal_xi_31_clone1->GetXaxis()->SetRange(1, 146);
        h_marginal_xi_31_clone1->Draw("][histsame");
        TH1D *h_marginal_xi_31_clone2 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone2->SetFillColor(kOrange + 8);
        h_marginal_xi_31_clone2->SetFillStyle(3002);
        h_marginal_xi_31_clone2->SetLineColor(kOrange + 8);
        h_marginal_xi_31_clone2->SetLineWidth(2);
    //    h_marginal_xi_31_clone2->GetXaxis()->SetRangeUser(0.864839, 1.10726);
        h_marginal_xi_31_clone2->GetXaxis()->SetRange(146, 182);
        h_marginal_xi_31_clone2->Draw("][histsame");
        TH1D *h_marginal_xi_31_clone3 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone3->SetFillColor(kBlue - 4);
        h_marginal_xi_31_clone3->SetFillStyle(3002);
        h_marginal_xi_31_clone3->SetLineColor(kBlue - 4);
        h_marginal_xi_31_clone3->SetLineWidth(2);
    //    h_marginal_xi_31_clone3->GetXaxis()->SetRangeUser(1.10726, 1.2279);
        h_marginal_xi_31_clone3->GetXaxis()->SetRange(182, 200);
        h_marginal_xi_31_clone3->Draw("][histsame");
        TH1D *h_marginal_xi_31_clone4 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone4->SetFillColor(kViolet - 3);
        h_marginal_xi_31_clone4->SetFillStyle(3002);
        h_marginal_xi_31_clone4->SetLineColor(kViolet - 3);
        h_marginal_xi_31_clone4->SetLineWidth(2);
    //    h_marginal_xi_31_clone4->GetXaxis()->SetRangeUser(1.2279, 1.47319);
        h_marginal_xi_31_clone4->GetXaxis()->SetRange(200, 237);
        h_marginal_xi_31_clone4->Draw("][histsame");
        TH1D *h_marginal_xi_31_clone5 = (TH1D*)h_marginal_xi_31->Clone();
        h_marginal_xi_31_clone5->SetFillColor(kBlack);
        h_marginal_xi_31_clone5->SetFillStyle(3002);
        h_marginal_xi_31_clone5->SetLineColor(kBlack);
        h_marginal_xi_31_clone5->SetLineWidth(2);
    //    h_marginal_xi_31_clone5->GetXaxis()->SetRangeUser(1.47319, 1.9);
        h_marginal_xi_31_clone5->GetXaxis()->SetRange(237, 301);
        h_marginal_xi_31_clone5->Draw("][histsame");
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
    */
        c_marginal_xi_31->SaveAs("c_marginal_xi_31.eps");
        c_marginal_xi_31->SaveAs("c_marginal_xi_31.pdf");
        c_marginal_xi_31->SaveAs("c_marginal_xi_31.png");


        TH1D *h_cumulative_xi_31 = new TH1D("h_cumulative_xi_31", "h_cumulative_xi_31",
                                      h_marginal_xi_31->GetNbinsX(),
                                      h_marginal_xi_31->GetXaxis()->GetBinLowEdge(1),
                                      h_marginal_xi_31->GetXaxis()->GetBinUpEdge(h_marginal_xi_31->GetNbinsX()));
        h_cumulative_xi_31->SetTitle(0);
        h_cumulative_xi_31->SetStats(0);

        sum = 0.0;
        for(Int_t i = 1; i <= h_marginal_xi_31->GetNbinsX(); ++ i)
        {
            sum += h_marginal_xi_31->GetBinContent(i) * bin_width_x;
            h_cumulative_xi_31->SetBinContent(i, sum);
        }

        TCanvas *c_cumulative_xi_31 = new TCanvas("c_cumulative_xi_31", "c_cumulative_xi_31");
        h_cumulative_xi_31->Draw("hist");

        TGraph *g_cumulative_xi_31_reverse = new TGraph(h_cumulative_xi_31->GetNbinsX());
        g_cumulative_xi_31_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_xi_31_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_xi_31->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_xi_31->GetBinContent(i + 1);
            g_cumulative_xi_31_reverse->SetPoint(i, y, x);
        }

        TCanvas *c_cumulative_xi_31_reverse = new TCanvas("c_cumulative_xi_31_reverse", "c_cumulative_xi_31_reverse");
        g_cumulative_xi_31_reverse->Draw("al");

        std::cout << "1 sigma: " << g_cumulative_xi_31_reverse->Eval(0.683) << std::endl;
        std::cout << "90 %: " << g_cumulative_xi_31_reverse->Eval(0.90) << std::endl;
        std::cout << "95 %: " << g_cumulative_xi_31_reverse->Eval(0.95) << std::endl;
        std::cout << "99 %: " << g_cumulative_xi_31_reverse->Eval(0.99) << std::endl;

        std::cout << "*** bin index ***" << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_reverse->Eval(0.683)) << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_reverse->Eval(0.90)) << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_reverse->Eval(0.95)) << std::endl;
        std::cout << h_marginal_xi_31->GetXaxis()->FindBin(g_cumulative_xi_31_reverse->Eval(0.99)) << std::endl;

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
    }

    ///////////////////////////////////////////////////////////////////////////
    // marginal T12
    ///////////////////////////////////////////////////////////////////////////

    TH1D *h_marginal_T12 = new TH1D("h_marginal_T12", "h_marginal_T12",
                                h_prob->GetNbinsY(),
                                h_prob->GetYaxis()->GetBinLowEdge(1),
                                h_prob->GetYaxis()->GetBinUpEdge(h_prob->GetNbinsY()));

    {
        const double bin_width_x = (1.9 + 0.1) / 301.0;
        const double bin_width_y = (1.45 - 0.95) / 301.0;

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
        h_marginal_T12->SetXTitle("^{150}Nd Amplitude Scale Factor (A_{0} = 3.45#times10^{-4} Bq)");
        h_marginal_T12->SetYTitle("Probability Density Function p(A_{^{150}Nd})");
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
//        h_marginal_T12_clone->SetMaximum(1.7);
//        h_marginal_T12_clone->SetMinimum(0.0);
//        h_marginal_T12_clone->Draw("axis");
        h_marginal_T12_clone->Draw("hist");

        {
            double sum = 0.0;
            const double width = h_marginal_T12_clone->GetXaxis()->GetBinWidth(1);
            int low = h_marginal_T12_clone->GetXaxis()->FindBin(1.19003 - 0.0694279);
            int high = h_marginal_T12_clone->GetXaxis()->FindBin(1.19003 + 0.0694279);
            for(Int_t i = low; i <= high; ++ i)
            {
                if((i == low) || (i == high))
                {
                    sum += 0.5 * h_marginal_T12_clone->GetBinContent(i) * width;
                }
                else
                {
                    sum += 0.5 * h_marginal_T12_clone->GetBinContent(i) * width;
                }
            }
            std::cout << "Integration between 1.19003 +- 0.0694279: " << sum << std::endl;
        }

        std::cout << "The estimate of sigma from the height is: " << (4.0 * std::atan(1.0)) / std::pow(peak_y, 2.0) << std::endl;

        TH1D *h_marginal_T12_estimated = new TH1D("h_marginal_T12_estimated", "h_marginal_T12_estimated",
                                                  h_marginal_T12->GetNbinsX(),
                                                  h_marginal_T12->GetXaxis()->GetBinLowEdge(1),
                                                  h_marginal_T12->GetXaxis()->GetBinUpEdge(h_marginal_T12->GetNbinsX()));
        for(Int_t i = 1; i <= h_marginal_T12_estimated->GetNbinsX(); ++ i)
        {
            double x = h_marginal_T12_estimated->GetXaxis()->GetBinCenter(i);
//            double content = peak_y * std::exp(-std::pow((x - 1.19003) / (0.240864 * 0.4), 2.0));
//            double content = peak_y * std::exp(-std::pow((x - 1.19003) / 0.0580154, 2.0));
//            double content = peak_y * std::exp(-std::pow((x - 1.19003) / 0.0694279, 2.0));
            //double content = peak_y * std::exp(-std::pow((x - 1.19003) / (0.112576 * 0.9), 2.0));
            //double content = peak_y * std::exp(-std::pow((x - 1.19333) / (0.112576 * 0.88), 2.0));
            double content = peak_y * std::exp(-0.5 * std::pow((x - 1.19333) / (0.0694279), 2.0));

            h_marginal_T12_estimated->SetBinContent(i, content);
        }
        h_marginal_T12_estimated->SetLineColor(kRed);
        h_marginal_T12_estimated->Draw("histsame");

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

        TCanvas *c_cumulative_T12 = new TCanvas("c_cumulative_T12", "c_cumulative_T12");
        h_cumulative_T12->Draw("hist");

        TGraph *g_cumulative_T12_reverse = new TGraph(h_cumulative_T12->GetNbinsX());
        g_cumulative_T12_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_T12_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_T12->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_T12->GetBinContent(i + 1);
            g_cumulative_T12_reverse->SetPoint(i, y, x);
        }

        TCanvas *c_cumulative_T12_reverse = new TCanvas("c_cumulative_T12_reverse", "c_cumulative_T12_reverse");
        g_cumulative_T12_reverse->Draw("al");

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

        ///////////////////////////////////////////////////////////////////////
        // 2 sided limit
        ///////////////////////////////////////////////////////////////////////

        TH1D *h_cumulative_T12_low = nullptr;
        TH1D *h_cumulative_T12_high = nullptr;

        double x_low, x_mid, x_high;
        double x_range;
        Int_t n_bins;

        if(is_even(h_marginal_T12->GetNbinsX()))
        {
            std::cout << "EVEN" << std::endl;

            x_low = h_marginal_T12->GetXaxis()->GetBinLowEdge(1);
            x_mid = h_marginal_T12->GetXaxis()->GetBinUpEdge(1 + (h_marginal_T12->GetNbinsX() / 2));
            x_high = h_marginal_T12->GetXaxis()->GetBinUpEdge(h_marginal_T12->GetNbinsX());
            n_bins = h_marginal_T12->GetNbinsX() / 2;

            std::cout << "n_bins=" << n_bins << " x_low=" << x_low << " x_mid=" << x_mid << " x_high=" << x_high << std::endl;

            h_cumulative_T12_low = 
                new TH1D("h_cumulative_T12_low", "h_cumulative_T12_low",
                         n_bins, x_low, x_mid);

            h_cumulative_T12_high = 
                new TH1D("h_cumulative_T12_high", "h_cumulative_T12_high",
                         n_bins, x_mid, x_high);
        }
        else
        {
            std::cout << "ODD" << std::endl;

            x_low = h_marginal_T12->GetXaxis()->GetBinLowEdge(1);
            //x_mid = h_marginal_T12->GetXaxis()->GetBinCenter((h_marginal_T12->GetNbinsX() + 1) / 2);
            Int_t index_mid = peak_j;
            //Int_t index_mid = max_bin_y;
            std::cout << ">>>>>>>> " << "peak_j=" << peak_j << " max_bin_y=" << max_bin_y << std::endl;
            x_mid = h_marginal_T12->GetXaxis()->GetBinCenter(index_mid);
            x_high = h_marginal_T12->GetXaxis()->GetBinUpEdge(h_marginal_T12->GetNbinsX());
            n_bins = (h_marginal_T12->GetNbinsX() + 1) / 2;
            x_range = x_high - x_low;

            std::cout << "n_bins=" << n_bins << " x_low=" << x_low << " x_mid=" << x_mid << " x_high=" << x_high << std::endl;

            h_cumulative_T12_low = 
                new TH1D("h_cumulative_T12_low", "h_cumulative_T12_low",
                         h_marginal_T12->GetNbinsX(), x_mid - x_range, x_mid);
            h_cumulative_T12_low->SetTitle(0);
            h_cumulative_T12_low->SetStats(0);

            h_cumulative_T12_high = 
                new TH1D("h_cumulative_T12_high", "h_cumulative_T12_high",
                         h_marginal_T12->GetNbinsX(), x_mid, x_mid + x_range);
            h_cumulative_T12_high->SetTitle(0);
            h_cumulative_T12_high->SetStats(0);


            double sum_low = 0.0;
            double sum_high = 0.0;
            //Int_t index_mid = (h_marginal_T12->GetNbinsX() + 1) / 2;
            std::cout << "index_mid=" << index_mid << std::endl;
            sum_low += 0.5 * h_marginal_T12->GetBinContent(index_mid) * bin_width_y;
            sum_high += 0.5 * h_marginal_T12->GetBinContent(index_mid) * bin_width_y;
            for(Int_t i = 1; i <= h_marginal_T12->GetNbinsX(); ++ i)
            {
                Int_t index_low = index_mid - i;
                Int_t index_high = index_mid + i;

                //std::cout << "index: " << index_low << " " << index_high << std::endl;
                if(index_low >= 1)
                {
                    sum_low += h_marginal_T12->GetBinContent(index_low) * bin_width_y;
                }
                if(i <= h_cumulative_T12_low->GetNbinsX())
                {   
                    h_cumulative_T12_low->SetBinContent(h_cumulative_T12_low->GetNbinsX() - (i - 1), sum_low);
                }

                if(index_high <= h_marginal_T12->GetNbinsX())
                {
                    sum_high += h_marginal_T12->GetBinContent(index_high) * bin_width_y;
                }
                if(i <= h_cumulative_T12_high->GetNbinsX())
                {
                    h_cumulative_T12_high->SetBinContent(i, sum_high);
                }
            }
        }

        TCanvas *c_cumulative_T12_low = new TCanvas("c_cumulative_T12_low", "c_cumulative_T12_low");
        h_cumulative_T12_low->Draw("hist");

        TCanvas *c_cumulative_T12_high = new TCanvas("c_cumulative_T12_high", "c_cumulative_T12_high");
        h_cumulative_T12_high->Draw("hist");

        // cumulative reverse low

        TGraph *g_cumulative_T12_low_reverse = new TGraph(h_cumulative_T12->GetNbinsX());
        g_cumulative_T12_low_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_T12_low_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_T12_low->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_T12_low->GetBinContent(i + 1);
            g_cumulative_T12_low_reverse->SetPoint(i, y, x);
        }

        //TCanvas *c_cumulative_T12_low_reverse = new TCanvas("c_cumulative_T12_low_reverse", "c_cumulative_T12_low_reverse");
        TCanvas *c_cumulative_T12_highlow_reverse = new TCanvas("c_cumulative_T12_highlow_reverse", "c_cumulative_T12_highlow_reverse");
        g_cumulative_T12_low_reverse->Draw("al");

        std::cout << "1 sigma: " << g_cumulative_T12_low_reverse->Eval(0.5 * 0.683) << std::endl;
        std::cout << "90 %: " << g_cumulative_T12_low_reverse->Eval(0.5 * 0.90) << std::endl;
        std::cout << "95 %: " << g_cumulative_T12_low_reverse->Eval(0.5 * 0.95) << std::endl;
        std::cout << "99 %: " << g_cumulative_T12_low_reverse->Eval(0.5 * 0.99) << std::endl;

        std::cout << "*** bin index ***" << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_low_reverse->Eval(0.5 * 0.683)) << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_low_reverse->Eval(0.5 * 0.90)) << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_low_reverse->Eval(0.5 * 0.95)) << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_low_reverse->Eval(0.5 * 0.99)) << std::endl;


        // cumulative reverse high

        TGraph *g_cumulative_T12_high_reverse = new TGraph(h_cumulative_T12->GetNbinsX());
        g_cumulative_T12_high_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_T12_high_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_T12_high->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_T12_high->GetBinContent(i + 1);
            g_cumulative_T12_high_reverse->SetPoint(i, y, x);
        }

        //TCanvas *c_cumulative_T12_high_reverse = new TCanvas("c_cumulative_T12_high_reverse", "c_cumulative_T12_high_reverse");
        g_cumulative_T12_high_reverse->Draw("lsame");

        std::cout << "1 sigma: " << g_cumulative_T12_high_reverse->Eval(0.5 * 0.683) << std::endl;
        std::cout << "90 %: " << g_cumulative_T12_high_reverse->Eval(0.5 * 0.90) << std::endl;
        std::cout << "95 %: " << g_cumulative_T12_high_reverse->Eval(0.5 * 0.95) << std::endl;
        std::cout << "99 %: " << g_cumulative_T12_high_reverse->Eval(0.5 * 0.99) << std::endl;

        std::cout << "*** bin index ***" << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_high_reverse->Eval(0.5 * 0.683)) << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_high_reverse->Eval(0.5 * 0.90)) << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_high_reverse->Eval(0.5 * 0.95)) << std::endl;
        std::cout << h_marginal_T12->GetXaxis()->FindBin(g_cumulative_T12_high_reverse->Eval(0.5 * 0.99)) << std::endl;

        std::cout << "1 sigma: T12 = " << peak_x << " " 
                  << "+ " << (g_cumulative_T12_high_reverse->Eval(0.5 * 0.683) - peak_x) << " "
                  << "- " << (peak_x - g_cumulative_T12_low_reverse->Eval(0.5 * 0.683)) << std::endl;
        std::cout << "90 %: T12 = " << peak_x << " " 
                  << "+ " << (g_cumulative_T12_high_reverse->Eval(0.5 * 0.90) - peak_x) << " "
                  << "- " << (peak_x - g_cumulative_T12_low_reverse->Eval(0.5 * 0.90)) << std::endl;
        std::cout << "95 %: T12 = " << peak_x << " " 
                  << "+ " << (g_cumulative_T12_high_reverse->Eval(0.5 * 0.95) - peak_x) << " "
                  << "- " << (peak_x - g_cumulative_T12_low_reverse->Eval(0.5 * 0.95)) << std::endl;
        std::cout << "99 %: T12 = " << peak_x << " " 
                  << "+ " << (g_cumulative_T12_high_reverse->Eval(0.5 * 0.99) - peak_x) << " "
                  << "- " << (peak_x - g_cumulative_T12_low_reverse->Eval(0.5 * 0.99)) << std::endl;


        ///////////////////////////////////////////////////////////////////////
        // compute mean and variance of marginal T12
        ///////////////////////////////////////////////////////////////////////

        const double marginal_T12_mean = peak_x;
        std::cout << "Mean: " << marginal_T12_mean << std::endl;
        const double marginal_T12_var = variance_x_fixed_mean(h_marginal_T12, marginal_T12_mean);
        std::cout << "Var: " << marginal_T12_var << std::endl;
        std::cout << "Standard Error: " << std::sqrt(marginal_T12_var) << std::endl;
        std::cout << "Alternative mean: " << mean_x(h_marginal_T12) << std::endl;

    }


}



void draw_gA_method_contour()
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



    ///////////////////////////////////////////////////////////////////////////
    // completely new code to create gA - MGT-3 phase space using contours
    ///////////////////////////////////////////////////////////////////////////

    double gA4_MGT2_min[3] = {0.0, 0.0, 0.0};
    double gA4_MGT2_max[3] = {0.0, 0.0, 0.0};
    {

        // workaround for root bug

        TH2D *h_mps_clone = (TH2D*)h_mps->Clone("tmpclone");
        double min_fval_SYSALL = 0.0;
        if(g_mode_fake_data == false)
        {
            min_fval_SYSALL = min_point_data_SYSALL_fval;
        }
        else
        {
            min_fval_SYSALL = min_point_fake_SYSALL_fval;
        }
        double clevels[3] = {min_fval_SYSALL + 2.30, min_fval_SYSALL + 4.61, min_fval_SYSALL + 9.21};
        h_mps_clone->SetContour(3, clevels);
        TCanvas *ctmp = new TCanvas("ctmp", "ctmp");
        h_mps_clone->Draw("CONTLIST");
        ctmp->Update();


        // get contours

	    TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
        TList *contLevel = nullptr;
        TGraph *curv = nullptr;
        TGraph *gc = nullptr;
        Int_t nGraphs = 0;
        Int_t TotalConts = 0;
        double xval0, yval0;
        TGraph *graph[3] = {nullptr, nullptr, nullptr};

        if(conts != nullptr)
        {
            TotalConts = conts->GetSize();
            std::cout << "size: " << TotalConts << std::endl;
        }
        else
        {
            std::cout << "get contours failed" << std::endl;
        }

        
        // get graphs for each contour, seems to be 1 graph object per contour

        for(Int_t i = 0; i < TotalConts; ++ i)
        {
            contLevel = (TList*)conts->At(i);
            std::cout << "Contour " << i << " has " << contLevel->GetSize() << " graphs" << std::endl;
            nGraphs += contLevel->GetSize();
        }
        nGraphs = 0;

        
        // get graph for each contour and iterate over points

        for(Int_t i = 0; i < TotalConts; ++ i)
        {
            contLevel = (TList*)conts->At(i);
            std::cout << "chi2=" << clevels[i] << std::endl;

            curv = (TGraph*)contLevel->First();
            for(Int_t j = 0; j < contLevel->GetSize(); ++ j) // size is 1
            {
                std::cout << "Graph " << nGraphs << " has " << curv->GetN() << " points" << std::endl;
                ++ nGraphs;

                graph[i] = new TGraph(curv->GetN());
                
                Int_t kk = 0;
                for(Int_t k = 0; k < curv->GetN(); ++ k)
                {
                    curv->GetPoint(k, xval0, yval0);

                    double A_150Nd = yval0;
                    double xi_31 = xval0;

                    const double A0_150Nd = 3.45e-04;
                    const double NA = 6.022e+23;
                    const double m = 36.6;
                    const double sec_to_yr = 1.0 / 31557600.0;
                    double T12 = (std::log(2.0) * NA * m) / (150.0 * A0_150Nd * A_150Nd) * sec_to_yr;
                    double gA4_MGT2 = std::pow(xi_31, 2.0) / (T12 * (psiN0 + xi_31 * psiN2));

                    /*
                    double MGT3 = h_gA_MGT3->GetXaxis()->GetBinCenter(i);
                    double gA = std::pow(tmp_value / std::pow(MGT3, 2.0), 1.0 / 4.0);
                    */

                    /*
                    for(double matrix_element = 0.01; matrix_element <= 0.2; matrix_element += 0.01)
                    {
                        double gA = std::pow(tmp_value / std::pow(matrix_element, 2.0), 1.0 / 4.0);
                        if(kk < curv->GetN())
                        {
                            graph[i]->SetPoint(kk, matrix_element, gA);
                            ++ kk;
                        }
                    }
                    */

                    if(k == 0)
                    {
                        gA4_MGT2_min[i] = gA4_MGT2;
                        gA4_MGT2_max[i] = gA4_MGT2;
                    }
                    else
                    {
                        if(gA4_MGT2 < gA4_MGT2_min[i]) gA4_MGT2_min[i] = gA4_MGT2;
                        if(gA4_MGT2 > gA4_MGT2_max[i]) gA4_MGT2_max[i] = gA4_MGT2;
                    }
                }

                std::cout << "i=" << i << " gA4_MGT2_min=" << gA4_MGT2_min[i] << std::endl;
                std::cout << "i=" << i << " gA4_MGT2_max=" << gA4_MGT2_max[i] << std::endl;
            }

        }


        double A_150Nd = min_point_data_SYSALL[1];
        double xi_31 = min_point_data_SYSALL[0];

        const double A0_150Nd = 3.45e-04;
        const double NA = 6.022e+23;
        const double m = 36.6;
        const double sec_to_yr = 1.0 / 31557600.0;
        double T12 = (std::log(2.0) * NA * m) / (150.0 * A0_150Nd * A_150Nd) * sec_to_yr;
        double gA4_MGT2 = std::pow(xi_31, 2.0) / (T12 * (psiN0 + xi_31 * psiN2));

        std::cout << "minimum point: " << gA4_MGT2 << std::endl;

        /*
        for(int i = 0; i < 3; ++ i)
        {
            contLevel = (TList*)conts->At(i);
            curv = (TGraph*)contLevel->First();
            for(Int_t k = 0; k < curv->GetN(); ++ k)
            {
                double x, y;
                curv->GetPoint(k, x, y);
                if(k == 0)
                {
                    xi_31_min[i] = x;
                    xi_31_max[i] = x;
                }
                if(x < xi_31_min[i]) xi_31_min[i] = x;
                if(x > xi_31_max[i]) xi_31_max[i] = x;
            }
            //curv = (TGraph*)contLevel->After(curv);
        
            std::cout << "xi_31: " << xi_31_min[i] << " " << xi_31_max[i] << std::endl;
        }
        */
        /*
        contLevel = (TList*)conts->At(0);
        curv = (TGraph*)contLevel->First();
        */
    /*
        TCanvas *cc = new TCanvas("cc", "cc");
        //graph[0]->Draw();
        std::cout << "curv->GetN()=" << curv->GetN() << std::endl;
        curv->Draw("AL");
        cc->Show();

        std::ofstream oftmp("oftmp.txt");
        for(Int_t i = 0; i < curv->GetN(); ++ i)
        {
            double x, y;
            curv->GetPoint(i, x, y);
            oftmp << x << " " << y << std::endl;
        }
        oftmp.close();
    */



    }
}



void draw_gA_method_forward_mapping()
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
    // convert to some new phase space
    // forward mapping method
    ///////////////////////////////////////////////////////////////////////////

    double cmin = std::numeric_limits<double>::infinity();
    double cmax = 0.0;
    TH2D *h_new = new TH2D("h_new", "h_new", 50, 0.0, 0.20, 50, 0.0, 1.2);
    h_new->SetStats(0);
    h_new->SetTitle(0);
    h_new->SetXTitle("M_{GT-3}^{2#nu}");
    h_new->SetYTitle("g_{A}");
    h_new->GetXaxis()->SetTitleFont(43);
    h_new->GetXaxis()->SetTitleSize(20);
    h_new->GetYaxis()->SetTitleFont(43);
    h_new->GetYaxis()->SetTitleSize(20);
    h_new->GetXaxis()->SetLabelFont(43);
    h_new->GetXaxis()->SetLabelSize(15);
    h_new->GetYaxis()->SetLabelFont(43);
    h_new->GetYaxis()->SetLabelSize(15);
    h_new->GetXaxis()->CenterTitle(true);
    h_new->GetYaxis()->CenterTitle(true);

    for(Int_t j = 1; j <= h_mps->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_mps->GetNbinsX(); ++ i)
        {
            double xi_31 = h_mps->GetXaxis()->GetBinCenter(i);
            double A_150Nd = h_mps->GetYaxis()->GetBinCenter(j);
            double T12 = (std::log(2.0) * NA * m) / (150.0 * A0_150Nd * A_150Nd) * sec_to_yr;

            //std::cout << "xi_31=" << xi_31 << " A_150Nd=" << A_150Nd << " T12=" << T12 << std::endl;

            // this is gAeff ** 4.0 * MGT-3 ** 2.0
            double value = std::pow(xi_31, 2.0) / (T12 * (psiN0 + xi_31 * psiN2));

            for(Int_t ii = 1; ii <= h_new->GetNbinsX(); ++ ii)
            {
                double matrix_element = h_new->GetXaxis()->GetBinCenter(ii);
                double gAeff = std::pow(value / std::pow(matrix_element, 2.0), 0.25);
                //std::cout << "MGT=" << matrix_element << std::endl;
                //std::cout << "value=" << value << std::endl;
                //std::cout << "gAeff=" << gAeff << std::endl;

                Int_t jj = h_new->GetYaxis()->FindBin(gAeff);

                //std::cout << "jj=" << jj << std::endl;
                //std::cin.get();

                if((1 <= jj) && (jj <= 50))
                {
                    double chi2 = h_mps->GetBinContent(i, j);
                    double content = h_new->GetBinContent(ii, jj);
                    if(content == 0.0)
                    {
                        h_new->SetBinContent(ii, jj, chi2);
                        if(chi2 < cmin) cmin = chi2;
                        if(chi2 > cmax) cmax = chi2;
                    }
                    else if(chi2 < content)
                    {
                        h_new->SetBinContent(ii, jj, chi2);
                        if(chi2 < cmin) cmin = chi2;
                        if(chi2 > cmax) cmax = chi2;
                    }
                }
            }
        }
    }

//    h_new->SetMaximum(cmax);
    h_new->GetZaxis()->SetRangeUser(cmin, cmin + 10.0);
    h_new->GetZaxis()->SetLabelFont(43);
    h_new->GetZaxis()->SetLabelSize(15);
//    std::cout << cmax << std::endl;


    TCanvas *c_new = new TCanvas("c_new", "c_new");
    c_new->SetTicks(2, 2);
    c_new->SetRightMargin(0.15);
    c_new->SetBottomMargin(0.15);
    if(1)
    {
        h_new->Draw("colz");
        c_new->Update();

        TPaletteAxis *palette2 = (TPaletteAxis*)h_new->GetListOfFunctions()->FindObject("palette");
        palette2->SetX1NDC(0.88 + 0.03);
        palette2->SetX2NDC(0.92 + 0.03);
        palette2->SetY1NDC(0.15);
        palette2->SetY2NDC(0.9);
        palette2->Draw();
        gPad->Modified();
        gPad->Update();
        c_new->Modified();
    }

    TH2D *h_new_clone = (TH2D*)h_new->Clone();
    for(Int_t j = 1; j <= h_new_clone->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_new_clone->GetNbinsX(); ++ i)
        {
            double content = h_new_clone->GetBinContent(i, j);
            if(content == 0.0)
            {
                //h_new_clone->SetBinContent(i, j, std::numeric_limits<double>::infinity());
                h_new_clone->SetBinContent(i, j, 10.0 * cmax);
            }
        }
    }
    /*
    for(Int_t j = 1; j <= h_new_clone->GetNbinsY(); ++ j)
    {
        h_new_clone->SetBinContent(1, j, 10.0 * cmax);
        h_new_clone->SetBinContent(h_new_clone->GetNbinsX(), j, 10.0 * cmax);
    }
    for(Int_t i = 1; i <= h_new_clone->GetNbinsX(); ++ i)
    {
        h_new_clone->SetBinContent(i, 1, 10.0 * cmax);
        h_new_clone->SetBinContent(i, h_new_clone->GetNbinsY(), 10.0 * cmax);
    }
    */

//    TFile *fouttmp = new TFile("fouttmp.root", "RECREATE");
//    h_new_clone->Write();
//    fouttmp->Close();

//    h_new_clone->SaveAs("h_new_clone.C");
//    std::cin.get();

    const double usermax = 10.0;
//    std::cout << "SetRangeUser(" << cmin << ", " << usermax << ")" << std::endl;
//    h_new_clone->GetZaxis()->SetRangeUser(cmin, usermax);
    double clevels[2] = {cmin + 2.30, cmin + 4.61};
    //double clevels[2] = {0.0, 1.0};
    h_new_clone->SetContour(2, clevels);
    h_new_clone->SetLineColor(kBlack);
    h_new_clone->SetLineWidth(2.0);
    if(1)
    {
        h_new_clone->Draw("cont3same");
    }
    else
    {
//        h_new_clone->Draw("cont4list");
    }
//    h_new_clone->GetXaxis()->SetRangeUser(0.0, 0.2);
//    h_new_clone->GetYaxis()->SetRangeUser(0.0, 1.2);

/*
    c_new->Update();

    TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    std::cout << "conts->GetSize()=" << conts->GetSize() << std::endl;

    TList *contlevel = nullptr;
    TGraph *contcurve = nullptr;

    for(int i = 0; i < conts->GetSize(); ++ i)
    {
        contlevel = (TList*)conts->At(i);
        std::cout << "Contour " << i << " has " << contlevel->GetSize() << " graphs" << std::endl;

        for(int j = 0; j < contlevel->GetSize(); ++ j)
        {
            
        }
    }
*/

    c_new->SaveAs("gA_phase_space_method_forward_mapping.eps");
    c_new->SaveAs("gA_phase_space_method_forward_mapping.pdf");
    c_new->SaveAs("gA_phase_space_method_forward_mapping.png");

}


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



#endif // NEWLOGLIKFITTER_DRAW_GA_H
