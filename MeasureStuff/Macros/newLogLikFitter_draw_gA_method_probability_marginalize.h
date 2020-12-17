#ifndef NEWLOGLIKFITTER_DRAW_GA_METHOD_PROBABILITY_MARGINALIZE_H
#define NEWLOGLIKFITTER_DRAW_GA_METHOD_PROBABILITY_MARGINALIZE_H

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

    /*
        TCanvas *c_cumulative_xi_31 = new TCanvas("c_cumulative_xi_31", "c_cumulative_xi_31");
        h_cumulative_xi_31->Draw("hist");
    */
        TGraph *g_cumulative_xi_31_reverse = new TGraph(h_cumulative_xi_31->GetNbinsX());
        g_cumulative_xi_31_reverse->SetTitle(0);
        for(Int_t i = 0; i < g_cumulative_xi_31_reverse->GetN(); ++ i)
        {
            double x = h_cumulative_xi_31->GetXaxis()->GetBinCenter(i + 1);
            double y = h_cumulative_xi_31->GetBinContent(i + 1);
            g_cumulative_xi_31_reverse->SetPoint(i, y, x);
        }
    /*
        TCanvas *c_cumulative_xi_31_reverse = new TCanvas("c_cumulative_xi_31_reverse", "c_cumulative_xi_31_reverse");
    */
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
        h_marginal_T12->SetXTitle("^{150}Nd Amplitude Scale Factor (#times 3.45#times10^{-4} Bq)");
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

    /*
        TCanvas *c_cumulative_T12_low = new TCanvas("c_cumulative_T12_low", "c_cumulative_T12_low");
        h_cumulative_T12_low->Draw("hist");

        TCanvas *c_cumulative_T12_high = new TCanvas("c_cumulative_T12_high", "c_cumulative_T12_high");
        h_cumulative_T12_high->Draw("hist");
    */
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
        //TCanvas *c_cumulative_T12_highlow_reverse = new TCanvas("c_cumulative_T12_highlow_reverse", "c_cumulative_T12_highlow_reverse");
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

#endif // NEWLOGLIKFITTER_DRAW_GA_METHOD_PROBABILITY_MARGINALIZE_H
