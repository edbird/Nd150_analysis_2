#ifndef NEWLOGLIKFITTER_DRAW_GA_H
#define NEWLOGLIKFITTER_DRAW_GA_H



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
    // convert to probability space
    ///////////////////////////////////////////////////////////////////////////

    TH2D *h_prob = new TH2D("h_prob", "h_prob",
                            h_mps->GetNbinsX(),
                            h_mps->GetXaxis()->GetBinLowEdge(1),
                            h_mps->GetXaxis()->GetBinUpEdge(h_mps->GetNbinsX()),
                            h_mps->GetNbinsY(),
                            h_mps->GetYaxis()->GetBinLowEdge(1),
                            h_mps->GetYaxis()->GetBinUpEdge(h_mps->GetNbinsY()));

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
    h_prob->SetYTitle("^{150}Nd Amplitude Scale Factor (A_{0} = 3.45#times10^{-4})");
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

    // marginal
    TH1D *h_marginal = new TH1D("h_marginal", "h_marginal",
                                h_prob->GetNbinsX(),
                                h_prob->GetXaxis()->GetBinLowEdge(1),
                                h_prob->GetXaxis()->GetBinUpEdge(h_prob->GetNbinsX()));

    for(Int_t i = 1; i <= h_prob->GetNbinsX(); ++ i)
    {
        double sum = 0.0;
        for(Int_t j = 1; j <= h_prob->GetNbinsY(); ++ j)
        {
            sum += h_prob->GetBinContent(i, j);
        }
        h_marginal->SetBinContent(i, sum);
    }
    const double bin_width = 2.0 / 301.0;
    h_marginal->Scale(1.0 / (sum * bin_width));

    sum = 0.0;
    for(Int_t i = 1; i <= h_marginal->GetNbinsX(); ++ i)
    {
        sum += h_marginal->GetBinContent(i);
    }
    std::cout << "integral in 1D P-space: " << sum * bin_width << std::endl;

    TCanvas *c_marginal = new TCanvas("c_marginal", "c_marginal");
    c_marginal->SetTicks(1, 2);
    c_marginal->SetBottomMargin(0.15);
    c_marginal->SetRightMargin(0.1);
    c_marginal->SetLeftMargin(0.12);
    h_marginal->SetTitle("");
    h_marginal->SetStats(0);
    h_marginal->GetZaxis()->SetLabelOffset(0.005);
    h_marginal->GetXaxis()->SetLabelSize(17.0);
    h_marginal->GetXaxis()->SetLabelFont(43);
    h_marginal->GetYaxis()->SetLabelSize(17.0);
    h_marginal->GetYaxis()->SetLabelFont(43);
    h_marginal->GetZaxis()->SetLabelSize(17.0);
    h_marginal->GetZaxis()->SetLabelFont(43);
    h_marginal->GetXaxis()->SetTitleSize(18.0);
    h_marginal->GetXaxis()->SetTitleFont(43);
    h_marginal->GetYaxis()->SetTitleSize(18.0);
    h_marginal->GetYaxis()->SetTitleFont(43);
    h_marginal->GetXaxis()->SetTitleOffset(1.5);
    h_marginal->GetYaxis()->SetTitleOffset(1.5);
    h_marginal->GetXaxis()->SetLabelOffset(0.01);
    h_marginal->GetYaxis()->SetLabelOffset(0.01);
    h_marginal->SetXTitle("#xi^{2#nu#beta#beta}_{31}");
    h_marginal->SetYTitle("Maginalized Probability Density Function P(#xi^{2#nu#beta#beta}_{31})");
    std::cout << h_marginal->GetBinWidth(1) << std::endl;
    std::cout << "INT=" << h_marginal->Integral() * h_marginal->GetBinWidth(1) << std::endl;
    TH1D *h_marginal_clone = (TH1D*)h_marginal->Clone();
//    h_marginal_clone->Scale(1.0 / (2.0 / 301.0));
//    for(Int_t i = 0; i <= h_marginal->GetNbinsX() + 1; ++ i)
    for(Int_t i = 1; i <= h_marginal->GetNbinsX(); ++ i)
    {
        const double content = h_marginal_clone->GetBinContent(i);
//        h_marginal_clone->SetBinContent(i, 301.0 / 2.0 * content);
//        std::cout << "i=" << i << " " << h_marginal->GetBinContent(i) << std::endl;
    }
    h_marginal_clone->SetMaximum(1.7);
    h_marginal_clone->SetMinimum(0.0);
    h_marginal_clone->Draw("axis");
    TH1D *h_marginal_clone1 = (TH1D*)h_marginal->Clone();
    h_marginal_clone1->SetFillColor(kRed);
    h_marginal_clone1->SetFillStyle(3002);
    h_marginal_clone1->SetLineColor(kRed);
    h_marginal_clone1->SetLineWidth(2);
//    h_marginal_clone1->GetXaxis()->SetRangeUser(-0.1, 0.864839);
    h_marginal_clone1->GetXaxis()->SetRange(1, 146);
    h_marginal_clone1->Draw("][histsame");
    TH1D *h_marginal_clone2 = (TH1D*)h_marginal->Clone();
    h_marginal_clone2->SetFillColor(kOrange + 8);
    h_marginal_clone2->SetFillStyle(3002);
    h_marginal_clone2->SetLineColor(kOrange + 8);
    h_marginal_clone2->SetLineWidth(2);
//    h_marginal_clone2->GetXaxis()->SetRangeUser(0.864839, 1.10726);
    h_marginal_clone2->GetXaxis()->SetRange(146, 182);
    h_marginal_clone2->Draw("][histsame");
    TH1D *h_marginal_clone3 = (TH1D*)h_marginal->Clone();
    h_marginal_clone3->SetFillColor(kBlue - 4);
    h_marginal_clone3->SetFillStyle(3002);
    h_marginal_clone3->SetLineColor(kBlue - 4);
    h_marginal_clone3->SetLineWidth(2);
//    h_marginal_clone3->GetXaxis()->SetRangeUser(1.10726, 1.2279);
    h_marginal_clone3->GetXaxis()->SetRange(182, 200);
    h_marginal_clone3->Draw("][histsame");
    TH1D *h_marginal_clone4 = (TH1D*)h_marginal->Clone();
    h_marginal_clone4->SetFillColor(kViolet - 3);
    h_marginal_clone4->SetFillStyle(3002);
    h_marginal_clone4->SetLineColor(kViolet - 3);
    h_marginal_clone4->SetLineWidth(2);
//    h_marginal_clone4->GetXaxis()->SetRangeUser(1.2279, 1.47319);
    h_marginal_clone4->GetXaxis()->SetRange(200, 237);
    h_marginal_clone4->Draw("][histsame");
    TH1D *h_marginal_clone5 = (TH1D*)h_marginal->Clone();
    h_marginal_clone5->SetFillColor(kBlack);
    h_marginal_clone5->SetFillStyle(3002);
    h_marginal_clone5->SetLineColor(kBlack);
    h_marginal_clone5->SetLineWidth(2);
//    h_marginal_clone5->GetXaxis()->SetRangeUser(1.47319, 1.9);
    h_marginal_clone5->GetXaxis()->SetRange(237, 301);
    h_marginal_clone5->Draw("][histsame");
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

    c_marginal->SaveAs("c_marginal.eps");
    c_marginal->SaveAs("c_marginal.pdf");
    c_marginal->SaveAs("c_marginal.png");


    TH1D *h_cumulative = new TH1D("h_cumulative", "h_cumulative",
                                  h_marginal->GetNbinsX(),
                                  h_marginal->GetXaxis()->GetBinLowEdge(1),
                                  h_marginal->GetXaxis()->GetBinUpEdge(h_marginal->GetNbinsX()));

    sum = 0.0;
    for(Int_t i = 1; i <= h_marginal->GetNbinsX(); ++ i)
    {
        sum += h_marginal->GetBinContent(i) * bin_width;
        h_cumulative->SetBinContent(i, sum);
    }

    TCanvas *c_cumulative = new TCanvas("c_cumulative", "c_cumulative");
    h_cumulative->Draw("hist");

    TGraph *g_cumulative_reverse = new TGraph(h_cumulative->GetNbinsX());
    for(Int_t i = 0; i < g_cumulative_reverse->GetN(); ++ i)
    {
        double x = h_cumulative->GetXaxis()->GetBinCenter(i);
        double y = h_cumulative->GetBinContent(i);
        g_cumulative_reverse->SetPoint(i, y, x);
    }

    TCanvas *c_cumulative_reverse = new TCanvas("c_cumulative_reverse", "c_cumulative_reverse");
    g_cumulative_reverse->Draw("al");

    std::cout << "1 sigma: " << g_cumulative_reverse->Eval(0.683) << std::endl;
    std::cout << "90 %: " << g_cumulative_reverse->Eval(0.90) << std::endl;
    std::cout << "95 %: " << g_cumulative_reverse->Eval(0.95) << std::endl;
    std::cout << "99 %: " << g_cumulative_reverse->Eval(0.99) << std::endl;

    std::cout << "*** bin index ***" << std::endl;
    std::cout << h_marginal->GetXaxis()->FindBin(g_cumulative_reverse->Eval(0.683)) << std::endl;
    std::cout << h_marginal->GetXaxis()->FindBin(g_cumulative_reverse->Eval(0.90)) << std::endl;
    std::cout << h_marginal->GetXaxis()->FindBin(g_cumulative_reverse->Eval(0.95)) << std::endl;
    std::cout << h_marginal->GetXaxis()->FindBin(g_cumulative_reverse->Eval(0.99)) << std::endl;

    double sum0 = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;
    for(Int_t i = 1; i <= h_marginal->GetNbinsX(); ++ i)
    {
        const double xi_31 = h_marginal->GetXaxis()->GetBinCenter(i);
        const double content = h_marginal->GetBinContent(i);

        sum0 += content;
        if(xi_31 < g_cumulative_reverse->Eval(0.683))
        {
            sum1 += content;
        }
        if(xi_31 < g_cumulative_reverse->Eval(0.90))
        {
            sum2 += content;
        }
        if(xi_31 < g_cumulative_reverse->Eval(0.95))
        {
            sum3 += content;
        }
        if(xi_31 < g_cumulative_reverse->Eval(0.99))
        {
            sum4 += content;
        }
    }
    std::cout << "sum0=" << sum0 << " sum1=" << sum1 << " sum2=" << sum2 << " sum3=" << sum3 << " sum4=" << sum4 << std::endl;

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
