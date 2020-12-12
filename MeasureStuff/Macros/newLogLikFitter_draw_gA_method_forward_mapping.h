#ifndef NEWLOGLIKFITTER_DRAW_GA_METHOD_FORWARD_MAPPING_H
#define NEWLOGLIKFITTER_DRAW_GA_METHOD_FORWARD_MAPPING_H

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

#endif // NEWLOGLIKFITTER_DRAW_GA_METHOD_FORWARD_MAPPING_H
