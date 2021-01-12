#ifndef NEWLOGLIKFITTER_DRAW_GA_METHOD_UNTITLED_H
#define NEWLOGLIKFITTER_DRAW_GA_METHOD_UNTITLED_H

// could not think of a good name


void draw_gA_method_untitled()
{

    std::vector<double> gamgt3_nd150_arg_x;
    std::vector<double> gamgt3_nd150_arg_y;
    std::ifstream fs_gamgt3_nd150_arg("gamgt3_nd150_arg.dat");
    while(fs_gamgt3_nd150_arg.good())
    {
        double x, y;
        fs_gamgt3_nd150_arg >> x >> y;
        gamgt3_nd150_arg_x.push_back(x);
        gamgt3_nd150_arg_y.push_back(y);
    }
    fs_gamgt3_nd150_arg.close();

    std::vector<double> gamgt3_nd150_ssd_x;
    std::vector<double> gamgt3_nd150_ssd_y;
    std::ifstream fs_gamgt3_nd150_ssd("gamgt3_nd150_ssd.dat");
    while(fs_gamgt3_nd150_ssd.good())
    {
        double x, y;
        fs_gamgt3_nd150_ssd >> x >> y;
        gamgt3_nd150_ssd_x.push_back(x);
        gamgt3_nd150_ssd_y.push_back(y);
    }
    fs_gamgt3_nd150_ssd.close();

    TGraph *g_gamgt3_nd150_arg = new TGraph(gamgt3_nd150_arg_x.size(), gamgt3_nd150_arg_x.data(), gamgt3_nd150_arg_y.data());
    TGraph *g_gamgt3_nd150_ssd = new TGraph(gamgt3_nd150_ssd_x.size(), gamgt3_nd150_ssd_x.data(), gamgt3_nd150_ssd_y.data());
    g_gamgt3_nd150_arg->SetLineColor(kBlue);
    g_gamgt3_nd150_ssd->SetLineColor(kGreen);
    g_gamgt3_nd150_arg->SetLineStyle(1);
    g_gamgt3_nd150_ssd->SetLineStyle(1);
    g_gamgt3_nd150_arg->SetLineWidth(2);
    g_gamgt3_nd150_ssd->SetLineWidth(2);

    // need to check values
    
    const double A_CV = 1.20166;

    // TODO: need to go around contours to find the CL
    const double xi_31_90CL_h = 1.28758;
    const double xi_31_1sigma_h = 1.0755;
    const double xi_31_CV = 0.740532;
    const double xi_31_1sigma_l = 0.540747;
    const double xi_31_90CL_l = 0.404129;

    const int NPOINTS = 1000;
    TGraph *g_xi_31_90CL_h = new TGraph(NPOINTS);
    TGraph *g_xi_31_1sigma_h = new TGraph(NPOINTS);
    TGraph *g_xi_31_CV = new TGraph(NPOINTS);
    TGraph *g_xi_31_1sigma_l = new TGraph(NPOINTS);
    TGraph *g_xi_31_90CL_l = new TGraph(NPOINTS);

    TGraph *g_xi_31_1sigma = new TGraph(2 * NPOINTS);
    TGraph *g_xi_31_90CL = new TGraph(2 * NPOINTS);

    for(int i = 0; i < NPOINTS; ++ i)
    {
        const double G0 = G0_ps_integral_yrinv;
        const double G2 = G2_ps_integral_yrinv;

        const double A0_150ND = 3.45e-04;
        const double NA = 6.022e+23;
        const double MASS = 36.6;
        const double YR_TO_SEC = 31557600.0;
        const double SEC_TO_YR = 1.0 / YR_TO_SEC;
        const double LOG2 = std::log(2.0);
        const double N_ISO = 150.0;
        const double NA_SEC_TO_YR = NA * SEC_TO_YR;

        const double activity = A_CV * A0_150ND;

        double T12_CV = LOG2 * (MASS / N_ISO) * (NA_SEC_TO_YR / activity);
        //std::cout << "T12_CV = " << T12_CV << " Yr" << std::endl;

        //double xi_31 = xi_31_CV;
        double T12 = T12_CV;

        double value_90CL_h = std::pow(xi_31_90CL_h, 2.0) / (T12 * (G0 + xi_31_90CL_h * G2));
        double value_1sigma_h = std::pow(xi_31_1sigma_h, 2.0) / (T12 * (G0 + xi_31_1sigma_h * G2));
        double value_CV = std::pow(xi_31_CV, 2.0) / (T12 * (G0 + xi_31_CV * G2));
        double value_1sigma_l = std::pow(xi_31_1sigma_l, 2.0) / (T12 * (G0 + xi_31_1sigma_l * G2));
        double value_90CL_l = std::pow(xi_31_90CL_l, 2.0) / (T12 * (G0 + xi_31_90CL_l * G2));

        // TODO: some way to make this smooth
        double x = 0.01 + 0.1 * (i / (double)(NPOINTS));
        double MGT_3 = x;
        double gA = std::pow(value_CV / std::pow(MGT_3, 2.0), 0.25);
        double y = gA;
        g_xi_31_CV->SetPoint(i, x, y);

        g_xi_31_1sigma_h->SetPoint(i, MGT_3, std::pow(value_1sigma_h / std::pow(MGT_3, 2.0), 0.25));
        g_xi_31_90CL_h->SetPoint(i, MGT_3, std::pow(value_90CL_h / std::pow(MGT_3, 2.0), 0.25));
        g_xi_31_90CL_l->SetPoint(i, MGT_3, std::pow(value_90CL_l / std::pow(MGT_3, 2.0), 0.25));
        g_xi_31_1sigma_l->SetPoint(i, MGT_3, std::pow(value_1sigma_l / std::pow(MGT_3, 2.0), 0.25));
    }

    for(int i = 0; i < NPOINTS; ++ i)
    {
        double x, y;

        g_xi_31_1sigma_l->GetPoint(i, x, y);
        g_xi_31_1sigma->SetPoint(i, x, y);
        
        g_xi_31_90CL_l->GetPoint(i, x, y);
        g_xi_31_90CL->SetPoint(i, x, y);
    }
    for(int i = 0; i < NPOINTS; ++ i)
    {
        double x, y;

        g_xi_31_1sigma_h->GetPoint(NPOINTS - i - 1, x, y);
        g_xi_31_1sigma->SetPoint(NPOINTS + i, x, y);
        
        g_xi_31_90CL_h->GetPoint(NPOINTS - i - 1, x, y);
        g_xi_31_90CL->SetPoint(NPOINTS + i, x, y);
    }

    for(int i = 0; i < 2 * NPOINTS; ++ i)
    {
        double x, y;
        g_xi_31_90CL->GetPoint(i, x, y);
        std::cout << x << " " << y << std::endl;
    }



    TCanvas *canvas = new TCanvas("canvas", "canvas");
    canvas->SetTicks(2, 2);
    canvas->SetBottomMargin(0.12);

    g_xi_31_CV->SetTitle(0);
    g_xi_31_CV->SetLineWidth(2);

    g_xi_31_CV->SetMaximum(1.25);
    g_xi_31_CV->SetMinimum(0.2);
    g_xi_31_CV->GetXaxis()->SetRangeUser(0.0, 0.1);

    g_xi_31_CV->GetXaxis()->SetTitle("M_{GT-3}^{2#nu}");
    g_xi_31_CV->GetYaxis()->SetTitle("g_{A}^{eff}");
    g_xi_31_CV->GetXaxis()->SetTitleFont(43);
    g_xi_31_CV->GetXaxis()->SetTitleSize(20);
    g_xi_31_CV->GetXaxis()->SetLabelFont(63);
    g_xi_31_CV->GetXaxis()->SetLabelSize(18);
    g_xi_31_CV->GetXaxis()->SetTitleOffset(1.2);
    g_xi_31_CV->GetYaxis()->SetTitleFont(43);
    g_xi_31_CV->GetYaxis()->SetTitleSize(20);
    g_xi_31_CV->GetYaxis()->SetLabelFont(63);
    g_xi_31_CV->GetYaxis()->SetLabelSize(18);
    g_xi_31_CV->GetYaxis()->ChangeLabel(1, 0, 0);

    g_xi_31_90CL->SetFillColor(kBlue - 9);
    g_xi_31_1sigma->SetFillColor(kViolet - 9);


    g_xi_31_CV->Draw("al");

    g_xi_31_90CL->Draw("fsame");
    g_xi_31_1sigma->Draw("fsame");

    //g_xi_31_CV->Draw("lsame");

//    g_xi_31_90CL_h->Draw("lsame");
//    g_xi_31_1sigma_h->Draw("lsame");
//    g_xi_31_1sigma_l->Draw("lsame");
//    g_xi_31_90CL_l->Draw("lsame");

    g_gamgt3_nd150_arg->Draw("lsame");
    g_gamgt3_nd150_ssd->Draw("lsame");

    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(63);
    latex.SetTextSize(20);

    TString s90;
    s90.Form("90 %% CL");
    TString s1sigma;
    s1sigma.Form("68.3 %% CL");

    latex.SetTextColor(kBlue - 9);
    latex.DrawLatex(0.73, 0.82, s90);
    latex.SetTextColor(kViolet - 9);
    latex.DrawLatex(0.73, 0.75, s1sigma);

    double posx = 0.15;
    double posy = 0.15;
    double sizex = 0.3;
    double sizey = 0.15;
    TLegend *leg = new TLegend(posx, posy, posx + sizex, posy + sizey);
    leg->AddEntry(g_gamgt3_nd150_arg, "sQRPA (Argonne)", "L");
    leg->AddEntry(g_gamgt3_nd150_ssd, "SSD", "L");
    leg->SetTextFont(63);
    leg->SetTextSize(20);
    leg->SetBorderSize(5);
    leg->SetShadowColor(kGray);
    leg->SetMargin(0.2);
    leg->Draw("BL");

}


#endif // NEWLOGLIKFITTER_DRAW_GA_METHOD_UNTITLED_H
