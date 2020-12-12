#ifndef NEWLOGLIKFITTER_DRAW_GA_METHOD_CONTOUR_H
#define NEWLOGLIKFITTER_DRAW_GA_METHOD_CONTOUR_H

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

#endif // NEWLOGLIKFITTER_DRAW_GA_METHOD_CONTOUR_H
