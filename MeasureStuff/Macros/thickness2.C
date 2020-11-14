


void thickness2()
{

    gStyle->SetPaintTextFormat(".1f");

    TString fname_base = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150";
    TString fname_nominal = fname_base + "/nd150_nominalThick/Nd150_2eNg_output_postprocessed_small.root";
    TString fname_double = fname_base + "/nd150_doubleThick/Nd150_2eNg_output_postprocessed_small.root";
    TString fname_half = fname_base + "/nd150_halfThick/Nd150_2eNg_output_postprocessed_small.root";

    ///////////////////////////////////////////////////////////////////////////
    // nominal
    ///////////////////////////////////////////////////////////////////////////

    TFile *f_nominal = new TFile(fname_nominal);
    TTree *tree_nominal = (TTree*)f_nominal->Get("Nd150_2eNg/Nd150_2eNg");
    Double_t electronEnergy_nominal[2];
    tree_nominal->SetBranchAddress("electronEnergy", electronEnergy_nominal);

    TH2D *h_nominal = new TH2D("h_nominal", "h_nominal", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH1D *hSingleEnergy_nominal = new TH1D("hSingleEnergy_nominal", "hSingleEnergy_nominal", 50, 0.0, 5.0);

    for(Long64_t ix = 0; ix < tree_nominal->GetEntries(); ++ ix)
    {
        tree_nominal->GetEntry(ix);
        h_nominal->Fill(electronEnergy_nominal[0], electronEnergy_nominal[1]);
        hSingleEnergy_nominal->Fill(electronEnergy_nominal[0]);
        hSingleEnergy_nominal->Fill(electronEnergy_nominal[1]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // double
    ///////////////////////////////////////////////////////////////////////////

    TFile *f_double = new TFile(fname_double);
    TTree *tree_double = (TTree*)f_double->Get("Nd150_2eNg/Nd150_2eNg");
    Double_t electronEnergy_double[2];
    tree_double->SetBranchAddress("electronEnergy", electronEnergy_double);

    TH2D *h_double = new TH2D("h_double", "h_double", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH1D *hSingleEnergy_double = new TH1D("hSingleEnergy_double", "hSingleEnergy_double", 50, 0.0, 5.0);

    for(Long64_t ix = 0; ix < tree_double->GetEntries(); ++ ix)
    {
        tree_double->GetEntry(ix);
        h_double->Fill(electronEnergy_double[0], electronEnergy_double[1]);
        hSingleEnergy_double->Fill(electronEnergy_double[0]);
        hSingleEnergy_double->Fill(electronEnergy_double[1]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // half
    ///////////////////////////////////////////////////////////////////////////

    TFile *f_half = new TFile(fname_half);
    TTree *tree_half = (TTree*)f_half->Get("Nd150_2eNg/Nd150_2eNg");
    Double_t electronEnergy_half[2];
    tree_half->SetBranchAddress("electronEnergy", electronEnergy_half);

    TH2D *h_half = new TH2D("h_half", "h_half", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH1D *hSingleEnergy_half = new TH1D("hSingleEnergy_half", "hSingleEnergy_half", 50, 0.0, 5.0);

    for(Long64_t ix = 0; ix < tree_half->GetEntries(); ++ ix)
    {
        tree_half->GetEntry(ix);
        h_half->Fill(electronEnergy_half[0], electronEnergy_half[1]);
        hSingleEnergy_half->Fill(electronEnergy_half[0]);
        hSingleEnergy_half->Fill(electronEnergy_half[1]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // set errors
    ///////////////////////////////////////////////////////////////////////////

    for(Int_t j = 1; j <= h_nominal->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_nominal->GetNbinsX(); ++ i)
        {
            Double_t content_nominal = h_nominal->GetBinContent(i, j);
            Double_t content_double = h_double->GetBinContent(i, j);
            Double_t content_half = h_half->GetBinContent(i, j);
            h_nominal->SetBinError(i, j, std::sqrt(content_nominal));
            h_double->SetBinError(i, j, std::sqrt(content_double));
            h_half->SetBinError(i, j, std::sqrt(content_half));
        }
    }

/*
    TCanvas *c_nominal = new TCanvas("c_nominal", "c_nominal");
    h_nominal->SetContour(1000);
    h_nominal->Draw("colz");

    TCanvas *c_double = new TCanvas("c_double", "c_double");
    h_double->SetContour(1000);
    h_double->Draw("colz");

    TCanvas *c_half = new TCanvas("c_half", "c_half");
    h_half->SetContour(1000);
    h_half->Draw("colz");
*/
    TH2D *h_ratio_double_nominal = new TH2D("h_ratio_double_nominal", "h_ratio_double_nominal", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_half_nominal = new TH2D("h_ratio_half_nominal", "h_ratio_half_nominal", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_double_virtual = new TH2D("h_ratio_double_virtual", "h_ratio_double_virtual", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_half_virtual = new TH2D("h_ratio_half_virtual", "h_ratio_half_virtual", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_ratio_double_nominal->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio_double_nominal->GetNbinsX(); ++ i)
        {
            Double_t content_double = h_double->GetBinContent(i, j);
            Double_t content_half = h_half->GetBinContent(i, j);
            Double_t content_nominal = h_nominal->GetBinContent(i, j);
            Double_t content_virtual = 0.5 * (content_double + content_half);
            
            // nominal
            //if((content_nominal != 0.0) && (content_double != 0.0))
            if(content_nominal != 0.0)
            {
                h_ratio_double_nominal->SetBinContent(i, j, content_double / content_nominal);
            }
            else
            {
                h_ratio_double_nominal->SetBinContent(i, j, 1.0);
            }
            //if((content_nominal != 0.0) && (content_half != 0.0))
            if(content_nominal != 0.0)
            {
                h_ratio_half_nominal->SetBinContent(i, j, content_half / content_nominal);
            }
            else
            {
                h_ratio_half_nominal->SetBinContent(i, j, 1.0);
            }

            // virtual
            //if((content_virtual != 0.0) && (content_double != 0.0))
            if(content_virtual != 0.0)
            {
                h_ratio_double_virtual->SetBinContent(i, j, content_double / content_virtual);
            }
            else
            {
                h_ratio_double_virtual->SetBinContent(i, j, 1.0);
            }
            //if((content_virtual != 0.0) && (content_half != 0.0))
            if(content_virtual != 0.0)
            {
                h_ratio_half_virtual->SetBinContent(i, j, content_half / content_virtual);
            }
            else
            {
                h_ratio_half_virtual->SetBinContent(i, j, 1.0);
            }
        }
    }


    TCanvas *c_ratio_double_nominal = new TCanvas("c_ratio_double_nominal", "c_ratio_double_nominal");
    h_ratio_double_nominal->SetContour(1000);
    h_ratio_double_nominal->Draw("colz");

    TCanvas *c_ratio_half_nominal = new TCanvas("c_ratio_half_nominal", "c_ratio_half_nominal");
    h_ratio_half_nominal->SetContour(1000);
    h_ratio_half_nominal->Draw("colz");
    
    TCanvas *c_ratio_double_virtual = new TCanvas("c_ratio_double_virtual", "c_ratio_double_virtual");
    h_ratio_double_virtual->SetContour(1000);
    h_ratio_double_virtual->Draw("colz");
    
    TCanvas *c_ratio_half_virtual = new TCanvas("c_ratio_half_virtual", "c_ratio_half_virtual");
    h_ratio_half_virtual->SetContour(1000);
    h_ratio_half_virtual->Draw("colz");


/*
    TH2D *h_ratio = new TH2D("h_ratio", "h_ratio", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_error = new TH2D("h_ratio_error", "h_ratio_error", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_diff = new TH2D("h_diff", "h_diff", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_diff->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_diff->GetNbinsX(); ++ i)
        {
            //Double_t content_nominal = h_nominal->GetBinContent(i, j);
            Double_t content_double = h_double->GetBinContent(i, j);
            Double_t content_half = h_half->GetBinContent(i, j);
            //Double_t error_nominal = h_nominal->GetBinError(i, j);
            Double_t error_double = h_double->GetBinError(i, j);
            Double_t error_half = h_half->GetBinError(i, j);
            if((content_double != 0.0) && (content_half != 0.0))
            {
                // required to ensure systematic has same magnitude as foil thickness
                Double_t content_virtual = 0.5 * (content_double + content_half);
                Double_t error_virtual = std::sqrt(content_virtual);
                h_ratio->SetBinContent(i, j, content_double / content_virtual);
                Double_t error = std::sqrt(std::pow(error_double / content_virtual, 2.0) + std::pow(error_virtual * content_double / (content_virtual * content_virtual), 2.0));
                h_ratio->SetBinError(i, j, error);
                h_ratio_error->SetBinContent(i, j, error);
                
                h_diff->SetBinContent(i, j, content_double - content_virtual);
            }
        }
    }



    TCanvas *c_diff = new TCanvas("c_diff", "c_diff");
    h_diff->SetContour(1000);
    h_diff->Draw("colz");

    TCanvas *c_ratio = new TCanvas("c_ratio", "c_ratio");
    h_ratio->SetContour(1000);
    h_ratio->Draw("colz");

    TCanvas *c_ratio_error = new TCanvas("c_ratio_error", "c_ratio_error");
    h_ratio_error->SetContour(1000);
    h_ratio_error->Draw("colz");
*/

    // create systematic 2d histograms
    //TH2D *h_sys_h = new TH2D("h_sys_h", "h_sys_h", 50, 0.0, 5.0, 50, 0.0, 5.0);
    //TH2D *h_sys_l = new TH2D("h_sys_l", "h_sys_l", 50, 0.0, 5.0, 50, 0.0, 5.0);
    //TH2D *h_sys = new TH2D("h_systematic_foil_thickness", "h_systematic_foil_thickness", 50, 0.0, 5.0, 50, 0.0, 5.0);
    /*
    for(Int_t j = 1; j <= h_ratio->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio->GetNbinsX(); ++ i)
        {
            Double_t content_double = h_double->GetBinContent(i, j);
            Double_t content_half = h_half->GetBinContent(i, j);
            if((content_double == 0.0) && (content_half == 0.0)) continue;
            Double_t error_ratio = h_ratio_error->GetBinContent(i, j);
            Double_t error = error_ratio;
            Double_t h = h_ratio->GetBinContent(i, j);
            //Double_t h = 0.0 + error;
            //Double_t l = 0.0 - error;
            //h_sys_h->SetBinContent(i, j, h);
            //h_sys_l->SetBinContent(i, j, l);
            h_sys->SetBinContent(i, j, h);
        }
    }
    */

//    TCanvas *c_sys = new TCanvas("c_sys", "c_sys");
//    h_sys->SetContour(1000);
//    h_sys->Draw("colz");



    TH2D *h_sys_h = new TH2D("h_systematic_foil_thickness_h", "h_systematic_foil_thickness_h", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_sys_l = new TH2D("h_systematic_foil_thickness_l", "h_systematic_foil_thickness_l", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_ratio_double_nominal->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio_double_nominal->GetNbinsX(); ++ i)
        {
            Double_t h = h_ratio_double_nominal->GetBinContent(i, j);
            Double_t l = h_ratio_half_nominal->GetBinContent(i, j);
            h_sys_h->SetBinContent(i, j, h);
            h_sys_l->SetBinContent(i, j, l);
        }
    } 

    TFile *f = new TFile("systematic_foil_thickness_nominal.root", "RECREATE");
    h_sys_h->Write();
    h_sys_l->Write();
    f->Close();

    for(Int_t j = 1; j <= h_ratio_double_virtual->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio_double_virtual->GetNbinsX(); ++ i)
        {
            Double_t h = h_ratio_double_virtual->GetBinContent(i, j);
            Double_t l = h_ratio_half_virtual->GetBinContent(i, j);
            h_sys_h->SetBinContent(i, j, h);
            h_sys_l->SetBinContent(i, j, l);
        }
    } 

    f = new TFile("systematic_foil_thickness_virtual.root", "RECREATE");
    h_sys_h->Write();
    h_sys_l->Write();
    f->Close();


/*
    TFile *systematic_foil_thickness = new TFile("systematic_foil_thickness.root", "RECREATE");
    h_sys_h->Write();
    //h_sys_l->Write();
    systematic_foil_thickness->Close();
*/
#if 0
    TCanvas *c_ratio_double_nominal = new TCanvas("c_ratio_double_nominal", "c_ratio_double_nominal");
    h_ratio_double_nominal->SetContour(1000);
    h_ratio_double_nominal->SetMarkerSize(0.5);
    h_ratio_double_nominal->Draw("colz");

    TCanvas *c_ratio_half_nominal = new TCanvas("c_ratio_half_nominal", "c_ratio_half_nominal");
    h_ratio_half_nominal->SetContour(1000);
    h_ratio_half_nominal->SetMarkerSize(0.5);
    h_ratio_half_nominal->Draw("colz");

    TCanvas *c_ratio_error_double_nominal = new TCanvas("c_ratio_error_double_nominal", "c_ratio_error_double_nominal");
    h_ratio_error_double_nominal->SetContour(1000);
    h_ratio_error_double_nominal->SetMarkerSize(0.5);
    h_ratio_error_double_nominal->Draw("colz");

    TCanvas *c_ratio_error_half_nominal = new TCanvas("c_ratio_error_half_nominal", "c_ratio_error_half_nominal");
    h_ratio_error_half_nominal->SetContour(1000);
    h_ratio_error_half_nominal->SetMarkerSize(0.5);
    h_ratio_error_half_nominal->Draw("colz");

    TCanvas *c_diff_double_nominal = new TCanvas("c_diff_double_nominal", "c_diff_double_nominal");
    h_diff_double_nominal->SetContour(1000);
    h_diff_double_nominal->Draw("colz");

    TCanvas *c_diff_half_nominal = new TCanvas("c_diff_half_nominal", "c_diff_half_nominal");
    h_diff_half_nominal->SetContour(1000);
    h_diff_half_nominal->Draw("colz");

    TCanvas *c_sys_h = new TCanvas("c_sys_h", "c_sys_h");
    h_sys_h->SetContour(1000);
    h_sys_h->Draw("colz");

    TCanvas *c_sys_l = new TCanvas("c_sys_l", "c_sys_l");
    h_sys_l->SetContour(1000);
    h_sys_l->Draw("colz");

    for(Int_t i = 1; i <= hSingleEnergy_nominal->GetNbinsX(); ++ i)
    {
        Double_t c_nominal = hSingleEnergy_nominal->GetBinContent(i);
        Double_t c_double = hSingleEnergy_double->GetBinContent(i);
        Double_t c_half = hSingleEnergy_half->GetBinContent(i);

        hSingleEnergy_nominal->SetBinError(i, std::sqrt(c_nominal));
        hSingleEnergy_double->SetBinError(i, std::sqrt(c_double));
        hSingleEnergy_half->SetBinError(i, std::sqrt(c_half));
    }

    //hSingleEnergy_double->Divide(hSingleEnergy_nominal);
    //hSingleEnergy_half->Divide(hSingleEnergy_nominal);

    for(Int_t i = 1; i <= hSingleEnergy_nominal->GetNbinsX(); ++ i)
    {
        Double_t c_nominal = hSingleEnergy_nominal->GetBinContent(i);
        Double_t c_double = hSingleEnergy_double->GetBinContent(i);
        Double_t c_half = hSingleEnergy_half->GetBinContent(i);

        Double_t e_nominal = hSingleEnergy_nominal->GetBinError(i);
        Double_t e_double = hSingleEnergy_double->GetBinError(i);
        Double_t e_half = hSingleEnergy_half->GetBinError(i);

        if((c_double != 0.0) && (c_nominal != 0.0))
        {
            hSingleEnergy_double->SetBinContent(i, c_double / c_nominal);
            hSingleEnergy_double->SetBinError(i,
                std::sqrt(std::pow(e_double / c_nominal, 2.0) + std::pow(e_nominal * c_double / (c_nominal * c_nominal), 2.0))
            );
        }

        if((c_half != 0.0) && (c_nominal != 0.0))
        {
            hSingleEnergy_half->SetBinContent(i, c_half / c_nominal);
            hSingleEnergy_half->SetBinError(i,
                std::sqrt(std::pow(e_half / c_nominal, 2.0) + std::pow(e_nominal * c_half / (c_nominal * c_nominal), 2.0))
            );
        }
    }

    TCanvas *c_SingleEnergy = new TCanvas("c_SingleEnergy", "c_SingleEnergy");
    hSingleEnergy_nominal->SetLineColor(kBlue);
    hSingleEnergy_double->SetLineColor(kRed);
    hSingleEnergy_half->SetLineColor(kGreen);
//    hSingleEnergy_nominal->Draw("hist");
    hSingleEnergy_double->Draw("E");
    hSingleEnergy_half->Draw("Esame");

    TH1D *hSingleEnergy_pull_double = new TH1D("hSingleEnergy_pull_double", "hSingleEnergy_pull_double", 50, 0.0, 5.0);
    TH1D *hSingleEnergy_pull_half = new TH1D("hSingleEnergy_pull_half", "hSingleEnergy_pull_half", 50, 0.0, 5.0);
    for(Int_t i = 1; i <= hSingleEnergy_pull_double->GetNbinsX(); ++ i)
    {
        if(hSingleEnergy_double->GetBinError(i) != 0.0)
        {
            hSingleEnergy_pull_double->SetBinContent(i, (hSingleEnergy_double->GetBinContent(i) - 1.0) / hSingleEnergy_double->GetBinError(i));
        }
        if(hSingleEnergy_half->GetBinError(i) != 0.0)
        {
            hSingleEnergy_pull_half->SetBinContent(i, (hSingleEnergy_half->GetBinContent(i) - 1.0) / hSingleEnergy_half->GetBinError(i));
        }
    }
    TCanvas *cSingleEnergy_pull = new TCanvas("cSingleEnergy_pull", "cSingleEnergy_pull");
    hSingleEnergy_pull_double->GetYaxis()->SetRangeUser(-4.0, 2.0);
    hSingleEnergy_pull_double->Draw("hist");
    hSingleEnergy_pull_half->SetLineColor(kRed);
    hSingleEnergy_pull_half->Draw("histsame");

    TH1D *h_systematic_thickness_1D = new TH1D("h_systematic_thickness_1D", "h_systematic_thickness_1D", 50, 0.0, 5.0);
    for(Int_t i = 1; i <= h_systematic_thickness_1D->GetNbinsX(); ++ i)
    {
        Double_t content1 = hSingleEnergy_double->GetBinError(i);
        Double_t content2 = hSingleEnergy_half->GetBinError(i);
        Double_t content = std::sqrt(std::pow(content1, 2.0)+ std::pow(content2, 2.0));
        h_systematic_thickness_1D->SetBinContent(i, content);
    }
    TCanvas *c_systematic_thickness_1D = new TCanvas("c_systematic_thickness_1D", "c_systematic_thickness_1D");
    h_systematic_thickness_1D->Draw("hist");


    TH2D *h_systematic_thickness = new TH2D("h_systematic_thickness", "h_systematic_thickness", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_systematic_thickness->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_systematic_thickness->GetNbinsX(); ++ i)
        {
            Double_t content1 = h_ratio_double_nominal->GetBinError(i, j);
            Double_t content2 = h_ratio_half_nominal->GetBinError(i, j);
            Double_t content = std::sqrt(std::pow(content1, 2.0)+ std::pow(content2, 2.0));
            h_systematic_thickness->SetBinContent(i, j, content);
        }
    }

    TCanvas *c_systematic_thickness = new TCanvas("c_systematic_thickness", "c_systematic_thickness");
    h_systematic_thickness->SetContour(1000);
    h_systematic_thickness->SetMarkerSize(0.5);
    h_systematic_thickness->Draw("colz");

#endif



}
