

void brem()
{

    
    TString fname_base = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150";
    TString fname_nominal = fname_base + "/nd150_nominalThick/Nd150_2eNg_output_postprocessed_small.root";
    TString fname_Hbrem = fname_base + "/Hbrem_2vBB/Nd150_2eNg_output_postprocessed_small.root";
    TString fname_Lbrem = fname_base + "/Lbrem_2vBB/Nd150_2eNg_output_postprocessed_small.root";

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
        int lowE_index = -1;
        int highE_index = -1;
        if(electronEnergy_nominal[0] > electronEnergy_nominal[1])
        {
            highE_index = 0;
            lowE_index = 1;
        }
        else
        {
            highE_index = 1;
            lowE_index = 0;
        }
        h_nominal->Fill(electronEnergy_nominal[lowE_index], electronEnergy_nominal[highE_index]);
        hSingleEnergy_nominal->Fill(electronEnergy_nominal[lowE_index]);
        hSingleEnergy_nominal->Fill(electronEnergy_nominal[highE_index]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Hbrem
    ///////////////////////////////////////////////////////////////////////////

    TFile *f_Hbrem = new TFile(fname_Hbrem);
    TTree *tree_Hbrem = (TTree*)f_Hbrem->Get("Nd150_2eNg/Nd150_2eNg");
    Double_t electronEnergy_Hbrem[2];
    tree_Hbrem->SetBranchAddress("electronEnergy", electronEnergy_Hbrem);

    TH2D *h_Hbrem = new TH2D("h_Hbrem", "h_Hbrem", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH1D *hSingleEnergy_Hbrem = new TH1D("hSingleEnergy_Hbrem", "hSingleEnergy_Hbrem", 50, 0.0, 5.0);

    for(Long64_t ix = 0; ix < tree_Hbrem->GetEntries(); ++ ix)
    {
        tree_Hbrem->GetEntry(ix);
        int lowE_index = -1;
        int highE_index = -1;
        if(electronEnergy_Hbrem[0] > electronEnergy_Hbrem[1])
        {
            highE_index = 0;
            lowE_index = 1;
        }
        else
        {
            highE_index = 1;
            lowE_index = 0;
        }
        h_Hbrem->Fill(electronEnergy_Hbrem[lowE_index], electronEnergy_Hbrem[highE_index]);
        hSingleEnergy_Hbrem->Fill(electronEnergy_Hbrem[lowE_index]);
        hSingleEnergy_Hbrem->Fill(electronEnergy_Hbrem[highE_index]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Lbrem
    ///////////////////////////////////////////////////////////////////////////

    TFile *f_Lbrem = new TFile(fname_Lbrem);
    TTree *tree_Lbrem = (TTree*)f_Lbrem->Get("Nd150_2eNg/Nd150_2eNg");
    Double_t electronEnergy_Lbrem[2];
    tree_Lbrem->SetBranchAddress("electronEnergy", electronEnergy_Lbrem);

    TH2D *h_Lbrem = new TH2D("h_Lbrem", "h_Lbrem", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH1D *hSingleEnergy_Lbrem = new TH1D("hSingleEnergy_Lbrem", "hSingleEnergy_Lbrem", 50, 0.0, 5.0);

    for(Long64_t ix = 0; ix < tree_Lbrem->GetEntries(); ++ ix)
    {
        tree_Lbrem->GetEntry(ix);
        int lowE_index = -1;
        int highE_index = -1;
        if(electronEnergy_Lbrem[0] > electronEnergy_Lbrem[1])
        {
            highE_index = 0;
            lowE_index = 1;
        }
        else
        {
            highE_index = 1;
            lowE_index = 0;
        }
        h_Lbrem->Fill(electronEnergy_Lbrem[lowE_index], electronEnergy_Lbrem[highE_index]);
        hSingleEnergy_Lbrem->Fill(electronEnergy_Lbrem[lowE_index]);
        hSingleEnergy_Lbrem->Fill(electronEnergy_Lbrem[highE_index]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // set errors
    ///////////////////////////////////////////////////////////////////////////

    for(Int_t j = 1; j <= h_Hbrem->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_Hbrem->GetNbinsX(); ++ i)
        {
            //Double_t content_nominal = h_nominal->GetBinContent(i, j);
            Double_t content_Hbrem = h_Hbrem->GetBinContent(i, j);
            Double_t content_Lbrem = h_Lbrem->GetBinContent(i, j);
            //h_nominal->SetBinError(i, j, std::sqrt(content_nominal));
            h_Hbrem->SetBinError(i, j, std::sqrt(content_Hbrem));
            h_Lbrem->SetBinError(i, j, std::sqrt(content_Lbrem));
        }
    }

/*
    TCanvas *c_Hbrem = new TCanvas("c_Hbrem", "c_Hbrem");
    h_Hbrem->SetContour(1000);
    h_Hbrem->Draw("colz");

    TCanvas *c_Lbrem = new TCanvas("c_Lbrem", "c_Lbrem");
    h_Lbrem->SetContour(1000);
    h_Lbrem->Draw("colz");
*/



    TH2D *h_ratio_Hbrem_nominal = new TH2D("h_ratio_Hbrem_nominal", "h_ratio_Hbrem_nominal", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_Lbrem_nominal = new TH2D("h_ratio_Lbrem_nominal", "h_ratio_Lbrem_nominal", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_Hbrem_virtual = new TH2D("h_ratio_Hbrem_virtual", "h_ratio_Hbrem_virtual", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_Lbrem_virtual = new TH2D("h_ratio_Lbrem_virtual", "h_ratio_Lbrem_virtual", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_ratio_Hbrem_nominal->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio_Hbrem_nominal->GetNbinsX(); ++ i)
        {
            Double_t content_Hbrem = h_Hbrem->GetBinContent(i, j);
            Double_t content_Lbrem = h_Lbrem->GetBinContent(i, j);
            Double_t content_nominal = h_nominal->GetBinContent(i, j);
            Double_t content_virtual = 0.5 * (content_Hbrem + content_Lbrem);
            
            // nominal
            //if((content_nominal != 0.0) && (content_Hbrem != 0.0))
            if(content_nominal != 0.0)
            {
                h_ratio_Hbrem_nominal->SetBinContent(i, j, content_Hbrem / (5.0 * content_nominal));
            }
            else
            {
                h_ratio_Hbrem_nominal->SetBinContent(i, j, 1.0);
            }
            //if((content_nominal != 0.0) && (content_Lbrem != 0.0))
            if(content_nominal != 0.0)
            {
                h_ratio_Lbrem_nominal->SetBinContent(i, j, content_Lbrem / (5.0 * content_nominal));
            }
            else
            {
                h_ratio_Lbrem_nominal->SetBinContent(i, j, 1.0);
            }

            // virtual
            //if((content_virtual != 0.0) && (content_Hbrem != 0.0))
            if(content_virtual != 0.0)
            {
                h_ratio_Hbrem_virtual->SetBinContent(i, j, content_Hbrem / content_virtual);
            }
            else
            {
                h_ratio_Hbrem_virtual->SetBinContent(i, j, 1.0);
            }
            //if((content_virtual != 0.0) && (content_Lbrem != 0.0))
            if(content_virtual != 0.0)
            {
                h_ratio_Lbrem_virtual->SetBinContent(i, j, content_Lbrem / content_virtual);
            }
            else
            {
                h_ratio_Lbrem_virtual->SetBinContent(i, j, 1.0);
            }
        }
    }



    TH2D *h_sys_h = new TH2D("h_systematic_brem_h", "h_systematic_brem_h", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_sys_l = new TH2D("h_systematic_brem_l", "h_systematic_brem_l", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_ratio_Hbrem_nominal->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio_Hbrem_nominal->GetNbinsX(); ++ i)
        {
            Double_t h = h_ratio_Hbrem_nominal->GetBinContent(i, j);
            Double_t l = h_ratio_Lbrem_nominal->GetBinContent(i, j);
            h_sys_h->SetBinContent(i, j, h);
            h_sys_l->SetBinContent(i, j, l);
        }
    } 

    TFile *f = new TFile("systematic_brem_nominal.root", "RECREATE");
    h_sys_h->Write();
    h_sys_l->Write();
    f->Close();

    for(Int_t j = 1; j <= h_ratio_Hbrem_virtual->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio_Hbrem_virtual->GetNbinsX(); ++ i)
        {
            Double_t h = h_ratio_Hbrem_virtual->GetBinContent(i, j);
            Double_t l = h_ratio_Lbrem_virtual->GetBinContent(i, j);
            h_sys_h->SetBinContent(i, j, h);
            h_sys_l->SetBinContent(i, j, l);
        }
    } 

    f = new TFile("systematic_brem_virtual.root", "RECREATE");
    h_sys_h->Write();
    h_sys_l->Write();
    f->Close();









/*
    TH2D *h_ratio = new TH2D("h_ratio", "h_ratio", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_error = new TH2D("h_ratio_error", "h_ratio_error", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_diff = new TH2D("h_diff", "h_diff", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_sys_h = new TH2D("h_systematic_brem_h", "h_systematic_brem_h", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_sys_l = new TH2D("h_systematic_brem_l", "h_systematic_brem_l", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_diff->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_diff->GetNbinsX(); ++ i)
        {
            //Double_t content_nominal = h_nominal->GetBinContent(i, j);
            Double_t content_Hbrem = h_Hbrem->GetBinContent(i, j);
            Double_t content_Lbrem = h_Lbrem->GetBinContent(i, j);
            //Double_t error_nominal = h_nominal->GetBinError(i, j);
            Double_t error_Hbrem = h_Hbrem->GetBinError(i, j);
            Double_t error_Lbrem = h_Lbrem->GetBinError(i, j);
            if((content_Hbrem != 0.0) && (content_Lbrem != 0.0))
            {
                // required to ensure systematic has same magnitude as foil thickness
                Double_t content_virtual = 0.5 * (content_Hbrem + content_Lbrem);
                Double_t error_virtual = std::sqrt(content_virtual);
                h_ratio->SetBinContent(i, j, content_Hbrem / content_virtual);
                Double_t error = std::sqrt(std::pow(error_Hbrem / content_virtual, 2.0) + std::pow(error_virtual * content_Hbrem / (content_virtual * content_virtual), 2.0));
                h_ratio->SetBinError(i, j, error);
                h_ratio_error->SetBinContent(i, j, error);
                
                h_diff->SetBinContent(i, j, content_Hbrem - content_virtual);

                h_sys_h->SetBinContent(i, j, content_Hbrem / content_virtual);
                h_sys_l->SetBinContent(i, j, content_Lbrem / content_virtual);
            }
            else
            {
                h_sys_h->SetBinContent(i, j, 1.0);
                h_sys_l->SetBinContent(i, j, 1.0);
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

/* 
    // create systematic 2d histograms
    //TH2D *h_sys_h = new TH2D("h_sys_h", "h_sys_h", 50, 0.0, 5.0, 50, 0.0, 5.0);
    //TH2D *h_sys_l = new TH2D("h_sys_l", "h_sys_l", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_sys = new TH2D("h_systematic_brem", "h_systematic_brem", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_ratio->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio->GetNbinsX(); ++ i)
        {
            Double_t content_Hbrem = h_Hbrem->GetBinContent(i, j);
            Double_t content_Lbrem = h_Lbrem->GetBinContent(i, j);
            if((content_Hbrem == 0.0) && (content_Lbrem == 0.0)) continue;
            Double_t error_ratio = h_ratio_error->GetBinContent(i, j);
            Double_t error = error_ratio;
            //Double_t h = 0.0 + error;
            //Double_t l = 0.0 - error;
            Double_t h = h_ratio->GetBinContent(i, j);
            //h_sys_h->SetBinContent(i, j, h);
            //h_sys_l->SetBinContent(i, j, l);
            h_sys->SetBinContent(i, j, h);
        }
    } 

    TCanvas *c_sys = new TCanvas("c_sys", "c_sys");
    h_sys->SetContour(1000);
    h_sys->Draw("colz");
*/
/*
    TFile *f = new TFile("systematic_brem.root", "RECREATE");
    //h_sys->Write();
    h_sys_h->Write();
    h_sys_l->Write();
    f->Close();
*/

}
