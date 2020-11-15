

void dEdX()
{

    
    TString fname_base = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150";
    TString fname_nominal = fname_base + "/nd150_nominalThick/Nd150_2eNg_output_postprocessed_small.root";
    TString fname_HdEdX = fname_base + "/HdEdX_2vBB/Nd150_2eNg_output_postprocessed_small.root";
    TString fname_LdEdX = fname_base + "/LdEdX_2vBB/Nd150_2eNg_output_postprocessed_small.root";

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
    // HdEdX
    ///////////////////////////////////////////////////////////////////////////

    TFile *f_HdEdX = new TFile(fname_HdEdX);
    TTree *tree_HdEdX = (TTree*)f_HdEdX->Get("Nd150_2eNg/Nd150_2eNg");
    Double_t electronEnergy_HdEdX[2];
    tree_HdEdX->SetBranchAddress("electronEnergy", electronEnergy_HdEdX);

    TH2D *h_HdEdX = new TH2D("h_HdEdX", "h_HdEdX", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH1D *hSingleEnergy_HdEdX = new TH1D("hSingleEnergy_HdEdX", "hSingleEnergy_HdEdX", 50, 0.0, 5.0);

    for(Long64_t ix = 0; ix < tree_HdEdX->GetEntries(); ++ ix)
    {
        tree_HdEdX->GetEntry(ix);
        int lowE_index = -1;
        int highE_index = -1;
        if(electronEnergy_HdEdX[0] > electronEnergy_HdEdX[1])
        {
            highE_index = 0;
            lowE_index = 1;
        }
        else
        {
            highE_index = 1;
            lowE_index = 0;
        }
        h_HdEdX->Fill(electronEnergy_HdEdX[lowE_index], electronEnergy_HdEdX[highE_index]);
        hSingleEnergy_HdEdX->Fill(electronEnergy_HdEdX[lowE_index]);
        hSingleEnergy_HdEdX->Fill(electronEnergy_HdEdX[highE_index]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // LdEdX
    ///////////////////////////////////////////////////////////////////////////

    TFile *f_LdEdX = new TFile(fname_LdEdX);
    TTree *tree_LdEdX = (TTree*)f_LdEdX->Get("Nd150_2eNg/Nd150_2eNg");
    Double_t electronEnergy_LdEdX[2];
    tree_LdEdX->SetBranchAddress("electronEnergy", electronEnergy_LdEdX);

    TH2D *h_LdEdX = new TH2D("h_LdEdX", "h_LdEdX", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH1D *hSingleEnergy_LdEdX = new TH1D("hSingleEnergy_LdEdX", "hSingleEnergy_LdEdX", 50, 0.0, 5.0);

    for(Long64_t ix = 0; ix < tree_LdEdX->GetEntries(); ++ ix)
    {
        tree_LdEdX->GetEntry(ix);
        int lowE_index = -1;
        int highE_index = -1;
        if(electronEnergy_LdEdX[0] > electronEnergy_LdEdX[1])
        {
            highE_index = 0;
            lowE_index = 1;
        }
        else
        {
            highE_index = 1;
            lowE_index = 0;
        }
        h_LdEdX->Fill(electronEnergy_LdEdX[lowE_index], electronEnergy_LdEdX[highE_index]);
        hSingleEnergy_LdEdX->Fill(electronEnergy_LdEdX[lowE_index]);
        hSingleEnergy_LdEdX->Fill(electronEnergy_LdEdX[highE_index]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // set errors
    ///////////////////////////////////////////////////////////////////////////

    for(Int_t j = 1; j <= h_HdEdX->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_HdEdX->GetNbinsX(); ++ i)
        {
            //Double_t content_nominal = h_nominal->GetBinContent(i, j);
            Double_t content_HdEdX = h_HdEdX->GetBinContent(i, j);
            Double_t content_LdEdX = h_LdEdX->GetBinContent(i, j);
            //h_nominal->SetBinError(i, j, std::sqrt(content_nominal));
            h_HdEdX->SetBinError(i, j, std::sqrt(content_HdEdX));
            h_LdEdX->SetBinError(i, j, std::sqrt(content_LdEdX));
        }
    }

/*
    TCanvas *c_HdEdX = new TCanvas("c_HdEdX", "c_HdEdX");
    h_HdEdX->SetContour(1000);
    h_HdEdX->Draw("colz");

    TCanvas *c_LdEdX = new TCanvas("c_LdEdX", "c_LdEdX");
    h_LdEdX->SetContour(1000);
    h_LdEdX->Draw("colz");
*/


    TH2D *h_ratio_HdEdX_nominal = new TH2D("h_ratio_HdEdX_nominal", "h_ratio_HdEdX_nominal", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_LdEdX_nominal = new TH2D("h_ratio_LdEdX_nominal", "h_ratio_LdEdX_nominal", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_HdEdX_virtual = new TH2D("h_ratio_HdEdX_virtual", "h_ratio_HdEdX_virtual", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_ratio_LdEdX_virtual = new TH2D("h_ratio_LdEdX_virtual", "h_ratio_LdEdX_virtual", 50, 0.0, 5.0, 50, 0.0, 5.0);
//    TH2D *h_special = new TH2D("h_special", "special", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_ratio_HdEdX_nominal->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio_HdEdX_nominal->GetNbinsX(); ++ i)
        {
            Double_t content_HdEdX = h_HdEdX->GetBinContent(i, j);
            Double_t content_LdEdX = h_LdEdX->GetBinContent(i, j);
            Double_t content_nominal = h_nominal->GetBinContent(i, j);
            Double_t content_virtual = 0.5 * (content_HdEdX + content_LdEdX);
            
            // nominal
            //if((content_nominal != 0.0) && (content_HdEdX != 0.0))
            if(content_nominal != 0.0)
            {
                h_ratio_HdEdX_nominal->SetBinContent(i, j, content_HdEdX / (5.0 * content_nominal));
            }
            else
            {
                h_ratio_HdEdX_nominal->SetBinContent(i, j, 1.0);
            }
            //if((content_nominal != 0.0) && (content_LdEdX != 0.0))
            if(content_nominal != 0.0)
            {
                h_ratio_LdEdX_nominal->SetBinContent(i, j, content_LdEdX / (5.0 * content_nominal));
            }
            else
            {
                h_ratio_LdEdX_nominal->SetBinContent(i, j, 1.0);
            }

            // virtual
            //if((content_virtual != 0.0) && (content_HdEdX != 0.0))
            if(content_virtual != 0.0)
            {
                h_ratio_HdEdX_virtual->SetBinContent(i, j, content_HdEdX / content_virtual);
            }
            else
            {
                h_ratio_HdEdX_virtual->SetBinContent(i, j, 1.0);
            }
            //if((content_virtual != 0.0) && (content_LdEdX != 0.0))
            if(content_virtual != 0.0)
            {
                h_ratio_LdEdX_virtual->SetBinContent(i, j, content_LdEdX / content_virtual);
            }
            else
            {
                h_ratio_LdEdX_virtual->SetBinContent(i, j, 1.0);
            }

//            h_special->SetBinContent(i, j, content_HdEdX / content_virtual * 
        }
    }

    TCanvas *c_ratio_HdEdX_nominal = new TCanvas("c_ratio_HdEdX_nominal", "c_ratio_HdEdX_nominal");
    h_ratio_HdEdX_nominal->SetContour(1000);
    h_ratio_HdEdX_nominal->Draw("colz");

    TCanvas *c_ratio_LdEdX_nominal = new TCanvas("c_ratio_LdEdX_nominal", "c_ratio_LdEdX_nominal");
    h_ratio_LdEdX_nominal->SetContour(1000);
    h_ratio_LdEdX_nominal->Draw("colz");
    
    TCanvas *c_ratio_HdEdX_virtual = new TCanvas("c_ratio_HdEdX_virtual", "c_ratio_HdEdX_virtual");
    h_ratio_HdEdX_virtual->SetContour(1000);
    h_ratio_HdEdX_virtual->Draw("colz");
    
    TCanvas *c_ratio_LdEdX_virtual = new TCanvas("c_ratio_LdEdX_virtual", "c_ratio_LdEdX_virtual");
    h_ratio_LdEdX_virtual->SetContour(1000);
    h_ratio_LdEdX_virtual->Draw("colz");


/*
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
    TH2D *h_sys = new TH2D("h_systematic_dEdX", "h_systematic_dEdX", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_ratio->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio->GetNbinsX(); ++ i)
        {
            Double_t content_HdEdX = h_HdEdX->GetBinContent(i, j);
            Double_t content_LdEdX = h_LdEdX->GetBinContent(i, j);
            if((content_HdEdX == 0.0) && (content_LdEdX == 0.0)) continue;
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
    TFile *f = new TFile("systematic_dEdX.root", "RECREATE");
    //h_sys->Write();
    h_sys_h->Write();
    h_sys_l->Write();
    f->Close();
*/

    TH2D *h_sys_h = new TH2D("h_systematic_dEdX_h", "h_systematic_dEdX_h", 50, 0.0, 5.0, 50, 0.0, 5.0);
    TH2D *h_sys_l = new TH2D("h_systematic_dEdX_l", "h_systematic_dEdX_l", 50, 0.0, 5.0, 50, 0.0, 5.0);
    for(Int_t j = 1; j <= h_ratio_HdEdX_nominal->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio_HdEdX_nominal->GetNbinsX(); ++ i)
        {
            Double_t h = h_ratio_HdEdX_nominal->GetBinContent(i, j);
            Double_t l = h_ratio_LdEdX_nominal->GetBinContent(i, j);
            h_sys_h->SetBinContent(i, j, h);
            h_sys_l->SetBinContent(i, j, l);
        }
    } 

    TFile *f = new TFile("systematic_dEdX_nominal.root", "RECREATE");
    h_sys_h->Write();
    h_sys_l->Write();
    f->Close();

    for(Int_t j = 1; j <= h_ratio_HdEdX_virtual->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_ratio_HdEdX_virtual->GetNbinsX(); ++ i)
        {
            Double_t h = h_ratio_HdEdX_virtual->GetBinContent(i, j);
            Double_t l = h_ratio_LdEdX_virtual->GetBinContent(i, j);
            h_sys_h->SetBinContent(i, j, h);
            h_sys_l->SetBinContent(i, j, l);
        }
    } 

    f = new TFile("systematic_dEdX_virtual.root", "RECREATE");
    h_sys_h->Write();
    h_sys_l->Write();
    f->Close();

}
