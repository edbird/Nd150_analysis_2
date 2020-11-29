


TH1D *h_singleEnergyTrue_after = nullptr;
TH1D *h_singleEnergyTrue_before = nullptr;
TH1D *h_totalEnergyTrue_after = nullptr;
TH1D *h_totalEnergyTrue_before = nullptr;






void GetTrueEnergy(
    Double_t* trueElectronEnergy,
    Float_t* Pxntu, Float_t* Pyntu, Float_t* Pzntu)
{
    const Double_t electron_rest_mass = 1.0e-3 * 0.51099895; // GeV
    const Double_t m = electron_rest_mass;
    const Double_t m2 = electron_rest_mass * electron_rest_mass;
    Double_t p0_2 = Pxntu[0] * Pxntu[0] + Pyntu[0] * Pyntu[0] + Pzntu[0] * Pzntu[0]; // GeV
    Double_t trueElectronEnergy_0 = std::sqrt(p0_2 + m2) - m; // GeV
    Double_t p1_2 = Pxntu[1] * Pxntu[1] + Pyntu[1] * Pyntu[1] + Pzntu[1] * Pzntu[1];
    Double_t trueElectronEnergy_1 = std::sqrt(p1_2 + m2) - m;

    trueElectronEnergy[0] = 1.0e+03 * trueElectronEnergy_0;
    trueElectronEnergy[1] = 1.0e+03 * trueElectronEnergy_1;
}





void efficiencyAfter()
{


    TString filePath = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/";
    TString thePath = "nd150/";
    TString sampleName = "nd150_rot_2n2b_m4/";

    TFile *inputFile = new TFile(filePath + thePath + sampleName + "Nd150_2eNg_output_postprocessed_small.root");
    
    if(!inputFile->IsOpen())
    {
        std::cout << "failed to open" << std::endl;
    }
    TTree *inputTree = (TTree*)inputFile->Get("Nd150_2eNg/Nd150_2eNg");

    int run;
    double eventTime;
    double radonWeight;
    double bi210Weight;
    double electronEnergy[2];
    double trueElectronEnergy[2];

    inputTree->SetBranchAddress("Run", &run);
    inputTree->SetBranchAddress("eventTime", &eventTime);
    inputTree->SetBranchAddress("radonWeight", &radonWeight);
    inputTree->SetBranchAddress("bi210Weight", &bi210Weight);
    inputTree->SetBranchAddress("electronEnergy", electronEnergy);
    inputTree->SetBranchAddress("trueElectronEnergy", trueElectronEnergy);

    Long_t events = (Long_t)inputTree->GetEntries();
    for(Long_t event_i = 0; event_i < events; ++ event_i)
    {
    
        inputTree->GetEvent(event_i);


        int highE_index = -1;
        int lowE_index = -1;
        if(electronEnergy[0] > electronEnergy[1])
        {
            highE_index = 0;
            lowE_index = 1;
        }
        else
        {
            highE_index = 1;
            lowE_index = 0;
        }

        const Double_t trueT0{trueElectronEnergy[highE_index]};
        const Double_t trueT1{trueElectronEnergy[lowE_index]};

        // set phase
        // 1869 - 3395 is P1
        // 3396 - 9186 is P2
        bool cut = true;
        Int_t thePhase = -1;
        if((1869 <= run) && (run <= 3395))
        {
            cut = false;
            thePhase = 0;
        }
        else if((3396 <= run) && (run <= 9186))
        {
            cut = false;
            thePhase = 1;
        }
        if(cut == true) continue;


        Double_t el_energy_0{electronEnergy[highE_index]};
        Double_t el_energy_1{electronEnergy[lowE_index]};

        ///////////////////////////////////////////////////////////////
        // standard cuts not implemented in preprocessing
        ///////////////////////////////////////////////////////////////

        if((el_energy_0 < 0.3) || (el_energy_1 < 0.3))
        {
            continue;
        }

        h_singleEnergyTrue_after->Fill(trueT0, 1.0);
        h_singleEnergyTrue_after->Fill(trueT1, 1.0);

        h_totalEnergyTrue_after->Fill(trueT0 + trueT1, 1.0);
    }

}


void efficiencyBefore()
{
 
    const int NUM_NAMES = 20;
    TString names[NUM_NAMES];
    names[0] = "foils_nd150_61_rot_nd150_101_01.root";
    names[1] = "foils_nd150_61_rot_nd150_101_02.root";
    names[2] = "foils_nd150_61_rot_nd150_102_01.root";
    names[3] = "foils_nd150_61_rot_nd150_102_02.root";
    names[4] = "foils_nd150_61_rot_nd150_103_01.root";
    names[5] = "foils_nd150_61_rot_nd150_103_02.root";
    names[6] = "foils_nd150_61_rot_nd150_104_01.root";
    names[7] = "foils_nd150_61_rot_nd150_104_02.root";
    names[8] = "foils_nd150_61_rot_nd150_105_01.root";
    names[9] = "foils_nd150_61_rot_nd150_105_02.root";
    names[10] = "foils_nd150_61_rot_nd150_106_01.root";
    names[11] = "foils_nd150_61_rot_nd150_106_02.root";
    names[12] = "foils_nd150_61_rot_nd150_107_01.root";
    names[13] = "foils_nd150_61_rot_nd150_107_02.root";
    names[14] = "foils_nd150_61_rot_nd150_108_01.root";
    names[15] = "foils_nd150_61_rot_nd150_108_02.root";
    names[16] = "foils_nd150_61_rot_nd150_109_01.root";
    names[17] = "foils_nd150_61_rot_nd150_109_02.root";
    names[18] = "foils_nd150_61_rot_nd150_110_01.root";
    names[19] = "foils_nd150_61_rot_nd150_110_02.root";

    TString filePath = "/mnt/ecb/unix/nemo2/reco/reco_20130206_fix/mc_8.0/rot/";

    TChain *tchain = new TChain("h10", "h10");
    for(int i = 0; i < NUM_NAMES; ++ i)
    {
        std::cout << "add: i=" << i << " name[i]=" << names[i] << std::endl; 
        tchain->Add(filePath + names[i]);
    }

    Int_t Nsc;
    Float_t Sc[2000][12];
    Int_t run = 0;
    Int_t Nvntu = 0;
    Float_t Xvntu = 0.0;
    Float_t Yvntu = 0.0;
    Float_t Zvntu = 0.0;
    Int_t Ntntu = 0;
    Float_t Pxntu[31];
    Float_t Pyntu[31];
    Float_t Pzntu[31];

    tchain->SetBranchAddress("Nsc", &Nsc);
    tchain->SetBranchAddress("Sc", Sc);
    tchain->SetBranchAddress("run", &run);
    tchain->SetBranchAddress("Nvntu", &Nvntu);
    tchain->SetBranchAddress("Xvntu", &Xvntu);
    tchain->SetBranchAddress("Yvntu", &Yvntu);
    tchain->SetBranchAddress("Zvntu", &Zvntu);
    tchain->SetBranchAddress("Ntntu", &Ntntu);
    tchain->SetBranchAddress("Pxntu", Pxntu);
    tchain->SetBranchAddress("Pyntu", Pyntu);
    tchain->SetBranchAddress("Pzntu", Pzntu);


    Double_t trueElectronEnergy[2];

    std::cout << "tchain->GetEntries()=" << tchain->GetEntries() << std::endl;
    std::cout << "number of generated MC was probably 1.0e+07" << std::endl;

    for(Long64_t ix{0}; ix < tchain->GetEntries(); ++ ix)
    {
        if(ix % 100000 == 0)
        {
            std::cout << "ix=" << (double)ix / (double)tchain->GetEntries() * 100.0 << " %" << std::endl;
        }

        tchain->GetEntry(ix);

        GetTrueEnergy(trueElectronEnergy, Pxntu, Pyntu, Pzntu);

        h_singleEnergyTrue_before->Fill(trueElectronEnergy[0], 1.0);
        h_singleEnergyTrue_before->Fill(trueElectronEnergy[1], 1.0);

        h_totalEnergyTrue_before->Fill(trueElectronEnergy[0] + trueElectronEnergy[1], 1.0);
    }

}



int efficiency()
{

    h_singleEnergyTrue_after = new TH1D("h_singleEnergyTrue_after", "h_singleEnergyTrue_after",
                                        50, 0.0, 5.0);
    h_singleEnergyTrue_before = new TH1D("h_singleEnergyTrue_before", "h_singleEnergyTrue_before",
                                        50, 0.0, 5.0);

    h_totalEnergyTrue_after = new TH1D("h_totalEnergyTrue_after", "h_totalEnergyTrue_after",
                                        50, 0.0, 5.0);
    h_totalEnergyTrue_before = new TH1D("h_totalEnergyTrue_before", "h_totalEnergyTrue_before",
                                        50, 0.0, 5.0);

    efficiencyAfter();

    TCanvas *c_singleEnergyTrue_after = new TCanvas("c_singleEnergyTrue_after", "c_singleEnergyTrue_after");
    h_singleEnergyTrue_after->Draw("hist");

    TCanvas *c_totalEnergyTrue_after = new TCanvas("c_totalEnergyTrue_after", "c_totalEnergyTrue_after");
    h_totalEnergyTrue_after->Draw("hist");

    efficiencyBefore();

    TCanvas *c_singleEnergyTrue_before = new TCanvas("c_singleEnergyTrue_before", "c_singleEnergyTrue_before");
    h_singleEnergyTrue_before->Draw("hist");

    TCanvas *c_totalEnergyTrue_before = new TCanvas("c_totalEnergyTrue_before", "c_totalEnergyTrue_before");
    h_totalEnergyTrue_before->Draw("hist");

    TH1D *h_singleElectronEfficiencyTrue = new TH1D("h_singleElectronEfficiencyTrue", "h_singleElectronEfficiencyTrue", 50, 0.0, 5.0);
    for(Int_t i = 1; i <= h_singleEnergyTrue_before->GetNbinsX(); ++ i)
    {
        double after = h_singleEnergyTrue_after->GetBinContent(i);
        double before = h_singleEnergyTrue_before->GetBinContent(i);
        if(before != 0.0)
        {
            h_singleElectronEfficiencyTrue->SetBinContent(i, after / before);
        }
    }

    TH1D *h_totalElectronEfficiencyTrue = new TH1D("h_totalElectronEfficiencyTrue", "h_totalElectronEfficiencyTrue", 50, 0.0, 5.0);
    for(Int_t i = 1; i <= h_totalEnergyTrue_before->GetNbinsX(); ++ i)
    {
        double after = h_totalEnergyTrue_after->GetBinContent(i);
        double before = h_totalEnergyTrue_before->GetBinContent(i);
        if(before != 0.0)
        {
            h_totalElectronEfficiencyTrue->SetBinContent(i, after / before);
        }
    }



    TCanvas *c_singleElectronEfficiencyTrue = new TCanvas("c_singleElectronEfficiencyTrue", "c_singleElectronEfficiencyTrue");
    h_singleElectronEfficiencyTrue->Draw("hist");

    TCanvas *c_totalElectronEfficiencyTrue = new TCanvas("c_totalElectronEfficiencyTrue", "c_totalElectronEfficiencyTrue");
    h_totalElectronEfficiencyTrue->Draw("hist");


    TFile *fout = new TFile("efficiency.root", "recreate");
    h_singleEnergyTrue_after->Write();
    h_singleEnergyTrue_before->Write();
    h_totalEnergyTrue_after->Write();
    h_totalEnergyTrue_before->Write();
    h_singleElectronEfficiencyTrue->Write();
    h_totalElectronEfficiencyTrue->Write();
    fout->Close();

    return 0;
}



void plot()
{
    TFile *fin = new TFile("efficiency.root");
    TH1D *h_singleElectronEfficiency = (TH1D*)fin->Get("h_singleElectronEfficiencyTrue");
    TH1D *h_totalElectronEfficiency = (TH1D*)fin->Get("h_totalElectronEfficiencyTrue");

    ///////////////////////////////////////////////////////////////////////////
    // h_singleElectronEfficiency
    ///////////////////////////////////////////////////////////////////////////

    {
        h_singleElectronEfficiency->SetStats(0);
        h_singleElectronEfficiency->SetTitle(0);
        h_singleElectronEfficiency->SetXTitle("True Single Electron Energy (T_{e}) [MeV]");
        h_singleElectronEfficiency->SetYTitle("Efficiency (#varepsilon)");
        h_singleElectronEfficiency->GetXaxis()->SetLabelFont(43);
        h_singleElectronEfficiency->GetXaxis()->SetLabelSize(15);
        h_singleElectronEfficiency->GetXaxis()->SetTitleFont(43);
        h_singleElectronEfficiency->GetXaxis()->SetTitleSize(20);
        h_singleElectronEfficiency->GetYaxis()->SetLabelFont(43);
        h_singleElectronEfficiency->GetYaxis()->SetLabelSize(15);
        h_singleElectronEfficiency->GetYaxis()->SetTitleFont(43);
        h_singleElectronEfficiency->GetYaxis()->SetTitleSize(20);

        h_singleElectronEfficiency->SetLineWidth(2);
        h_singleElectronEfficiency->SetLineColor(kBlue);
        h_singleElectronEfficiency->SetFillColor(kBlue);
        h_singleElectronEfficiency->SetFillStyle(3002);

        h_singleElectronEfficiency->GetYaxis()->SetRangeUser(0.0, 0.09);

        TCanvas *c = new TCanvas("c", "c");
        c->SetTicks(2, 2);
        h_singleElectronEfficiency->Draw("hist");

        double integral = h_singleElectronEfficiency->Integral();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(43);
        latex.SetTextSize(20);

        TString latexstring_ratio_1;
        TString latexstring_ratio_1b;
        TString latexstring_ratio_2;
        TString latexstring_ratio_2b;
        TString latexstring_ratio_3;
        latexstring_ratio_1.Form("#int");
        latexstring_ratio_1b.Form("dT_{e}");
        latexstring_ratio_2.Form("#int");
        latexstring_ratio_2b.Form("#varepsilon dT_{e}");
        latexstring_ratio_3.Form("= %.2f", integral / 28.0);

        const double latex_x = 0.55;
        const double latex_y = 0.5;

        latex.DrawLatex(latex_x, latex_y, latexstring_ratio_1);
        latex.DrawLatex(latex_x + 0.1, latex_y, latexstring_ratio_1b);
        latex.DrawLatex(latex_x + 0.025, latex_y - 0.03, "#scale[0.7]{0 MeV}");
        latex.DrawLatex(latex_x + 0.025, latex_y + 0.05, "#scale[0.7]{5 MeV}");
    //    latex.DrawLatex(0.3, 0.8, "#int _{0 MeV}^{0.3 MeV} ^{150}Nd(#xi_{31}^{2#nu}=0.7)");
        latex.DrawLatex(latex_x, latex_y + 0.13, latexstring_ratio_2);
        latex.DrawLatex(latex_x + 0.1, latex_y + 0.13, latexstring_ratio_2b);
        latex.DrawLatex(latex_x + 0.025, latex_y + 0.13 - 0.03, "#scale[0.7]{0 MeV}");
        latex.DrawLatex(latex_x + 0.025, latex_y + 0.13 + 0.05, "#scale[0.7]{5 MeV}");
        latex.DrawLatex(latex_x + 0.18, latex_y + 0.742 - 0.67, latexstring_ratio_3);

        TLine *ldiv = new TLine(latex_x, latex_y + 0.08, latex_x + 0.165, latex_y + 0.08);
        ldiv->SetLineWidth(2.0);
        ldiv->SetNDC();
        ldiv->Draw();

        TLegend *leg = new TLegend(0.55, 0.15, 0.85, 0.25);
        leg->AddEntry(h_singleElectronEfficiency, "^{150}Nd Efficiency", "F");
        leg->SetBorderSize(5);
        leg->SetShadowColor(kGray + 2);
        leg->SetTextFont(63);
        leg->SetTextSize(20);
        leg->Draw("BR");

        c->SaveAs("150Nd_singleefficiency.eps");
        c->SaveAs("150Nd_singleefficiency.pdf");
        c->SaveAs("150Nd_singleefficiency.png");
    
    }



    ///////////////////////////////////////////////////////////////////////////
    // h_totalElectronEfficiency
    ///////////////////////////////////////////////////////////////////////////

    {
        h_totalElectronEfficiency->SetStats(0);
        h_totalElectronEfficiency->SetTitle(0);
        h_totalElectronEfficiency->SetXTitle("True Total Electron Energy (T_{e}) [MeV]");
        h_totalElectronEfficiency->SetYTitle("Efficiency (#varepsilon)");
        h_totalElectronEfficiency->GetXaxis()->SetLabelFont(43);
        h_totalElectronEfficiency->GetXaxis()->SetLabelSize(15);
        h_totalElectronEfficiency->GetXaxis()->SetTitleFont(43);
        h_totalElectronEfficiency->GetXaxis()->SetTitleSize(20);
        h_totalElectronEfficiency->GetYaxis()->SetLabelFont(43);
        h_totalElectronEfficiency->GetYaxis()->SetLabelSize(15);
        h_totalElectronEfficiency->GetYaxis()->SetTitleFont(43);
        h_totalElectronEfficiency->GetYaxis()->SetTitleSize(20);

        h_totalElectronEfficiency->SetLineWidth(2);
        h_totalElectronEfficiency->SetLineColor(kBlue);
        h_totalElectronEfficiency->SetFillColor(kBlue);
        h_totalElectronEfficiency->SetFillStyle(3002);

        h_totalElectronEfficiency->GetYaxis()->SetRangeUser(0.0, 0.14);

        TCanvas *c = new TCanvas("c2", "c2");
        c->SetTicks(2, 2);
        h_totalElectronEfficiency->Draw("hist");

        double integral = h_totalElectronEfficiency->Integral();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(43);
        latex.SetTextSize(20);

        TString latexstring_ratio_1;
        TString latexstring_ratio_1b;
        TString latexstring_ratio_2;
        TString latexstring_ratio_2b;
        TString latexstring_ratio_3;
        latexstring_ratio_1.Form("#int");
        latexstring_ratio_1b.Form("dT_{ee}");
        latexstring_ratio_2.Form("#int");
        latexstring_ratio_2b.Form("#varepsilon dT_{e}");
        latexstring_ratio_3.Form("= %.2f", integral / 28.0);

        const double latex_x = 0.15;
        const double latex_y = 0.65;

        latex.DrawLatex(latex_x, latex_y, latexstring_ratio_1);
        latex.DrawLatex(latex_x + 0.1, latex_y, latexstring_ratio_1b);
        latex.DrawLatex(latex_x + 0.025, latex_y - 0.03, "#scale[0.7]{0 MeV}");
        latex.DrawLatex(latex_x + 0.025, latex_y + 0.05, "#scale[0.7]{5 MeV}");
    //    latex.DrawLatex(0.3, 0.8, "#int _{0 MeV}^{0.3 MeV} ^{150}Nd(#xi_{31}^{2#nu}=0.7)");
        latex.DrawLatex(latex_x, latex_y + 0.13, latexstring_ratio_2);
        latex.DrawLatex(latex_x + 0.1, latex_y + 0.13, latexstring_ratio_2b);
        latex.DrawLatex(latex_x + 0.025, latex_y + 0.13 - 0.03, "#scale[0.7]{0 MeV}");
        latex.DrawLatex(latex_x + 0.025, latex_y + 0.13 + 0.05, "#scale[0.7]{5 MeV}");
        latex.DrawLatex(latex_x + 0.18, latex_y + 0.742 - 0.67, latexstring_ratio_3);

        TLine *ldiv = new TLine(latex_x, latex_y + 0.08, latex_x + 0.165, latex_y + 0.08);
        ldiv->SetLineWidth(2.0);
        ldiv->SetNDC();
        ldiv->Draw();

        TLegend *leg = new TLegend(0.55, 0.15, 0.85, 0.25);
        leg->AddEntry(h_totalElectronEfficiency, "^{150}Nd Efficiency", "F");
        leg->SetBorderSize(5);
        leg->SetShadowColor(kGray + 2);
        leg->SetTextFont(63);
        leg->SetTextSize(20);
        leg->Draw("BR");

        c->SaveAs("150Nd_totalefficiency.eps");
        c->SaveAs("150Nd_totalefficiency.pdf");
        c->SaveAs("150Nd_totalefficiency.png");

    }


}


