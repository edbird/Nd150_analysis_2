void testPalette()
{
    TCanvas *c1 = new TCanvas("c1","test");
    gStyle->SetPalette(1,0);
    TH2D *h = new TH2D("h","h",100,0,100,100,0,100);
    h->Fill(20,20,20);
    h->Draw("colz");
    c1->Update();
    TPaletteAxis *palette =
    (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.2);
    palette->SetX2NDC(0.12);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.5);
    c1->Modified();
}
