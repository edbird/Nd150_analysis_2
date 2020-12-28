#ifndef SYSTEMATICS_H
#define SYSTEMATICS_H



class Systematics
{


    public:


    Systematics()
    {

        reset();

        aux_data_is_loaded = false;

        h_systematic_foil_thickness_virtual_h = nullptr;
        h_systematic_foil_thickness_virtual_l = nullptr;
        h_systematic_dEdX_virtual_h = nullptr;
        h_systematic_dEdX_virtual_l = nullptr;
        h_systematic_brem_virtual_h = nullptr;
        h_systematic_brem_virtual_l = nullptr;
        h_systematic_foil_thickness_nominal_h = nullptr;
        h_systematic_foil_thickness_nominal_l = nullptr;
        h_systematic_dEdX_nominal_h = nullptr;
        h_systematic_dEdX_nominal_l = nullptr;
        h_systematic_brem_nominal_h = nullptr;
        h_systematic_brem_nominal_l = nullptr;

        load_aux_data();
    }



    void reset()
    { 
        systematic_energy_offset = 0.0;
        systematic_energy_scale = 0.0;
        systematic_efficiency = 0.0;
        systematic_enrichment = 0.0;
        systematic_energy_offsetsmall = 0.0;
        
        systematic_foil_thickness_virtual = 0.0;
        systematic_dEdX_virtual = 0.0;
        systematic_brem_virtual = 0.0;
        
        systematic_foil_thickness_nominal = 0.0;
        systematic_dEdX_nominal = 0.0;
        systematic_brem_nominal = 0.0;

        systematic_energy_scale_small = 0.0;
        systematic_energy_gaussian_smear = 0.0;

        systematic_bkg_bi214_int = 0.0;
        systematic_bkg_bi207_int = 0.0;
        systematic_bkg_tl208_int = 0.0;
        systematic_bkg_eu152_int = 0.0;
        systematic_bkg_k40_int = 0.0;
        systematic_bkg_mylar = 0.0;
        systematic_bkg_bi214_sfoil_swire = 0.0;
        systematic_bkg_bi214_air = 0.0;
        systematic_bkg_tl208_air = 0.0;
        systematic_bkg_external = 0.0;
        systematic_bkg_neighbour = 0.0;

        systematic_optical_correction = 0.0;
    }


    void load_aux_data()
    {
        if(aux_data_is_loaded == false)
        {
            {
                TFile *f1 = new TFile("systematic_foil_thickness_virtual.root");
                h_systematic_foil_thickness_virtual_h = (TH2D*)f1->Get("h_systematic_foil_thickness_h");
                h_systematic_foil_thickness_virtual_l = (TH2D*)f1->Get("h_systematic_foil_thickness_l");
                if(h_systematic_foil_thickness_virtual_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_foil_thickness_virtual.root" << std::endl;
                }
                if(h_systematic_foil_thickness_virtual_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_foil_thickness_virtual.root" << std::endl;
                }
            }

            {
                TFile *f2 = new TFile("systematic_dEdX_virtual.root");
                h_systematic_dEdX_virtual_h = (TH2D*)f2->Get("h_systematic_dEdX_h");
                h_systematic_dEdX_virtual_l = (TH2D*)f2->Get("h_systematic_dEdX_l");
                if(h_systematic_dEdX_virtual_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_dEdX_virtual.root" << std::endl;
                }
                if(h_systematic_dEdX_virtual_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_dEdX_virtual.root" << std::endl;
                }
            }

            {
                TFile *f3 = new TFile("systematic_brem_virtual.root");
                h_systematic_brem_virtual_h = (TH2D*)f3->Get("h_systematic_brem_h");
                h_systematic_brem_virtual_l = (TH2D*)f3->Get("h_systematic_brem_l");
                if(h_systematic_brem_virtual_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_brem_virtual.root" << std::endl;
                }
                if(h_systematic_brem_virtual_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_brem_virtual.root" << std::endl;
                }
            }

            {
                TFile *f4 = new TFile("systematic_foil_thickness_nominal.root");
                h_systematic_foil_thickness_nominal_h = (TH2D*)f4->Get("h_systematic_foil_thickness_h");
                h_systematic_foil_thickness_nominal_l = (TH2D*)f4->Get("h_systematic_foil_thickness_l");
                if(h_systematic_foil_thickness_nominal_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_foil_thickness_nominal.root" << std::endl;
                }
                if(h_systematic_foil_thickness_nominal_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_foil_thickness_nominal.root" << std::endl;
                }
            }

            {
                TFile *f5 = new TFile("systematic_dEdX_nominal.root");
                h_systematic_dEdX_nominal_h = (TH2D*)f5->Get("h_systematic_dEdX_h");
                h_systematic_dEdX_nominal_l = (TH2D*)f5->Get("h_systematic_dEdX_l");
                if(h_systematic_dEdX_nominal_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_dEdX_nominal.root" << std::endl;
                }
                if(h_systematic_dEdX_nominal_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_dEdX_nominal.root" << std::endl;
                }
            }

            {
                TFile *f6 = new TFile("systematic_brem_nominal.root");
                h_systematic_brem_nominal_h = (TH2D*)f6->Get("h_systematic_brem_h");
                h_systematic_brem_nominal_l = (TH2D*)f6->Get("h_systematic_brem_l");
                if(h_systematic_brem_nominal_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_brem_nominal.root" << std::endl;
                }
                if(h_systematic_brem_nominal_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_brem_nominal.root" << std::endl;
                }
            }

            // n3-nl-model
            {
                std::ifstream ifs("n3-nl-model.txt");
                if(ifs.is_open())
                {
                    std::vector<double> energy;
                    std::vector<double> correction_nominal;
                    std::vector<double> correction_high;
                    std::vector<double> correction_low;
                    
                    double e, cl, cn, ch;
                    while(ifs.good())
                    {
                        ifs >> e >> cl >> cn >> ch;
                        energy.push_back(e);
                        correction_low.push_back(cl);
                        correction_nominal.push_back(cn);
                        correction_high.push_back(ch);
                    }

                    ifs.close();

                    int np = energy.size();
                    g_systematic_optical_correction_l = new TGraph(np);
                    g_systematic_optical_correction_n = new TGraph(np);
                    g_systematic_optical_correction_h = new TGraph(np);

                    for(int i = 0; i < np; ++ i)
                    {
                        double x = energy.at(i);
                        double cn = correction_nominal.at(i);
                        double cl = correction_low.at(i);
                        double ch = correction_high.at(i);
                        g_systematic_optical_correction_l->SetPoint(i, x, cl);
                        g_systematic_optical_correction_n->SetPoint(i, x, cn);
                        g_systematic_optical_correction_h->SetPoint(i, x, ch);
                    }

                    // debug
                    new TCanvas;
                    g_systematic_optical_correction_n->SetMinimum(1.0 - 0.05);
                    g_systematic_optical_correction_n->SetMaximum(1.0 + 0.05);
                    g_systematic_optical_correction_n->Draw("al");
                    g_systematic_optical_correction_l->Draw("al");
                    g_systematic_optical_correction_h->Draw("al");
                }
                else
                {
                    std::cout << "ERROR: Could not open file n3-nl-model.txt" << std::endl;
                    throw "n3-nl-model";
                }
            }

            aux_data_is_loaded = true;
        }
        
    }


    // initial systematics
    double systematic_energy_offset;                // 0.1 MeV energy offset
    double systematic_energy_scale;                 // 1.2 % energy scale
    double systematic_efficiency;                   // 5.55 % efficiency
    double systematic_enrichment;                   // 0.5 % enrichment
    double systematic_energy_offsetsmall;           // 0.3 keV energy offset
    double systematic_foil_thickness_virtual;       // FT V
    double systematic_dEdX_virtual;                 // DE V
    double systematic_brem_virtual;                 // B V
    double systematic_foil_thickness_nominal;       // FT N
    double systematic_dEdX_nominal;                 // DE N
    double systematic_brem_nominal;                 // B N

    double systematic_energy_scale_small;           // 0.2 % energy scale
    double systematic_energy_gaussian_smear;        // 1.0 % Gaussian smear

    // background systematics
    double systematic_bkg_bi214_int;                // 214Bi / 214Pb INT
    double systematic_bkg_bi207_int;                // 207 Bi INT
    double systematic_bkg_tl208_int;                // 208 Tl INT
    double systematic_bkg_eu152_int;                // 152/154 Eu INT
    double systematic_bkg_k40_int;                  // 40 K INT
    double systematic_bkg_mylar;                    // 214Bi / 214Pb mylar
    double systematic_bkg_bi214_sfoil_swire;        // 214 Bi / 214 Pb sfoil/swire
    double systematic_bkg_bi214_air;                // 214 Bi AIR
    double systematic_bkg_tl208_air;                // 208 Tl AIR
    double systematic_bkg_external;                 // external (all)
    double systematic_bkg_neighbour;                // neighbour (all)

    // optical correction systematic
    double systematic_optical_correction;           // optical correction systematic



    bool aux_data_is_loaded;

    TH2D *h_systematic_foil_thickness_virtual_h;
    TH2D *h_systematic_foil_thickness_virtual_l;
    TH2D *h_systematic_dEdX_virtual_h;
    TH2D *h_systematic_dEdX_virtual_l;
    TH2D *h_systematic_brem_virtual_h;
    TH2D *h_systematic_brem_virtual_l;
    TH2D *h_systematic_foil_thickness_nominal_h;
    TH2D *h_systematic_foil_thickness_nominal_l;
    TH2D *h_systematic_dEdX_nominal_h;
    TH2D *h_systematic_dEdX_nominal_l;
    TH2D *h_systematic_brem_nominal_h;
    TH2D *h_systematic_brem_nominal_l;

    TGraph *g_systematic_optical_correction_l;
    TGraph *g_systematic_optical_correction_n;
    TGraph *g_systematic_optical_correction_h;

}gSystematics;



void gSystematics_print()
{
        std::cout << "SYSTEMATICS: CONSTANT OFFSET: " << gSystematics.systematic_energy_offset << " MeV" << std::endl;
        std::cout << "SYSTEMATICS: CONSTANT SCALE: " << gSystematics.systematic_energy_scale << " %" << std::endl;
        std::cout << "SYSTEMATICS: EFFICIENCY: " << gSystematics.systematic_efficiency << " %" << std::endl;
        std::cout << "SYSTEMATICS: ENRICHMENT: " << gSystematics.systematic_enrichment << " %" << std::endl;
        std::cout << "SYSTEMATICS: OFFSETSMALL: " << gSystematics.systematic_energy_offsetsmall << " MeV" << std::endl;

        std::cout << "SYSTEMATICS: FOIL THICKNESS (V): " << gSystematics.systematic_foil_thickness_virtual << " " << std::endl;
        std::cout << "SYSTEMATICS: ENERGY LOSS (V): " << gSystematics.systematic_dEdX_virtual << " " << std::endl;
        std::cout << "SYSTEMATICS: BREMSSTRAHLUNG (V): " << gSystematics.systematic_brem_virtual << " " << std::endl;
        
        std::cout << "SYSTEMATICS: FOIL THICKNESS (N): " << gSystematics.systematic_foil_thickness_nominal << " " << std::endl;
        std::cout << "SYSTEMATICS: ENERGY LOSS (N): " << gSystematics.systematic_dEdX_nominal << " " << std::endl;
        std::cout << "SYSTEMATICS: BREMSSTRAHLUNG (N): " << gSystematics.systematic_brem_nominal << " " << std::endl;

        std::cout << "SYSTEMATICS: CONSTANT SCALE SMALL: " << gSystematics.systematic_energy_scale_small << " %" << std::endl;
        std::cout << "SYSTEMATICS: ENERGY GAUSSIAN SMEAR: " << gSystematics.systematic_energy_gaussian_smear << " %" << std::endl;

        std::cout << "SYSTEMATICS: BKG 214Bi/214Pb INT: " << gSystematics.systematic_bkg_bi214_int << " " << std::endl;
        std::cout << "SYSTEMATICS: BKG 207Pb INT: " << gSystematics.systematic_bkg_bi207_int << " " << std::endl;
        std::cout << "SYSTEMATICS: BKG 228Ac/212Bi/208Tl INT: " << gSystematics.systematic_bkg_tl208_int << " " << std::endl;
        std::cout << "SYSTEMATICS: BKG 152Eu/154Eu INT: " << gSystematics.systematic_bkg_eu152_int << " " << std::endl;
        std::cout << "SYSTEMATICS: BKG 40k/234mPa INT: " << gSystematics.systematic_bkg_k40_int << " " << std::endl;
        std::cout << "SYSTEMATICS: BKG 214Bi/214Pb MYLAR: " << gSystematics.systematic_bkg_mylar << " " << std::endl;
        std::cout << "SYSTEMATICS: BKG 214Bi/214Pb SFOIL/SWIRE: " << gSystematics.systematic_bkg_bi214_sfoil_swire << " " << std::endl;
        std::cout << "SYSTEMATICS: BKG 214Bi/214Pb AIR: " << gSystematics.systematic_bkg_bi214_air << " " << std::endl;
        std::cout << "SYSTEMATICS: BKG 208Tl AIR: " << gSystematics.systematic_bkg_tl208_air << " " << std::endl;
        std::cout << "SYSTEMATICS: BKG EXTERNAL: " << gSystematics.systematic_bkg_external << " " << std::endl;
        std::cout << "SYSTEMATICS: BKG NEIGHBOUR FOIL: " << gSystematics.systematic_bkg_neighbour << " " << std::endl;
        
        std::cout << "SYSTEMATICS: OPTICAL CORRECTION ??? : " << gSystematics.systematic_optical_correction << " ??? " << std::endl;

}


#endif // SYSTEMATICS_H
