#ifndef NEWLOGLIKFITTER_H
#define NEWLOGLIKFITTER_H

#include "TH1.h"
#include "TFractionFitter.h"
#include "TArray.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "THStack.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include <string>
#include "TLegend.h"
#include "TColor.h"
#include "TVector.h"
#include <vector>
#include <map>

//#include "InputFiles.h"

#include "../include/InputNumberDef.h"
#include "../include/InputFileDef.h"
#include "../include/InputNameDef.h"
#include "../include/InputColorDef.h"


const Int_t BACKGROUND_MODE_A = 0; // FREE (constrained) backgrounds. Do not add as systematics.
const Int_t BACKGROUND_MODE_B = 1; // FIXED backgrounds. Add in a systematics
const Int_t BACKGROUND_MODE = BACKGROUND_MODE_B;

const Int_t MODE_PARALLEL = 0;

bool g_mode_fake_data = false;

bool V_ENABLE_SYSALL = false;

bool OPTICAL_CORRECTION_ENABLE = false;


const int N_SYSTEMATICS = 25;
bool V_ENABLE_SYSn[N_SYSTEMATICS] =
{
/* 01 */    false,  // false    // constant 1.0 MeV shift
/* 02 */    false,  // false    // scale factor: m = 1 % + 0.2 %
/* 03 */    true,   // true     // +- 5.55 % efficiency
/* 04 */    true,   // true     // +- 0.50 % enrichment
/* 05 */    false,  // true     // +- 3 keV
/* 06 */    true,   // true     // foil thickness
/* 07 */    true,   // true     // dE/dX
/* 08 */    false,  // false    // brem
/* 09 */    false,  // false    // foil thickness (nominal)
/* 10 */    false,  // false    // dE/dX (nominal)
/* 11 */    true,   // true     // brem (nominal)
/* 12 */    true,   // true     // 0.2 % energy scale factor m = 0.2 %
/* 13 */    true,   // true     // 1.0 % Gaussian smear
/* 14 */    true,   // true     // BKG: 214Bi/214Pb (int)
/* 15 */    true,   // true     // BKG: 207Bi (int)
/* 16 */    true,   // true     // BKG: 228Ac/212Bi/208Tl (int)
/* 17 */    true,   // true     // BKG: 152Eu/154Eu (int)
/* 18 */    true,   // true     // BKG: 40K/234mPa (int)
/* 19 */    true,   // true     // BKG: 214Bi/214Pb (mylar)
/* 20 */    true,   // true     // BKG: 214Bi/214Pb (sfoil, swire)
/* 21 */    true,   // true     // BKG: 214Bi/214Pb (air)
/* 22 */    true,   // true     // BKG: 208Tl (air)
/* 23 */    true,   // true     // BKG: EXTERNAL
/* 24 */    true,   // true     // BKG: NEIGHBOUR
/* 25 */    true,   // true     // optical correction related systematic
};

// TODO: HOW CAN THE MIN POINT HAVE MOVED IF OPTICAL CORRECTION DISABLED?
// answer: min point depends only on shifted pseduodata not V SYS matrix method
// however I thought new phase space scans included optical correction?

bool DRAWSYS_ENABLE_SYSn[N_SYSTEMATICS] =
{
/* 01 */    false,  // false    // constant 1.0 MeV shift
/* 02 */    false,  // false    // scale factor: m = 1 % + 0.2 %
/* 03 */    true,   // true     // +- 5.55 % efficiency
/* 04 */    true,   // true     // +- 0.50 % enrichment
/* 05 */    false,  // true     // +- 3 keV
/* 06 */    true,   // true     // foil thickness
/* 07 */    true,   // true     // dE/dX
/* 08 */    false,  // false    // brem
/* 09 */    false,  // false    // foil thickness (nominal)
/* 10 */    false,  // false    // dE/dX (nominal)
/* 11 */    true,   // true     // brem (nominal)
/* 12 */    true,   // true     // 0.2 % energy scale factor m = 0.2 %
/* 13 */    true,   // true     // 1.0 % Gaussian smear
/* 14 */    true,   // true     // BKG: 214Bi/214Pb (int)
/* 15 */    true,   // true     // BKG: 207Bi (int)
/* 16 */    true,   // true     // BKG: 228Ac/212Bi/208Tl (int)
/* 17 */    true,   // true     // BKG: 152Eu/154Eu (int)
/* 18 */    true,   // true     // BKG: 40K/234mPa (int)
/* 19 */    true,   // true     // BKG: 214Bi/214Pb (mylar)
/* 20 */    true,   // true     // BKG: 214Bi/214Pb (sfoil, swire)
/* 21 */    true,   // true     // BKG: 214Bi/214Pb (air)
/* 22 */    true,   // true     // BKG: 208Tl (air)
/* 23 */    true,   // true     // BKG: EXTERNAL
/* 24 */    true,   // true     // BKG: NEIGHBOUR
/* 25 */    true,   // true     // optical correction related systematic
};



///////////////////////////////////////////////////////////////////////////////
// add random variance to 150Nd MC
//
// this is used to test the SYSALL contours for fake data
// they are much larger in real data - this may be due to variance
// NOTE: removed. issue was due to fitting data with systematics enabled
// but not using different fval in drawing code.


// used to simulate fake_data for min point of (0.7, 1.15) which is close
// to the fit point for real data
double ND150_FAKEDATA_SCALE_FACTOR = 1.0;

// the phase, as string either "1" for Phase 1 or "2" for Phase 2
//Int_t thePhase = -1;
//TString Phase = "BOTH"; // remove these now, fit both phases simultaniously
//const Int_t gMODE_thePhase = -1; // -1=both, 0=P1, 1=P2
//const TString gMODE_Phase = "BOTH"; // BOTH/P1/P2
bool gEnablePhase1;
bool gEnablePhase2;
//std::map<int, int> map_1d_channel_to_phase;
//std::map<int, int> map_2d_channel_to_phase;

parameter_group g_pg;

//bool g_mode_fake_data = false;
//bool g_mode_fake_data = true;

// true = fit to fake data
// false = fit to real data
// fit is done using MinimizeFCNAxialVector
// the exact fit mode / choice of algorithm is set elsewhere
//bool g_mode_fake_data = false;
//bool g_mode_fake_data = false;
// this might change during program execution
// for example, it is necessary to set this to true to fit the fake data
// to measure systematic effects, before switching back to false to fit
// real data for the rest of the algorithm

//std::string g_datetimestamp_string;

// globals required in logLikelihood function but cannot be passed as
// parameters easily
double xi_31_baseline; // remove this? TODO
TH2D *h_nEqNull;
TH2D *h_nEqTwo;
double psiN0;
double psiN2;
double bb_Q;

const double G0_ps_integral_MeV = 0.759721E-45; //0.744684E-45;
const double G0_ps_integral_yrinv = 0.364244E-16; //0.357034E-16;
const double G2_ps_integral_MeV = 0.429791E-45; //0.420438E-45;
const double G2_ps_integral_yrinv = 0.206061E-16; //0.201577E-16;

// TODO: set numberParams dynamically
// using input file
//static const int numberParams = 40; // NOTE: moved above

//const int numberPhases = 2;
//const int number1DHists_perphase = 6;
//const int number2DHists_perphase = 1;
const int number1DHists = 6; //number1DHists_perphase * numberPhases;
const int number2DHists = 1; //number2DHists_perphase * numberPhases;


///////////////////////////////////////////////////////////////////////////////
// channel enable flags

int channel_enable_1D[number1DHists] =
{
0, // ch  0 = hTotalE        (P1&2)
1, // ch  1 = hSingleEnergy  (P1&2)
0, // ch  2 = hHighEnergy    (P1&2)
0, // ch  3 = hLowEnergy     (P1&2)
0, // ch  4 = hEnergySum     (P1&2)
0  // ch  5 = hEnergyDiff    (P1&2)
};

int channel_enable_2D[number2DHists] =
{
0  // ch  0 = hHighLowEnergy    (P1&2)
};


///////////////////////////////////////////////////////////////////////////////
// number of free params
// this is used in the draw code
// must be updated manually
// it is a combination of two things
// in the gParameterGroup, params can be HARD/SOFT constrained or FREE
// they are only free parameters if this mode is FREE
// then some code blocks call Fix() on MnUserParameterState which effectively
// sets the mode to HARD
// i could find a way to do this automatically but it isn't worth the time
// at this late stage
// 2020-11-03
int gNumberFreeParams = 0;



///////////////////////////////////////////////////////////////////////////////
// histogram names

const std::string channel_histname_1D[number1DHists] =
{
    "hTotalE_",
    "hSingleEnergy_",
    "hHighEnergy_",
    "hLowEnergy_",
    "hEnergySum_",
    "hEnergyDiff_"
};

const std::string channel_histname_2D[number2DHists] =
{
    "hHighLowEnergy_"
};


///////////////////////////////////////////////////////////////////////////////
// draw flags

const bool channel_enable_draw_1D[number1DHists] =
{
    true,
    true,
    false,
    false,
    false,
    false
};

const bool channel_enable_draw_2D[number2DHists] =
{
    false
};

// global variable to hold chisquare result
// there is probably a better way to do this using Eval
// but I am too lazy to implement it properly
//Double_t global_chisquare; // TODO: removed, no longer used? 

TObjArray *allMCSamples1D[number1DHists];
TObjArray *allMCSamples2D[number2DHists];
TObjArray *allDataSamples1D;
TObjArray *allDataSamples2D;
TObjArray *allFakeDataSamples1D;
TObjArray *allFakeDataSamples2D;

// stop printing unneccessary debug statements
bool rebuild_fake_data_first_run = true;
bool rebuild_background_MC_first_run = true;

// select between data / fake data in transparent way
//TObjArray *allDataOrFakeDataSamples1D;
//TObjArray *allDataOrFakeDataSamples2D;


std::vector<std::pair<double,double>> ll_walk;
std::vector<std::pair<double,double>> ll_walk_save;


///////////////////////////////////////////////////////////////////////////////
// Fix all background channels (forced all background channels to be fit
// as HARD parameters)
///////////////////////////////////////////////////////////////////////////////

bool FORCE_BKG_HARD = false;



///////////////////////////////////////////////////////////////////////////////
// state machine
///////////////////////////////////////////////////////////////////////////////

bool V_ENABLE_STAT = true; // leave on
// TODO: might need to update code which enables/disables bins
// now that SYS are included

// enable/disable systematics
//bool V_ENABLE_SYSALL = false;

/*
const int N_BKG_SYSTEMATICS = 43;
bool V_ENABLE_SYSBKG = true;

*/
// TODO: min points



// set to true to recalculate V_PHYS_MATHMORE
// ... and maybe some other stuff I got confused about it
bool recalculate_V_PHYS_MATHMORE = true;

// set to true after changing which systematics are enabled
// related to functions inside MinimizeFCNAxialVector::operator()
bool recalculate_V_PHYS_SYS = true;


const bool ENABLE_MIN_POINT_SYSn[N_SYSTEMATICS + 1] =
{
/* 01 */    false,  // false    // constant 1.0 MeV shift
/* 02 */    false,  // false    // scale factor: m = 1 % + 0.2 %
/* 03 */    true,   // true     // +- 5.55 % efficiency
/* 04 */    false,  // false    // +- 0.50 % enrichment
/* 05 */    false,  // false    // +- 3 keV
/* 06 */    true,   // true     // foil thickness
/* 07 */    true,   // true     // dE/dX
/* 08 */    false,  // false    // brem
/* 09 */    false,  // false    // foil thickness (nominal)
/* 10 */    false,  // false    // dE/dX (nominal)
/* 11 */    true,   // true     // brem (nominal)
/* 12 */    true,   // true     // 0.2 % energy scale factor m = 0.2 %
/* 13 */    false,  // true     // 1.0 % Gaussian smear
/* 14 */    false,  // true     // BKG: 214Bi/214Pb (int)
/* 15 */    false,  // true     // BKG: 207Bi (int)
/* 16 */    false,  // true     // BKG: 228Ac/212Bi/208Tl (int)
/* 17 */    false,  // true     // BKG: 152Eu/154Eu (int)
/* 18 */    false,  // true     // BKG: 40K/234mPa (int)
/* 19 */    false,  // true     // BKG: 214Bi/214Pb (mylar)
/* 20 */    false,  // true     // BKG: 214Bi/214Pb (sfoil, swire)
/* 21 */    false,  // true     // BKG: 214Bi/214Pb (air)
/* 22 */    false,  // true     // BKG: 208Tl (air)
/* 23 */    false,  // true     // BKG: EXTERNAL
/* 24 */    false,  // true     // BKG: NEIGHBOUR
/* 25 */    true,   // true     // BKG: ALL
/* 26 */    true    // true     // optical correction related systematic
};
// enable/disable drawing of minimum points for each systematic fit
//const bool ENABLE_MIN_POINT_SYS1 = true; // +- 0.1 MeV
//const bool ENABLE_MIN_POINT_SYS2 = true; // +- 1.2 % scale
//const bool ENABLE_MIN_POINT_SYS3 = true; // +- 5.55 % efficiency
//const bool ENABLE_MIN_POINT_SYS4 = true; // +- 0.50 % enrichment
//const bool ENABLE_MIN_POINT_SYS5 = true; // +- 3 keV





// TODO: get rid of the stack
std::vector<bool> V_ENABLE_SYSALL_stack;
std::vector<bool> V_ENABLE_SYSn_stack[N_SYSTEMATICS]; //TODO anything here to init?
//std::vector<bool> V_ENABLE_SYS1_stack;
//std::vector<bool> V_ENABLE_SYS2_stack;
//std::vector<bool> V_ENABLE_SYS3_stack;
//std::vector<bool> V_ENABLE_SYS4_stack;
//std::vector<bool> V_ENABLE_SYS5_stack;
/*std::vector<bool> V_ENABLE_SYSBKG_stack;*/


///////////////////////////////////////////////////////////////////////////////
// state machine support functions

void check_V_ENABLE_SYS_stack()
{

    std::size_t V_ENABLE_SYSALL_stack_size = V_ENABLE_SYSALL_stack.size();

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        std::size_t stack_size = V_ENABLE_SYSn_stack[i].size();
        if(V_ENABLE_SYSALL_stack_size != stack_size)
        {
            std::cout << "ERROR: " << __func__ << " i=" << i << std::endl;
            std::cout << "V_ENABLE_SYSALL_stack_size=" << V_ENABLE_SYSALL_stack_size << std::endl;
            std::cout << "stack_size=" << stack_size << std::endl;
            throw "V_ENABLE_SYSx_stack size error";
        }
    }

    /*
    std::size_t stack_size = V_ENABLE_SYSBKG_stack.size();
    if(V_ENABLE_SYSALL_stack_size != stack_size)
    {
        std::cout << "ERROR: " << __func__ << std::endl;
        std::cout << "V_ENABLE_SYSALL_stack_size=" << V_ENABLE_SYSALL_stack_size << std::endl;
        std::cout << "stack_size=" << stack_size << std::endl;
        throw "V_ENABLE_SYSBKG_stack size error";
    }
    */

    return;
}


void V_ENABLE_SYS_stack_push()
{
    check_V_ENABLE_SYS_stack();

    V_ENABLE_SYSALL_stack.push_back(V_ENABLE_SYSALL);

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        V_ENABLE_SYSn_stack[i].push_back(V_ENABLE_SYSn[i]);
    }

    //V_ENABLE_SYSBKG_stack.push_back(V_ENABLE_SYSBKG);

    check_V_ENABLE_SYS_stack();
}


void V_ENABLE_SYS_stack_pop()
{
    check_V_ENABLE_SYS_stack();

    if(V_ENABLE_SYSALL_stack.size() > 0)
    {
        V_ENABLE_SYSALL = V_ENABLE_SYSALL_stack.back();
        V_ENABLE_SYSALL_stack.pop_back();

        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            V_ENABLE_SYSn[i] = V_ENABLE_SYSn_stack[i].back();
            V_ENABLE_SYSn_stack[i].pop_back();
        }

        //V_ENABLE_SYSBKG = V_ENABLE_SYSBKG_stack.back();
        //V_ENABLE_SYSBKG_stack.pop_back();
    }
    else
    {
        std::cout <<  "V_ENABLE_SYS_stack_pop ERROR" << std::endl;
        throw "V_ENABLE_SYS_stack_pop ERROR";
    }

    check_V_ENABLE_SYS_stack();
}




///////////////////////////////////////////////////////////////////////////////
// minimization points
///////////////////////////////////////////////////////////////////////////////

// systematic: SYS1
// h = high
// l = low

struct MinPoint
{
    MinPoint()
    {
        A = 0.0;
        A_err = 0.0;
        xi_31 = 0.0;
        xi_31_err = 0.0;
        eff = 0.0;
        eff_err = 0.0;
        fval = 0.0;
    }

    double A;
    double A_err;
    double xi_31;
    double xi_31_err;
    double eff;
    double eff_err;
    double fval;
}
min_point_data,
min_point_fake,
min_point_data_SYSALL,
min_point_fake_SYSALL,
min_point_data_HSD,
min_point_data_HSD_CH0,
min_point_data_SSD;

// min point data, xi free
//double min_point_data[2] = {0.0, 0.0};
//double min_point_data_fval = 0.0;
//double min_point_data_err[2] = {0.0, 0.0};

// min point fake xi free
//double min_point_fake[2] = {0.0, 0.0};
//double min_point_fake_fval = 0.0;
//double min_point_fake_err[2] = {0.0, 0.0};

// min point data with systematics, xi free
//double min_point_data_SYSALL[2] = {0.0, 0.0};
//double min_point_data_SYSALL_fval = 0.0;
//double min_point_data_SYSALL_err[2] = {0.0, 0.0};

// min point fake with systematics, xi free
//double min_point_fake_SYSALL[2] = {0.0, 0.0};
//double min_point_fake_SYSALL_fval = 0.0;
//double min_point_fake_SYSALL_err[2] = {0.0, 0.0};

// min point data, xi=HSD, total electron energy
//double min_point_data_HSD_CH0[2] = {0.0, 0.0};
//double min_point_data_HSD_CH0_fval = 0.0;
//double min_point_data_HSD_CH0_err[2] = {0.0, 0.0};

// min point data, xi=HSD
//double min_point_data_HSD[2] = {0.0, 0.0};
//double min_point_data_HSD_fval = 0.0;
//double min_point_data_HSD_err[2] = {0.0, 0.0};

// min point data, xi=SSD
//double min_point_data_SSD[2] = {0.0, 0.0};
//double min_point_data_SSD_fval = 0.0;
//double min_point_data_SSD_err[2] = {0.0, 0.0};

// +- 0.1 MeV
// +- 1.2 % scale
// +- 5.55 % efficiency
// +- 0.50 % enrichment
// +- 3 keV
// foil thickness (V)
// dEdX (V)
// brem (V)
// foil thickness (N)
// dEdX (N)
// brem (N)
// +- 0.2 % scale
// 1.0 % Gaussian smear
// BKG: 214Bi/214Pb (int)
// BKG: 207Bi (int)
// BKG: 228Ac/212Bi/208Tl (int)
// BKG: 152Eu/154Eu (int)
// BKG: 40K/234mPa (int)
// BKG: 214Bi/214Pb (mylar)
// BKG: 214Bi/214Pb (sfoil, swire)
// BKG: 214Bi/214Pb (air)
// BKG: 208Tl (air)
// BKG: EXTERNAL
// BKG: NEIGHBOUR

MinPoint min_point_fake_sysn_h[N_SYSTEMATICS + 1];
MinPoint min_point_fake_sysn_l[N_SYSTEMATICS + 1];

/*
double min_point_fake_sysn_h[N_SYSTEMATICS + 1][2] =
{
    {0.0, 0.0},
    {0.0, 0.0}, // 2
    {0.0, 0.0},
    {0.0, 0.0}, // 4
    {0.0, 0.0},
    {0.0, 0.0}, // 6
    {0.0, 0.0},
    {0.0, 0.0}, // 8
    {0.0, 0.0},
    {0.0, 0.0}, // 10
    {0.0, 0.0},
    {0.0, 0.0}, // 12
    {0.0, 0.0},
    {0.0, 0.0}, // 14
    {0.0, 0.0},
    {0.0, 0.0}, // 16
    {0.0, 0.0},
    {0.0, 0.0}, // 18
    {0.0, 0.0},
    {0.0, 0.0}, // 20
    {0.0, 0.0},
    {0.0, 0.0}, // 22
    {0.0, 0.0},
    {0.0, 0.0}, // 24
    {0.0, 0.0},
    {0.0, 0.0}  // 26
};
double min_point_fake_sysn_h_err[N_SYSTEMATICS + 1][2] =
{
    {0.0, 0.0},
    {0.0, 0.0}, // 2
    {0.0, 0.0},
    {0.0, 0.0}, // 4
    {0.0, 0.0},
    {0.0, 0.0}, // 6
    {0.0, 0.0},
    {0.0, 0.0}, // 8
    {0.0, 0.0},
    {0.0, 0.0}, // 10
    {0.0, 0.0},
    {0.0, 0.0}, // 12
    {0.0, 0.0},
    {0.0, 0.0}, // 14
    {0.0, 0.0},
    {0.0, 0.0}, // 16
    {0.0, 0.0},
    {0.0, 0.0}, // 18
    {0.0, 0.0},
    {0.0, 0.0}, // 20
    {0.0, 0.0},
    {0.0, 0.0}, // 22
    {0.0, 0.0},
    {0.0, 0.0}, // 24
    {0.0, 0.0},
    {0.0, 0.0}  // 26
};

double min_point_fake_sysn_h_fval[N_SYSTEMATICS + 1] =
{
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, // 8
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, // 16
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, // 24
    0.0, 0.0
};

double min_point_fake_sysn_l[N_SYSTEMATICS + 1][2] =
{
    {0.0, 0.0},
    {0.0, 0.0}, // 2
    {0.0, 0.0},
    {0.0, 0.0}, // 4
    {0.0, 0.0},
    {0.0, 0.0}, // 6
    {0.0, 0.0},
    {0.0, 0.0}, // 8
    {0.0, 0.0},
    {0.0, 0.0}, // 10
    {0.0, 0.0},
    {0.0, 0.0}, // 12
    {0.0, 0.0},
    {0.0, 0.0}, // 14
    {0.0, 0.0},
    {0.0, 0.0}, // 16
    {0.0, 0.0},
    {0.0, 0.0}, // 18
    {0.0, 0.0},
    {0.0, 0.0}, // 20
    {0.0, 0.0},
    {0.0, 0.0}, // 22
    {0.0, 0.0},
    {0.0, 0.0}, // 24
    {0.0, 0.0},
    {0.0, 0.0}  // 26
};
double min_point_fake_sysn_l_err[N_SYSTEMATICS + 1][2] =
{
    {0.0, 0.0},
    {0.0, 0.0}, // 2
    {0.0, 0.0},
    {0.0, 0.0}, // 4
    {0.0, 0.0},
    {0.0, 0.0}, // 6
    {0.0, 0.0},
    {0.0, 0.0}, // 8
    {0.0, 0.0},
    {0.0, 0.0}, // 10
    {0.0, 0.0},
    {0.0, 0.0}, // 12
    {0.0, 0.0},
    {0.0, 0.0}, // 14
    {0.0, 0.0},
    {0.0, 0.0}, // 16
    {0.0, 0.0},
    {0.0, 0.0}, // 18
    {0.0, 0.0},
    {0.0, 0.0}, // 20
    {0.0, 0.0},
    {0.0, 0.0}, // 22
    {0.0, 0.0},
    {0.0, 0.0}, // 24
    {0.0, 0.0},
    {0.0, 0.0}  // 26
};

double min_point_fake_sysn_l_fval[N_SYSTEMATICS + 1] =
{
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, // 8
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, // 16
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, // 24
    0.0, 0.0
};
*/


void min_point_save(
    const std::string &min_point_fname,
    MinPoint &min_point)
{
    std::ofstream ofs(min_point_fname);
    ofs << min_point.xi_31 << " " << min_point.A << std::endl;
    ofs << min_point.fval << std::endl;
    ofs << min_point.xi_31_err << " " << min_point.A_err << std::endl;
    ofs << min_point.eff << " " << min_point.eff_err << std::endl;

    const double N_ISOTOPE = 150.0;
    const double N_AVOGADRO = 6.022e+23;
    const double MASS = 36.55;

    const double ACTIVITY0 = 3.45e-4;
    const double ACTIVITY = ACTIVITY0 * min_point.A;
    const double ACTIVITY_ERR = ACTIVITY0 * min_point.A_err;

    const double YEARTOSEC = 1.0 / 31557600.0;
    const double N_AVOGADRO_YEARTOSEC = N_AVOGADRO * YEARTOSEC;

    const double T12 = std::log(2.0) * (MASS  / N_ISOTOPE) * (N_AVOGADRO_YEARTOSEC / ACTIVITY);
    const double T12_ERR = T12 * ACTIVITY_ERR / ACTIVITY; // all positive

    ofs << "T_{\\frac{1}{2}} = " << T12 << " +- " << T12_ERR << std::endl;
    ofs << "#xi_{31} = " << min_point.xi_31 << " +- " << min_point.xi_31_err << std::endl;
    ofs << "#chi^{2} = " << min_point.fval << std::endl;
    ofs << "#varepsilon = " << min_point.eff << " +- " << min_point.eff_err << std::endl;

    ofs.close();
}



int min_point_load(
    const std::string &min_point_fname,
    MinPoint &min_point)
{
    int ret = 0;
    std::ifstream ifs(min_point_fname);
    if(ifs.is_open())
    {
        ifs >> min_point.xi_31 >> min_point.A;
        ifs >> min_point.fval;
        ifs >> min_point.xi_31_err >> min_point.A_err;
        ifs >> min_point.eff >> min_point.eff_err;
        ifs.close();
    }
    else
    {
        ret = -1;
    }
    return ret;
}
    

/*
void min_point_save(
    const std::string &min_point_fname,
    double *const min_point,
    double *const min_point_err,
    const double &min_point_fval)
{
    std::ofstream ofs(min_point_fname);
    ofs << min_point[0] << " " << min_point[1] << std::endl;
    ofs << min_point_fval << std::endl;
    ofs << min_point_err[0] << " " << min_point_err[1] << std::endl;

    const double N_ISOTOPE = 150.0;
    const double N_AVOGADRO = 6.022e+23;
    const double MASS = 36.55;

    const double ACTIVITY0 = 3.45e-4;
    const double ACTIVITY = ACTIVITY0 * min_point[1];
    const double ACTIVITY_ERR = ACTIVITY0 * min_point_err[1];

    const double YEARTOSEC = 1.0 / 31557600.0;
    const double N_AVOGADRO_YEARTOSEC = N_AVOGADRO * YEARTOSEC;

    const double T12 = std::log(2.0) * (MASS  / N_ISOTOPE) * (N_AVOGADRO_YEARTOSEC / ACTIVITY);
    const double T12_ERR = T12 * ACTIVITY_ERR / ACTIVITY; // all positive

    ofs << "T_{\\frac{1}{2}} = " << T12 << " +- " << T12_ERR << std::endl;
    ofs << "#xi_{31} = " << min_point[0] << " +- " << min_point_err[0] << std::endl;
    ofs << "#chi^{2} = " << min_point_fval << std::endl;

    ofs.close();
}


int min_point_load(
    const std::string &min_point_fname,
    double *const min_point,
    double *const min_point_err,
    double &min_point_fval)
{
    int ret = 0;
    std::ifstream ifs(min_point_fname);
    if(ifs.is_open())
    {
        ifs >> min_point[0] >> min_point[1];
        ifs >> min_point_fval;
        ifs >> min_point_err[0] >> min_point_err[1];
        ifs.close();
    }
    else
    {
        ret = -1;
    }
    return ret;
}
*/

///////////////////////////////////////////////////////////////////////////////
// systematic objects - Phase 1
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// nominal

std::vector<double> *systematic_nominal_1D_P1[number1DHists] =
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *systematic_nominal_2D_P1[number2DHists] =
{
    nullptr
};

std::vector<double> *systematic_n_low_1D_P1[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 4
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 8
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 12
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 16
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 20
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 24
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *systematic_n_low_2D_P1[N_SYSTEMATICS][number2DHists] =
{
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 4
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 8
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 12
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 16
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 20
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 24
    {nullptr}
};

std::vector<double> *systematic_n_high_1D_P1[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 4
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 8
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 12
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 16
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 20
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 24
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *systematic_n_high_2D_P1[N_SYSTEMATICS][number2DHists] =
{
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 4
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 8
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 12
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 16
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 20
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 24
    {nullptr}
};

/* not implemented
std::vector<double> *sys_V_MATRIX_coeff_1D_P1[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};
*/

std::vector<double> *systematic_n_V_MATRIX_coeff_1D_P1[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 4
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 8
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 12
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 16
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 20
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 24
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

/* not implemented
std::vector<double> *sys_V_MATRIX_coeff_2D_P1[number2DHists] =
{
    nullptr
};
*/

std::vector<double> *systematic_n_V_MATRIX_coeff_2D_P1[N_SYSTEMATICS][number2DHists] =
{
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 4
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 8
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 12
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 16
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 20
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 24
    {nullptr}
};


///////////////////////////////////////////////////////////////////////////////
// systematic objects - Phase 2
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// nominal

std::vector<double> *systematic_nominal_1D_P2[number1DHists] =
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *systematic_nominal_2D_P2[number2DHists] =
{
    nullptr
};

std::vector<double> *systematic_n_low_1D_P2[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 4
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 8
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 12
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 16
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 20
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 24
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *systematic_n_low_2D_P2[N_SYSTEMATICS][number2DHists] =
{
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 4
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 8
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 12
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 16
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 20
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 24
    {nullptr}
};

std::vector<double> *systematic_n_high_1D_P2[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 4
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 8
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 12
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 16
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 20
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 24
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *systematic_n_high_2D_P2[N_SYSTEMATICS][number2DHists] =
{
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 4
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 8
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 12
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 16
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 20
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 24
    {nullptr}
};

/* not implemented
std::vector<double> *sys_V_MATRIX_coeff_1D_P2[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};
*/

std::vector<double> *systematic_n_V_MATRIX_coeff_1D_P2[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 4
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 8
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 12
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 16
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 20
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 24
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

// TODO
/* not implemented
std::vector<double> *sys_V_MATRIX_coeff_2D_P2[number2DHists] =
{
    nullptr
};
*/

std::vector<double> *systematic_n_V_MATRIX_coeff_2D_P2[N_SYSTEMATICS][number2DHists] =
{
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 4
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 8
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 12
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 16
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 20
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 24
    {nullptr}
};




///////////////////////////////////////////////////////////////////////////////
// chi2 objects for V MATRIX method, Phase 1
///////////////////////////////////////////////////////////////////////////////

TMatrixD *V_PHYS_1D_P1_MATHMORE[number1DHists] =
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

TMatrixD *V_PHYS_2D_P1_MATHMORE[number2DHists] =
{
    nullptr
};

std::vector<bool> *V_ENABLE_BIN_1D_P1[number1DHists] =
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<bool> *V_ENABLE_BIN_2D_P1[number2DHists] =
{
    nullptr
};

// not used
/*
std::vector<double> *V_PHYS_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};
*/

// changes depending on M
std::vector<double> *V_PHYS_STAT_1D_P1_data[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *V_PHYS_STAT_2D_P1_data[number2DHists] = 
{
    nullptr
};

std::vector<double> *V_PHYS_SYSn_1D_P1_data[N_SYSTEMATICS][number1DHists] =
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 4
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 8
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 12
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 16
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 20
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 24
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *V_PHYS_SYSn_2D_P1_data[N_SYSTEMATICS][number2DHists] =
{
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 4
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 8
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 12
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 16
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 20
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 24
    {nullptr}
};


std::vector<double> *V_PHYS_SYSALL_1D_P1_data[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *V_PHYS_SYSALL_2D_P1_data[number2DHists] = 
{
    nullptr
};

std::vector<double> *D_1D_P1_data[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *D_2D_P1_data[number2DHists] = 
{
    nullptr
};

std::vector<double> *D_minus_M_1D_P1_data[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *D_minus_M_2D_P1_data[number2DHists] = 
{
    nullptr
};

std::vector<double> *M_1D_P1_data[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *M_2D_P1_data[number2DHists] = 
{
    nullptr
};

///////////////////////////////////////////////////////////////////////////////
// chi2 objects for V MATRIX method, Phase 2
///////////////////////////////////////////////////////////////////////////////

TMatrixD *V_PHYS_1D_P2_MATHMORE[number1DHists] =
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

TMatrixD *V_PHYS_2D_P2_MATHMORE[number2DHists] =
{
    nullptr
};

std::vector<bool> *V_ENABLE_BIN_1D_P2[number1DHists] =
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<bool> *V_ENABLE_BIN_2D_P2[number2DHists] =
{
    nullptr
};

/*
std::vector<double> *V_PHYS_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};
*/

std::vector<double> *V_PHYS_STAT_1D_P2_data[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *V_PHYS_STAT_2D_P2_data[number2DHists] = 
{
    nullptr
};

std::vector<double> *V_PHYS_SYSn_1D_P2_data[N_SYSTEMATICS][number1DHists] =
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 4
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 8
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 12
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 16
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 20
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}, // 24
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *V_PHYS_SYSn_2D_P2_data[N_SYSTEMATICS][number2DHists] =
{
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 4
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 8
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 12
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 16
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 20
    {nullptr},
    {nullptr},
    {nullptr},
    {nullptr}, // 24
    {nullptr}
};

std::vector<double> *V_PHYS_SYSALL_1D_P2_data[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *V_PHYS_SYSALL_2D_P2_data[number2DHists] = 
{
    nullptr
};


std::vector<double> *D_1D_P2_data[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *D_2D_P2_data[number2DHists] = 
{
    nullptr
};

std::vector<double> *D_minus_M_1D_P2_data[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *D_minus_M_2D_P2_data[number2DHists] = 
{
    nullptr
};

std::vector<double> *M_1D_P2_data[number1DHists] = 
{
    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
};

std::vector<double> *M_2D_P2_data[number2DHists] = 
{
    nullptr
};



void
draw_V_PHYS_init_helper(TH2D* histo, const TString &latexstring, const TString& xtitle, const TString& ytitle)
{
    //histo->SetTitle(0);
    histo->SetTitle(latexstring);
    histo->SetStats(0);
    histo->SetContour(1000);

    histo->GetXaxis()->SetTitle(xtitle);
    histo->GetYaxis()->SetTitle(ytitle);

    histo->GetZaxis()->SetLabelOffset(0.005);
    histo->GetXaxis()->SetLabelSize(17.0);
    histo->GetXaxis()->SetLabelFont(43);
    histo->GetYaxis()->SetLabelSize(17.0);
    histo->GetYaxis()->SetLabelFont(43);
    histo->GetZaxis()->SetLabelSize(17.0);
    histo->GetZaxis()->SetLabelFont(43);
    histo->GetXaxis()->SetTitleSize(18.0);
    histo->GetXaxis()->SetTitleFont(43);
    histo->GetYaxis()->SetTitleSize(18.0);
    histo->GetYaxis()->SetTitleFont(43);
    histo->GetXaxis()->SetTitleOffset(1.5);
    histo->GetYaxis()->SetTitleOffset(1.2);
    histo->GetXaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetLabelOffset(0.01);

}


void 
draw_V_PHYS_canvas_helper(TCanvas *canvas, TH2D* histo, Color_t textcolor)
{
    canvas->cd();
    canvas->SetTicks(2, 2);
    canvas->SetRightMargin(0.15);
    canvas->SetBottomMargin(0.15);

    // TODO: try something like
    // https://root-forum.cern.ch/t/drawing-histograms-with-colz-how-to-move-the-z-axis-color-legend/41915

    canvas->Update();
    TPaletteAxis *palette = (TPaletteAxis*)histo->GetListOfFunctions()->FindObject("palette");
    if(palette != nullptr)
    {
        palette->SetX1NDC(0.88 + 0.03);
        palette->SetX2NDC(0.92 + 0.03);
        palette->SetY1NDC(0.15);
        palette->SetY2NDC(0.9);
        palette->Draw();
    }
    gPad->Modified();
    gPad->Update();
    canvas->Modified();

    TString latexstring = histo->GetTitle();
    histo->SetTitle(0);

    // latex top right label
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(43);
    latex.SetTextSize(20);
    latex.SetTextColor(textcolor);
    latex.DrawLatex(0.2, 0.8, latexstring);
}



bool DRAW_V_PHYS_ = false;

void draw_V_PHYS_(
    const TString &dir
    )
{

    const int channel = 1;

    ///////////////////////////////////////////////////////////////////////////
    // Phase 1 - initialize pointers: set to null

    // stat
    TH2D *draw_V_PHYS_STAT_P1 = nullptr;

    // sys N
    TH2D *draw_V_PHYS_SYSn_P1[N_SYSTEMATICS];
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        draw_V_PHYS_SYSn_P1[i] = nullptr;
    }

    // sys all
    TH2D *draw_V_PHYS_SYSALL_P1 = nullptr;
    
    // PHYS
    TH2D *draw_V_PHYS_P1 = nullptr;
    TH2D *draw_V_PHYSINV_P1 = nullptr;
    TH2D *draw_V_DmMPHYSINVDmM_P1 = nullptr;
    
    ///////////////////////////////////////////////////////////////////////////
    // Phase 2 - initialize pointers: set to null

    // stat
    TH2D *draw_V_PHYS_STAT_P2 = nullptr;

    // sys N
    TH2D *draw_V_PHYS_SYSn_P2[N_SYSTEMATICS]; 
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        draw_V_PHYS_SYSn_P2[i] = nullptr;
    }

    // sys all
    TH2D *draw_V_PHYS_SYSALL_P2 = nullptr;

    // PHYS
    TH2D *draw_V_PHYS_P2 = nullptr;
    TH2D *draw_V_PHYSINV_P2 = nullptr;
    TH2D *draw_V_DmMPHYSINVDmM_P2 = nullptr;


    ///////////////////////////////////////////////////////////////////////////
    // initialize pointers: allocate memory

    ///////////////////////////////////////////////////////
    // PHASE 1
    ///////////////////////////////////////////////////////

    draw_V_PHYS_STAT_P1 =
        new TH2D("draw_V_PHYS_STAT_P1", "draw_V_PHYS_STAT_P1",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYS_STAT_P1,
        "V^{STAT} (Phase 1)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        TString name;
        name.Form("draw_V_PHYS_SYS%d_P1", i);
        draw_V_PHYS_SYSn_P1[i] = new TH2D(name, name, 50, 0.0, 5.0, 50, 0.0, 5.0);
        TString title;
        title.Form("V^{SYS%d} (Phase 1)", i);
        draw_V_PHYS_init_helper(
            draw_V_PHYS_SYSn_P1[i], title,
            "V_{ij} bin i Energy (MeV)",
            "V_{ij} bin j Energy (MeV)");
    }

    draw_V_PHYS_SYSALL_P1
        = new TH2D("draw_V_PHYS_SYSALL_P1", "draw_V_PHYS_SYSALL_P1",
                   50, 0.0, 5.0,
                   50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYS_SYSALL_P1,
        "V^{SYS} (Phase 1)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    draw_V_PHYS_P1 =
        new TH2D("draw_V_PHYS_P1", "draw_V_PHYS_P1",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYS_P1,
        "V (Phase 1)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    draw_V_PHYSINV_P1 = 
        new TH2D("draw_V_PHYSINV_P1", "draw_V_PHYSINV_P1",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYSINV_P1,
        "V^{-1} (Phase 1)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");
    
    draw_V_DmMPHYSINVDmM_P1 = 
        new TH2D("draw_V_DmMPHYSINVDmM_P1", "draw_V_DmMPHYSINVDmM_P1",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_DmMPHYSINVDmM_P1,
        "V^{-1}_{RES^{2}} (Phase 1)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");
    

    ///////////////////////////////////////////////////////
    // PHASE 2
    ///////////////////////////////////////////////////////
    
    draw_V_PHYS_STAT_P2 = new TH2D(
        "draw_V_PHYS_STAT_P2", "draw_V_PHYS_STAT_P2",
        50, 0.0, 5.0,
        50, 0.0, 5.0);
    draw_V_PHYS_init_helper(draw_V_PHYS_STAT_P2,
        "V^{STAT} (Phase 2)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        TString name;
        name.Form("draw_V_PHYS_SYS%d_P2", i);
        draw_V_PHYS_SYSn_P2[i] =
            new TH2D(name, name, 50, 0.0, 5.0, 50, 0.0, 5.0);

        TString title;
        title.Form("V^{SYS%d} (Phase 2)", i);
        draw_V_PHYS_init_helper(
            draw_V_PHYS_SYSn_P2[i],
            title,
            "V_{ij} bin i Energy (MeV)",
            "V_{ij} bin j Energy (MeV)");
    }

    draw_V_PHYS_SYSALL_P2 =
        new TH2D("draw_V_PHYS_SYSALL_P2", "draw_V_PHYS_SYSALL_P2",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYS_SYSALL_P2,
        "V^{SYS}_{ij} (Phase 2)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    draw_V_PHYS_P2 =
        new TH2D("draw_V_PHYS_P2", "draw_V_PHYS_P2",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYS_P2,
        "V (Phase 2)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    draw_V_PHYSINV_P2 =
        new TH2D("draw_V_PHYSINV_P2", "draw_V_PHYSINV_P2",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYSINV_P2,
        "V^{-1} (Phase 2)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    draw_V_DmMPHYSINVDmM_P2 =
        new TH2D("draw_V_DmMPHYSINVDmM_P2", "draw_V_DmMPHYSINVDmM_P2",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_DmMPHYSINVDmM_P2,
        "V^{-1}_{RES^{2}} (Phase 2)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    //std::cout << "ALLOC " << __func__ << std::endl;
                        

    ///////////////////////////////////////////////////////////////////////////
    // set contents of TH2D objects


    double fval_fake = 0.0;

    ///////////////////////////////////////////////////////////
    // Phase 1

    int j_counter = 0;
    for(Int_t j = 0; j < 50; ++ j)
    {
        if(V_ENABLE_BIN_1D_P1[channel]->at(j) == true)
        {
            // do nothing
        }
        else
        {
            continue;
        }

        int i_counter = 0;
        for(Int_t i = 0; i < 50; ++ i)
        {
            if(V_ENABLE_BIN_1D_P1[channel]->at(i) == true)
            {
                // do nothing
            }
            else
            {
                continue;
            }

            double cstat = 0.0;
            double csys = 0.0;
            double sysn[N_SYSTEMATICS];
            double content_MATHMORE = 0.0;
            double D_minus_M_content_1 = 0.0;
            double D_minus_M_content_2 = 0.0;

            if(V_ENABLE_STAT == true)
            {
                if(i == j)
                {
                    cstat = V_PHYS_STAT_1D_P1_data[channel]->at(i + j * 50);
                }
            }
            if(V_ENABLE_SYSALL == true)
            {
                csys = V_PHYS_SYSALL_1D_P1_data[channel]->at(i + j * 50);
            }

            for(int n = 0; n < N_SYSTEMATICS; ++ n)
            {
                sysn[n] = V_PHYS_SYSn_1D_P1_data[n][channel]->at(i + j * 50);
            }
            
            // assumes inversion done - which it should be if det != 0.0
            D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->at(i);
            content_MATHMORE = V_PHYS_1D_P1_MATHMORE[channel]->operator[](j_counter).operator[](i_counter);
            D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->at(j);

            // set

            draw_V_PHYS_STAT_P1->SetBinContent(i + 1, j + 1, cstat);

            for(int n = 0; n < N_SYSTEMATICS; ++ n)
            {
                draw_V_PHYS_SYSn_P1[n]->SetBinContent(i + 1, j + 1, sysn[n]);
            }

            draw_V_PHYS_SYSALL_P1->SetBinContent(i + 1, j + 1, csys);
            draw_V_PHYS_P1->SetBinContent(i + 1, j + 1, cstat + csys);
            draw_V_PHYSINV_P1->SetBinContent(i + 1, j + 1, content_MATHMORE);
            draw_V_DmMPHYSINVDmM_P1->SetBinContent(i + 1, j + 1, D_minus_M_content_1 * content_MATHMORE * D_minus_M_content_2);

            fval_fake += D_minus_M_content_1 * content_MATHMORE * D_minus_M_content_2;

            ++ i_counter;
        }

        ++ j_counter;
    }

    ///////////////////////////////////////////////////////////
    // Phase 2

    j_counter = 0;
    for(Int_t j = 0; j < 50; ++ j)
    {
        if(V_ENABLE_BIN_1D_P2[channel]->at(j) == true)
        {
            // do nothing
        }
        else
        {
            continue;
        }

        int i_counter = 0;
        for(Int_t i = 0; i < 50; ++ i)
        {
            if(V_ENABLE_BIN_1D_P2[channel]->at(i) == true)
            {
                // do nothing
            }
            else
            {
                continue;
            }

            double cstat = 0.0;
            double csys = 0.0;
            double sysn[N_SYSTEMATICS];
            double content_MATHMORE = 0.0;
            double D_minus_M_content_1 = 0.0;
            double D_minus_M_content_2 = 0.0;

            if(V_ENABLE_STAT == true)
            {
                if(i == j)
                {
                    cstat = V_PHYS_STAT_1D_P2_data[channel]->at(i + j * 50);
                }
            }
            if(V_ENABLE_SYSALL == true)
            {
                csys = V_PHYS_SYSALL_1D_P2_data[channel]->at(i + j * 50);
            }

            for(int n = 0; n < N_SYSTEMATICS; ++ n)
            {
                sysn[n] = V_PHYS_SYSn_1D_P2_data[n][channel]->at(i + j * 50);
                csys += sysn[n];
            }
            
            // assumes inversion done - which it should be if det != 0.0
            D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->at(i);
            content_MATHMORE = V_PHYS_1D_P2_MATHMORE[channel]->operator[](j_counter).operator[](i_counter);
            D_minus_M_content_2 = D_minus_M_1D_P2_data[channel]->at(j);
 
            // set

            draw_V_PHYS_STAT_P2->SetBinContent(i + 1, j + 1, cstat);
            for(int n = 0; n < N_SYSTEMATICS; ++ n)
            {
                draw_V_PHYS_SYSn_P2[n]->SetBinContent(i + 1, j + 1, sysn[n]);
            }

            draw_V_PHYS_SYSALL_P2->SetBinContent(i + 1, j + 1, csys);
            draw_V_PHYS_P2->SetBinContent(i + 1, j + 1, cstat + csys);
            draw_V_PHYSINV_P2->SetBinContent(i + 1, j + 1, content_MATHMORE);
            draw_V_DmMPHYSINVDmM_P2->SetBinContent(i + 1, j + 1, D_minus_M_content_1 * content_MATHMORE * D_minus_M_content_2);

            fval_fake += D_minus_M_content_1 * content_MATHMORE * D_minus_M_content_2;

            ++ i_counter;
        }

        ++ j_counter;
    }





    ///////////////////////////////////////////////////////////////////////////
    // draw canvas

    ///////////////////////////////////////////////////////
    // PHASE 1
    ///////////////////////////////////////////////////////

    TCanvas *canvas_V_PHYS_STAT_P1 = new TCanvas("canvas_V_PHYS_STAT_P1", "canvas_V_PHYS_STAT_P1");
    draw_V_PHYS_STAT_P1->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_STAT_P1, draw_V_PHYS_STAT_P1, kBlack);
    canvas_V_PHYS_STAT_P1->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_STAT_P1.png");

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        TString name;
        name.Form("canvas_V_PHYS_SYS%d_P1", i);
        TString oname;
        oname = name + ".png";
        TCanvas *canvas_V_PHYS_SYSn_P1 = new TCanvas(name, name);
        draw_V_PHYS_SYSn_P1[i]->Draw("colz");
        if(i == 0 || i == 4)
        {
            draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSn_P1, draw_V_PHYS_SYSn_P1[i], kWhite);
        }
        else
        {
            draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSn_P1, draw_V_PHYS_SYSn_P1[i], kBlack);
        }
        canvas_V_PHYS_SYSn_P1->SaveAs(TString("./") + dir + "/" + oname);
    }

    TCanvas *canvas_V_PHYS_SYSALL_P1 = new TCanvas("canvas_V_PHYS_SYSALL_P1", "canvas_V_PHYS_SYSALL_P1");
    draw_V_PHYS_SYSALL_P1->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSALL_P1, draw_V_PHYS_SYSALL_P1, kWhite);
    canvas_V_PHYS_SYSALL_P1->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_SYSALL_P1.png");

    TCanvas *canvas_V_PHYS_P1 = new TCanvas("canvas_V_PHYS_P1", "canvas_V_PHYS_P1");
    draw_V_PHYS_P1->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_P1, draw_V_PHYS_P1, kWhite);
    canvas_V_PHYS_P1->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_P1.png");

    TCanvas *canvas_V_PHYSINV_P1 = new TCanvas("canvas_V_PHYSINV_P1", "canvas_V_PHYSINV_P1");
    draw_V_PHYSINV_P1->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYSINV_P1, draw_V_PHYSINV_P1, kWhite);
    canvas_V_PHYSINV_P1->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYSINV_P1.png");

    TCanvas *canvas_V_DmMPHYSINVDmM_P1 = new TCanvas("canvas_V_DmMPHYSINVDmM_P1", "canvas_V_DmMPHYSINVDmM_P1");
    draw_V_DmMPHYSINVDmM_P1->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_DmMPHYSINVDmM_P1, draw_V_DmMPHYSINVDmM_P1, kWhite);
    canvas_V_DmMPHYSINVDmM_P1->SaveAs(TString("./") + dir + "/" + "canvas_V_DmMPHYSINVDmM_P1.png");

    ///////////////////////////////////////////////////////
    // PHASE 2
    ///////////////////////////////////////////////////////

    TCanvas *canvas_V_PHYS_STAT_P2 = new TCanvas("canvas_V_PHYS_STAT_P2", "canvas_V_PHYS_STAT_P2");
    draw_V_PHYS_STAT_P2->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_STAT_P2, draw_V_PHYS_STAT_P2, kBlack);
    canvas_V_PHYS_STAT_P2->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_STAT_P2.png");

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        TString name;
        name.Form("canvas_V_PHYS_SYS%d_P2", i);
        TString oname;
        oname = name + ".png";
        TCanvas *canvas_V_PHYS_SYSn_P2 = new TCanvas(name, name);
        draw_V_PHYS_SYSn_P2[i]->Draw("colz");
        if(i == 0 || i == 4)
        {
            draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSn_P2, draw_V_PHYS_SYSn_P2[i], kWhite);
        }
        else
        {
            draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSn_P2, draw_V_PHYS_SYSn_P2[i], kBlack);
        }
        canvas_V_PHYS_SYSn_P2->SaveAs(TString("./") + dir + "/" + oname);
    }

    TCanvas *canvas_V_PHYS_SYSALL_P2 = new TCanvas("canvas_V_PHYS_SYSALL_P2", "canvas_V_PHYS_SYSALL_P2");
    draw_V_PHYS_SYSALL_P2->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSALL_P2, draw_V_PHYS_SYSALL_P2, kWhite);
    canvas_V_PHYS_SYSALL_P2->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_SYSALL_P2.png");

    TCanvas *canvas_V_PHYS_P2 = new TCanvas("canvas_V_PHYS_P2", "canvas_V_PHYS_P2");
    draw_V_PHYS_P2->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_P2, draw_V_PHYS_P2, kWhite);
    canvas_V_PHYS_P2->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_P2.png");

    TCanvas *canvas_V_PHYSINV_P2 = new TCanvas("canvas_V_PHYSINV_P2", "canvas_V_PHYSINV_P2");
    draw_V_PHYSINV_P2->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYSINV_P2, draw_V_PHYSINV_P2, kWhite);
    canvas_V_PHYSINV_P2->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYSINV_P2.png");

    TCanvas *canvas_V_DmMPHYSINVDmM_P2 = new TCanvas("canvas_V_DmMPHYSINVDmM_P2", "canvas_V_DmMPHYSINVDmM_P2");
    draw_V_DmMPHYSINVDmM_P2->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_DmMPHYSINVDmM_P2, draw_V_DmMPHYSINVDmM_P2, kWhite);
    canvas_V_DmMPHYSINVDmM_P2->SaveAs(TString("./") + dir + "/" + "canvas_V_DmMPHYSINVDmM_P2.png");

    std::cout << "printed all canvas" << std::endl;
}

#endif
