#include <iostream>
#include "mtr0.h"
#include "TFile.h"
#include "ramC/rambo2.h"
#include <tclap/CmdLine.h>
#include "TVectorD.h"
using namespace std;

rambo2 ramCM;
EvtVector4R k1, k2, P, k3;
double mc, s, T, U;
double NG, NJ0, NJ1, NJ2;
double matr[nMatr];
TH2D * hMatr[nMatr];

double gen_event(bool debug) {
    ramCM.next();
    //    double wt = 1;
    P = *ramCM.getV(0);
    k3 = *ramCM.getV(1);

    // calculate scalar products
    T = -(k1 - k3).mass2();
    U = -(k2 - k3).mass2();
    mc = Mcc / 2;

    // debug print
    if (debug) {
        //        cout << "wt=" << wt << endl;
        cout << " k1=" << k1 << " k2=" << k2 << endl;
        cout << "P=" << P << " Mpsi=" << P.mass() << " k3=" << k3 << endl;
        cout << " s=" << s << " T=" << T << " U=" << U << " s+t+u=" << s - T - U << endl;
    };
    return 1.;
}

void fill_matr() {
    NG = 1. / sqrt(s * T * U);
    NJ0 = 2. / sqrt(s * T * U);
    NJ1 = 1. / sqrt(s * T * U);
    NJ2 = 1. / Mcc;
    matr[0] = calc_matr0_0();
    matr[1] = calc_matr0_1();
    matr[2] = calc_matr0_2();
    matr[3] = calc_matr0_3();
    matr[4] = calc_matr0_4();
    matr[5] = calc_matr0_5();
    matr[6] = calc_matr0_6();
    matr[7] = calc_matr0_7();
    matr[8] = calc_matr0_8();
    matr[9] = calc_matr0_9();
    matr[10] = calc_matr0_10();
    matr[11] = calc_matr0_11();
    matr[12] = calc_matr0_12();
    matr[13] = calc_matr0_13();
    matr[14] = calc_matr0_14();
    matr[15] = calc_matr0_15();
    matr[16] = calc_matr0_16();
    matr[17] = calc_matr0_17();
    matr[18] = calc_matr0_18();
    matr[19] = calc_matr0_19();
    matr[20] = calc_matr0_20();
    matr[21] = calc_matr0_21();
    matr[22] = calc_matr0_22();
    matr[23] = calc_matr0_23();
    matr[24] = calc_matr0_24();
    matr[25] = calc_matr0_25();

};

double sMin, sMax, alpha;
int nEv;
string output_name;
int nTBin = 10, nsBin = 10;

bool init_commandline_args(int argc, char **argv) {
    try {
        TCLAP::CmdLine cmd("Calculates  matrix elements for gg->psi g process in delta approximation", ' ', "0.9");
        //        TCLAP::ValueArg<float> s_arg("s","s","s",false,10,"float",cmd);
        TCLAP::ValueArg<int> n_arg("n", "n", "log_10(nEv)", false, 4, "float", cmd);
        TCLAP::ValueArg<string> out_arg("o", "o", "name for output root file", false, "matr0.root", "string", cmd);
        TCLAP::ValueArg<float> sMin_arg("m", "sMin", "minimum s", false, -1, "float", cmd);
        TCLAP::ValueArg<float> sMax_arg("M", "sMax", "maximum s", false, 0, "float", cmd);
        TCLAP::ValueArg<float> alpha_arg("a", "alpha", "parameter of smart s sampling", false, 3, "float", cmd);
        TCLAP::ValueArg<int> tb_agr("t", "tb", "number of t bins", false, 50, "int", cmd);
        TCLAP::ValueArg<int> sb_arg("s", "sb", "number of s bins", false, 100, "int", cmd);

        cmd.parse(argc, argv);
        sMin = sMin_arg.getValue();
        if (sMin < Mcc * Mcc) sMin = Mcc * Mcc;
        sMax = sMax_arg.getValue();
        if (sMax < sMin) sMax = 1.1 * sMin;
        alpha = alpha_arg.getValue();
        nTBin = tb_agr.getValue();
        nsBin = sb_arg.getValue();
        if (sMax < sMin) sMax = sMin;
        if (n_arg.getValue() > 10) {
            cout << " log_10(nEv)=" << n_arg.getValue() << " is too lagre!" << endl;
            return false;
        };
        nEv = pow(10, n_arg.getValue());
        output_name = out_arg.getValue();
    } catch (TCLAP::ArgException e) {
        cout << " error " << e.error() << " for arg " << e.argId() << endl;
        return false;
    }
    return true;
}

int main(int argc, char **argv) {
    cout << "delta" << endl;
    if (!init_commandline_args(argc, argv)) return -1;
    cout<<" sMin="<<sMin<<" sMax="<<sMax<<" alpha="<<alpha<<" nEv="<<nEv<<endl;

    TFile out_file(output_name.c_str(), "RECREATE");
    TVectorD dd(4);
    dd[0] = sMin;
    dd[1] = sMax;
    dd[2] = alpha;
    dd[3] = 0;
    dd.Write("stats");
    
    // init histograms
    for (int iH = 0; iH < nMatr; ++iH) {
        string name = "matr_" + i_to_string(iH);
        name = "h" + i_to_string(iH);
        hMatr[iH] = new TH2D(name.c_str(), name.c_str(), nTBin, 0, 1, nsBin, 0, 1);
        hMatr[iH]->Sumw2();
    }


    ramCM.setMass(0, Mcc);
    ramCM.setMass(1, 0.);

    for (int iEv = 0; iEv < nEv; ++iEv) {
        double xs = ramCM.random_generator->rand(0, 1);
        s = sMin + (sMax - sMin) * pow(xs, alpha);
        double ecm = sqrt(s);
        ramCM.setECM(ecm);
        k1.set(ecm / 2, 0, 0, ecm / 2);
        k2.set(ecm / 2, 0, 0, -ecm / 2);
        bool debug = (iEv < 3);
        if (iEv % (nEv / 10) == 0) cout << iEv << endl;
        if (debug) cout << "--------------- Debug print at iEv=" << iEv << "-------------" << endl;
        double weight = gen_event(debug);
        double nT = T / (s - Mcc * Mcc);
        fill_matr();
        for (int iH = 0; iH < nMatr; ++iH)
            hMatr[iH]->Fill(nT, xs, matr[iH]);


    };


    for (int iH = 0; iH < nMatr; ++iH) {
        hMatr[iH]->Scale(1. * nTBin * nsBin / hMatr[iH]->GetEntries());
        hMatr[iH]->Write();
    };
    out_file.Save();



    return 0;
}
