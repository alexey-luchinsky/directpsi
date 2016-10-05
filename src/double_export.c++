#include <iostream>
#include "TFile.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "ramC/Random.h"
#include <tclap/CmdLine.h>
#include "utils.h"

using namespace std;
using namespace TCLAP;



string in_fileName, out_fileName;
int nTBin, nsBin;
double delta;

bool init_command_args(int argc, char **argv) {
    try {
        TCLAP::CmdLine cmd("Create (nT,s) tables from obtained using calcInts.exe program data", ' ', "0.9");
        TCLAP::ValueArg<string> in_agr("i", "in", "input file name", false, "matr.root", "string", cmd);
        TCLAP::ValueArg<string> out_arg("o", "out", "output file name", false, "interpolation.root", "string", cmd);
        TCLAP::ValueArg<int> tb_agr("t", "tb", "number of t bins", false, 50, "int", cmd);
        TCLAP::ValueArg<int> sb_arg("s", "sb", "number of s bins", false, 100, "int", cmd);
        TCLAP::ValueArg<float> delta_arg("d", "delta", "wave function parameter", false, 0.4, "float", cmd);
        cmd.parse(argc, argv);
        in_fileName = in_agr.getValue();
        out_fileName = out_arg.getValue();
        nTBin = tb_agr.getValue();
        nsBin = sb_arg.getValue();
        delta = delta_arg.getValue();

        return true;
    } catch (TCLAP::ArgException e) {
        cout << " error " << e.error() << " for arg " << e.argId() << endl;
        return false;
    };
}

int main(int argc, char **argv) {
    if (!init_command_args(argc, argv))
        return -1;

    // determine normalization coefficient
    double sum = 0;
    Random random_generator;
    for (int i = 0; i < 1e6; ++i) {
        double Q = random_generator.rand(0, Mcc / 2);
        double PI = acos(-1.);
        double mc = sqrt(Mcc * Mcc / 4 - Q * Q);
        double mtr = (Q * Q + 3 * pow(Mcc + 2 * mc, 2)) / (Mcc + 2 * mc) / pow(mc, 3. / 2);
        sum += 4 * PI * Q * Q * wave_function(Q*Q, delta) * mtr * Mcc / 2;
    };
    double Gee = 92.9e-6 * 5.791e-2;
    double ec = 2. / 3;
    double alphaQED = 1. / 137;
    sum = sum / 1e6;
    double norm = sqrt(12*PI*Gee/Mcc)/(4*PI*alphaQED*ec*sum/Mcc);
    cout << " normalization coefficient : " << norm << endl;

    TFile in_file(in_fileName.c_str(), "READ");
    if (!in_file.IsOpen()) {
        cout << " Cannot open file " << in_fileName << endl;
        return -1;
    };
    TVectorD stats = *((TVectorD*) in_file.Get("stats"));


    TNtuple *tup = (TNtuple*) in_file.Get("tup");
    double sMin = stats(0), sMax = stats(1), alpha = stats(2);
    int nEv = tup->GetEntries();
    cout << " nEv=" << nEv << " sMin=" << sMin << " sMax=" << sMax << " alpha=" << alpha << " delta=" << delta << endl;

    TFile out_file(out_fileName.c_str(), "RECREATE");
    out_file.cd();
    TVectorD dd(4);
    dd[0] = sMin;
    dd[1] = sMax;
    dd[2] = alpha;
    dd[3] = delta;
    dd.Write("stats");


    float matr[nMatr];
    TH2D * hMatr[nMatr];
    for (int iH = 0; iH < nMatr; ++iH) {
        string name = "matr_" + i_to_string(iH);
        tup->SetBranchAddress(name.c_str(), &matr[iH]);
        name = "h" + i_to_string(iH);
        hMatr[iH] = new TH2D(name.c_str(), name.c_str(), nTBin, 0, 1, nsBin, 0, 1);
        hMatr[iH]->Sumw2();
    }

    float weight, wf, cosPsi, q2, xs;
    tup->SetBranchAddress("wt", &weight);
    tup->SetBranchAddress("cosPsi", &cosPsi);
    tup->SetBranchAddress("xs", &xs);
    tup->SetBranchAddress("q2", &q2);

    for (int iEv = 0; iEv < nEv; ++iEv) {
        tup->GetEntry(iEv);

        double wt = 4 * PI * q2 * Mcc / 2; // d^3q = 
        double WF = wave_function(q2, delta);
        double nT = (1 + cosPsi) / 2;
        if (iEv < 10) cout << "matr0=" << matr[0] << " q2=" << q2 << " wt=" << wt << " WF=" << WF << " nT=" << nT << " xs=" << xs << endl;
        for (int iH = 0; iH < nMatr; ++iH) {
            if(iH>0) matr[iH] *= norm;
            hMatr[iH]->Fill(nT, xs, matr[iH] * wt * WF);
        };
    };
    for (int iH = 0; iH < nMatr; ++iH) {
        hMatr[iH]->Scale(1. * nTBin * nsBin / hMatr[iH]->GetEntries());
        hMatr[iH]->Write();
    };
    out_file.Save();

    return 0;
}
