#include <iostream>
#include "TFile.h"
#include "TNtuple.h"
#include "TH2D.h"

using namespace std;

string f_to_string(double v) {
    char c[30];
    sprintf(c, "%f", v);
    return string(c);
}

string i_to_string(int v) {
    char c[30];
    sprintf(c, "%d", v);
    return string(c);
}

double wave_function(double q2, double delta) {
    const double PI = acos(-1.);
    return exp(-q2 / delta / delta) * pow(delta, 3) / sqrt(PI);
}

const int nMatr = 26;
TH2D * hMatr[nMatr];

void load_integrals(TFile *in_file) {
    TNtuple *tup = (TNtuple*) in_file->Get("tup");

    double sMin = tup->GetMinimum("s"), sMax = tup->GetMaximum("s");
    int nEv = tup->GetEntries();
    cout << " nEv=" << nEv << " sMin=" << sMin << " sMax=" << sMax << endl;
    float matr[nMatr];

    for (int iH = 0; iH < nMatr; ++iH) {
        string name = "matr_" + i_to_string(iH);
        tup->SetBranchAddress(name.c_str(), &matr[iH]);
        name = "h" + i_to_string(iH);
        hMatr[iH] = new TH2D(name.c_str(), name.c_str(), 50, -1, 1, 10, sMin, sMax);
        hMatr[iH]->Sumw2();
    }
    float weight, wf, cosPsi, q2, s;
    tup->SetBranchAddress("wt", &weight);
    tup->SetBranchAddress("cosPsi", &cosPsi);
    tup->SetBranchAddress("s", &s);

    for (int i = 0; i < nEv; ++i) {
        tup->GetEntry(i);
        double delta = 0.4;
        double WF = wave_function(q2, delta);
        for (int iH = 0; iH < nMatr; ++iH)
            hMatr[iH]->Fill(cosPsi, s, matr[iH] * weight * WF);
    }
}

int main(void) {
    cout << "partonic" << endl;
    string name = "matr.root";
    cout << "Loading integrals from " << name << endl;
    TFile *in_file = new TFile(name.c_str());
    load_integrals(in_file);

    TH2D *h20 = hMatr[20];
    TAxis *xa = h20->GetXaxis();
    cout << h20->GetEntries() << " " << xa->GetBinCenter(1) << " " << xa->GetBinCenter(xa->GetNbins()) << endl;
    cout << h20->Interpolate(.1, 13) << endl;
    in_file->Close();
    return 0;
}
