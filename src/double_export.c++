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
    const double PI=acos(-1.);
    return exp(-q2/delta/delta)*pow(delta,3)/sqrt(PI);
}



int main(void) {
    cout << "double_export" << endl;

    TFile in_file("matr.root", "READ");
    TNtuple *tup = (TNtuple*) in_file.Get("tup");
    double sMin = tup->GetMinimum("s"), sMax = tup->GetMaximum("s");
    int nEv = tup->GetEntries();
    cout << " nEv=" << nEv << " sMin=" << sMin << " sMax=" << sMax << endl;

    TFile out_file("interpolate.root","RECREATE");
    out_file.cd();

    const int nMatr = 26;
    float matr[nMatr];
    const int nTBin = 50;
    TH2D * hMatr[nMatr];
    for (int iH = 0; iH < nMatr; ++iH) {
        string name = "matr_" + i_to_string(iH);
        tup->SetBranchAddress(name.c_str(), &matr[iH]);
        name = "h" + i_to_string(iH);
        hMatr[iH] = new TH2D(name.c_str(), name.c_str(), nTBin, 0, 1, 10, sMin, sMax);
        hMatr[iH]->Sumw2();
    }
    float weight, wf, cosPsi, q2, s;
    tup->SetBranchAddress("wt", &weight);
    tup->SetBranchAddress("cosPsi", &cosPsi);
    tup->SetBranchAddress("s", &s);

    for (int iEv = 0; iEv < nEv; ++iEv) {
        tup->GetEntry(iEv);
        double delta = 0.4;
        double WF = wave_function(q2, delta);
        double nT = (1 + cosPsi) / 2;
        //        if(iEv<10 || nT<0 || nT>1.) cout<<"nT="<<nT<<endl;
        for (int iH = 0; iH < nMatr; ++iH)
            hMatr[iH]->Fill(nT, s, matr[iH] * weight * WF);
    }
    for (int iH = 0; iH < nMatr; ++iH) {
        hMatr[iH]->Scale(1. / hMatr[iH]->GetEntries());
        hMatr[iH]->Write();
    };
    out_file.Save();

    return 0;
}
