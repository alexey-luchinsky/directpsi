#include <iostream>
#include "TFile.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "ramC/Random.h"

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
TH2D * hMatr[nMatr]; // hMatr[nT, s], nT=-(k1-k3)^2/(s-Mcc)^2 in [0,1])
const double Mcc = 3.1;
const int nTBin = 50;

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
    for (int iH = 0; iH < nMatr; ++iH)
        hMatr[iH]->Scale(1. / hMatr[iH]->GetEntries());
}

void saveHST(TH1D *hist, TString name, bool print = false) {
    if (print) cout << " Saving " << name << endl;
    FILE *file = fopen(name.Data(), "w");
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        fprintf(file, "%e %e %e\n", hist->GetBinCenter(i), hist->GetBinContent(i) / hist->GetBinWidth(i), hist->GetBinError(i) / hist->GetBinWidth(i));
    };
    if (print) cout << "\t Histogram sum=" << hist->GetSum() << endl;
    fclose(file);
};

TH1D* genPT(double s) {
    double pTmax = (s - Mcc * Mcc) / (2 * sqrt(s));
    cout << " pTmax=" << pTmax << endl;
    int nBin = 20;
    TH1D* hist = new TH1D("pT2", "pT2", nBin, 0, pTmax * pTmax);
    Random random_generator;
    int nEv = 1e5;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double nT = random_generator.rand(0, 1), t = (Mcc * Mcc - s) * nT, u = Mcc * Mcc - s - t,
                pT2 = t * u / s;
        if (iEv < 10) cout << nT << " " << t << " " << u << " " << pT2 << " " << pT2 / pTmax / pTmax << endl;
        for(int iH=1; iH<nMatr; ++iH) {
            double mtr=hMatr[iH]->Interpolate(nT,s);
            hist->Fill(pT2,pow(mtr,2));
        }
    }

    return hist;
}

int main(void) {
    cout << "partonic" << endl;
    string name = "matr.root";
    cout << "Loading integrals from " << name << endl;
    TFile *in_file = new TFile(name.c_str());
    load_integrals(in_file);

    TH1D *h = genPT(14.9);
    saveHST(h, "hPt2.hst");

    //    TH2D *h20 = hMatr[20];
    //    TAxis *xa = h20->GetXaxis();
    //    cout << h20->GetEntries() << " " << xa->GetBinCenter(1) << " " << xa->GetBinCenter(xa->GetNbins()) << endl;
    //    cout << h20->Interpolate(.1, 13) << endl;
    //    in_file->Close();
    return 0;
}
