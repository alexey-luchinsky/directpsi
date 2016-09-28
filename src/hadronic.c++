#include <iostream>
#include <vector>
#include "ramC/Random.h"
#include "LHAPDF/LHAPDF.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH2D.h"
#include <tclap/CmdLine.h>

using namespace std;
using namespace LHAPDF;

const int nMatr = 26;
TH2D * hMatr[nMatr]; // hMatr[nT, s], nT=-(k1-k3)^2/(s-Mcc)^2 in [0,1])
TFile *hist_file;

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


void saveHST(TH1D *hist, TString name, bool print = false) {
    if (print) cout << " Saving " << name << endl;
    FILE *file = fopen(name.Data(), "w");
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        fprintf(file, "%e %e %e\n", hist->GetBinCenter(i), hist->GetBinContent(i) / hist->GetBinWidth(i), hist->GetBinError(i) / hist->GetBinWidth(i));
    };
    if (print) cout << "\t Histogram sum=" << hist->GetSum() << endl;
    fclose(file);
};

string in_fileName, out_fileName;
double S;
int nEv;

bool load_integrals() {
    hist_file = new TFile(in_fileName.c_str(), "READ");
    if(!hist_file->IsOpen()) return false;
    for (int iH = 0; iH < nMatr; ++iH) {
        string name = "h" + i_to_string(iH);
        hMatr[iH] = (TH2D*) hist_file->Get(name.c_str());
    }
    return true;
}

bool init_commandline_args(int argc, char **argv) {
    try {
        TCLAP::CmdLine cmd("Calculate hadronic process using double_export.exe program data", ' ', "0.9");
        TCLAP::ValueArg<string> in_agr("i", "in", "input file name", false, "interpolation.root", "string", cmd);
        TCLAP::ValueArg<string> out_arg("o", "out", "output file name", false, "PP.root", "string", cmd);
        TCLAP::ValueArg<float> S_arg("s", "s", "squared energy of hadronic reaction", false, 500, "float", cmd);
        TCLAP::ValueArg<float> n_arg("n","n","log_10(nEv)",false,6,"float",cmd);
        cmd.parse(argc, argv);
        in_fileName = in_agr.getValue();
        out_fileName = out_arg.getValue();
        S = S_arg.getValue();
        nEv=pow(10,n_arg.getValue());
        if(nEv>1e10) {
            cout<<" nEv="<<nEv<<" is too large!"<<endl;
            return false;
        }
        return true;
    } catch (TCLAP::ArgException e) {
        cout << " error " << e.error() << " for arg " << e.argId() << endl;
        return false;
    }
}

void save_pdf(const LHAPDF::PDF *pdf, double q2) {
    ofstream pdf_file;
    pdf_file.open("pdf.txt");
    double dx=1e-4;
    for(double x=dx; x<1.-dx; x+=dx) {
        pdf_file << x << " " << pdf->xfxQ2(0, x, q2) / x<< endl;
    }
    pdf_file.close();
}

int main(int argc, char **argv) {
    init_commandline_args(argc, argv);
    const string setname = "CT10";
    const int imem = 0;
    const PDF *pdf = mkPDF(setname);

    double Mcc = 3.1, Mcc2 = Mcc*Mcc, q2 = Mcc2;
    save_pdf(pdf, q2);

    if(!load_integrals()) {
        cout<<"Cannot open file "<<in_fileName<<endl;
        return -1;
    }

    Random random_generator;
    TFile out_file(out_fileName.c_str(), "RECREATE");
    TNtuple tup("tup", "tup", "hatS:pT2:nT:x1:x2:y:mtr2:pdf1:pdf2:wt");

    double sMin = hMatr[0]->GetYaxis()->GetBinLowEdge(1);
    if (sMin < Mcc2) {
        cout << " root file sMin=" << sMin << " lower than Mcc2=" << Mcc2 << ". Setting sMin=Mcc2" << endl;
        sMin = Mcc2;
    };
    double sMax = hMatr[0]->GetYaxis()->GetBinUpEdge(hMatr[0]->GetYaxis()->GetNbins());
    if (sMax > S) {
        cout << " sMax=" << sMax << " higher than S=" << S << " setting sMax=S" << endl;
        sMax = S;
    };
    cout << "S=" << S << " sMin=" << sMin << " sMax=" << sMax << endl;
    
    TH1D *hh=new TH1D("hh","hh",30, sMin, sMax);
    hh->Sumw2();
    
    for (int iEv = 0; iEv < nEv; ++iEv) {
        bool debug = (iEv < 0);
        if (debug) cout << "----- Debug print at i=" << iEv << "---------" << endl;
        double wt = 1;
        double s = random_generator.rand(sMin, sMax);
        if (debug) cout << "\t s=" << s << endl;

        wt *= (sMax - sMin) / S; // hat s
        // y
        double yMax = log(S / s) / 2;
        double y = random_generator.rand(-yMax, yMax);
        if (debug) cout << "\t y=" << y << endl;
        wt *= 2 * yMax; // y
        
        double x1 = sqrt(s / S) * exp(y), pdf1 = pdf->xfxQ2(0, x1, q2) / x1;
        if (debug) cout << "\t x1=" << x1 << " pdf1=" << pdf1 << endl;
        double x2 = sqrt(s / S) * exp(-y), pdf2 = pdf->xfxQ2(0, x2, q2) / x2;

        double nT = random_generator.rand(0, 1), t = (Mcc2 - s) * nT, u = Mcc2 - s - t,
                pT2 = t * u / s;
        wt *= (s-Mcc2);

        double mtr2 = 0, mtr;
        if (debug)
            cout << " nT=" << nT << " s=" << s << " wt=" << wt << endl;
        for (int iH = 1; iH < nMatr; ++iH) {
            mtr = hMatr[iH]->Interpolate(nT, s);
            mtr2 += pow(mtr, 2);
        };
        //    TNtuple tup("tup","tup","hatS:pT2:nT:x1:x2:y:mtr2:pdf1:pdf2:wt");
        tup.Fill(s, pT2, nT, x1, x2, y, mtr2, pdf1, pdf2, wt);
        
        // conversion from matr2 to dsdt
        double PI=acos(-1.);
        wt *= 1./(64*PI*s)*4/s;
        int nTbin=hMatr[0]->GetXaxis()->GetNbins(), nsbin=hMatr[0]->GetYaxis()->GetNbins();
        mtr = hMatr[0]->Interpolate(nT,s)*nTbin*nsbin;
        mtr2=pow(mtr,2);
        if(iEv<10)
            cout<<"mtr2="<<mtr2<<endl;
        hh->Fill(s,pdf1*pdf2*wt);
    };
    tup.Write();
    out_file.Save();
    hh->Scale(1./nEv);
    saveHST(hh, "dSigmaDs_1.hst");
    
    return 0;
}
