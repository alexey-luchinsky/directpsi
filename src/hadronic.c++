// This is me
#include <iostream>
#include <vector>
#include "ramC/Random.h"
#include "LHAPDF/LHAPDF.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TVectorD.h"
#include <tclap/CmdLine.h>
#include "utils.h"

using namespace std;
using namespace LHAPDF;

TH2D * hMatr[nMatr]; // hMatr[nT, xs], nT=-(k1-k3)^2/(s-Mcc)^2 in [0,1])
TFile *hist_file;




string in_fileName, out_fileName;
double delta, S;
double sMin, sMax, alpha;
int nEv;

bool load_integrals() {
    hist_file = new TFile(in_fileName.c_str(), "READ");
    if (!hist_file->IsOpen()) return false;

    // read statistics
    TVectorD stats = *((TVectorD*) hist_file->Get("stats"));
    sMin = stats(0);
    sMax = stats(1);
    alpha = stats(2);
    delta = stats(3);

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
        TCLAP::ValueArg<float> n_arg("n", "n", "log_10(nEv)", false, 6, "float", cmd);
        cmd.parse(argc, argv);
        in_fileName = in_agr.getValue();
        out_fileName = out_arg.getValue();
        S = S_arg.getValue();
        nEv = pow(10, n_arg.getValue());
        if (nEv > 1e10) {
            cout << " nEv=" << nEv << " is too large!" << endl;
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
    double dx = 1e-4;
    for (double x = dx; x < 1. - dx; x += dx) {
        pdf_file << x << " " << pdf->xfxQ2(0, x, q2) / x << endl;
    }
    pdf_file.close();
}

int main(int argc, char **argv) {
    init_commandline_args(argc, argv);
    const string setname = "CT10";
    const int imem = 0;
    const PDF *pdf = mkPDF(setname);

    double Mcc = 3.1, Mcc2 = Mcc*Mcc, scale2 = Mcc2;
    save_pdf(pdf, scale2);

    if (!load_integrals()) {
        cout << "Cannot open file " << in_fileName << endl;
        return -1;
    }
    cout << " sMin=" << sMin << " sMax=" << sMax << " alpha=" << alpha << " delta=" << delta << endl;


    Random random_generator;
    TFile out_file(out_fileName.c_str(), "RECREATE");
    TNtuple tup("tup", "tup", "hatS:pT2:xF:nT:x1:x2:y:mtr2:mtr20:pdf1:pdf2:wt");

    //    double sMin = hMatr[0]->GetYaxis()->GetBinLowEdge(1);
    if (sMin < Mcc2) {
        cout << " root file sMin=" << sMin << " lower than Mcc2=" << Mcc2 << ". Setting sMin=Mcc2" << endl;
        sMin = Mcc2;
    };
    //    double sMax = hMatr[0]->GetYaxis()->GetBinUpEdge(hMatr[0]->GetYaxis()->GetNbins());
    if (sMax > S) {
        cout << " sMax=" << sMax << " higher than S=" << S << " setting sMax=S" << endl;
        sMax = S;
    };
    cout << "S=" << S << " sMin=" << sMin << " sMax=" << sMax << endl;
    cout << "delta=" << delta << endl;
    cout << " alpha=" << alpha << endl;

    TH1D *h0 = new TH1D("h0", "h0", 30, sMin, sMax);
    h0->Sumw2();
    TH1D *hAll = new TH1D("hAll", "hAll", 30, sqrt(sMin), sqrt(sMax));
    hAll->Sumw2();

    for (int iEv = 0; iEv < nEv; ++iEv) {
        bool debug = (iEv < 0);
        if (debug) cout << "----- Debug print at i=" << iEv << "---------" << endl;
        double wt = 1;
        double xs = random_generator.rand(0, 1);
        double s = sMin + (sMax - sMin) * pow(xs, alpha);
        wt *= alpha * (sMax - sMin) * pow(xs, alpha - 1) / S;

        if (debug) cout << "\t s=" << s << endl;

        // y
        double yMax = log(S / s) / 2;
        double y = random_generator.rand(-yMax, yMax);
        if (debug) cout << "\t y=" << y << endl;
        wt *= 2 * yMax; // y


        double x1 = sqrt(s / S) * exp(y), pdf1 = pdf->xfxQ2(0, x1, scale2) / x1;
        if (debug) cout << "\t x1=" << x1 << " pdf1=" << pdf1 << endl;
        double x2 = sqrt(s / S) * exp(-y), pdf2 = pdf->xfxQ2(0, x2, scale2) / x2;

        double nT = random_generator.rand(0, 1), t = (Mcc2 - s) * nT, u = Mcc2 - s - t,
                pT2 = t * u / s;
        double cosPsi=2*nT-1;
        double xF = (s+Mcc2)/(2*s)*(x1-x2)+(s-Mcc2)/(2*s)*(x1+x2)*cosPsi;
        wt *= (s - Mcc2);
        // conversion to dsdt
        wt *= 1. / (64 * PI * s)*4 / s;
        // symmetry, etc
        wt *= 1./(2*2*8*2*8);
        // tranfer to nb
        wt *= 0.389e6;

        double mtr2 = 0, mtr;
        if (debug)
            cout << " nT=" << nT << " s=" << s << " wt=" << wt << endl;
        for (int iH = 1; iH < nMatr; ++iH) {
            // alpha_s
            double PI = acos(-1.);
            double alphaQCD = pdf->alphasQ2(scale2);
            mtr = hMatr[iH]->Interpolate(nT, xs);
            mtr2 += pow(4 * PI*alphaQCD, 3) * pow(mtr, 2);
        };
        double mtr0 = hMatr[0]->Interpolate(nT, xs), mtr20 = pow(mtr0, 2);
//    TNtuple tup("tup", "tup", "hatS:pT2:xF:nT:x1:x2:y:mtr2:mtr20:pdf1:pdf2:wt");
        tup.Fill(s, pT2, xF, nT, x1, x2, y, mtr2, mtr20, pdf1, pdf2, wt);

    };
    tup.Write();
    out_file.Save();

    tup.Project("h0", "hatS", "mtr20*pdf1*pdf2*wt");
    h0->Scale(1. / nEv);
    saveHST(h0, ("dSigmaDs_" + f_to_string(S) + "_" + f_to_string(delta) + "_1.hst").c_str());

    tup.Project("hAll", "sqrt(hatS)", "mtr2*pdf1*pdf2*wt");
    hAll->Scale(1. / nEv);
    saveHST(hAll, ("dSigmaDs_" + f_to_string(S) + "_" + f_to_string(delta) + "_All.hst").c_str());
    
    TH1D *hPt=new TH1D("hPt","hPt",20,0.5,2.5); hPt->Sumw2();
    tup.Project("hPt", "sqrt(pT2)", "mtr2*pdf1*pdf2*wt");
    hPt->Scale(1./nEv);
    saveHST(hPt,("hPt_"+f_to_string(S)+"_"+f_to_string(delta)+".hst").c_str());
    
    TH1D *hXF=new TH1D("hXF","XF",20,-0.4,0.4); hXF->Sumw2();
    tup.Project("hXF", "xF", "mtr2*pdf1*pdf2*wt");
    hXF->Scale(1./nEv);
    saveHST(hXF,("hXF_"+f_to_string(S)+"_"+f_to_string(delta)+".hst").c_str());


    return 0;
}
