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
#include "ramC/rambo2.h"

using namespace std;
using namespace LHAPDF;

TH2D * hMatr[nMatr]; // hMatr[nT, xs], nT=-(k1-k3)^2/(s-Mcc)^2 in [0,1])
TFile *hist_file;




string in_fileName, out_fileName;
double delta, S;
double sMin, sMax, alpha;
int nEv;
string pdfName;
int nBins;

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

string prefix;
int nDebug;

bool init_commandline_args(int argc, char **argv) {
    try {
        TCLAP::CmdLine cmd("Calculate hadronic process using double_export.exe program data", ' ', "0.9");
        TCLAP::ValueArg<string> in_agr("i", "in", "input file name", false, "interpolation.root", "string", cmd);
        TCLAP::ValueArg<string> out_arg("o", "out", "output file name", false, "PP.root", "string", cmd);
        TCLAP::ValueArg<float> e_arg("e", "e", "energy of hadronic reaction", false, 500, "float", cmd);
        TCLAP::ValueArg<float> n_arg("n", "n", "log_10(nEv)", false, 6, "float", cmd);
        TCLAP::ValueArg<string> pdfName_arg("p", "pdf", "pdf set name", false, "CT10", "string", cmd);
        TCLAP::ValueArg<string> prefix_arg("","prefix","prefix for saved files",false,"","string",cmd);
        TCLAP::ValueArg<int> nBins_arg("","nBins","number of histogram bins",false,10,"int",cmd);
        TCLAP::ValueArg<int> nDebug_arg("","debud","number of debug events",false,0,"int",cmd);
        cmd.parse(argc, argv);
        in_fileName = in_agr.getValue();
        out_fileName = out_arg.getValue();
        S = pow(e_arg.getValue(),2);
        nEv = pow(10, n_arg.getValue());
        if (nEv > 1e10) {
            cout << " nEv=" << nEv << " is too large!" << endl;
            return false;
        }
        pdfName = pdfName_arg.getValue();
        prefix=prefix_arg.getValue();
        nBins=nBins_arg.getValue();
        nDebug=nDebug_arg.getValue();
        return true;
    } catch (TCLAP::ArgException e) {
        cout << " error " << e.error() << " for arg " << e.argId() << endl;
        return false;
    }
}

void save_pdf(const LHAPDF::PDF *pdf, double q2) {
    ofstream pdf_file;
    pdf_file.open((prefix+"pdf.txt").c_str());
    double dx = 1e-4;
    for (double x = dx; x < 1. - dx; x += dx) {
        pdf_file << x << " " << pdf->xfxQ2(0, x, q2) / x << endl;
    }
    pdf_file.close();
}

int main(int argc, char **argv) {
    init_commandline_args(argc, argv);

    const int imem = 0;
    const PDF *pdf = mkPDF(pdfName);

    double Mcc2 = Mcc*Mcc, scale2 = Mcc2;
    save_pdf(pdf, scale2);

    if (!load_integrals()) {
        cout << "Cannot open file " << in_fileName << endl;
        return -1;
    }
    cout << " sMin=" << sMin << " sMax=" << sMax << " alpha=" << alpha << " delta=" << delta << endl;


    Random random_generator;
    rambo2 ram;
    ram.setMass(0, Mcc);
    ram.setMass(1, 0.);

    TFile out_file((prefix+out_fileName).c_str(), "RECREATE");
    TNtuple tup("tup", "tup", "hatS:pT2:xF:nT:x1:x2:y:yPsi:mtr2:mtr20:pdf1:pdf2:wt");
    // initialize  histograms
    TH1D *h_mFinal=new TH1D("mFinal","mFinal",nBins, sqrt(sMin), 5); h_mFinal->Sumw2();
    TH1D *h_pT2=new TH1D("pT2","pT2",nBins,0,(S-Mcc2)/(2*sqrt(S))); h_pT2->Sumw2();
    TH1D *h_xF=new TH1D("xF","xF",nBins,-2,2); h_xF->Sumw2();
    TH1D *h_yPsi=new TH1D("yPsi","yPsi",nBins,-1,1); h_yPsi->Sumw2();
    
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


    EvtVector4R k1, k2, P, k3;

    for (int iEv = 0; iEv < nEv; ++iEv) {
        bool debug = (iEv < nDebug);
//        if( iEv % ())
        if (debug) cout << "----- Debug print at i=" << iEv << "---------" << endl;
        double wt = 1;
        double xs = random_generator.rand(0, 1);
        double s = sMin + (sMax - sMin) * pow(xs, alpha);
        wt *= alpha * (sMax - sMin) * pow(xs, alpha - 1) / S;
        double ecm = sqrt(s);
        if(debug) {
            cout<<" ecm="<<ecm<<";"<<endl;
        }
        ram.setECM(ecm);
        if (!ram.next()) continue;
        P = *ram.getV(0);   k3 = *ram.getV(1);
        if(debug) {
            cout<<" before boost:"<<endl;
            cout<<"\t P="<<P<<"; k3="<<k3<<";"<<endl;
        }

        // y
        double yMax = log(S / s) / 2;
        double y = random_generator.rand(-yMax, yMax);
        if (debug) cout << "\t y=" << y << endl;
        wt *= 2 * yMax; // y


        double x1 = sqrt(s / S) * exp(y), pdf1 = pdf->xfxQ2(0, x1, scale2) / x1;
        if (debug) cout << "\t x1=" << x1 << "; pdf1=" << pdf1 <<";"<< endl;
        double x2 = sqrt(s / S) * exp(-y), pdf2 = pdf->xfxQ2(0, x2, scale2) / x2;

        k1.set(sqrt(S)*x1/2, 0, 0, sqrt(S)*x1/2);
        k2.set(sqrt(S)*x2/2, 0, 0, -sqrt(S)*x2/2);
        double gamma=(x1+x2)/(2*sqrt(x1*x2)), beta=(x1-x2)/(x1+x2);
        P.set( gamma*(P.get(0)+beta*P.get(3)), P.get(1), P.get(2), gamma*(P.get(3)+beta*P.get(0)));
        k3.set( gamma*(k3.get(0)+beta*k3.get(3)), k3.get(1), k3.get(2), gamma*(k3.get(3)+beta*k3.get(0)));
        
        
        double t=(k1-P).mass2(), nT=t/(Mcc2-s), u=(k1-k3).mass2();
        double pT2=get_pT2(P);
        double xF=2*P.get(3)/ecm;
        if (debug) {
            cout << "\t s=" << s << ";"<<endl;
            cout<<" x1="<<x1<<"; x2="<<x2<<";"<<endl;
            cout<<" gamma="<<gamma<<"; beta="<<beta<<";"<<endl;
            cout << "\t k1=" << k1 << "; k2=" << k2 << ";"<<endl;
            cout << "\t aP=" << P <<"; P^2="<<P.mass2()<< "; ak3=" << k3 << "; k3^2="<<k3.mass2()<<";"<<endl;
            cout<<" Ptot="<<k1+k2<<"="<<P+k3<<endl;
            cout << "\t t="<<t<<"; nT="<<nT<<" u="<<u<<"; s+t+u="<<s+t+u<<endl;
        };

        
        wt *= (s - Mcc2);
        // conversion to dsdt
        wt *= 1. / (64 * PI * s)*4 / s;
        // symmetry, etc
        wt *= 1. / (2 * 2 * 8 * 2 * 8);
        // tranfer to nb
        wt *= 0.389*1e6;

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
        tup.Fill(s, pT2, xF, nT, x1, x2, y, getRapidity(P), mtr2, mtr20, pdf1, pdf2, wt);
        
        // fill histograms
        h_mFinal->Fill(sqrt(s),mtr2*pdf1*pdf2*wt);
        h_pT2->Fill(get_pT2(P), mtr2*pdf1*pdf2*wt);
        h_yPsi->Fill(getRapidity(P),mtr2*pdf1*pdf2*wt);
        h_xF->Fill(xF,mtr2*pdf1*pdf2*wt);

    };
    tup.Write();
    
    string hist_name="_e"+f_to_string(sqrt(S))+"_d"+f_to_string(delta)+"_a"+f_to_string(alpha)+"_"+pdfName+".hst";
    
    h_mFinal->Scale(1./nEv); h_mFinal->Write(); saveHST(h_mFinal,prefix+"m"+hist_name);
    h_pT2->Scale(1./nEv); h_pT2->Write(); saveHST(h_pT2,prefix+"pT2"+hist_name);
    h_yPsi->Scale(1./nEv); h_yPsi->Write(); saveHST(h_yPsi,prefix+"yPsi"+hist_name);
    h_xF->Scale(1./nEv); h_xF->Write(); saveHST(h_xF,prefix+"xF"+hist_name);
    out_file.Save();

    return 0;
}
