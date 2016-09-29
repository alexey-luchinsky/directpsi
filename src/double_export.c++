#include <iostream>
#include "TFile.h"
#include "TNtuple.h"
#include "TH2D.h"
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
        TCLAP::ValueArg<string> in_agr("i","in","input file name",false,"matr.root","string",cmd);
        TCLAP::ValueArg<string> out_arg("o","out","output file name",false,"interpolation.root","string",cmd);
        TCLAP::ValueArg<int> tb_agr("t","tb","number of t bins",false,50,"int",cmd);
        TCLAP::ValueArg<int> sb_arg("s","sb","number of s bins",false,100,"int",cmd);
        TCLAP::ValueArg<float> delta_arg("d","delta","wave function parameter",false,0.4,"float",cmd);
        cmd.parse(argc, argv);
        in_fileName=in_agr.getValue();
        out_fileName=out_arg.getValue();
        nTBin=tb_agr.getValue();
        nsBin=sb_arg.getValue();
        delta=delta_arg.getValue();
        
        return true;
    }    catch (TCLAP::ArgException e) {
        cout << " error " << e.error() << " for arg " << e.argId() << endl;
        return false;
    };
}

int main(int argc, char **argv) {
    if(!init_command_args(argc, argv))
        return -1;
    
    // determine WF normalization coefficient
    double sum=0, Mcc=3.1;
    Random random_generator;
    for(int i=0; i<1e6; ++i) {
        double Q=random_generator.rand(0,Mcc/2);
        double PI=acos(-1.);
        sum += 4*PI*Q*Q*wave_function(Q*Q,delta)*Mcc/2;
    };
    sum=sum/1e6;
    cout<<" WF normalization : "<<sum<<endl;

    TFile in_file(in_fileName.c_str(), "READ");
    if(!in_file.IsOpen()) {
        cout<<" Cannot open file "<<in_fileName<<endl;
        return -1;
    };
    TNtuple *tup = (TNtuple*) in_file.Get("tup");
    double sMin = tup->GetMinimum("s"), sMax = tup->GetMaximum("s");
    int nEv = tup->GetEntries();
    cout << " nEv=" << nEv << " sMin=" << sMin << " sMax=" << sMax << " delta="<<delta<<endl;

    TFile out_file(out_fileName.c_str(), "RECREATE");
    out_file.cd();

    const int nMatr = 26;
    float matr[nMatr];
    TH2D * hMatr[nMatr];
    for (int iH = 0; iH < nMatr; ++iH) {
        string name = "matr_" + i_to_string(iH);
        tup->SetBranchAddress(name.c_str(), &matr[iH]);
        name = "h" + i_to_string(iH);
        hMatr[iH] = new TH2D(name.c_str(), name.c_str(), nTBin, 0, 1, nsBin, sMin, sMax);
        hMatr[iH]->Sumw2();
    }
    float weight, wf, cosPsi, q2, s;
    tup->SetBranchAddress("wt", &weight);
    tup->SetBranchAddress("cosPsi", &cosPsi);
    tup->SetBranchAddress("s", &s);

    for (int iEv = 0; iEv < nEv; ++iEv) {
        tup->GetEntry(iEv);
        double WF = wave_function(q2, delta);
        double nT = (1 + cosPsi) / 2;
        for (int iH = 0; iH < nMatr; ++iH)
            hMatr[iH]->Fill(nT, s, matr[iH] * weight * WF);
    };
    for (int iH = 0; iH < nMatr; ++iH) {
        hMatr[iH]->Scale(1. / hMatr[iH]->GetEntries());
        hMatr[iH]->Write();
    };
    out_file.Save();

    return 0;
}
