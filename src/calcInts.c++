#include <iostream>
#include <fstream>
//#include "math.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"

#include <stdlib.h>
#include <stdio.h>

using namespace std;

#include "ramC/Random.h"
#include "ramC/rambo2.h"

double calc_matr0(); // 1
double calc_matr1(); // -1 -1 -> 0 -1
double calc_matr2(); // -1 -1 -> 1 -1
double calc_matr3(); // -1 -1 -> 2 -1
double calc_matr4(); // -1 -1 -> 0 1
double calc_matr5(); // -1 -1 -> 1 1
double calc_matr6(); // -1 -1 -> 2 1
double calc_matr7(); // -1 1 -> 0 -1
double calc_matr8(); // -1 1 -> 1 -1
double calc_matr9(); // -1 1 -> 2 -1
double calc_matr10(); // -1 1 -> 0 1
double calc_matr11(); // -1 1 -> 1 1
double calc_matr12(); // -1 1 -> 2 1
double calc_matr13(); // 1 -1 -> 0 -1
double calc_matr14(); // 1 -1 -> 1 -1
double calc_matr15(); // 1 -1 -> 2 -1
double calc_matr16(); // 1 -1 -> 0 1
double calc_matr17(); // 1 -1 -> 1 1
double calc_matr18(); // 1 -1 -> 2 1
double calc_matr19(); // 1 1 -> 0 -1
double calc_matr20(); // 1 1 -> 1 -1
double calc_matr21(); // 1 1 -> 2 -1
double calc_matr22(); // 1 1 -> 0 1
double calc_matr23(); // 1 1 -> 1 1
double calc_matr24(); // 1 1 -> 2 1
double calc_matr25(); // 1 0 -> 1 1


Random random_generator(0);
rambo2 ramCM;
EvtVector4R k1, k2, P, k3, q;
double Mcc, s;
double T, U, k1q, Q, k2q, qq, mc;

const int nMatr = 26;
double matr[nMatr];
TH1D *hMatr[nMatr];
TNtuple *tup;

// psi wave function
double delta;

double WaveFunction(double Q, double del) {
    return exp(-pow(Q / del, 2));
    //return 1;
};

double gen_event(bool debug) {
    ramCM.next();
    double wt = 1;
    P = *ramCM.getV(0);
    k3 = *ramCM.getV(1);

    // throw relative quarks momentum
    Q = random_generator.rand(0, Mcc / 2);
    wt *= pow(Q, 2);
    double cosQ = random_generator.rand(-1., 1.), sinQ = sqrt(1. - cosQ * cosQ);
    wt *= 2;
    const double PI = acos(-1);
    double phi = random_generator.rand(0, 2 * PI);
    wt *= 2 * PI;
    q.set(0, Q * sinQ * cos(phi), Q * sinQ * sin(phi), Q * cosQ);
    q.applyBoostTo(P);

    // calculate scalar products
    T = -(k1 - k3).mass2();
    U = -(k2 - k3).mass2();
    k1q = k1*q, k2q = k2*q, qq = q.mass2();
    mc = sqrt(Mcc * Mcc / 4 - Q * Q);

    // debug print
    if (debug) {
        cout << "wt=" << wt << endl;
        cout << " k1=" << k1 << " k2=" << k2 << endl;
        cout << "P=" << P << " Mpsi=" << P.mass() << " k3=" << k3 << endl;
        cout << " s=" << s << " T=" << T << " U=" << U << " s+t+u=" << s - T - U << endl;
        cout << " Q^2=" << Q * Q << endl;
        cout << "q=" << q << " q^2=" << q.mass2() << " Pq=" << P * q << endl;
    };


    return wt;
}

void fill_matr() {
    matr[0] = calc_matr0();
    matr[1] = calc_matr1();
    matr[2] = calc_matr2();
    matr[3] = calc_matr3();
    matr[4] = calc_matr4();
    matr[5] = calc_matr5();
    matr[6] = calc_matr6();
    matr[7] = calc_matr7();
    matr[8] = calc_matr8();
    matr[9] = calc_matr9();
    matr[10] = calc_matr10();
    matr[11] = calc_matr11();
    matr[12] = calc_matr12();
    matr[13] = calc_matr13();
    matr[14] = calc_matr14();
    matr[15] = calc_matr15();
    matr[16] = calc_matr16();
    matr[17] = calc_matr17();
    matr[18] = calc_matr18();
    matr[19] = calc_matr19();
    matr[20] = calc_matr20();
    matr[21] = calc_matr21();
    matr[22] = calc_matr22();
    matr[23] = calc_matr23();
    matr[24] = calc_matr24();
    matr[25] = calc_matr25();
}

void calc_integrals(double s_, int nEv) {
    Float_t values[nMatr+6];
    
    // clear histograms
    s = s_;
    for (int i = 0; i < nMatr; ++i)
        hMatr[i]->Reset();

    double ecm = sqrt(s);
    k1.set(ecm / 2, 0, 0, ecm / 2);
    k2.set(ecm / 2, 0, 0, -ecm / 2);

    
    // g g -> psi g random space
    ramCM.setECM(ecm);
    ramCM.setMass(0, Mcc);
    ramCM.setMass(1, 0.);

    for (int iEv = 0; iEv < nEv; ++iEv) {
        bool debug = (iEv < 3);
        if (iEv % (nEv / 10) == 0) cout << iEv << endl;
        if (debug) cout << "--------------- Debug print at iEv=" << iEv << "-------------" << endl;
        double weight = gen_event(debug);
        fill_matr();

        double cosPsi = P.get(3) / P.d3mag();
        double wave_function = WaveFunction(Q, delta);

        for(int i=0; i<nMatr; ++i) {
            hMatr[i]->Fill(cosPsi,matr[i]*weight*wave_function);
            values[i]=matr[i];
        };
        values[nMatr+1]=s;
        values[nMatr+2]=cosPsi;
        values[nMatr+3]=weight;
        values[nMatr+4]=wave_function;
        values[nMatr+5]=Q*Q;
        tup->Fill(values);
        if (debug) {
            cout<<" weight="<<weight<<endl;
            cout<<" delta="<<delta<<endl;
            cout << " WF=" << wave_function << endl;
            cout<<" values[nMatr+2]="<<values[nMatr+2]<<endl;
        };

    };

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

string f_to_string(double v) {
  char c[20];
  sprintf(c,"%f",v);
  return string(c);
}

string i_to_string(int v) {
  char c[20];
  sprintf(c,"%d",v);
  return string(c);
}

bool file_exists(string name)
{
    if (FILE * file = fopen(name.c_str(), "r"))
    {
        fclose(file);
        return true;
    }
    return false;
}

int main(int argc, char **argv) {
    cout<<"calcInts [s=10] [delta=0.4] [n=1e4]"<<endl;
    s=10.; delta=0.4;
    int nEv=1e4;
    if(argc>1) s=atof(argv[1]);
    if(argc>2) delta=atof(argv[2]);
    if(argc>3) nEv=atof(argv[3]);


    
    cout<<" s="<<s<<" delta="<<delta<<" nEv="<<nEv<<endl;
    cout<<" seed="<<random_generator.get_seed()<<endl;
    
    
    TFile file("matr.root", "RECREATE");
    // create tuple
    string field_names="";      
    for(int i=0; i<=nMatr; ++i) {
        field_names +="matr_"+i_to_string(i)+":";
    };
    field_names+="s:";
    field_names+="cosPsi:";
    field_names+="wt:";
    field_names+="wf:";
    field_names+="q2";
    cout<<field_names<<endl;
    tup=new TNtuple("tup","tup",field_names.c_str());
    
    
    // init histograms
    for (int i = 0; i < nMatr; ++i) {
        string name = "hMatr_" + i_to_string(i);
        hMatr[i] = new TH1D(name.c_str(), name.c_str(), 20, -1, 1);
        hMatr[i]->Sumw2();
    }

    // initial gluon momenta

    Mcc = 3.1;
    cout<<"s="<<s<<" delta="<<delta<<endl;
    string data_path="dat/"+f_to_string(s)+"_"+f_to_string(delta)+"/";
    // check if file exists
    if(file_exists(data_path+"hMatr0.hst")) {
        cout<<"file exists already, exiting"<<endl;
        return 0;
    }
    system(("mkdir -p "+data_path).c_str());
    calc_integrals(s, nEv);


    // normalize histograms
    double int0 = hMatr[0]->Integral();
    for (int i = 0; i < nMatr; ++i)
        hMatr[i]->Scale(1. / int0); // normalize integrals to intQ(1)
    
    // save histograms
    for (int i = 0; i < nMatr; ++i) {
//        hMatr[i]->Write();
        string fileName=data_path+"hMatr"+i_to_string(i)+".hst";
        saveHST(hMatr[i],fileName.c_str());
    };
    tup->Write();
    file.Save();



    return 0;
}
