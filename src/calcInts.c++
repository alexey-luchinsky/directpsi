#include <iostream>
#include <fstream>
//#include "math.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TVectorD.h"

#include <tclap/CmdLine.h>

#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

using namespace std;

#include "ramC/Random.h"
#include "ramC/rambo2.h"
#include "mtr.h"

rambo2 ramCM;
EvtVector4R k1, k2, P, k3, q;
double s;
double T, U, k1q, Q, k2q, qq, mc;
double NG, NJ0, NJ1, NJ2;

int nEv;
double sMin, sMax, alpha;
string output_name;


double matr[nMatr];
TNtuple *tup;




double gen_event(bool debug) {
    ramCM.next();
//    double wt = 1;
    P = *ramCM.getV(0);
    k3 = *ramCM.getV(1);

    // throw relative quarks momentum
    Q = ramCM.random_generator->rand(0, Mcc / 2);
//    wt *= pow(Q, 2);
    double cosQ = ramCM.random_generator->rand(-1., 1.), sinQ = sqrt(1. - cosQ * cosQ);
//    wt *= 2;
    const double PI = acos(-1);
    double phi = ramCM.random_generator->rand(0, 2 * PI);
//    wt *= 2 * PI;
    q.set(0, Q * sinQ * cos(phi), Q * sinQ * sin(phi), Q * cosQ);
    q.applyBoostTo(P);

    // calculate scalar products
    T = -(k1 - k3).mass2();
    U = -(k2 - k3).mass2();
    k1q = k1*q, k2q = k2*q, qq = q.mass2();
    mc = sqrt(Mcc * Mcc / 4 - Q * Q);

    // debug print
    if (debug) {
//        cout << "wt=" << wt << endl;
        cout << " k1=" << k1 << " k2=" << k2 << endl;
        cout << "P=" << P << " Mpsi=" << P.mass() << " k3=" << k3 << endl;
        cout << " s=" << s << " T=" << T << " U=" << U << " s+t+u=" << s - T - U << endl;
        cout << " Q^2=" << Q * Q << endl;
        cout << "q=" << q << " q^2=" << q.mass2() << " Pq=" << P * q << endl;
    };


    return 1.;
}

void fill_matr() {
    NG = 1. / sqrt(s * T * U);
    NJ0 = 2. / sqrt(s * T * U);
    NJ1 = 1. / sqrt(s * T * U);
    NJ2 = 1. / Mcc;

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

void calc_integrals(int nEv) {
    Float_t values[nMatr + 5];



    // g g -> psi g random space
    ramCM.setMass(0, Mcc);
    ramCM.setMass(1, 0.);

    for (int iEv = 0; iEv < nEv; ++iEv) {
        if( (iEv % (nEv/10))==0)
            cout<<"======== "<<(int)(100.*iEv/nEv)<<"% =========="<<endl;
        double xs=ramCM.random_generator->rand(0,1);
        s = sMin+(sMax-sMin)*pow(xs,alpha);

        double ecm = sqrt(s);
        ramCM.setECM(ecm);
        k1.set(ecm / 2, 0, 0, ecm / 2);
        k2.set(ecm / 2, 0, 0, -ecm / 2);

        bool debug = (iEv < 3);
        if (debug) cout << "--------------- Debug print at iEv=" << iEv << "-------------" << endl;
        double weight = gen_event(debug);
        fill_matr();

        double cosPsi = P.get(3) / P.d3mag();

        for (int i = 0; i < nMatr; ++i) {
            if(i>0) matr[i] = matr[i]/(Mcc+2*mc)/pow(mc,3./2);
            values[i] = matr[i];
        };
        values[nMatr + 1] = xs;
        values[nMatr + 2] = cosPsi;
        values[nMatr + 3] = weight;
        //        values[nMatr+4]=0;
        values[nMatr + 4] = Q*Q;
        //        values[nMatr+6]=0;
        tup->Fill(values);
        if (debug) {
            cout << " weight=" << weight << endl;
            cout << " values[nMatr+2]=" << values[nMatr + 2] << endl;
        };

    };

}


bool file_exists(string name) {
    if (FILE * file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    }
    return false;
}

bool init_commandline_args(int argc, char **argv) {
    try {
        TCLAP::CmdLine cmd("Calculates integrated over q helicity matrix elements for gg->psi g process", ' ', "0.9");
        //        TCLAP::ValueArg<float> s_arg("s","s","s",false,10,"float",cmd);
        TCLAP::ValueArg<int> n_arg("n", "n", "log_10(nEv)", false, 4, "float", cmd);
        TCLAP::ValueArg<string> out_arg("o", "o", "name for output root file", false, "matr.root", "string", cmd);
        TCLAP::ValueArg<float> sMin_arg("m", "sMin", "minimum s", false, -1, "float", cmd);
        TCLAP::ValueArg<float> sMax_arg("M", "sMax", "maximum s", false, 0, "float", cmd);
        TCLAP::ValueArg<float> alpha_arg("a","alpha","parameter of smart s sampling",false,3,"float",cmd);
        cmd.parse(argc, argv);
        sMin = sMin_arg.getValue(); if(sMin<Mcc*Mcc) sMin=Mcc*Mcc;
        sMax = sMax_arg.getValue(); if(sMax<sMin) sMax=1.1*sMin;
        alpha = alpha_arg.getValue();
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
    if (!init_commandline_args(argc, argv)) return -1;
    cout << " sMin=" << sMin << " sMax=" << sMax << endl;
    cout<<" alpha="<<alpha<<endl;
    cout << " nEv=" << nEv << endl;
    cout << " output_name=" << output_name << endl;



//    cout << " s=" << s << " nEv=" << nEv << endl;
    cout << " seed=" << ramCM.random_generator->get_seed() << endl;


    //    string tuple_file_name="matr_"+f_to_string(s)+".root";
    TFile file(output_name.c_str(), "RECREATE");
    TVectorD stats(3);
    stats[0]=sMin; stats[1]=sMax; stats[2]=alpha;
    stats.Write("stats");
    
    
    // create tuple
    string field_names = "";
    for (int i = 0; i <= nMatr; ++i) {
        field_names += "matr_" + i_to_string(i) + ":";
    };
    field_names += "xs:";
    field_names += "cosPsi:";
    field_names += "wt:";
    field_names += "q2";
    cout << field_names << endl;
    tup = new TNtuple("tup", "tup", field_names.c_str());



    calc_integrals(nEv);



    tup->Write();
    file.Save();



    return 0;
}
