#include "utils.h"

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

void saveHST(TH1D *hist, TString name, bool print) {
    if (print) cout << " Saving " << name << endl;
    FILE *file = fopen(name.Data(), "w");
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        fprintf(file, "%e %e %e\n", hist->GetBinCenter(i), hist->GetBinContent(i) / hist->GetBinWidth(i), hist->GetBinError(i) / hist->GetBinWidth(i));
    };
    if (print) cout << "\t Histogram sum=" << hist->GetSum() << endl;
    fclose(file);
};

double wave_function(double q2, double delta) {
    return exp(-q2 / delta / delta);
}

double get_pT2(EvtVector4R P) {
    return pow(P.get(1),2) + pow(P.get(2),2);
}

double getRapidity(EvtVector4R P) {
    return log((P.get(3)+P.get(0))/(P.get(3)-P.get(0)))/2;
}