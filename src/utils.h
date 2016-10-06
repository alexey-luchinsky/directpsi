/* 
 * File:   newfile.h
 * Author: luchinsky
 *
 * Created on September 29, 2016, 11:26 AM
 */

#ifndef NEWFILE_H
#define	NEWFILE_H

#include <iostream>
#include "TH2D.h"
#include "EvtGenBase/EvtVector4R.hh"
using namespace std;

string f_to_string(double v);
string i_to_string(int v);
void saveHST(TH1D *hist, TString name, bool print = false);
double wave_function(double q2, double delta);
const int nMatr = 26;
const double PI=acos(-1.);
double const Mcc=3.1;
double get_pT2(EvtVector4R P);
#endif	/* NEWFILE_H */
