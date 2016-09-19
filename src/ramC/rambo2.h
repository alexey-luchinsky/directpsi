/* 
 * File:   rambo2.h
 * Author: luchinsky
 *
 * Created on March 21, 2016, 10:37 AM
 * 
 * 2D LIPS without RAMBO
 * RPP agreements
 */

#ifndef RAMBO2_H
#define	RAMBO2_H

#include "../EvtGenBase/EvtVector4R.hh"
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "Random.h"
using namespace std;

class rambo2 {
public:
    rambo2(double ecm_=10, double m1_=0, double m2_=0);
//    rambo2(const rambo2& orig);
//    virtual ~rambo2();
    void setECM(double ecm_);
    double getECM();
    void setMass(int i, double m);
    double getMass(int i);
    EvtVector4R *getV(int i);
    double next();
    double getWT();
    double rnd(double min=0, double max=1);
    Random *random_generator;
private:
    double PI;
    double ecm, m1, m2;
    EvtVector4R k1, k2;
    double WT;

};

#endif	/* RAMBO2_H */

