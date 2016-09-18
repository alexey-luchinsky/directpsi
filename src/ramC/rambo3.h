/* 
 * File:   rambo3.h
 * Author: luchinsky
 *
 * Created on March 21, 2016, 11:34 AM
 * 3D LIPS q k3 -> k1 k2 k3
 */

#ifndef RAMBO3_H
#define	RAMBO3_H

#include "../EvtGenBase/EvtVector4R.hh"
#include <iostream>
#include <stdlib.h>
#include "Random.h"


using namespace std;
class rambo3 {
public:
    rambo3(double ecm_=10, double m1_=0, double m2_=0, double m3_=0);
    void setECM(double ecm_);
    
    double next();
    
    EvtVector4R *getV(int i);

    double ecm, m1, m2,m3;
    EvtVector4R k1,k2,k3,q;
    double q2;
    double min_cos, max_cos, min_q2, max_q2;
    double WT;
    double PI;
};

#endif	/* RAMBO3_H */

