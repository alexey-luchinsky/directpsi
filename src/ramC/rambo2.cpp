/* 
 * File:   rambo2.cpp
 * Author: luchinsky
 * 
 * Created on March 21, 2016, 10:37 AM
 */

#include "rambo2.h"

rambo2::rambo2(double ecm_, double m1_, double m2_) {
    ecm=ecm_; m1=m1_; m2=m2_;
    PI=acos(-1);
//    std::cout<<rnd()<<std::endl;
}

double rambo2::getECM() {
    return ecm;
}

void rambo2::setECM(double ecm_) {
    ecm=ecm_;
}

double rambo2::getMass(int i) {
    if(i==0) return m1;
    else if(i==1) return m2;
    else {
        std::cout<<"wrong mass number "<<i<<" in rambo2::getMass()"<<std::endl;
        return -1;
    }
}

void rambo2::setMass(int i, double m) {
    if(i==0) m1=m;
    else if(i==1) m2=m;
    else 
        std::cout<<"wrong mass number "<<i<<" in rambo2::setMass()"<<std::endl;
}

EvtVector4R* rambo2::getV(int i) {
    if(i==0) return &k1;
    else if(i==1) return &k2;
    else {
        std::cout<<"wrong mass number "<<i<<" in rambo2::getV()"<<std::endl;
        return NULL;
    };
}

double rambo2::getWT() {
    return WT;
}

double rambo2::rnd(double min, double max) {
    return min+(max-min)*rand()/RAND_MAX;
}

double rambo2::next() {
    double e1=(ecm*ecm+m1*m1-m2*m2)/(2*ecm);
    double e2=(ecm*ecm+m2*m2-m1*m1)/(2*ecm);
    double mom=sqrt(e1*e1-m1*m1);
    double cosT=rnd(-1,1), sinT=sqrt(1.-cosT*cosT);
    double phi=rnd(0., 2*PI), cosPhi=cos(phi), sinPhi=sin(phi);
    k1.set(e1,mom*sinT*cosPhi, mom*sinT*sinPhi, mom*cosT);
    k2.set(e2,-mom*sinT*cosPhi, -mom*sinT*sinPhi, -mom*cosT);
    WT=mom/(ecm*4*pow(2*PI,6))*2*(2*PI);
    return WT;
}