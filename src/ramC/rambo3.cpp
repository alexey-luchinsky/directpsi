/* 
 * File:   rambo3.cpp
 * Author: luchinsky
 * 
 * Created on March 21, 2016, 11:34 AM
 */

#include "rambo3.h"

// 3D LIPS q k3 -> k1 k2 k3

rambo3::rambo3(double ecm_, double m1_, double m2_, double m3_) {
    m1 = m1_;
    m2 = m2_;
    m3 = m3_;
    min_cos = -1, max_cos = 1;
    PI = acos(-1);
    setECM(ecm_);
    random_generator = new Random(0);
}

void rambo3::setECM(double ecm_) {
    ecm = ecm_;
    min_q2 = pow(m1 + m2, 2);
    max_q2 = pow(ecm - m3, 2);
}

double rambo3::next() {
    do {
        WT = 1;
        // q k3 phase space
        q2 = random_generator->rand(min_q2, max_q2);
        WT *= (max_q2 - min_q2);
        double cosT = random_generator->rand(min_cos, max_cos), sinT = sqrt(1. - cosT * cosT);
        WT *= (max_cos - min_cos);
        double phi = random_generator->rand(0, 2 * PI), sinPhi = sin(phi), cosPhi = cos(phi);
        WT *= 2 * PI;
        double eQ = (ecm * ecm + q2 - m3 * m3) / (2 * ecm);
        double e3 = (ecm * ecm + m3 * m3 - q2) / (2 * ecm);
        double mom = sqrt(eQ * eQ - q2);
        q.set(eQ, mom * sinT*cosPhi, mom * sinT*sinPhi, mom * cosT);
        k3.set(e3, -mom * sinT*cosPhi, -mom * sinT*sinPhi, -mom * cosT);
        WT *= mom / (ecm * 4 * pow(2 * PI, 6));

        // q -> k1 k2 phase space
        double e1 = (q2 + m1 * m1 - m2 * m2) / (2 * sqrt(q2));
        double e2 = (q2 + m2 * m2 - m1 * m1) / (2 * sqrt(q2));
        mom = sqrt(e1 * e1 - m1 * m1);
        cosT = random_generator->rand(-1, 1);
        sinT = sqrt(1. - cosT * cosT);
        WT *= 2;
        phi = random_generator->rand(0, 2 * PI);
        sinPhi = sin(phi);
        cosPhi = cos(phi);
        WT *= 2 * PI;
        k1.set(e1, mom * sinT*cosPhi, mom * sinT*sinPhi, mom * cosT);
        k2.set(e2, -mom * sinT*cosPhi, -mom * sinT*sinPhi, -mom * cosT);
        WT *= mom / (sqrt(q2)*4 * pow(2 * PI, 6));

        // boost
        k1.applyBoostTo(q);
        k2.applyBoostTo(q);

        WT *= pow(2 * PI, 3);
    } while (!(WT * WT > 0));
    return WT;

}

EvtVector4R *rambo3::getV(int i) {
    if (i == 0) return &k1;
    else if (i == 1) return &k2;
    else if (i == 2) return &k3;
    else {
        std::cout << "Wrong particle number " << i << " in rambo2::getV" << std::endl;
        return NULL;
    }
}