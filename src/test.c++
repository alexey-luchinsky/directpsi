#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

double s;
double T, U, k1q, Q, k2q, qq, mc;
double NG, NJ0, NJ1, NJ2;

#include "ramC/rambo2.h"
TEST_CASE("ram2","[rambo]") {
    rambo2 ram;
    double ecm=10, m1=0.1, m2=0.2;
    ram.setECM(ecm); ram.setMass(0,m1); ram.setMass(1,m2);
    REQUIRE(Approx(ram.getECM())==ecm);
    REQUIRE(Approx(ram.getMass(0))==m1);
    REQUIRE(Approx(ram.getMass(1))==m2);
    double wt=ram.next();
    EvtVector4R k1=*ram.getV(0), k2=*ram.getV(1);
    REQUIRE(Approx(k1.mass2())==m1*m1);
    REQUIRE(Approx(k2.mass2())==m2*m2);
    double const PI=acos(-1);
    ram.setMass(0,m1); ram.setMass(1,m1);
    REQUIRE(Approx(pow(2*PI,4)*wt).epsilon(1e-4)==sqrt(1.-4*m1*m1/ecm/ecm)/(8*PI));
};

#include "ramC/rambo3.h"
TEST_CASE("ram3","[rambo]") {
    double ecm=10;
    rambo3 ram(ecm);
    double m1=0.1, m2=0.2, m3=0.3;
    ram.m1=m1; ram.m2=m2; ram.m3=m3;
    double wt=ram.next();
    EvtVector4R k1=*ram.getV(0), k2=*ram.getV(1), k3=*ram.getV(2);
    REQUIRE(Approx(k1.mass2())==m1*m1);
    REQUIRE(Approx(k2.mass2())==m2*m2);
    REQUIRE(Approx(k3.mass2())==m3*m3);
    REQUIRE(Approx((k1+k2+k3).mass2())==ecm*ecm);

};
