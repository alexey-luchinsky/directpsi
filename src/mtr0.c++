
#include "mtr0.h"

double calc_matr0_0() { return 1;}

double calc_matr0_1() {
// -1 -1 -> 0 -1
	 double mtr=0;
	 mtr *= pow(NG,3)*NJ0;
	 return mtr;
}

double calc_matr0_2() {
// -1 -1 -> 1 -1
	 double mtr=0;
	 mtr *= pow(NG,3)*NJ1;
	 return mtr;
}

double calc_matr0_3() {
// -1 -1 -> 2 -1
	 double mtr=0;
	 mtr *= pow(NG,3)*NJ2;
	 return mtr;
}

double calc_matr0_4() {
// -1 -1 -> 0 1
	 double mtr=(-32*pow(Mcc,3)*pow(s,3)*pow(T,2)*pow(U,2))/((s - T)*(s - U));
	 mtr *= pow(NG,3)*NJ0;
	 return mtr;
}

double calc_matr0_5() {
// -1 -1 -> 1 1
	 double mtr=(64*pow(Mcc,3)*pow(s,3)*pow(T,2)*pow(U,2))/((s - T)*(s - U));
	 mtr *= pow(NG,3)*NJ1;
	 return mtr;
}

double calc_matr0_6() {
// -1 -1 -> 2 1
	 double mtr=0;
	 mtr *= pow(NG,3)*NJ2;
	 return mtr;
}

double calc_matr0_7() {
// -1 1 -> 0 -1
	 double mtr=(-32*pow(Mcc,3)*pow(s,2)*pow(T,3)*pow(U,2))/((s - T)*(T + U));
	 mtr *= pow(NG,3)*NJ0;
	 return mtr;
}

double calc_matr0_8() {
// -1 1 -> 1 -1
	 double mtr=(64*Mcc*pow(s,2)*pow(T,3)*pow(U,2)*(pow(s,2)*(T - U) - U*pow(T + U,2) + s*(-pow(T,2) + T*U + 2*pow(U,2))))/((s - T)*(s - U)*pow(T + U,2));
	 mtr *= pow(NG,3)*NJ1;
	 return mtr;
}

double calc_matr0_9() {
// -1 1 -> 2 -1
	 double mtr=(128*pow(Mcc,5)*pow(s,2)*pow(T,3)*pow(U,2))/((s - T)*(s - U)*pow(T + U,2));
	 mtr *= pow(NG,3)*NJ2;
	 return mtr;
}

double calc_matr0_10() {
// -1 1 -> 0 1
	 double mtr=(32*pow(Mcc,3)*pow(s,2)*pow(T,2)*pow(U,3))/((s - U)*(T + U));
	 mtr *= pow(NG,3)*NJ0;
	 return mtr;
}

double calc_matr0_11() {
// -1 1 -> 1 1
	 double mtr=(-64*Mcc*pow(s,2)*pow(T,2)*pow(U,3)*(pow(s,2)*(T - U) + T*pow(T + U,2) + s*(-2*pow(T,2) - T*U + pow(U,2))))/((s - T)*(s - U)*pow(T + U,2));
	 mtr *= pow(NG,3)*NJ1;
	 return mtr;
}

double calc_matr0_12() {
// -1 1 -> 2 1
	 double mtr=(-128*pow(Mcc,5)*pow(s,2)*pow(T,2)*pow(U,3))/((s - T)*(s - U)*pow(T + U,2));
	 mtr *= pow(NG,3)*NJ2;
	 return mtr;
}

double calc_matr0_13() {
// 1 -1 -> 0 -1
	 double mtr=(-32*pow(Mcc,3)*pow(s,2)*pow(T,2)*pow(U,3))/((s - U)*(T + U));
	 mtr *= pow(NG,3)*NJ0;
	 return mtr;
}

double calc_matr0_14() {
// 1 -1 -> 1 -1
	 double mtr=(-64*Mcc*pow(s,2)*pow(T,2)*pow(U,3)*(pow(s,2)*(T - U) + T*pow(T + U,2) + s*(-2*pow(T,2) - T*U + pow(U,2))))/((s - T)*(s - U)*pow(T + U,2));
	 mtr *= pow(NG,3)*NJ1;
	 return mtr;
}

double calc_matr0_15() {
// 1 -1 -> 2 -1
	 double mtr=(-128*pow(Mcc,5)*pow(s,2)*pow(T,2)*pow(U,3))/((s - T)*(s - U)*pow(T + U,2));
	 mtr *= pow(NG,3)*NJ2;
	 return mtr;
}

double calc_matr0_16() {
// 1 -1 -> 0 1
	 double mtr=(32*pow(Mcc,3)*pow(s,2)*pow(T,3)*pow(U,2))/((s - T)*(T + U));
	 mtr *= pow(NG,3)*NJ0;
	 return mtr;
}

double calc_matr0_17() {
// 1 -1 -> 1 1
	 double mtr=(64*Mcc*pow(s,2)*pow(T,3)*pow(U,2)*(pow(s,2)*(T - U) - U*pow(T + U,2) + s*(-pow(T,2) + T*U + 2*pow(U,2))))/((s - T)*(s - U)*pow(T + U,2));
	 mtr *= pow(NG,3)*NJ1;
	 return mtr;
}

double calc_matr0_18() {
// 1 -1 -> 2 1
	 double mtr=(128*pow(Mcc,5)*pow(s,2)*pow(T,3)*pow(U,2))/((s - T)*(s - U)*pow(T + U,2));
	 mtr *= pow(NG,3)*NJ2;
	 return mtr;
}

double calc_matr0_19() {
// 1 1 -> 0 -1
	 double mtr=(32*pow(Mcc,3)*pow(s,3)*pow(T,2)*pow(U,2))/((s - T)*(s - U));
	 mtr *= pow(NG,3)*NJ0;
	 return mtr;
}

double calc_matr0_20() {
// 1 1 -> 1 -1
	 double mtr=(64*pow(Mcc,3)*pow(s,3)*pow(T,2)*pow(U,2))/((s - T)*(s - U));
	 mtr *= pow(NG,3)*NJ1;
	 return mtr;
}

double calc_matr0_21() {
// 1 1 -> 2 -1
	 double mtr=0;
	 mtr *= pow(NG,3)*NJ2;
	 return mtr;
}

double calc_matr0_22() {
// 1 1 -> 0 1
	 double mtr=0;
	 mtr *= pow(NG,3)*NJ0;
	 return mtr;
}

double calc_matr0_23() {
// 1 1 -> 1 1
	 double mtr=0;
	 mtr *= pow(NG,3)*NJ1;
	 return mtr;
}

double calc_matr0_24() {
// 1 1 -> 2 1
	 double mtr=0;
	 mtr *= pow(NG,3)*NJ2;
	 return mtr;
}

double calc_matr0_25() {
// 1 0 -> 1 1
	 double mtr=0;
	 mtr *= pow(NG,3)*NJ1;
	 return mtr;
}

