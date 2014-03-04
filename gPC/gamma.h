//
//  gamma.h
//  gPC
//
//  Created by 马 征 on 14-2-26.
//  Copyright (c) 2014年 马 征. All rights reserved.
//

#ifndef gPC_gamma_h
#define gPC_gamma_h

#include "dynamicVector.h"

double gammln(const double xx) {                           //return ln(gamma(xx)) at xx, when xx>0
    double x,tmp,y,ser;
    static const double cof[14]={57.1562356658629235,-59.5979603554754912, 14.1360979747417471,-0.491913816097620199,.339946499848118887e-4, .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3, -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3, .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5}; if (xx <= 0) throw("bad arg in gammln");
    y=x=xx;
    tmp = x+5.24218750000000000;                           //Rational 671/128.
    tmp = (x+0.5)*log(tmp)-tmp;
    ser = 0.999999999999997092;
    for (int j=0;j<14;j++) ser += cof[j]/++y;
    return tmp+log(2.5066282746310005*ser/x);
}

double factrl(const int n) {                                //Returns the value n! as a floating-point number.
//    static vector<double, 171> a;
    static dynamicVector<double> a(171, 0.0);
    static bool init=true;
    if (init) {
        init = false;
        a(0) = 1.;
        for (int i=1;i<171;i++) a(i) = i*a[i-1];
    }
    if (n < 0 || n > 170) throw("factrl out of range");
    return a[n];
}

double factln(const int n) {                                 //Returns ln.n!.
    static const int NTOP=2000;
//    static vector<double, NTOP> a;
    static dynamicVector<double> a(NTOP, 0.0);
    static bool init=true;
    if (init) {
        init = false;
        for (int i=0;i<NTOP;i++) a(i) = gammln(i+1.); }
    if (n < 0) throw("negative arg in factln");
    if (n < NTOP) return a[n];
    return gammln(n+1.);                                       //Out of range of table.
}

double bico(const int n, const int k) {                //Returns the binomial coefficientas C(n,k) as a floating-point number.
    if (n<0 || k<0 || k>n) throw("bad args in bico");
    if (n<171) return floor(0.5+factrl(n)/(factrl(k)*factrl(n-k)));
    return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
    
}

double beta(const double z, const double w) {//Returns the value of the beta function B(z,w).
    return exp(gammln(z)+gammln(w)-gammln(z+w));
}



#endif
