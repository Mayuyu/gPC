//
//  main.cpp
//  gPC
//
//  Created by 马 征 on 14-2-7.
//  Copyright (c) 2014年 马 征. All rights reserved.
//

#include <iostream>
#include "dynamicMatrix.h"

using namespace std;

int M=21;
int N=21;
double leftx=-1.0;
double bottomv=-1.0;
double dt=0.05;
double dx=0.1;
double dv=0.1;

double C1(double x, double v) {
    return v;
}

double C2(double x, double v) {
    return -1.0;
}

double initial(double x, double v) {
    if (x+v<1.0) {
        return 1.0;
    }
    else {
        return 0.0;
    }
}

template <class T>
void march(dynamicMatrix<T>& u) {
    dynamicMatrix<T> tmp=u;
    for (int i=1; i<u.height()-1; i++) {
        for (int j=1; j<u.width()-1; j++) {
            double x,v,c1,c2,ux,uv;
            x=leftx+i*dx;
            v=bottomv+j*dv;
            c1=C1(x,v);
            c2=C2(x,v);
            if (c1>0.0) {
                ux=(tmp(i,j,"read")-tmp(i-1,j,"read"))/dx;
            }
            else {
                ux=(tmp(i+1,j,"read")-tmp(i,j,"read"))/dx;
            }
            if (c2>0.0) {
                uv=(tmp(i,j,"read")-tmp(i,j-1,"read"))/dv;
            }
            else {
                uv=(tmp(i,j+1,"read")-tmp(i,j,"read"))/dv;
            }
            u(i,j)=tmp(i,j)-dt*(c1*ux+c2*uv);
        }
    }
}


int main(int argc, const char * argv[]) {
    dynamicMatrix<double> u(M,N,0.0);
    for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
            u(i,j)=initial(leftx+i*dx, bottomv+j*dv);
        }
    }
//    for (int t=0; t<1; t++) {
//        march(u);
//    }
    print(u);
    return 0;
}

