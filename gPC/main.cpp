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

int M=101;
int N=101;
int timesteps=100;
double leftx=-1.5;
double bottomv=-1.5;
double dt=0.01;
double dx=0.03;
double dv=0.03;

double C1(double x, double v) {
    return v;
}

double C2(double x, double v) {
    return -0.2;
}

double initial(double x, double v) {
    if (x*x+v*v<1.0 && x>=0 && v<0) {
        return 1.0;
    }
    if (x*x+v*v<1.0 && x<=0 && v>0) {
        return 1.0;
    }
    else {
        return 0.0;
    }
}

template <class T>
void march(dynamicMatrix<T>& u, const double& r) {
    dynamicMatrix<T> tmp(u);
    for (int t=0; t<timesteps; t++) {
        for (int i=1; i<u.height()-1; i++) {
            for (int j=1; j<u.width()-1; j++) {
                double x,v,c1,c2,ux,uv;
                x=leftx+i*dx;
                v=bottomv+j*dv;
                c1=C1(x,v);
                c2=C2(x,v)+r;
                if (c1>0.0) {
                    ux=(u(i,j,"read")-u(i-1,j,"read"))/dx;
                }
                else {
                    ux=(u(i+1,j,"read")-u(i,j,"read"))/dx;
                }
                if (c2>0.0) {
                    uv=(u(i,j,"read")-u(i,j-1,"read"))/dv;
                }
                else {
                    uv=(u(i,j+1,"read")-u(i,j,"read"))/dv;
                }
                tmp(i,j)=u(i,j,"read")-dt*(c1*ux+c2*uv);
            }
        }
        u=tmp;
    }
}

template <class T>
T abs(const T& a) {
    if (a>0) {
        return a;
    }
    else {
        return -a;
    }
}

template <class T>
T norm1(const dynamicMatrix<T>& u) {
    T sum=0.0;
    for (int i=0; i<u.height(); i++) {
        for (int j=0; j<u.width(); j++) {
            sum+=abs(u(i,j,"read"));
        }
    }
    return sum/(u.height()*u.width());
}



int main(int argc, const char * argv[]) {
    dynamicMatrix<double> u(M,N,0.0),f(M,N,0.0);
    for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
            double x,v;
            x=leftx+i*dx;
            v=bottomv+j*dv;
            u(i,j)=initial(x,v);
            if ((x-v-0.1)*(x-v-0.1)+(v+0.2)*(v+0.2)<1.0 && (x-v)>=0.1 && v<-0.2) {
                f(i,j)=1.0;
            }
            if ((x-v-0.1)*(x-v-0.1)+(v+0.2)*(v+0.2)<1.0 && (x-v)<=0.1 && v>-0.2) {
                f(i,j)=1.0;
            }
            else
                f(i,j)=0.0;
        }
    }
    march(u,0.0);
    cout << norm1(f-u) << endl;
    return 0;
}

