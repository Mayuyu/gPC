//
//  main.cpp
//  gPC
//
//  Created by 马 征 on 14-2-7.
//  Copyright (c) 2014年 马 征. All rights reserved.
//

#include <iostream>
#include <ctime>
#include "dynamicMatrix.h"

using namespace std;

const double epsilon=0.00000001;

int M=102;
int N=102;
int timesteps=100;
double leftx=-1.515;
double bottomv=-1.515;
double dt=0.01;
double dx=0.03;
double dv=0.03;

double ran() {
    return (double)(rand()/(double)(RAND_MAX/2))-1.0;
}

double C1(double x, double v) {
    return v;
}

double VR(double x) {
    if (x<-epsilon) {
        return 0.2;
    }
    else {
        return 0.0;
    }
}

double VL(double x) {
    if (x>epsilon) {
        return 0.0;
    }
    else {
        return 0.2;
    }
}

int locate(double v) {
    for (int j=0; j<N; j++) {
        if (bottomv+j*dv-v<=epsilon && bottomv+(j+1)*dv-v>epsilon) {
            return j;
        }
    }
    return -1;
}

double initial(double x, double v) {
    if (x*x+v*v<1.0 && x>=0.0 && v<0.0) {
        return 1.0;
    }
    else if (x*x+v*v<1.0 && x<=0.0 && v>0.0) {
        return 1.0;
    }
    else {
        return 0.0;
    }
}

template <class T>
void initial(dynamicMatrix<T>& u) {
    for (int i=0; i<u.height(); i++) {
        for (int j=0; j<u.width(); j++) {
            u(i,j)=initial(leftx+i*dx,bottomv+j*dv);
        }
    }
}

template <class T>
void initial(dynamicMatrix<dynamicVector<T> >& u) {
    for (int i=0; i<u.height(); i++) {
        for (int j=0; j<u.width(); j++) {
            u(i,j)(0)=initial(leftx+i*dx,bottomv+j*dv);
        }
    }
}

template <class T>
void config(dynamicMatrix<T>& L) {
    for (int i=0; i<L.height(); i++) {
        if (i!=0) {
            L(i,i-1)=-(double)i/sqrt(4.0*i*i-1);
        }
        if (i!=L.height()-1) {
            L(i,i+1)=-(double)(i+1)/sqrt((2.0*i+1)*(2.0*i+3));
        }
    }
}

template <class T>
void march(dynamicMatrix<T>& u, const double& r) {
    dynamicMatrix<T> tmp(u);
    for (int i=1; i<u.height()-1; i++) {
        for (int j=1; j<u.width()-1; j++) {
            double x,v,c1,c2,Vl,Vr,u1,u2,ux,uv;
            v=bottomv+j*dv;
            x=leftx+i*dx;
            c1=C1(x,v);
            if (c1>0) {
                u1=u(i,j,"read");
                Vl=VL(x-0.5*dx);
                Vr=VR(x-0.5*dx);
                if ((Vl-Vr)>epsilon) {
                    if ((v*v+2.0*(Vr-Vl))>epsilon) {
                        double vl=sqrt(v*v+2.0*(Vr-Vl));
                        int k=locate(vl);
                        if (k==-1) {
                            u2=0.0;
                        }
                        else {
                            u2=(bottomv+(k+1)*dv-vl)*u(i-1,k,"read")/dv+(vl-bottomv-k*dv)*u(i-1,k+1,"read")/dv;
                        }
                    }
                    else {
                        u2=u(i,locate(-v),"read");
                    }
                }
                else if ((Vl-Vr)<-epsilon) {
                    double vl=sqrt(v*v+2.0*(Vr-Vl));
                    int k=locate(vl);
                    if (k==-1) {
                        u2=0.0;
                    }
                    else {
                        u2=(bottomv+(k+1)*dv-vl)*u(i-1,k,"read")/dv+(vl-bottomv-k*dv)*u(i-1,k+1,"read")/dv;
                    }
                }
                else {
                    u2=u(i-1,j,"read");
                }
            }
            else if (c1<0) {
                u2=u(i,j,"read");
                Vl=VL(x+0.5*dx);
                Vr=VR(x+0.5*dx);
                if ((Vl-Vr)<-epsilon) {
                    if ((v*v+2.0*(Vl-Vr))>epsilon) {
                        double vr=-sqrt(v*v+2.0*(Vl-Vr));
                        int k=locate(vr);
                        if (k==-1) {
                            u1=0.0;
                        }
                        else {
                            u1=(bottomv+(k+1)*dv-vr)*u(i+1,k,"read")/dv+(vr-bottomv-k*dv)*u(i+1,k+1,"read")/dv;
                        }
                    }
                    else {
                        u1=u(i,locate(-v),"read");
                    }
                }
                else if ((Vl-Vr)>epsilon) {
                    double vr=-sqrt(v*v+2.0*(Vl-Vr));
                    int k=locate(vr);
                    if (k==-1) {
                        u1=0.0;
                    }
                    else {
                        u1=(bottomv+(k+1)*dv-vr)*u(i+1,k,"read")/dv+(vr-bottomv-k*dv)*u(i+1,k+1,"read")/dv;
                    }
                }
                else {
                    u1=u(i+1,j,"read");
                }
            }
            else {
                u1=u2=0.0;
            }
            ux=(u1-u2)/dx;
            c2=(VR(x-0.5*dx)-VL(x+0.5*dx))/dx-0.1*r;
            if (c2>epsilon) {
                uv=(u(i,j,"read")-u(i,j-1,"read"))/dv;
            }
            else if (c2<-epsilon) {
                uv=(u(i,j+1,"read")-u(i,j,"read"))/dv;
            }
            else {
                uv=0.0;
            }
            tmp(i,j)=u(i,j,"read")-dt*(c1*ux+c2*uv);
        }
    }
    u=tmp;
}


template <class T>
void solve(dynamicMatrix<T>& u, const double& r) {
    initial(u);
    for (int t=0; t<timesteps; t++) {
        march(u, r);
    }
}

template <class T>
void march(dynamicMatrix<dynamicVector<T> >& u) {
    int order=u(0,0,"read").dim();
    dynamicMatrix<dynamicVector<T> > tmp(u);
    dynamicMatrix<T> L(order,order,0.0), I(order,order,0.0);
    for (int i=0; i<I.width(); i++) {
        I(i,i)=1.0;
    }
    config(L);
    initial(u);
    for (int t=0; t<timesteps; t++) {
        for (int i=1; i<u.height()-1; i++) {
            for (int j=1; j<u.width()-1; j++) {
                double x,v,c1,c2,Vl,Vr;
                dynamicVector<T> u1(order,0.0),u2(order,0.0),ux(order,0.0),uv(order,0.0);
                v=bottomv+j*dv;
                x=leftx+i*dx;
                c1=C1(x,v);
                if (c1>0) {
                    u1=u(i,j,"read");
                    Vl=VL(x-0.5*dx);
                    Vr=VR(x-0.5*dx);
                    if ((Vl-Vr)>epsilon) {
                        if ((v*v+2.0*(Vr-Vl))>epsilon) {
                            double vl=sqrt(v*v+2.0*(Vr-Vl));
                            int k=locate(vl);
                            if (k==-1) {
                                u2=0.0;
                            }
                            else {
                                u2=(bottomv+(k+1)*dv-vl)*u(i-1,k,"read")/dv+(vl-bottomv-k*dv)*u(i-1,k+1,"read")/dv;
                            }
                        }
                        else {
                            u2=u(i,locate(-v),"read");
                        }
                    }
                    else if ((Vl-Vr)<-epsilon) {
                        double vl=sqrt(v*v+2.0*(Vr-Vl));
                        int k=locate(vl);
                        if (k==-1) {
                            u2=0.0;
                        }
                        else {
                            u2=(bottomv+(k+1)*dv-vl)*u(i-1,k,"read")/dv+(vl-bottomv-k*dv)*u(i-1,k+1,"read")/dv;
                        }
                    }
                    else {
                        u2=u(i-1,j,"read");
                    }
                }
                else if (c1<0) {
                    u2=u(i,j,"read");
                    Vl=VL(x+0.5*dx);
                    Vr=VR(x+0.5*dx);
                    if ((Vl-Vr)<-epsilon) {
                        if ((v*v+2.0*(Vl-Vr))>epsilon) {
                            double vr=-sqrt(v*v+2.0*(Vl-Vr));
                            int k=locate(vr);
                            if (k==-1) {
                                u1=0.0;
                            }
                            else {
                                u1=(bottomv+(k+1)*dv-vr)*u(i+1,k,"read")/dv+(vr-bottomv-k*dv)*u(i+1,k+1,"read")/dv;
                            }
                        }
                        else {
                            u1=u(i,locate(-v),"read");
                        }
                    }
                    else if ((Vl-Vr)>epsilon) {
                        double vr=-sqrt(v*v+2.0*(Vl-Vr));
                        int k=locate(vr);
                        if (k==-1) {
                            u1=0.0;
                        }
                        else {
                            u1=(bottomv+(k+1)*dv-vr)*u(i+1,k,"read")/dv+(vr-bottomv-k*dv)*u(i+1,k+1,"read")/dv;
                        }
                    }
                    else {
                        u1=u(i+1,j,"read");
                    }
                }
                else {
                    u1=0;;
                    u2=0;
                }
                ux=(u1-u2)/dx;
                c2=(VR(x-0.5*dx)-VL(x+0.5*dx))/dx;
                if (c2>0.0) {
                    uv=(u(i,j,"read")-u(i,j-1,"read"))/dv;
                }
                else {
                    uv=(u(i,j+1,"read")-u(i,j,"read"))/dv;
                }
                tmp(i,j)=u(i,j,"read")-dt*(c1*I*ux+c2*I*uv+0.1*L*uv);
            }
        }
        u=tmp;
    }
}


template <class T>
void MC(dynamicMatrix<T>& u, const int m) {
    dynamicMatrix<T> tmp(u);
    for (int i=0; i<m; i++) {
        march(tmp, ran());
        u+=tmp;
    }
    u=((double)(1.0/m))*u;
}

template <class T>
T norm1(const dynamicMatrix<T>& u) {
    T sum=0.0;
    int count=0;
    for (int i=0; i<u.height(); i++) {
        for (int j=0; j<u.width(); j++) {
            if (i<=(u.height()/2-1) && bottomv+j*dv+0.5*dv<-sqrt(bottomv*bottomv-0.4)) {
                sum+=0.0;
            }
            else {
                sum+=fabs(u(i,j,"read"));
                count++;
            }
        }
    }
    return sum/count;
}

double exact(double x, double v) {
    if (x>=0 && v<sqrt(0.4) && x<v) {
        return 1.0;
    }
    else if (x>=0 && x<1 && v<0 && v>0.5*(x-sqrt(2.0-x*x)) ) {
        return 1.0;
    }
    else if (x<=0 && v<x && v>-sqrt(0.6) && x<(1.0-sqrt(0.6-v*v)/sqrt(v*v+0.4))*v) {
        return 1.0;
    }
    else if (x<=0 && v>0 && x>-1 && v<0.5*(x+sqrt(2.0-x*x))) {
        return 1.0;
    }
    else if (x>=0 && v>sqrt(0.4) && v>x && v<sqrt(1.4) && x>(1.0-sqrt(1.4-v*v)/sqrt(v*v-0.4))*v) {
        return 1.0;
    }
    else {
        return 0.0;
    }
}

template <class T>
void exact(dynamicMatrix<T>& u,const T& r) {
    for (int i=0; i<u.height(); i++) {
        for (int j=0; j<u.width(); j++) {
            double x=leftx+i*dx;
            double v=bottomv+j*dv;
            u(i,j)=exact(x, v);
        }
    }
}


int main(int argc, const char * argv[]) {
    srand((unsigned)time(NULL));
    dynamicMatrix<double> u(M,N,0.0),uL0(M,N,0.0),test(M,N,0.0),f(M,N,0.0);
//    for (int k=3; k<13; k++) {
//        dynamicVector<double> I(k,0.0);
//        dynamicMatrix<dynamicVector<double> > uL(M,N,I);
//        march(uL);
//        for (int i=0; i<uL.height(); i++) {
//            for (int j=0; j<uL.width(); j++) {
//                uL0(i,j)=uL(i,j,"read")[0];
//            }
//        }
//        cout << norm1(uL0-test) << endl;
//        test=uL0;
//        
//    }
    dynamicVector<double> I(5,0.0);
    dynamicMatrix<dynamicVector<double> > uL(M,N,I);
//    march(uL);
//    for (int i=0; i<uL.height(); i++) {
//        for (int j=0; j<uL.width(); j++) {
//            uL0(i,j)=uL(i,j,"read")[0];
//        }
//    }
//
//    for (int k=1; k<30; k++) {
//        MC(u, 10000);
//        cout << norm1(u-uL0) << endl;
//    }
//    return 0;
    solve(u, 0.0);
    exact(test, 0.0);
    cout << norm1(u-test) << endl;
    return 0;
}

























/*******************************
 
    smooth potential
 
 *****************************/

//int M=501;
//int N=501;
//int timesteps=500;
//double leftx=-1.5;
//double bottomv=-1.5;
//double dt=0.002;
//double dx=0.006;
//double dv=0.006;
//
//double ran() {
//    return (double)(rand()/(double)(RAND_MAX/2))-1.0;
//}
//
//double C1(double x, double v) {
//    return v;
//}
//
//double C2(double x, double v) {
//    return -0.2;
//}
//
//double initial(double x, double v) {
//    if (x*x+v*v<1.0 && x>=0.0 && v<0.0) {
//        return 1.0;
//    }
//    if (x*x+v*v<1.0 && x<=0.0 && v>0.0) {
//        return 1.0;
//    }
//    else {
//        return 0.0;
//    }
//}
//
//template <class T>
//void initial(dynamicMatrix<T>& u) {
//    for (int i=0; i<u.height(); i++) {
//        for (int j=0; j<u.width(); j++) {
//            u(i,j)=initial(leftx+i*dx,bottomv+j*dv);
//        }
//    }
//}
//
//template <class T>
//void initial(dynamicMatrix<dynamicVector<T> >& u) {
//    for (int i=0; i<u.height(); i++) {
//        for (int j=0; j<u.width(); j++) {
//            u(i,j)(0)=initial(leftx+i*dx,bottomv+j*dv);
//        }
//    }
//}
//
//template <class T>
//void config(dynamicMatrix<T>& L) {
//    for (int i=0; i<L.height(); i++) {
//        if (i!=0) {
//            L(i,i-1)=-(double)i/sqrt(4.0*i*i-1);
//        }
//        if (i!=L.height()-1) {
//            L(i,i+1)=-(double)(i+1)/sqrt((2.0*i+1)*(2.0*i+3));
//        }
//    }
//}
//
//template <class T>
//void march(dynamicMatrix<T>& u, const double& r) {
//    dynamicMatrix<T> tmp(u);
//    initial(u);
//    for (int t=0; t<timesteps; t++) {
//        for (int i=1; i<u.height()-1; i++) {
//            for (int j=1; j<u.width()-1; j++) {
//                double x,v,c1,c2,ux,uv;
//                x=leftx+i*dx;
//                v=bottomv+j*dv;
//                c1=C1(x,v);
//                c2=C2(x,v)-0.1*r;
//                if (c1>0.0) {
//                    ux=(u(i,j,"read")-u(i-1,j,"read"))/dx;
//                }
//                else {
//                    ux=(u(i+1,j,"read")-u(i,j,"read"))/dx;
//                }
//                if (c2>0.0) {
//                    uv=(u(i,j,"read")-u(i,j-1,"read"))/dv;
//                }
//                else {
//                    uv=(u(i,j+1,"read")-u(i,j,"read"))/dv;
//                }
//                tmp(i,j)=u(i,j,"read")-dt*(c1*ux+c2*uv);
//            }
//        }
//        u=tmp;
//    }
//}
//
//template <class T>
//void march(dynamicMatrix<dynamicVector<T> >& u) {
//    int order=u(0,0,"read").dim();
//    dynamicMatrix<dynamicVector<T> > tmp(u);
//    dynamicMatrix<T> L(order,order,0.0), I(order,order,0.0);
//    for (int i=0; i<I.width(); i++) {
//        I(i,i)=1.0;
//    }
//    config(L);
//    initial(u);
//    for (int t=0; t<timesteps; t++) {
//        for (int i=1; i<u.height()-1; i++) {
//            for (int j=1; j<u.width()-1; j++) {
//                double x,v,c1,c2;
//                dynamicVector<T> ux,uv;
//                x=leftx+i*dx;
//                v=bottomv+j*dv;
//                c1=C1(x,v);
//                c2=C2(x,v);
//                if (c1>0.0) {
//                    ux=(u(i,j,"read")-u(i-1,j,"read"))/dx;
//                }
//                else {
//                    ux=(u(i+1,j,"read")-u(i,j,"read"))/dx;
//                }
//                if (c2>0.0) {
//                    uv=(u(i,j,"read")-u(i,j-1,"read"))/dv;
//                }
//                else {
//                    uv=(u(i,j+1,"read")-u(i,j,"read"))/dv;
//                }
//                tmp(i,j)=u(i,j,"read")-dt*(c1*I*ux+c2*I*uv+0.1*L*uv);
//            }
//        }
//        u=tmp;
//    }
//}
//
//
//template <class T>
//void MC(dynamicMatrix<T>& u, const int m) {
//    dynamicMatrix<T> tmp(u);
//    for (int i=0; i<m; i++) {
//        march(tmp, ran());
//        u+=tmp;
//    }
//    u=((double)(1.0/m))*u;
//}
//
//template <class T>
//T abs(const T& a) {
//    if (a>0) {
//        return a;
//    }
//    else {
//        return -a;
//    }
//}
//
//template <class T>
//T norm1(const dynamicMatrix<T>& u) {
//    T sum=0.0;
//    for (int i=0; i<u.height(); i++) {
//        for (int j=0; j<u.width(); j++) {
//            sum+=fabs(u(i,j,"read"));
//        }
//    }
//    return sum/(u.height()*u.width());
//}
//
//template <class T>
//void exact(dynamicMatrix<T>& u,const T& r) {
//    for (int i=0; i<u.height(); i++) {
//        for (int j=0; j<u.width(); j++) {
//            double x=leftx+i*dx;
//            double v=bottomv+j*dv;
//            u(i,j)=initial(x-v-0.1-0.05*r, v+0.2+0.1*r);
//        }
//    }
//}
//
//
//int main(int argc, const char * argv[]) {
//    srand((unsigned)time(NULL));
//    dynamicVector<double> I(4,0.0);
//    dynamicMatrix<double> u(M,N,0.0),uL0(M,N,0.0),test(4,4,0.0),f(M,N,0.0);
//    dynamicMatrix<dynamicVector<double> > uL(M,N,I);
//    march(u, 0.0);
//    exact(f, 0.0);
//    cout << norm1(u-f) << endl;
////    march(uL);
////    for (int i=0; i<uL.height(); i++) {
////        for (int j=0; j<uL.width(); j++) {
////            uL0(i,j)=uL(i,j,"read")[0];
////        }
////    }
////    MC(u, 100);
////    cout << norm1(uL0-u) << endl;
//    return 0;
//}

