//
//  main.cpp
//  gPC
//
//  Created by 马 征 on 14-2-7.
//  Copyright (c) 2014年 马 征. All rights reserved.
//

#include <iostream>
#include "dynamicMatrix.h"
#include "list.h"

template<class T> class difference2:public list<dynamicMatrix<T> >{
public:
    difference2(int = 0,int = 0,const T& = 0,const T& = 0,const T& = 0,const T& = 0,
                const T& = 0,const T& = 0,const T& = 0,const T& = 0,const T& = 0);
    const difference2<T>& operator+=(const difference2<T>&);
    const difference2<T>& operator-=(const difference2<T>&);
    const difference2& operator*=(const T&);
    T& operator()(int i,int j,int k,int l){
        return (*this->item[(k-i+1)*3+l-j+1])(i,j);
    }  //  (i,j,k,l)th element
    const T& operator()(int i,int j,int k,int l, const char*) const{
        return (*this->item[(k-i+1)*3+l-j+1])(i,j,"read");
    }  //  (i,j,k,l)th element (read only)
    int width() const{ return this->item[0]->width(); }
    int height() const{ return this->item[0]->height(); }
};
template<class T>
difference2<T>::difference2(int m,int n,const T&a,const T&b,
                            const T&c,const T&d,const T&e,const T&f,const T&g,
                            const T&h,const T&i){
    this->number = 9;
    this->item = new dynamicMatrix<T>*[9];
    this->item[0] = new dynamicMatrix<T>(m,n,a);
    this->item[1] = new dynamicMatrix<T>(m,n,b);
    this->item[2] = new dynamicMatrix<T>(m,n,c);
    this->item[3] = new dynamicMatrix<T>(m,n,d);
    this->item[4] = new dynamicMatrix<T>(m,n,e);
    this->item[5] = new dynamicMatrix<T>(m,n,f);
    this->item[6] = new dynamicMatrix<T>(m,n,g);
    this->item[7] = new dynamicMatrix<T>(m,n,h);
    this->item[8] = new dynamicMatrix<T>(m,n,i);
}

template<class T>
const difference2<T>& difference2<T>::operator+=(const difference2<T>&d){
    for(int i=0; i<this->number; i++)
        *this->item[i] += d[i];
    return *this;
}  //  adding another difference2 to the current one

template<class T>
const difference2<T>& difference2<T>::operator-=(const difference2<T>&d){
    for(int i=0; i<this->number; i++)
        *this->item[i] -= d[i];
    return *this;
}  //  subtracting another difference2 from the current one

template<class T>
const difference2<T>& difference2<T>::operator*=(const T&t){
    for(int i=0; i<this->number; i++)
        *this->item[i] *= t;
    return *this;
}  //  multiplying the difference2 operator by a scalar T

template<class T>
const difference2<T> operator+(const difference2<T>&d1,
                               const difference2<T>&d2){
    return difference2<T>(d1) += d2;
}  //  addition of two difference2s

template<class T>
const difference2<T> operator-(const difference2<T>&d1,
                               const difference2<T>&d2){
    return difference2<T>(d1) -= d2;
}  //  subtraction of two difference2s

template<class T>
const difference2<T> operator*(const T&t, const difference2<T>&d){
    return difference2<T>(d) *= t;
}  //  scalar times difference2

template<class T>
const difference2<T> operator*(const difference2<T>&d, const T&t){
    return difference2<T>(d) *= t;
}  //  difference2 times scalar

template<class T>
const dynamicMatrix<T> operator*(const difference2<T>&d, const dynamicMatrix<T>&v){
    dynamicMatrix<T> dv(v.height(),v.width(),0);
    for(int i=0; i<v.height(); i++)
        for(int j=0; j<v.width(); j++)
            for(int k=max(0,i-1); k<=min(v.height()-1,i+1); k++)
                for(int l=max(0,j-1); l<=min(v.width()-1,j+1); l++)
                    dv(i,j) += d(i,j,k,l,"read")*v(k,l,"read");
    return dv;
}  //  difference2 times dynamicMatrix

template<class T>
const dynamicMatrix<T> operator/(const dynamicMatrix<T>&f, const difference2<T>&d){
    dynamicMatrix<T> x(f);
    for(int iteration=0; iteration < 100; iteration++)
        for(int i=0; i<f.height(); i++)
            for(int j=0; j<f.width(); j++){
                double residual = f(i,j,"read");
                for(int k=max(0,i-1); k<=min(f.height()-1,i+1); k++)
                    for(int l=max(0,j-1); l<=min(f.width()-1,j+1); l++)
                        residual -= d(i,j,k,l,"read")*x(k,l,"read");
                x(i,j) += residual/d(i,j,i,j,"read");
            }
    return x;
}  //  solving d*x=f approximately by 10 GS iterations

template<class T> class xytGrid:public list<dynamicMatrix<T> >{
public:
    xytGrid(int = 0,int = 0,int = 0,const T& = 0);
    int timeSteps() const{ return this->size(); }  // number of time steps
    int width() const{ return this->item[0]->width(); }  // width of grid
    int height() const{ return this->item[0]->height(); }  // height of grid
    dynamicMatrix<T>& operator()(int i){if(this->item[i])return *this->item[i];else assert(false);}//ith plane
    T& operator()(int i,int j,int k){return (*this->item[i])(j,k); }// (i,j,k)th comp.
};

template<class T>
xytGrid<T>::xytGrid(int m,int n,int l,const T&a){
    this->number = m;
    this->item = m ? new dynamicMatrix<T>*[m] : 0;
    for(int i=0; i<m; i++)
        this->item[i] = new dynamicMatrix<T>(n,l,a);
}  //  constructor

double F(double, double, double){return 0.;}
double C1(double, double, double){return 0.;}
double C2(double, double, double){return 0.;}
double Alpha(double, double, double){return 0.;}
double G(double, double, double){return 0.;}
double Initial(double x, double y){return (1.-x*x)*(1.-y*y);}
const double Epsilon=1.;

template<class T>
void convDif(difference2<T>&d,dynamicMatrix<T>&f,double hx,
             double hy,double deltaT,double t){
    for(int k=0; k<d.height(); k++)
        for(int j=0; j<d.width(); j++){
            if(t>deltaT/2)f(k,j)=F(j*hx,k*hy,t-deltaT/2);
            double c1=C1(j*hx,k*hy,t)/hx;
            if(c1>0.){
                d(k,j,k,j)=c1;
                d(k,j,k,j-1)=-c1;
            }
            else{
                d(k,j,k,j)=-c1;
                d(k,j,k,j+1)=c1;
            }
            double c2=C2(j*hx,k*hy,t)/hy;
            if(c2>0.){
                d(k,j,k,j)+=c2;
                d(k,j,k-1,j)=-c2;
            }
            else{
                d(k,j,k,j)-=c2;
                d(k,j,k+1,j)=c2;
            }
        }
    d += Epsilon * difference2<T>(d.height(),d.width(),
                                  0,-1/hy/hy,0,-1/hx/hx,2/hx/hx+2/hy/hy,-1/hx/hx,0,-1/hy/hy,0);
    for(int k=0; k<d.height(); k++){
        d(k,0,k,0) += d(k,0,k,-1);
        d(k,0,k,0) -= d(k,0,k,-1) * hx * Alpha(0.,k*hy,t);
        if(t>deltaT/2){
            f(k,0) -= d(k,0,k,-1) * hx * G(0,k*hy,t-deltaT/2);
            f(k,d.width()-1) -= d(k,d.width()-1,k,d.width())
            * G(d.width()*hx,k*hy,t-deltaT/2);
        }
    }
    for(int j=0; j<d.width(); j++){
        d(0,j,0,j) += d(0,j,-1,j);
        d(0,j,0,j) -= d(0,j,-1,j) * hy * Alpha(j*hx,0.,t);
        if(t>deltaT/2){
            f(0,j) -= d(0,j,-1,j) * hy * G(j*hx,0.,t-deltaT/2);
            f(d.height()-1,j) -= d(d.height()-1,j,d.height(),j)
            * G(j*hx,d.height()*hy,t-deltaT/2);
        }
    }
}  //  set the convection-diffusion matrix and right-hand side

template<class T>
void march(xytGrid<T>&g, double hx, double hy, double deltaT){
    difference2<T> I(g.height(),g.width(),0,0,0,0,1,0,0,0,0);
    for(int k=0; k<g.height(); k++)
        for(int j=0; j<g.width(); j++)
            g(0,k,j) = Initial(j*hx,k*hy);
    dynamicMatrix<T> f(g.height(),g.width());
    difference2<T> previous(g.height(),g.width());
    convDif(previous,f,hx,hy,deltaT,0);
    for(int time=1; time < g.timeSteps(); time++){
        difference2<T> current(g.height(),g.width());
        convDif(current,f,hx,hy,deltaT,time*deltaT);
        g(time) = ((I - 0.5 * deltaT * previous) * g[time-1] + deltaT * f)
        / (I + 0.5 * deltaT * current);
        previous = current;
    }
    print(g[g.timeSteps()-1]);
}  //  semi-implicit time marching

class domain2{
    xytGrid<double> g;
    double Time;
    double Width;
    double Height;
public:
    domain2(double T, double Lx, double Ly, double accuracy):
    g((int)(T/accuracy)+1,(int)(Ly/accuracy)+1,(int)(Lx/accuracy)+1),
    Time(T),Width(Lx),Height(Ly){}  //  constructor
    void solveConvDif(){march(g,Width/g.width(),Height/g.height(),
                              Time/g.timeSteps());}  //  solve convection-diffusion equation
};

int main(int argc, const char * argv[])
{
    domain2 D2(10.,.3,1.,.1);
    D2.solveConvDif();
    return 0;
}

