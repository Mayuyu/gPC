//
//  dynamicVector.h
//  gPC
//
//  Created by 马 征 on 14-2-7.
//  Copyright (c) 2014年 马 征. All rights reserved.
//

#ifndef gPC_dynamicVector_h
#define gPC_dynamicVector_h

#include<stdio.h>
#include<math.h>

void print(double d){
    printf("%f  ",d);
}

int power(int basis, unsigned exp){
    return exp ? basis * power(basis,exp-1) : 1;
}  //  "basis" to the "exp"

int max(int a, int b){return a>b ? a : b;}
int min(int a, int b){return a<b ? a : b;}

template<class T> class dynamicVector{
protected:
    int dimension;
    T* component;
public:
    dynamicVector(int = 0, const T& = 0);
    dynamicVector(const dynamicVector&);
    const dynamicVector& operator=(const dynamicVector&);
    const dynamicVector& operator=(const T&);
    ~dynamicVector(){delete [] component;}  //  destructor
    int dim() const{ return dimension; }  //  return the dimension
    T& operator()(int i){ return component[i]; }  //read/write ith component
    const T& operator[](int i) const{ return component[i]; }  //read only
    const dynamicVector& operator+=(const dynamicVector&);
    const dynamicVector& operator-=(const dynamicVector&);
    const dynamicVector& operator*=(const T&);
    const dynamicVector& operator/=(const T&);
};

template<class T>
dynamicVector<T>::dynamicVector(int dim, const T& a) : dimension(dim),
component(dim ? new T[dim] : 0){
    for(int i = 0; i < dim; i++)
        component[i] = a;
}  //  constructor

template<class T>
dynamicVector<T>::dynamicVector(const dynamicVector<T>& v) : dimension(v.dimension),
component(v.dimension ? new T[v.dimension] : 0){
    for(int i = 0; i < v.dimension; i++)
        component[i] = v.component[i];
}  //  copy constructor

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator=(const dynamicVector<T>& v){
    if(this != &v){
        if(dimension > v.dimension)
            delete [] (component + v.dimension);
        if(dimension < v.dimension){
            delete [] component;
            component = new T[v.dimension];
        }
        for(int i = 0; i < v.dimension; i++)
            component[i] = v.component[i];
        dimension = v.dimension;
    }
    return *this;
}  //  assignment operator

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator=(const T& a){
    for(int i = 0; i < dimension; i++)
        component[i] = a;
    return *this;
}  //  assignment operator with a scalar argument

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator+=( const dynamicVector<T>&v){
    for(int i = 0; i < dimension; i++)
        component[i] += v[i];
    return *this;
}  //  adding a dynamicVector to the current dynamicVector

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator-=( const dynamicVector<T>&v){
    for(int i = 0; i < dimension; i++)
        component[i] -= v[i];
    return *this;
}  //  subtracting a dynamicVector from the current dynamicVector

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator*=(const T& a){
    for(int i = 0; i < dimension; i++)
        component[i] *= a;
    return *this;
}  //  multiplying the current dynamicVector by a scalar

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator/=(const T& a){
    for(int i = 0; i < dimension; i++)
        component[i] /= a;
    return *this;
}  //  dividing the current dynamicVector by a scalar

template<class T>
const dynamicVector<T> operator+(const dynamicVector<T>&u, const dynamicVector<T>&v){
    return dynamicVector<T>(u) += v;
}  //  dynamicVector plus dynamicVector

template<class T>
const dynamicVector<T> operator-(const dynamicVector<T>&u, const dynamicVector<T>&v){
    return dynamicVector<T>(u) -= v;
}  //  dynamicVector minus dynamicVector

template<class T>
const dynamicVector<T> operator*(const dynamicVector<T>&u, const T& a){
    return dynamicVector<T>(u) *= a;
}  //  dynamicVector times scalar

template<class T>
const dynamicVector<T> operator*(const T& a, const dynamicVector<T>&u){
    return dynamicVector<T>(u) *= a;
}  //  T times dynamicVector

template<class T>
const dynamicVector<T> operator/(const dynamicVector<T>&u, const T& a){
    return dynamicVector<T>(u) /= a;
}  //  dynamicVector divided by scalar

template<class T>
const dynamicVector<T> operator-(const dynamicVector<T>&u){
    return dynamicVector<T>(u) *= -1.;
}  //  negative of a dynamicVector

template<class T>
T operator*(const dynamicVector<T>&u, const dynamicVector<T>&v){
    T sum = 0;
    for(int i = 0; i < u.dim(); i++)
        sum += u[i] * +v[i];
    return sum;
}  //  dynamicVector times dynamicVector (inner product)

template<class T>
void print(const dynamicVector<T>&v){
    for(int i = 0;i < v.dim(); i++)
        printf("v[%d]=%f;  ",i,(double)v[i]);
    printf("\n");
}  //  printing a dynamicVector

#endif
