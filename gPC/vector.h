//
//  vector.h
//  gPC
//
//  Created by 马 征 on 14-2-26.
//  Copyright (c) 2014年 马 征. All rights reserved.
//

#ifndef gPC_vector_h
#define gPC_vector_h

#include <iostream>

using namespace std;

//int power(int basis, unsigned exp){
//    return exp ? basis * power(basis,exp-1) : 1;
//}  //  "basis" to the "exp"

//int max(int a, int b){return a>b ? a : b;}
//int min(int a, int b){return a<b ? a : b;}

template<class T, int N> class vector{
    T component[N];
public:
    vector(const T& = 0);
    vector(const T&a,const T&b){
        component[0] = a; component[1] = b;
    }  // constructor for 2-d vectors
    vector(const T&a,const T&b,const T&c){
        component[0] = a; component[1] = b; component[2] = c;
    }  // constructor for 3-d vectors
    vector(const vector&);
    const vector& operator=(const vector&);
    const vector& operator=(const T&);
    ~vector(){}  //  destructor
    const T& operator[](int i) const{ return component[i]; }  //ith component
    T& operator()(int i){ return component[i];
    }  //read/write ith component
    void set(int i,const T& a){ component[i] = a; }  //  change ith component
    const vector& operator+=(const vector&);
    const vector& operator-=(const vector&);
    const vector& operator*=(const T&);
    const vector& operator/=(const T&);
};


template<class T, int N>
vector<T,N>::vector(const T& a){
    for(int i = 0; i < N; i++)
        component[i] = a;
}  //  constructor

template<class T, int N>
vector<T,N>::vector(const vector<T,N>& v){
    for(int i = 0; i < N; i++)
        component[i] = v.component[i];
}  //  copy constructor

template<class T, int N>
const vector<T,N>& vector<T,N>::operator=(const vector<T,N>& v){
    if(this != &v)
        for(int i = 0; i < N; i++)
            component[i] = v.component[i];
    return *this;
}  //  assignment operator

template<class T, int N>
const vector<T,N>& vector<T,N>::operator=(const T& a){
    for(int i = 0; i < N; i++)
        component[i] = a;
    return *this;
}  //  assignment operator with a scalar argument

template<class T, int N>
const vector<T,N>& vector<T,N>::operator+=(const vector<T,N>&v){
    for(int i = 0; i < N; i++)
        component[i] += v[i];
    return *this;
}  //  adding a vector to the current vector

template<class T, int N>
const vector<T,N>& vector<T,N>::operator-=(const vector<T,N>&v){
    for(int i = 0; i < N; i++)
        component[i] -= v[i];
    return *this;
}  //  subtracting a vector from the current vector

template<class T, int N>
const vector<T,N>& vector<T,N>::operator*=(const T& a){
    for(int i = 0; i < N; i++)
        component[i] *= a;
    return *this;
}  //  multiplying the current vector by a scalar

template<class T, int N>
const vector<T,N>& vector<T,N>::operator/=(const T& a){
    for(int i = 0; i < N; i++)
        component[i] /= a;
    return *this;
}  //  multiplying the current vector by a scalar

template<class T, int N>
const vector<T,N> operator+(const vector<T,N>&u, const vector<T,N>&v){
    return vector<T,N>(u) += v;
}  //  vector plus vector

template<class T, int N>
const vector<T,N> operator-(const vector<T,N>&u, const vector<T,N>&v){
    return vector<T,N>(u) -= v;
}  //  vector minus vector

template<class T, int N>
const vector<T,N> operator*(const vector<T,N>&u, const T& a){
    return vector<T,N>(u) *= a;
}  //  vector times scalar

template<class T, int N>
const vector<T,N> operator*(const T& a, const vector<T,N>&u){
    return vector<T,N>(u) *= a;
}  //  'T' times vector

template<class T, int N>
const vector<T,N> operator/(const vector<T,N>&u, const T& a){
    return vector<T,N>(u) /= a;
}  //  vector times scalar

template<class T, int N>
const vector<T,N>& operator+(const vector<T,N>&u){
    return u;
}  //  negative of a vector

template<class T, int N>
const vector<T,N> operator-(const vector<T,N>&u){
    return vector<T,N>(u) *= -1;
}  //  negative of a vector

template<class T, int N>
const T operator*(const vector<T,N>&u, const vector<T,N>&v){
    T sum = 0;
    for(int i = 0; i < N; i++)
        sum += u[i] * +v[i];
    return sum;
}  //  vector times vector (inner product)

template<class T, int N>
const T squaredNorm(const vector<T,N>&u){
    return u*u;
}  //  squared l2 norm

template<class T, int N>
const T l2norm(const vector<T,N>&u){
    return sqrt(u*u);
}  //  l2 norm

template<class T, int N>
ostream& operator <<(ostream& out, const vector<T,N>&v){
    out << "(" ;
    for (int i = 0; i< N; i++){
        if(i != N-1)
            out << "v[" << i << "]=" << v[i] << ", ";
        else
            out << "v[" << i << "]=" << v[i] << ")";
    }
    return out;
}  //  printing a vector


#endif
