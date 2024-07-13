#ifndef ELLIPTICSOLVER_H
#define ELLIPTICSOLVER_H

#include <iostream>
#include <fstream>
#include "math_matriz.h"

using namespace std;

template <class T>
class EllipticSolver
{
private:
    /* data */
    void buildMatrixA();
    void buildVectorB(double (*f)(double, double), double (*g1)(double), double (*g2)(double), double (*g3)(double), double (*g4)(double));
public:
    Matriz<T> A;
    T xl, xr, yb, yt, dx, dy;
    int M, N;
    EllipticSolver(T xl, T xr, T yb, T yt, int M, int N); 
    ~EllipticSolver();
    void solve();
};

template <class T>
EllipticSolver<T>::EllipticSolver(T xl, T xr, T yb, T yt, int M, int N)
{
    this->xl = xl;
    this->xr = xr;
    this->yb = yb;
    this->yt = yt;
    this->M = M;
    this->N = N;
    this->dx = (xr - xl)/M;
    this->dy = (yt - yb)/N;
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;

    int m = M + 1;
    int n = N + 1;
    vector<T> x(m);
    vector<T> y(n);
    
    for(int i = 0; i < m; i++) {
        x[i] = xl + i * dx;
        //cout << x[i] << " ";
    }
    //cout << endl;
    for(int j = 0; j < n; j++) {
        y[j] = yb + j * dy;
        //cout << y[j] << " ";
    }
    //cout << endl;
}

template <class T>
EllipticSolver<T>::~EllipticSolver()
{
}

template <class T>
void EllipticSolver<T>::buildMatrixA() {
    int m = M + 1;
    int n = N + 1;
    A = Matriz<T>(m*n, m*n); // ojo, n representa y, m representa x
    cout << A << endl;
}

template <class T>
void EllipticSolver<T>::solve() {
    buildMatrixA();
}

#endif // ELLIPTICSOLVER_H