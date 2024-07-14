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
    Matriz<T> L;
    Matriz<T> D;
    Matriz<T> U;
    void buildSystem(double (*f)(double, double), double (*g1)(double), double (*g2)(double), double (*g3)(double), double (*g4)(double));
    void SOR();
public:
    Matriz<T> A;
    vector<T> b;
    T xl, xr, yb, yt, dx, dy;
    int M, N;
    EllipticSolver(T xl, T xr, T yb, T yt, T dx, T dy); 
    ~EllipticSolver();
    void solve(double (*f)(double, double), double (*g1)(double), double (*g2)(double), double (*g3)(double), double (*g4)(double));
};

template <class T>
EllipticSolver<T>::EllipticSolver(T xl, T xr, T yb, T yt, T dx, T dy)
{
    this->xl = xl;
    this->xr = xr;
    this->yb = yb;
    this->yt = yt;
    this->dx = dx;
    this->dy = dy;
    this->M = (xr - xl) / dx;
    this->N = (yt - yb) / dy;
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;
    cout << "M = " << M << endl;
    cout << "N = " << N << endl;
}

template <class T>
EllipticSolver<T>::~EllipticSolver()
{
}

template <class T>
void EllipticSolver<T>::buildSystem(double (*f)(double, double), double (*g1)(double), double (*g2)(double), double (*g3)(double), double (*g4)(double)) {
    int m = M + 1;
    int n = N + 1;
    A = Matriz<T>(m*n, m*n); // ojo, n representa y, m representa x
    b = vector<T>(m*n);
    T dx2 = dx * dx;
    T dy2 = dy * dy;

    vector<T> x(m);
    vector<T> y(n);
    
    for(int i = 0; i < m; i++) { x[i] = xl + i * dx; }
    for(int j = 0; j < n; j++) { y[j] = yb + j * dy; }
    
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            if(j == 1) {
                A(i + j*m, i + j*m) = 1.0;
                b[i + j*m] = g1(x[i]);
            } else if(j == n - 2) {
                A(i + j*m, i + j*m) = 1.0;
                b[i + j*m] = g2(x[i]);
            } else if(i == 0) {
                A(i + j*m, i + j*m) = 1.0;
                b[i + j*m] = g3(y[j]);
            } else if(i == m - 1) {
                A(i + j*m, i + j*m) = 1.0;
                b[i + j*m] = g4(y[j]);
            } else {
                A(i + j*m, i + j*m) = -2/dx2 - 2/dy2;
                A(i + 1 + j*m, i + j*m) = 1/dx2;
                A(i - 1 + j*m, i + j*m) = 1/dx2;
                A(i + j*m, i + 1 + j*m) = 1/dy2;
                A(i + j*m, i - 1 + j*m) = 1/dy2;
                b[i + j*m] = f(x[i], y[j]);
            }
        }
    }


    cout << A << endl;

    for(int i = 0; i < m*n; i++) {
        cout << b[i] << " ";
    }
    cout << endl;
}

template <class T>
void EllipticSolver<T>::SOR() {
    // Construir matriz D, L y U
    D = Matriz<T>(A.nrow(), A.ncol());
    L = Matriz<T>(A.nrow(), A.ncol());
    U = Matriz<T>(A.nrow(), A.ncol());
    for(int i = 0; i < A.nrow(); i++) {
        for(int j = 0; j < A.ncol(); j++) {
            if(i == j) { D(i, j) = A(i, j); } else { D(i, j) = 0.0; }
            // Fila > columna es inferior
            if(i > j) { L(i, j) = A(i, j); } else { L(i, j) = 0.0; }
            // Fila < columna es superior
            if(i < j) { U(i, j) = A(i, j); } else { U(i, j) = 0.0; }
        }
    }
}

template <class T>
void EllipticSolver<T>::solve(double (*f)(double, double), double (*g1)(double), double (*g2)(double), double (*g3)(double), double (*g4)(double)) {
    buildSystem(f, g1, g2, g3, g4);
    SOR();
}

#endif // ELLIPTICSOLVER_H