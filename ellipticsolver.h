#ifndef ELLIPTICSOLVER_H
#define ELLIPTICSOLVER_H

#include <iostream>
#include <fstream>
#include "math_matriz.h"

using namespace std;

template <class T>
T distance(T x1, T y1, T x2, T y2) {
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

template <class T>
class EllipticSolver
{
private:
    Matriz<T> L;
    Matriz<T> D;
    Matriz<T> U;
    vector<T> x;
    vector<T> y;
    void buildSystem(T (*f)(T, T), T (*g1)(T), T (*g2)(T), T (*g3)(T), T (*g4)(T), T px, T py, T r);
    void SOR(T omega);
public:
    EllipticSolver(T xl, T xr, T yb, T yt, T dx, T dy); 
    ~EllipticSolver();
    Matriz<T> A;
    vector<T> b;
    Matriz<T> phi; // Respuesta
    T xl, xr, yb, yt, dx, dy;
    int M, N; //linspace(a, b, M) <--> arange(a, b, dx) [a, b)
    int m, n;    
    void solve(T (*f)(T, T), T (*g1)(T), T (*g2)(T), T (*g3)(T), T (*g4)(T), T omega, T px, T py, T r);
    void set_empty_disk(T px, T py, T r);
    void save(const string path);
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
    this->m = M + 1;
    this->n = N + 1;
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;
    cout << "M = " << M << endl;
    cout << "N = " << N << endl;
    cout << "m = " << m << endl;
    cout << "n = " << n << endl;
    x = vector<T>(m);
    y = vector<T>(n);
    for(int i = 0; i < m; i++) { x[i] = xl + i * dx; }
    for(int j = 0; j < n; j++) { y[j] = yb + j * dy; }
}

template <class T>
EllipticSolver<T>::~EllipticSolver()
{
}

template <class T>
void EllipticSolver<T>::set_empty_disk(T px, T py, T r) {    
    T d; // distancia
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            d = distance(x[i], y[j], px, py);
            if(d <= r) {
                cout << x[i] << " " << y[j] << endl;
                cout << phi(i, j) << endl;
            }
        }
    }
}

template <class T>
void EllipticSolver<T>::buildSystem(T (*f)(T, T), T (*g1)(T), T (*g2)(T), T (*g3)(T), T (*g4)(T), T px, T py, T r) {
    A = Matriz<T>(m*n, m*n); // ojo, n representa y, m representa x
    b = vector<T>(m*n);
    T dx2 = dx * dx;
    T dy2 = dy * dy;


    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    
    // ptos interiores
    for(int i = 2; i < m; i++) {
        for(int j = 2; j < n; j++) {
            A(i + (j-1)*m - 1, i + (j-1)*m - 1) = -2/dx2 - 2/dy2;
            A(i + (j-1)*m - 1, i - 1 + (j-1)*m - 1) = 1/dx2;
            A(i + (j-1)*m - 1, i + 1 + (j-1)*m - 1) = 1/dx2;
            A(i + (j-1)*m - 1, i + (j - 2)*m - 1) = 1/dy2;
            A(i + (j-1)*m - 1, i + j*m - 1) = 1/dy2;
            b[i + (j-1)*m - 1] = f(x[i], y[j]);
        }
    }

    // condiciones de borde
    for(int i = 0; i < m; i++) {
        int j = 0;
        A(i + j*m, i + j*m) = 1.0;
        b[i + j*m] = g1(x[i]);
        j = n - 1;
        A(i + j*m, i + j*m) = 1.0;
        b[i + j*m] = g2(x[i]);
    }

    for(int j = 0; j < n; j++) {
        int i = 0;
        A(i + j*m, i + j*m) = 1.0;
        b[i + j*m] = g3(y[j]);
        i = m - 1;
        A(i + j*m, i + j*m) = 1.0;
        b[i + j*m] = g4(y[j]);
    }


    // Force Phi(x,y) = 0 at some points
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            double d = distance(x[i], y[j], px, py);
            if(d <= r) {
                // cautela!
                A(i + j*m, i + j*m) = -1000000;//;
                //cout << A(i + j*m, i + j*m) << endl; //= 1.0;
                b[i + j*m] = 1.0;
            }
        }
    }

}

template <class T>
void EllipticSolver<T>::SOR(T w) {
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

    Matriz<T> wL = w * L;
    Matriz<T> inv = inversa(D + wL);
    vector<T> xsol(b.size(), 0.0), xprev(b.size(), 0.0); //xsol(m*n)

    cout << "Begin SOR" << endl;
    int iter = 0;
    int max_iter = 500;
    while(iter < max_iter) {
        xprev = xsol;
        xsol = inv*((1 - w)*D*xprev - (w*U)*xprev) + (w*inv)*b;
        if(norm2(xprev - xsol) < 1e-14) { cout << "TOL reached" << endl; break; }
        iter++;
    }
    cout << "End SOR" << endl;
    
    phi = Matriz<T>(m, n);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            phi(i, j) = xsol[i + j*m];
        }
    }
}

template <class T>
void EllipticSolver<T>::solve(T (*f)(T, T), T (*g1)(T), T (*g2)(T), T (*g3)(T), T (*g4)(T), T omega, T px, T py, T r) {
    buildSystem(f, g1, g2, g3, g4, px, py ,r);
    SOR(omega);
}

template <class T>
void EllipticSolver<T>::save(const string path) {
    phi.save(path);

    ofstream myfile;
    myfile.open(path + ".params");
    myfile << "{" << endl;
    myfile << "    \"xl\": " << xl << "," << endl;
    myfile << "    \"xr\": " << xr << "," << endl;
    myfile << "    \"yb\": " << yb << "," << endl;
    myfile << "    \"yt\": " << yt << "," << endl;
    myfile << "    \"M\": " << M << "," << endl;
    myfile << "    \"N\": " << N <<  endl;
    myfile << "}" << endl;
    myfile.close();
}

#endif // ELLIPTICSOLVER_H