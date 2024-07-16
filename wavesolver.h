#ifndef WAVESOLVER_H
#define WAVESOLVER_H

#include <string>
#include "math_matriz.h"

using namespace std;

template <class T>
class WaveSolver {
    public:
    WaveSolver(T a, T b, T t0, T tmax, T c, T M, T N);
    ~WaveSolver();
    Matriz<T> U;
    Matriz<T> A; // parametrizacion del espacio (linealización). U(0), U(1), U(2)
    T a;
    T b;
    T t0;
    T tmax;
    T c;
    T M;
    T N;
    void solve_finitediff(T (*f)(T), T (*g)(T), T (*l)(T), T(*r)(T), T (*s)(T, T));
    void save(const string path);
    private:
    void gen_matriz_A(const int, const T);
    void solve();
};


template <class T>
WaveSolver<T>::WaveSolver(T a, T b, T t0, T tmax, T c,T M, T N) {
    this->a = a;
    this->b = b;
    this->t0 = t0;
    this->tmax = tmax;
    this->c = c;
    this->M = M;
    this->N = N;
}

template <class T>
WaveSolver<T>::~WaveSolver(){}

/*
    Construir la matriz A de ecuación de onda
*/
template <class T>
void WaveSolver<T>::gen_matriz_A(const int M, const T sigma) {
    A = Matriz<T>(M, M);
    for(int i = 0; i < M; i ++) {
        for(int j = 0; j < M; j ++) {
            // diagonal principal
            if(i == j) {
                A(i, j) = 2 - 2 * sigma * sigma;
            }
            // diagonal inferior
            else if(i - 1 == j) {
                A(i, j) = sigma * sigma;
            }
            // diagonal superior
            else if(i + 1 == j) {
                A(i, j) = sigma * sigma;
            }
        }
    }
}

template <class T>
void WaveSolver<T>::solve_finitediff(T (*f)(T), T (*g)(T), T (*l)(T), T(*r)(T), T (*s)(T, T)) {
    vector<T> xx(M, 0); // x = (x0, x1, x2, ..., x[m-1])
    vector<T> tt(N, 0); //
    
    U = Matriz<T>(M, N);

    T dx = (b - a) / (M - 1);
    T dt = (tmax - t0) / (N - 1);
    T sigma = c * dt / dx;
    T sigma2 = sigma * sigma;

    cout << "c*dt/dx = " << sigma << endl;
    
    for(int i = 0; i < M; i++) { xx[i] = a + i * dx; }
    for(int j = 0; j < N; j++) { tt[j] = t0 + j * dt; }

    gen_matriz_A(M, sigma);
    cout << "Matriz A generada" << endl;

    // Condición inicial
    for(int i = 0; i < M; i++) {
        U(i, 0) = f(xx[i]);
    }

    /*  U(0,0) U(1,0), U(2,0) ... U(M-1, 0) // t=0   f(x, 0)
        U(0,1) U(1,1), U(2,1) ... U(M-1, 1) // t=1
        ...
        U(0,25) U(1,25), U(2,25) ... U(M-1, 25) // t=25
    */
    

    // Condición de borde izquierda y derecha
    for(int j = 1; j < N; j++) {
        U(0, j) = l(tt[j]);
        U(M-1, j) = r(tt[j]);
    }
    
    cout << U.nrow() << " x " << U.ncol() << endl;

    // 1er time-step
    vector<double> A_times_U0 = (0.5 * A) * U.col(0);

    for(int i = 1; i < M-1; i++) {
        U(i, 1) = A_times_U0[i] + dt * g(xx[i]) + dt*dt * s(xx[i], tt[0]);
        if(i == 1) {
            U(i, 1) += 0.5 * sigma2 * U(0, 0);
        } else if(i == M - 2) {
            U(i, 1) += 0.5 * sigma2 * U(M-1, 0);
        }
    }

    cout << "1er time-step" << endl;

    // Siguientes time-steps
    vector<double> A_times_Uj;
    for(int j = 2; j < N; j++) { // t = 2 en adelante
        A_times_Uj = A * U.col(j - 1);
        for(int i = 1; i < M-1; i++) {
            U(i, j) = A_times_Uj[i] - U(i, j - 2) + dt*dt * s(xx[i], tt[j-1]);
            if(i == 1) {
                U(i, j) += sigma2 * U(0, j - 1);
            } else if(i == M - 2) {
                U(i, j) += sigma2 * U(M - 1, j - 1);
            }
        }
    }

    cout << "Fin" << endl;
}

template <class T>
void WaveSolver<T>::save(const string path) {
    
    // Guardar resultado y parametros a disco
    U.save(path);
    ofstream myfile;
    myfile.open(path + ".params");
    myfile << "{" << endl;
    myfile << "    \"a\": " << a << "," << endl;
    myfile << "    \"b\": " << b << "," << endl;
    myfile << "    \"t0\": " << t0 << "," << endl;
    myfile << "    \"tmax\": " << tmax << "," << endl;
    myfile << "    \"c\": " << c << "," << endl;
    myfile << "    \"M\": " << M << "," << endl;
    myfile << "    \"N\": " << N << endl;
    myfile << "}" << endl;
    myfile.close();
}




#endif //WAVESOLVER_H