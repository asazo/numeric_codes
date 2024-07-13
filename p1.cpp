#include <cmath>
#include "wavesolver.h"

using namespace std;

// Condición inicial U(x, 0)
double f(double x) {
    return 0.0;//sin(3.141692 * x);
}

// Condición de velocidad dU/dt(x, 0)
double g(double x) {
    return 0.0;
}

// Condición de borde izquierda l(a, t) = l(t) 
double l(double t) {
    return 0.0;
}

// Condición de borde derecha  r(b, t) = r(t)
double r(double t) {
    return 0.0;
}

// Source
double s(double x, double t) {
    return sin(t);
}

int main() {
    double a = 0;
    double b = 1;
    double t0 = 0;
    double tmax = 25.0;
    double c = 0.1;
    double M = 30;
    double N = 80; 
    WaveSolver<double> wavesolver = WaveSolver<double>(a, b, t0, tmax, c, M, N);
    wavesolver.solve_finitediff(f, g, l, r, s);
    wavesolver.save("./p1.txt");
    return 0;
}