#include "ellipticsolver.h"

int main() {
    double xl = 0;
    double xr = 3;
    double yb = 0;
    double yt = 2;
    int M = 8;
    int N = 5;
    EllipticSolver ellipticsolver = EllipticSolver<double>(xl, xr, yb, yt, M, N);
    ellipticsolver.solve();
    return 0;
}