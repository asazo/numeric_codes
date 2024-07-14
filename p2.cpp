#include <cmath>
#include "ellipticsolver.h"


//phi(x, 0) = g1(x)
double g1(double x) {
    return log(x*x + 1);
}

//phi(x, 2) = g2(x)
double g2(double x) {
    return log(x*x + 4);
}

//phi(0, y) = g3(y)
double g3(double y) {
    return 2*log(y);
}

//phi(3, y) = g4(y)
double g4(double y) {
    return log(y*y + 1);
}

double rho_0(double x, double y) {
    return pow((x-1)*(x-1) + (y - 0.5)*(y - 0.5), 2) - pow((x - 1), 3) - pow((y- 0.5), 3);
}

double rho(double x, double y) {
    double rho0 = rho_0(x, y);
    return 200 * rho0 * (rho0 <= 0);
}

double f(double x, double y) {
    return 0.0; //-rho(x, y);
}

int main() {
    double xl = 0;
    double xr = 1;
    double yb = 1;
    double yt = 2;
    double dx = 0.25;
    double dy = 0.25;
    double omega = 1.0;
    EllipticSolver ellipticsolver = EllipticSolver<double>(xl, xr, yb, yt, dx, dy);
    ellipticsolver.solve(f, g1, g2 ,g3, g4, omega);
    cout << ellipticsolver.phi << endl;
    return 0;
}