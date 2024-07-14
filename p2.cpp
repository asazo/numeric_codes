#include <cmath>
#include "ellipticsolver.h"


//phi(x, 0) = g1(x)
double g1(double x) {
    return -1.0;//sin(M_PI * x);//log(x*x + 1);
}

//phi(x, 2) = g2(x)
double g2(double x) {
    return -1.0;//sin(M_PI * x);//log(x*x + 4);
}

//phi(0, y) = g3(y)
double g3(double y) {
    return 1.0;
}

//phi(3, y) = g4(y)
double g4(double y) {
    return 1.0;
}

double rho_0(double x, double y) {
    return pow((x-1)*(x-1) + (y - 0.5)*(y - 0.5), 2) - pow((x - 1), 3) - pow((y- 0.5), 3);
}

double rho(double x, double y) {
    double rho0 = rho_0(x, y);
    return 200 * rho0 * (rho0 <= 0);
}

double f(double x, double y) {
    return -rho(x, y);
}

int main() {
    double xl = 0.0;
    double xr = 3.0;
    double yb = 0.0;
    double yt = 2.0;
    double dx = 0.15;//0.15;
    double dy = 0.1;//0.10;
    double omega = 1.8;
    double px = 2.5;
    double py = 1.5;
    double r = 0.3;
    EllipticSolver ellipticsolver = EllipticSolver<double>(xl, xr, yb, yt, dx, dy);
    ellipticsolver.solve(f, g1, g2 ,g3, g4, omega, px, py, r);
    ellipticsolver.set_empty_disk(px, py ,r);
    ellipticsolver.save("./p2.txt");
    return 0;
}