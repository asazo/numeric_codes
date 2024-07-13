#include "solver.h"
#include <cmath>


using namespace std;


double f1(double t, double y1, double y2) {
    return y2;
}

double f2(double t, double y1, double y2) {
    return y1*y1*y1 - 5*y1 - exp(-t*t)*y1;
}


int main() {
    double ta = 0;
    double tb = 4;
    double ya = 1.0;
    double yb = 0.0;
    double h = 1e-4;
    // Pendiente (slope) A
    double Sa = -1.0;
    // Pendiente (slope) B
    double Sb = 2.0;
    //elegir metodo
    string method = "rk4";

    BVPSolver bvpsolver = BVPSolver(ta, tb, ya, yb, h);
    bvpsolver.solve_shooting(f1, f2, Sa, Sb, method);
    bvpsolver.save("p3.txt");
    
    /*double y0 = 1.0;
    double t0 = 0.0;
    double tmax = 1.0;
    double h = 0.1;
    cout << "Resolviendo P3 ..." << endl;
    
    odesolver.rk4(f, t0, tmax, y0, h);
    odesolver.save("p3.txt");
 
    cout << "Archivo guardado en disco." << endl;*/
}