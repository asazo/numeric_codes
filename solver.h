#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

int sign(double n) {
    return (n > 0) - (n < 0);
}

class ODESolver {
    public:
    ODESolver(string method);
    string method;
    vector<double> t;
    vector<vector<double>> ysol;
    void euler_system(double (*f1)(double, double, double), double(*f2)(double, double, double), double t0, double tmax, double y1_0, double y2_0, double h);
    void rk4(double (*f)(double, double), double t0, double tmax, double y0, double h);
    void rk4_system(double (*f1)(double, double, double), double(*f2)(double, double, double), double t0, double tmax, double y1_0, double y2_0, double h);
    void solve(double (*f1)(double, double, double), double(*f2)(double, double, double), double t0, double tmax, double y1_0, double y2_0, double h);
    void save(const string path);
};

ODESolver::ODESolver(string method) {
    this->method = method;
}

void ODESolver::solve(double (*f1)(double, double, double), double(*f2)(double, double, double), double t0, double tmax, double y1_0, double y2_0, double h) {
    if(method == "euler") {
        euler_system(f1, f2, t0, tmax, y1_0, y2_0, h);
    } else if(method == "rk4") {
        rk4_system(f1, f2, t0, tmax, y1_0, y2_0, h);
    }
}

/*
    Guardar resultado del solver a disco
*/
void ODESolver::save(const string path) {
    ofstream output_file;
    output_file.open(path);
    for(int i = 0; i < t.size(); i++) {
        output_file <<std::fixed << std::setw(11) << std::setprecision(6) << t[i];
        if(i < t.size() - 1) { output_file << ","; }
    }
    output_file << endl;
    // Guardar por cada solucion y1, y2, ... yn
    for(int i = 0; i < ysol.size(); i++) {
        for(int j = 0; j < ysol[i].size(); j++) {
            output_file <<std::fixed << std::setw(11) << std::setprecision(6) << ysol[i][j];
            if(j < ysol[i].size() - 1) { output_file << ","; }
        }
        output_file << endl;
    }
    output_file.close();
}

/*
  Implementación Runge Kutta 4
  Args:
    f(t, y):  lado derecho de IVP
    t0, tmax: intervalo de resolución
    y0:       condicion inicial
    h:        step size
*/
void ODESolver::rk4(double (*f)(double, double), double t0, double tmax, double y0, double h) {
    double ti = t0;
    double yi = y0;
    double s1, s2, s3, s4;

    ysol.push_back(vector<double>());
    t.push_back(t0);
    ysol[0].push_back(y0);
    
    for(int i = 0; ti < tmax; i++) {
        s1 = f(ti, yi);
        s2 = f(ti + 0.5 * h, yi + 0.5 * s1);
        s3 = f(ti + 0.5 * h, yi + 0.5 * s2);
        s4 = f(ti + h, yi + h * s3);
        yi = yi + h * (s1 + 2 * s2 + 2 * s3 + s4) / 6.0;
        ysol[0].push_back(yi);
        ti = t0 + (i + 1) * h;
        t.push_back(ti);
    }
}

void ODESolver::euler_system(double (*f1)(double, double, double), double(*f2)(double, double, double), double t0, double tmax, double y1_0, double y2_0, double h) {
    double ti = t0;
    double y1_i = y1_0;
    double y2_i = y2_0;

    ysol.push_back(vector<double>());
    ysol.push_back(vector<double>());
    ysol[0].push_back(y1_0);
    ysol[1].push_back(y2_0);
    t.push_back(ti);
    
    for(int i = 0; ti < tmax; i++) {
        y1_i = y1_i + h*f1(ti, y1_i, y2_i);
        y2_i = y2_i + h*f2(ti, y1_i, y2_i);
    
        ysol[0].push_back(y1_i);
        ysol[1].push_back(y2_i);
        ti = t0 + (i + 1) * h;
        t.push_back(ti);
    }
}


void ODESolver::rk4_system(double (*f1)(double, double, double), double(*f2)(double, double, double), double t0, double tmax, double y1_0, double y2_0, double h) {
    double ti = t0;
    double y1_i = y1_0;
    double y2_i = y2_0;
    double s1_1, s2_1, s3_1, s4_1;
    double s1_2, s2_2, s3_2, s4_2;

    // Por cada function, un vector solucion
    ysol.push_back(vector<double>());
    ysol.push_back(vector<double>());
    ysol[0].push_back(y1_0);
    ysol[1].push_back(y2_0);
    t.push_back(ti);

    for(int i = 0; ti < tmax; i++) {
        s1_1 = f1(ti, y1_i, y2_i);
        s1_2 = f2(ti, y1_i, y2_i);

        s2_1 = f1(ti + 0.5 * h, y1_i + 0.5 * s1_1, y2_i + 0.5 * s1_2);
        s2_2 = f2(ti + 0.5 * h, y1_i + 0.5 * s1_1, y2_i + 0.5 * s1_2);

        s3_1 = f1(ti + 0.5 * h, y1_i + 0.5 * s2_1, y2_i + 0.5 * s2_2);
        s3_2 = f2(ti + 0.5 * h, y1_i + 0.5 * s2_1, y2_i + 0.5 * s2_2);

        s4_1 = f1(ti + h, y1_i + h * s3_1, y2_i + h * s3_2);
        s4_2 = f2(ti + h, y1_i + h * s3_1, y2_i + h * s3_2);
        
        y1_i = y1_i + h * (s1_1 + 2 * s2_1 + 2 * s3_1 + s4_1) / 6.0;
        y2_i = y2_i + h * (s1_2 + 2 * s2_2 + 2 * s3_2 + s4_2) / 6.0;

        ysol[0].push_back(y1_i);
        ysol[1].push_back(y2_i);
        ti = t0 + (i + 1) * h;
        t.push_back(ti);
    }
}


class BVPSolver  {
    public:
    BVPSolver(double ta, double tbd, double ya, double yb, double h);
    double ta;
    double tb;
    double ya;
    double yb;
    double h;
    vector<double> t;
    vector<vector<double>> ysol;
    void solve_shooting(double (*f1)(double, double, double), double (*f2)(double, double, double), double Sa, double Sb, string method);
    void save(const string path);
};

BVPSolver::BVPSolver(double ta, double tb, double ya, double yb, double h) {
    this->ta = ta;
    this->tb = tb;
    this->ya = ya;
    this->yb = yb;
    this->h = h;
}

void BVPSolver::save(const string path) {
    ofstream output_file;
    output_file.open(path);
    for(int i = 0; i < t.size(); i++) {
        output_file <<std::fixed << std::setw(11) << std::setprecision(6) << t[i];
        if(i < t.size() - 1) { output_file << ","; }
    }
    output_file << endl;
    // Guardar por cada solucion y1, y2, ... yn
    for(int i = 0; i < ysol.size(); i++) {
        for(int j = 0; j < ysol[i].size(); j++) {
            output_file <<std::fixed << std::setw(11) << std::setprecision(6) << ysol[i][j];
            if(j < ysol[i].size() - 1) { output_file << ","; }
        }
        output_file << endl;
    }
    output_file.close();
}


/*
    Implementación método del disparo
    Args:
        f(t, y, y'):  lado derecho de BVP
        t0, tmax: intervalo de resolución
        y0:       condicion inicial
        h:        step size
*/
void BVPSolver::solve_shooting(double (*f1)(double, double, double), double (*f2)(double, double, double), double Sa, double Sb, string method) {
    // Funcion para metodo de biseccion F(s), buscamos F(s) = 0
    auto F = [](double (*f1)(double, double, double), double (*f2)(double, double, double), double ta, double tb, double ya,
     double yb, double h, double s, string method) {
        ODESolver odesolver = ODESolver(method);
        odesolver.solve(f1, f2, ta, tb, ya, s, h);
        return odesolver.ysol[0].back() - yb;
    };

    cout << "F(Sa) = " << F(f1, f2, ta, tb, ya, yb, h, Sa, method) << endl;
    cout << "F(Sb) = " << F(f1, f2, ta, tb, ya, yb, h, Sb, method) << endl;

    // Metodo de bisección
    int n = 1;
    int max_iter = 100;
    double tol = 1e-14;
    double c, Fc, FA;

    while (n < max_iter) {
        c = 0.5 * (Sa + Sb);
        //cout << c << endl;
        Fc = F(f1, f2, ta, tb, ya, yb, h, c, method); 
        if (Fc == 0) {
            cout << "Resultado en pendiente " << c << endl;
            break;
        } else if((Sb - Sa) *0.5 < tol) {
            cout << "Tolerancia alcanzada " << (Sb - Sa) *0.5 << endl; 
            cout << "F(c) = " << Fc << endl; 
            break;
        }
        n++;

        FA = F(f1, f2, ta, tb, ya, yb, h, Sa, method); 
        if (sign(Fc) == sign(FA)) {
            Sa = c;
        } else {
            Sb = c;
        }
    }

    ODESolver odesolver = ODESolver(method);
    odesolver.rk4_system(f1, f2, ta, tb, ya, c, h);
    ysol = odesolver.ysol;
    t = odesolver.t;
}

#endif //SOLVER_H