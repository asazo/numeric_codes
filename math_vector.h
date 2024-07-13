#ifndef VECTORES_H
#define VECTORES_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

/*Declaracion de la clase*/

class Vector
{
  private:  // variables miembros
   vector<double> v;
  
  public: // metodos, constructores y destructores
    Vector();                      //  constructor por defecto
    Vector(int, double = 0);       //  constructor parametrizado
    ~Vector();                     //  destructor
    Vector(vector <double>);       //  constructor a partir de un vector
    Vector(const Vector &);      //  constructor de copia (copia superficial)
    double & operator()(int);        //  lectura-escritura utilizando el operador ()
    double operator()(int) const;    //  lectura utilizando el operador ()
    double & operator[](int);        //  lectura-escritura utilizando el operador []
    double operator[](int) const;    //  lectura utilizando el operador []
    double modulo();                 //  modulo de un vector
    int size();                      //  tamaño del vector
};

/*Implementacion de la clase*/

Vector::Vector(){}

Vector::Vector(int n, double x)
{
  vector<double> aux(n,x);
  this-> v = aux;
}

Vector::~Vector(){}

Vector::Vector(const Vector & u)
{
  v = u.v;
}

Vector::Vector(vector<double> u)
{
  v = u;
}

double & Vector::operator()(int i)
{
  return v[i];
}

double Vector::operator()(int i) const
{
  return v[i];
}

double & Vector::operator[](int i)
{
  return v[i];
}

double Vector::operator[](int i) const
{
  return v[i];
}

double Vector::modulo()
{
  double u = 0;
  for(int i = 0; i < int(v.size()); i+=1)
  {
    u += v[i] * v[i];
  }
  return sqrt(u);
}

int Vector::size()
{
  return v.size();
}

// Sobrecarga del operador de insercion <<
ostream & operator << (ostream & os, Vector & u)
{
  int n = u.size();
  os << "(";
  for(int i = 0; i < n-1; i+=1)
    {
      os << u[i] << ", ";
    }
  os << u[n-1] << ")";
  return os;
}

// Operaciones aritmeticas
Vector operator + (Vector x, Vector y)
{
  vector<double> aux;
  for(int i = 0; i < x.size(); i+=1)
  {
    aux.push_back(x[i] + y[i]);
  }
  return Vector(aux);
}

Vector operator - (Vector x, Vector y)
{
  vector<double> aux;
  for(int i = 0; i < x.size(); i+=1)
  {
    aux.push_back(x[i] - y[i]);
  }
  return Vector(aux);
}

Vector operator * (double a, Vector x)
{
  vector<double> aux;
  for(int i = 0; i < x.size(); i+=1)
  {
    aux.push_back(x[i] * a);
  }
  return Vector(aux);
}


Vector operator * (Vector x, double a)
{
  vector<double> aux;
  for(int i = 0; i < x.size(); i+=1)
  {
    aux.push_back(x[i] * a);
  }
  return Vector(aux);
}

double operator * (Vector x, Vector y){
  double S = 0.0;
  for(int i = 0; i < x.size(); i+=1)
  {
    S = S + x[i] * y[i];
  }
  return S;
}

Vector operator / (Vector x, double a)
{
  vector<double> aux;
  for(int i = 0; i < x.size(); i+=1)
  {
    aux.push_back(a / x[i]);
  }
  return Vector(aux);
} 

#endif