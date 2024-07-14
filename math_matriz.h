//////////////////////////////
// 12-06-2023
// Julio Perez V.
// Plantilla de matrices nxm
//////////////////////////////

#ifndef MATRIZ_H
#define MATRIZ_H

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include <fstream>
#include <assert.h>

using namespace std;

                        /*Declaracion de la clase*/

template <class T>
class Matriz
{
    private: // variables miembro
    vector<vector<T>> matriz;

    public: // metodos, constructores y destructores
    Matriz();                                   // constructor por defecto
    ~Matriz();                                  // destructor
    Matriz(int, int);                           // constructor vacio
    Matriz(const Matriz &);                     // constructor de copia
    Matriz(vector<vector<T>>);                  // matriz creada a partir de un vector de vectores
    Matriz(int);                                // matriz identidad
    Matriz<T> operator = (const Matriz<T> &);   // operador de asignacion
    T operator()(int, int) const;               // lectura
    T & operator()(int, int);                   // lectura-escritura
    vector<vector<T>> get();                    // obtencion de elementos
    int nrow() const;                            // numero de filas
    int ncol() const;                            // numero de columnas
    vector<T> row(int);                         // obtener fila
    vector<T> col(int);                         // obtener columna 
    void save(string);                          // guardar a disco
    //Matriz<T> transpuesta(const Matriz<T> &);    // matriz transpuesta
    
    
};

template <class T>
T determinante(const Matriz<T> &);          // determinante
template <class T>
Matriz<T> inversa(const Matriz<T> &);       // matriz inversa


/*Implementacion de la clase*/

// Constructor por defecto
template <class T>
Matriz<T>::Matriz(){}

// Destructor
template <class T>
Matriz<T>::~Matriz(){}

// Constructor vacio
template <class T>
Matriz<T>::Matriz(int a, int b)
{
    vector<vector<T>> m;
    vector<T> m1;
    for (int i = 0; i < b; i+=1)
    {
        m1.push_back(0);
    }
    for (int j = 0; j < a; j+=1)
    {
        m.push_back(m1);
    }
    matriz = m;
}

// Constructor de copia
template <class T>
Matriz<T>::Matriz(const Matriz & a)
{
    matriz = a.matriz;
}

// Constructor a partir de un vector de vectores
template <class T>
Matriz<T>::Matriz(vector<vector<T>> a)
{
    matriz = a;
}

// Matriz identidad
template <class T>
Matriz<T>::Matriz(int n)
{
    vector<vector<T>> aux;
    for (int i = 0; i < n; i+=1)
    {
        vector<T> aux1(n, T(0));
        aux1[i] = T(1);          // asigna el valor '1' a los elementos de la diagonal de la matriz
        aux.push_back(aux1);
    }
    matriz = aux;
}

// Operador de asignacion
template <class T>
Matriz<T> Matriz<T>::operator = (const Matriz<T> & x)
{
    if (this != & x)
    {
        matriz = x.matriz;
    }
    return *this;
}

// Lectura 
template <class T>
T Matriz<T>::operator()(int i, int j) const
{
    return matriz[i][j];
}

// Lectura-escritura 
template <class T>
T & Matriz<T>::operator()(int i, int j)
{
    return matriz[i][j];
}

// Elementos de la matriz
template <class T>
vector<vector<T>> Matriz<T>::get()
{
    return matriz;    
}

// Filas de la matriz
template <class T>
int Matriz<T>::nrow() const
{
    return matriz.size();
}

// Columnas de la matriz
template <class T>
int Matriz<T>::ncol() const
{
    return matriz[0].size();
}

template<class T>
vector<T> Matriz<T>::row(int i) {
    return matriz[i];
}

template<class T>
vector<T> Matriz<T>::col(int j) {
    vector<T> res(matriz.size());
    for(int i = 0; i < matriz.size(); i++) {
        res[i] = matriz[i][j];
    }
    return res;
}

// Sobrecarga del operador <<
template <class T>
ostream & operator << (ostream & os, Matriz<T> a)
{
    vector<vector<T>> x = a.get();
    os << "[";
    for (int i = 0; i < x.size(); i+=1)
    {
        if (i == (x.size() - 1))
        {
            os << "[";
            for (int j = 0; j < x[i].size(); j+=1)
            {
                if (j == x[i].size() - 1)
                {
                    os << x[i][j];
                }
                else
                {
                    os << x[i][j] << ", ";
                }
            }
            os << "]]";
        }
        else
        {
            os << "[";
            for (int j = 0; j < x[i].size(); j+=1)
            {
                if (j == x[i].size() - 1)
                {
                    os << x[i][j];
                }
                else
                {
                    os << x[i][j] << ", ";
                }
            }
            os << "]," << endl;
        }
    }
    return os;
}

                       /*Operaciones aritmeticas*/

// Suma entre dos matrices
template <class T>
Matriz<T> operator + (Matriz<T> & a1, Matriz<T> & b1)
{
    Matriz<T> suma(a1.nrow(), a1.ncol());
    for (int i = 0; i < a1.nrow(); i+=1)
    {
        for (int j = 0; j < a1.ncol(); j+=1)
        {
            suma(i,j) = a1(i,j) + b1(i,j);
        }
    }
    return suma;
}

// Resta entre dos matrices
template <class T>
Matriz<T> operator - (Matriz<T> & a1, Matriz<T> & b1)
{
    Matriz<T> resta(a1.nrow(), a1.ncol());
    for (int i = 0; i < a1.nrow(); i+=1)
    {
        for (int j = 0; j < a1.ncol(); j+=1)
        {
            resta(i,j) = a1(i,j) - b1(i,j);
        }
    }
    return resta;
}

// Producto de escalar 
template <class T>
Matriz<T> operator * (T a1, Matriz<T> & b1)
{
    Matriz<T> producto(b1.nrow(), b1.ncol());
    for (int i = 0; i < b1.nrow(); i+=1)
    {
        for (int j = 0; j < b1.ncol(); j+=1)
        {
            producto(i,j) = a1*b1(i,j);
        }
    }
    return producto;
}

// Producto matriz vector
template <class T>
vector<T> operator * (const Matriz<T> & a, const vector<T> &x) 
{
    vector<T> resultado(a.nrow(), 0);
    for(int i = 0; i < a.nrow(); i++) {
        for(int j = 0; j < a.ncol(); j++) {
            resultado[i] += a(i, j) * x[j];
        }
    }
    return resultado;
}

// Producto entre dos matrices
template <class T>
Matriz<T> operator * (const Matriz<T> & a1, const Matriz<T> & b1)
{
    Matriz<T> producto(a1.nrow(), b1.ncol());
    for (int i = 0; i < a1.nrow(); i+=1)
    {
        for (int j = 0; j < b1.ncol(); j+=1)
        {
            for (int k = 0; k < a1.ncol(); k+=1)
            {
                producto(i,j) = producto(i,j) +  a1(i,k)*b1(k,j);
            }
        }
    }
    return producto;
}

// Sobrecarga del operador <<
template <class T>
ostream & operator << (ostream & os, vector<T> v) {
    os << "[";
    for(int i = 0; i < v.size(); i++) {
        os << v[i];
        if(i + 1 < v.size()) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

// Suma dos vectores
template <class T>
vector<T> operator + (const vector<T> & v1, const vector<T> v2) {
    assert(v1.size() == v2.size());
    vector<T> res(v1.size());
    for(int i = 0; i < v1.size(); i++) {
        res[i] = v1[i] + v2[i];
    }
    return res;
}

// Resta dos vectores
template <class T>
vector<T> operator - (const vector<T> & v1, const vector<T> v2) {
    assert(v1.size() == v2.size());
    vector<T> res(v1.size());
    for(int i = 0; i < v1.size(); i++) {
        res[i] = v1[i] - v2[i];
    }
    return res;
}

// Norma de vector
template <class T>
T norm2(const vector<T> & v) {
    T res = 0.0;
    for(int i = 0; i < v.size(); i++) { res += (v[i] * v[i]); }
    return sqrt(res);
}


///////////// Transpuesta /////////////
template <class T>
Matriz<T> transpuesta(const Matriz<T> & M)
{
    int filas = M.nrow();
    int columnas = M.ncol();

    Matriz<T> m_transpuesta(columnas, filas);

    for (int i = 0; i < filas; i+=1)
    {
        for (int j = 0; j < columnas; j+=1)
        {
            m_transpuesta(j,i) = M(i,j);
        }
    }
    return m_transpuesta;
}

///////////// Eliminacion Gaussiana /////////////
int contador = 0; // variable global
template <class T>
Matriz<T> gauss(const Matriz<T> & x)
{
    Matriz<T> m_gauss = x; // <--- matriz solucion
    
    // Reinicio de la variable contador
    contador = 0;

    for (int i = 0; i < m_gauss.nrow(); i+=1)
    {
        T pivote = m_gauss(i, i);
        int fila_pivote_maximo = i;
        for (int j = 0; j < m_gauss.nrow(); j+=1)
        {
            if (abs(m_gauss(i, j)) > abs(pivote))
            {
                pivote = m_gauss(i, j);
                fila_pivote_maximo = j;
            }
        }
        if (fila_pivote_maximo != i)
        {
            for (int j = 0; j < m_gauss.ncol(); j+=1)
            {
                double dato_i_j = m_gauss(i, j);
                m_gauss(i, j) = m_gauss(fila_pivote_maximo, j);
                m_gauss(fila_pivote_maximo, j) = dato_i_j;
            }
            contador += 1;
        }
        for (int j = 0; j < m_gauss.nrow(); j+=1)
        {
            if (i != j)
            {
                T factor = m_gauss(j, i)/m_gauss(i, i);
                for (int k = 0; k < m_gauss.ncol(); k+=1)
                {
                    m_gauss(j, k) = m_gauss(j, k) - factor*m_gauss(i, k);
                }
            }
        }
    }
    return m_gauss;
}

///////////// Determinante /////////////
template <class T>
T determinante(const Matriz<T> & x)
{
    // Matriz reducida
    Matriz<T> m_reducida = gauss(x);
    T determinante = 1.0;

    determinante = determinante*pow(-1, contador);

    for (int i = 0; i < x.nrow(); i+=1)
    {
        determinante = determinante * m_reducida(i, i);
    }
    return determinante;
}

///////////// Inversa /////////////
template <class T>
Matriz<T> inversa(const Matriz<T> & x)
{
    // Matriz extendida
    Matriz<T> m_extendida(x.nrow(), x.ncol()*2);

    // Matriz original
    for (int i = 0; i < x.nrow(); i+=1)
    {
        for (int j = 0; j < x.ncol(); j+=1)
        {
            m_extendida(i, j) = x(i, j);
        }
    }

    // Matriz identidad
    for (int i = 0; i < x.nrow(); i+=1)
    {
        m_extendida(i, i + x.ncol()) = 1;
    }
    
    // Eliminacion gaussiana a la matriz extendida
    m_extendida = gauss(m_extendida);

    for (int i = 0; i < x.nrow(); i+=1)
    {
        T factor = m_extendida(i, i);
        
        for (int j = 0; j < x.ncol(); j+=1)
        {
            m_extendida(i, j + x.ncol()) = m_extendida(i, j + x.ncol())/factor;

        }
    }

    // Matriz inversa
    Matriz<T> m_inversa(x.nrow(), x.ncol());

    // La matriz inversa estara formada por los elementos de la matriz identidad
    // luego de aplicarle Gauss
    for (int i = 0; i < x.nrow(); i+=1)
    {
        for (int j = 0; j < x.ncol(); j+=1)
        {
            m_inversa(i, j) = m_extendida(i, j + x.ncol());
        }
    }
    return m_inversa;
}

template <class T>
void Matriz<T>::save(string path) {
    ofstream myfile;
    myfile.open(path);
    for(int i = 0; i < matriz.size(); i++) {
        for (int j = 0; j < matriz[i].size(); j++) {
            myfile <<std::fixed << std::setw(11) << std::setprecision(6) << matriz[i][j];
            if(j <  matriz[i].size() - 1) {
                myfile << ",";
            }
        }
        myfile << endl;
    }
    myfile.close();
}

#endif