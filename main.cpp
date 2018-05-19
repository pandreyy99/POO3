#include <iostream>
#include <fstream>
#include "polinom.h"

using namespace std;

int main() {
    ifstream fin("C:\\Users\\Andrei\\Documents\\GitHub\\POO3\\polinom.txt");
    polinom<float> p, p1, p2;
    fin >> p1 >> p2;
    cout << "P1[ 3 ] = " << p1[3] << '\n';
    p = p1 + p2;
    cout << p << '\n';
    p = p1 * p2;
    cout << p << '\n';
    cout << p1 * (float) 5 << '\n';
    cout << (float) 5 * p1 << '\n';
    cout << "Aici incepe impartirea : " << '\n';
    p = p1 / p2;
    cout << p;
    return 0;
}
