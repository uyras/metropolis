#ifndef HAMILTONIANS_H
#define HAMILTONIANS_H

#include <cmath>
#include <string>
#include "misc.h"

inline double hamiltonian_dipolar(const Part& a, const Part& b){
    Vect rij = {a.p.x-b.p.x, a.p.y-b.p.y, a.p.z-b.p.z};
    double r2, r, r5,E;
    r2 = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
    r = sqrt(r2); 
    r5 = r2 * r2 * r; //радиус в пятой
    
    E = //энергия считается векторным методом, так как она не нужна для каждой оси
            (( (a.m.x * b.m.x + a.m.y * b.m.y + a.m.z * b.m.z) * r2)
                -
                (3 * (b.m.x * rij.x + b.m.y * rij.y + b.m.z * rij.z) * (a.m.x * rij.x + a.m.y * rij.y + a.m.z * rij.z)  )) / r5;
    return E;
}

inline hamiltonian_t hamiltonian_selector(string name){
    return hamiltonian_dipolar; //пока доступен только один гамильтониан
}

#endif //HAMILTONIANS_H