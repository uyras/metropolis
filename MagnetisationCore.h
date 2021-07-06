#ifndef MAGNETISATIONCORE_H
#define MAGNETISATIONCORE_H

#include <vector>
#include <map>
#include <gmpxx.h>
#include "PartArray.h"

class MagnetisationCore{

public:

    MagnetisationCore(const Vect & vector, const std::vector<unsigned> & spins);

    void init(PartArray & sys);
    long getMOFull();

    void method(unsigned spinId);

    Vect vector;
    std::vector<unsigned> spins;
    std::vector< double > magnetisationValues;

    double mOld;
    double dbg = false;
    mpf_class mv;
    mpf_class mv2;
    PartArray * _sys;
};

#endif //MAGNETISATIONCORE_H