#ifndef CORELLATIONCORE_H
#define CORELLATIONCORE_H

#include <vector>
#include <map>
#include <gmpxx.h>
#include "PartArray.h"

class CorrelationCore{

public:
    std::vector< std::forward_list < Part* > > correlationNeighbours;
    std::map< std::pair<unsigned, unsigned>, short > correlationValues;
    unsigned correlationPairsNum;

    CorrelationCore(double minRange, double maxRange, unsigned methodVar);

    /**
     * @brief Сaches the neighbours and energies for further calculations
     * 
     * @param sys system
     */
    void init(const PartArray & sys);
    long getCPFull(const PartArray & sys);

    char method(const Part* partA,const Part* partB);
    char fxor(const Part* partA,const Part* partB);
    char feSign(const Part* partA, const Part* partB);
    char fScalar(const Part* partA, const Part* partB);

    double _minRange;
    double _maxRange;
    unsigned _methodVar;
    long cpOld;
    double dbg = false;
    mpf_class cp;
    mpf_class cp2;
};

#endif //CORELLATIONCORE_H