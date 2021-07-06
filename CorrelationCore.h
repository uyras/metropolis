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

    CorrelationCore(
        double minRange, 
        double maxRange, 
        unsigned methodVar,
        const std::vector<unsigned> & spins);

    /**
     * @brief Сaches the neighbours and energies for further calculations
     * 
     * @param sys system
     */
    void init(const PartArray & sys);
    long getCPFull(const PartArray & sys);

    void method(const Part* partA);
    char fxor(const Part* partA,const Part* partB);
    char feSign(const Part* partA, const Part* partB);
    char fScalar(const Part* partA, const Part* partB);

    double _minRange;
    double _maxRange;
    unsigned _methodVar;
    std::vector<unsigned> spins;

    long cpOld;
    double dbg = false;
    mpf_class cp;
    mpf_class cp2;

private:
    void initMethod1(Part* partA, Part* partB);
    void initMethod2( Part* partA, Part* partB);
    void initMethod3( Part* partA, Part* partB);
    double minRange2;
    double maxRange2;
};

#endif //CORELLATIONCORE_H