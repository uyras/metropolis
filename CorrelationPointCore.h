#ifndef CORELLATIONPOINTCORE_H
#define CORELLATIONPOINTCORE_H

#include <vector>
#include <map>
#include <gmpxx.h>
#include "PartArray.h"

class CorrelationPointCore{

public:
    std::vector< std::forward_list < Part* > > correlationPointSpins;
    std::vector< std::forward_list < Part* > > correlationNeighbours;
    std::map< std::pair<unsigned, unsigned>, short > correlationValues;
    unsigned correlationPairsNum;
    float spinsInPoint;

    CorrelationPointCore(const std::vector<double> & X, 
        const std::vector<double> & Y, 
        double distance, double minRange, double maxRange);

    /**
     * @brief Сaches the neighbours and energies for further calculations
     * 
     * @param sys system
     */
    void init(const PartArray & sys);
    long getCPFull(const PartArray & sys);

    unsigned pointCount() const {return X.size();}

    void method(const Part* partA);

    double _minRange;
    double _maxRange;
    double _distance;
    long cpOld;
    double dbg = false;
    mpf_class cp;
    mpf_class cp2;
    std::vector<double> X;
    std::vector<double> Y;
};

#endif //CORELLATIONPOINTCORE_H