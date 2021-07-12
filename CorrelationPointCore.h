#ifndef CORELLATIONPOINTCORE_H
#define CORELLATIONPOINTCORE_H

#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
#include "PartArray.h"
#include "CalculationParameter.h"

class CorrelationPointCore: public CalculationParameter
{

public:
    std::vector< std::forward_list < Part* > > correlationPointSpins;
    std::vector< std::forward_list < Part* > > correlationNeighbours;
    std::map< std::pair<unsigned, unsigned>, short > correlationValues;
    unsigned correlationPairsNum;
    float spinsInPoint;

    CorrelationPointCore(
        const std::string & parameterId,
        PartArray * prototype,
        const std::vector<double> & X, 
        const std::vector<double> & Y, 
        double distance, double minRange, double maxRange);

    virtual bool check(unsigned N) const;
    virtual void printHeader(unsigned) const;
    virtual bool init(PartArray * sys);
    virtual bool unInit();

    virtual void iterate(unsigned id);
    virtual void incrementTotal();
    virtual double getTotalDouble(unsigned steps){ mpf_class tmp = this->cp / steps; return tmp.get_d();}
    virtual double getTotal2Double(unsigned steps){ mpf_class tmp = this->cp2 / steps; return tmp.get_d();}

    virtual CorrelationPointCore * copy() { return new CorrelationPointCore(*this); }

private:

    long getFullTotal(const PartArray * _sys) const;
    unsigned pointCount() const {return X.size();}

    short method(const Part* partA, const Part* partB) const;

    double _minRange;
    double _maxRange;
    double _distance;
    long cpOld;
    mpf_class cp;
    mpf_class cp2;
    std::vector<double> X;
    std::vector<double> Y;
};

#endif //CORELLATIONPOINTCORE_H