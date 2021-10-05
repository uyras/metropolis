#ifndef CORELLATIONCORE_H
#define CORELLATIONCORE_H

#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
#include "PartArray.h"
#include "CalculationParameter.h"

class CorrelationCore: public CalculationParameter
{

public:
    std::vector< std::forward_list < Part* > > correlationNeighbours;
    std::map< std::pair<unsigned, unsigned>, short > correlationValues;
    double correlationPairsNum;

    CorrelationCore(
        const std::string & parameterId,
        PartArray * prototype,
        double minRange, 
        double maxRange, 
        unsigned methodVar,
        const std::vector<uint64_t> & spins);

    virtual bool check(unsigned) const;
    virtual void printHeader(unsigned) const;
    virtual bool init(PartArray * sys);

    virtual void iterate(unsigned id);
    virtual void incrementTotal();
    virtual mpf_class getTotal(unsigned steps){ return this->cp / steps;}
    virtual mpf_class getTotal2(unsigned steps){ return this->cp2 / steps;}

    virtual CorrelationCore * copy() { return new CorrelationCore(*this); }


private:
    char method(const Part* partA, const Part* partB) const;
    char fxor(const Part* partA, const Part* partB) const;
    char feSign(const Part* partA, const Part* partB) const;
    char fScalar(const Part* partA, const Part* partB) const;


    void initMethod1(Part* partA, Part* partB);
    void initMethod2( Part* partA, Part* partB);
    void initMethod3( Part* partA, Part* partB);

    long getFullTotal(const PartArray * _sys) const;

    double minRange2;
    double maxRange2;

    double _minRange;
    double _maxRange;
    unsigned _methodVar;
    std::vector<uint64_t> spins;

    mpf_class cp;
    mpf_class cp2;
    long cpOld;

    bool areNeighbours(Part* partA, Part* partB)
    { 
        std::forward_list < Part* > &tmp = correlationNeighbours[partA->Id()];
        return std::find(tmp.begin(), tmp.end(), partB) != tmp.end();
    }
};

#endif //CORELLATIONCORE_H