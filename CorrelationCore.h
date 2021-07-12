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
    unsigned correlationPairsNum;

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
    virtual bool unInit();

    virtual void iterate(unsigned id);
    virtual void incrementTotal();
    virtual double getTotalDouble(unsigned steps){ mpf_class tmp = this->cp / steps; return tmp.get_d();}
    virtual double getTotal2Double(unsigned steps){ mpf_class tmp = this->cp2 / steps; return tmp.get_d();}

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
};

#endif //CORELLATIONCORE_H