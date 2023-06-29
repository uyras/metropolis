#ifndef CORELLATIONPOINTCORE_H
#define CORELLATIONPOINTCORE_H

#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
#include "PartArray.h"
#include "CalculationParameter.h"
#include <dos2.h>

class CorrelationPointCore: public CalculationParameter
{

public:
    std::vector< std::forward_list < Part* > > correlationPointSpins;
    std::vector< std::forward_list < Part* > > correlationNeighbours;
    std::vector< std::vector< std::forward_list < Part* > > > correlationNeighboursByPoint; //fills only when histogram enabled
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

    virtual void iterate(unsigned id);
    virtual void incrementTotal();
    virtual mpf_class getTotal(unsigned steps){ return this->cp / steps; }
    virtual mpf_class getTotal2(unsigned steps){ return this->cp2 / steps; }

    virtual CorrelationPointCore * copy() { return new CorrelationPointCore(*this); }


    void enableHistogram(std::string filename) { this->_histogramEnabled = true; this->_histogramFilename = filename; }
    void disableHistogram(){ this->_histogramEnabled = false; }
    inline bool histogramEnabled(){ return this->_histogramEnabled; }

    void save(unsigned num);

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

    bool _histogramEnabled;
    std::string _histogramFilename;
    Dos2<int> dos;
};

#endif //CORELLATIONPOINTCORE_H