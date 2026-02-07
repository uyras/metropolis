#ifndef CORELLATIONPOINTCORE_H
#define CORELLATIONPOINTCORE_H

#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
#include <forward_list>
#include "CalculationParameter.h"
#include "misc.h"
#include "Hamiltonians.h"
#include "dos2.h"

class CorrelationPointCore: public CalculationParameter
{

public:
    /**
     * @brief Список спинов, которые принадлежат каждой корреляционной точке
     */
    std::vector< std::forward_list < size_t > > correlationPointSpins;

    /**
     * @brief Список спинов-соседей для каждого спина 
     * (@todo переделать, т.к. спин может принадлежать разным точкам и там будут разные соседи)
     */
    std::vector< std::forward_list < size_t > > correlationNeighbours;

    /**
     * @brief Список спинов-соседей для каждого спина, с привязкой к каждой точке
     * Сперва идет нумерация по точкам, потом по спинам, и в конце перечисляются соседи
     * (@todo переделать, т.к. по второму измерению очень большой массив получается)
     */
    std::vector< std::vector< std::forward_list < size_t > > > correlationNeighboursByPoint; //fills only when histogram enabled
    std::map< std::pair<unsigned, unsigned>, short > correlationValues;
    unsigned correlationPairsNum;
    float spinsInPoint;

    CorrelationPointCore(
        const config_section_t &sect,
        const ConfigManager *conf);

    void findPointsArray();

    static string name(){ return "correlationpoint"; }

    virtual bool check(unsigned N) const;
    virtual void printHeader() const;
    virtual void init(const state_t &state);

    virtual void iterate(size_t id);
    virtual void incrementTotal();
    virtual mpf_class getTotal(unsigned steps){ return this->cp / steps; }
    virtual mpf_class getTotal2(unsigned steps){ return this->cp2 / steps; }

    virtual CorrelationPointCore * copy() { return new CorrelationPointCore(*this); }


    void enableHistogram(std::string filename) { this->_histogramEnabled = true; this->_histogramFilename = filename; }
    void disableHistogram(){ this->_histogramEnabled = false; }
    inline bool histogramEnabled(){ return this->_histogramEnabled; }

    void save(unsigned num);

private:

    long getFullTotal(const state_t &_state) const;
    unsigned pointCount() const {return X.size();}

    short method(const size_t idA, const size_t idB) const;

    double _minRange;
    double _maxRange;
    double _distance;
    long cpOld;
    mpf_class cp;
    mpf_class cp2;
    std::vector<double> X;
    std::vector<double> Y;

    /**
     * @brief В этой переменной сохраняется текущее состояние системы, для которой актуально значение cpOld.  
     */
    state_t currentState;

    bool _histogramEnabled;
    std::string _histogramFilename;
    Dos2<int> dos;
};

#endif //CORELLATIONPOINTCORE_H