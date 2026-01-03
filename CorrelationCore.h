#ifndef CORELLATIONCORE_H
#define CORELLATIONCORE_H

#include <inicpp/inicpp.h>
#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
#include <cstdint>
#include <forward_list>
#include "CalculationParameter.h"
#include "ConfigManager.h"
#include "misc.h"
#include "Hamiltonians.h"

class CorrelationCore: public CalculationParameter
{

public:
    CorrelationCore(
        const config_section_t &sect,
        const ConfigManager *conf);

    static string name(){ return "correlation"; }

    virtual bool check(unsigned) const;
    virtual void printHeader() const;
    virtual void init(const state_t &state);

    virtual void iterate(size_t id);
    virtual void incrementTotal();
    virtual mpf_class getTotal(unsigned steps){ return this->cp / steps;}
    virtual mpf_class getTotal2(unsigned steps){ return this->cp2 / steps;}

    virtual CorrelationCore * copy() { return new CorrelationCore(*this); }


private:
    char method(const size_t idA, const size_t idB) const;
    char fxor(const size_t idA, const size_t idB) const;
    char feSign(const size_t idA, const size_t idB) const;
    char fScalar(const size_t idA, const size_t idB) const;

    void initMethod1( size_t idA, size_t idB );
    void initMethod2( size_t idA, size_t idB );
    void initMethod3( size_t idA, size_t idB );

    long getFullTotal(const state_t &_state) const;

    /**
     * @brief Для каждого спина список соседних, с которыми считаются корреляции
     */
    std::vector< std::forward_list < size_t > > correlationNeighbours;

    /**
     * @brief 
     */
    std::map< std::pair<size_t, size_t>, signed short > correlationValues;
    double correlationPairsNum;

    double _minRange;
    double _maxRange;
    unsigned _methodVar;

    /**
     * @brief Список спинов, для которых считается CorrelationCore
     */
    std::vector<size_t> spins;

    /**
     * @brief В этой переменной сохраняется текущее состояние системы, для которой актуально значение cpOld.  
     */
    state_t currentState;

    mpf_class cp;
    mpf_class cp2;
    long cpOld;

    bool areNeighbours(size_t idA, size_t idB)
    { 
        std::forward_list < size_t > &tmp = correlationNeighbours[idA];
        return std::find(tmp.begin(), tmp.end(), idB) != tmp.end();
    }
};

#endif //CORELLATIONCORE_H