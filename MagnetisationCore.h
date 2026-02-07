#ifndef MAGNETISATIONCORE_H
#define MAGNETISATIONCORE_H

#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
#include <cstdint>
#include "CalculationParameter.h"
#include "misc.h"

class MagnetisationCore: public CalculationParameter
{

public:
    MagnetisationCore(
        const config_section_t &sect,
        const ConfigManager *conf);

    static string name(){ return "magnetisation"; }


    virtual bool check(unsigned) const;
    virtual void printHeader() const;
    virtual void init(const state_t &state);

    virtual void iterate(size_t id);
    virtual void incrementTotal();
    virtual mpf_class getTotal(unsigned steps){ return this->mv / steps; }
    virtual mpf_class getTotal2(unsigned steps){ return this->mv2 / steps; }

    virtual MagnetisationCore * copy() { return new MagnetisationCore(*this); }

private:
    double method(unsigned spinId, signed char spinState) const;

    double getFullTotal(const state_t &_state) const;

    Vect vector;
    std::vector<size_t> spins;
    std::vector< double > magnetisationValues;

    /**
     * @brief В этой переменной сохраняется текущее состояние системы, для которой актуально значение cpOld.  
     */
    state_t currentState;
    double mOld;
    mpf_class mv;
    mpf_class mv2;
    double _sumModule;
};

#endif //MAGNETISATIONCORE_H