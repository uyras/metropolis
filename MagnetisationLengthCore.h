#ifndef MAGNETISATIONLENGTHCORE_H
#define MAGNETISATIONLENGTHCORE_H

#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
#include <cstdint>
#include "CalculationParameter.h"
#include "misc.h"

class MagnetisationLengthCore: public CalculationParameter
{

public:
    MagnetisationLengthCore(
        const config_section_t &sect,
        const ConfigManager *conf);

    static string name(){ return "magnetisationlength"; }


    virtual bool check(unsigned) const;
    virtual void printHeader() const;
    virtual void init(const state_t &state);

    virtual void iterate(size_t id);
    virtual void incrementTotal();
    virtual mpf_class getTotal(unsigned steps){ return this->mv / steps; }
    virtual mpf_class getTotal2(unsigned steps){ return this->mv2 / steps; }

    virtual MagnetisationLengthCore * copy() { return new MagnetisationLengthCore(*this); }

private:
    Vect method(unsigned spinId, signed char spinState) const;

    Vect getFullTotal(const state_t &_state) const;

    std::vector<size_t> spins;
    state_t enabledSpins;

    /**
     * @brief В этой переменной сохраняется текущее состояние системы, 
     * для которой актуально значение mOld.  
     */
    state_t currentState;
    Vect mOld;
    mpf_class mv;
    mpf_class mv2;
};

#endif //MAGNETISATIONLENGTHCORE_H