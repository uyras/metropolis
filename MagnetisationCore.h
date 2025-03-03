#ifndef MAGNETISATIONCORE_H
#define MAGNETISATIONCORE_H

#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
#include <cstdint>
#include "CalculationParameter.h"

class MagnetisationCore: public CalculationParameter
{

public:

    MagnetisationCore(const std::string & parameterId,
        shared_ptr<MagneticSystem> prototype,
        const Vect & vector, 
        const std::vector<size_t> & spins);

    virtual bool check(unsigned) const;
    virtual void printHeader(unsigned) const;
    virtual bool init(state_t state);

    virtual void iterate(unsigned id);
    virtual void incrementTotal();
    virtual mpf_class getTotal(unsigned steps){ return this->mv / steps; }
    virtual mpf_class getTotal2(unsigned steps){ return this->mv2 / steps; }

    virtual MagnetisationCore * copy() { return new MagnetisationCore(*this); }

    void setModule(bool module);

private:
    double method(unsigned spinId, const PartArray * _sys) const;

    double getFullTotal(const PartArray * _sys) const;

    Vect vector;
    std::vector<size_t> spins;
    std::vector< double > magnetisationValues;

    double mOld;
    mpf_class mv;
    mpf_class mv2;
    double _sumModule;
};

#endif //MAGNETISATIONCORE_H