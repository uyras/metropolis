#ifndef MAGNETISATIONCORE_H
#define MAGNETISATIONCORE_H

#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
#include <cstdint>
#include "PartArray.h"
#include "CalculationParameter.h"

class MagnetisationCore: public CalculationParameter
{

public:

    MagnetisationCore(const std::string & parameterId,
        PartArray * prototype,
        const Vect & vector, 
        const std::vector<uint64_t> & spins);

    virtual bool check(unsigned) const;
    virtual void printHeader(unsigned) const;
    virtual bool init(PartArray * sys);

    virtual void iterate(unsigned id);
    virtual void incrementTotal();
    virtual mpf_class getTotal(unsigned steps){ return this->mv / steps; }
    virtual mpf_class getTotal2(unsigned steps){ return this->mv2 / steps; }
    virtual mpf_class getTotal4(unsigned steps){ return this->mv4 / steps; }

    virtual MagnetisationCore * copy() { return new MagnetisationCore(*this); }

    void setModule(bool module);

private:
    double method(unsigned spinId, const PartArray * _sys) const;

    double getFullTotal(const PartArray * _sys) const;

    Vect vector;
    std::vector<uint64_t> spins;
    std::vector< double > magnetisationValues;

    double mOld;
    mpf_class mv;
    mpf_class mv2;
    mpf_class mv4;
    double _sumModule;
};

#endif //MAGNETISATIONCORE_H