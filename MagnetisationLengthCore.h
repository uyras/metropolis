#ifndef MAGNETISATIONLENGTHCORE_H
#define MAGNETISATIONLENGTHCORE_H

#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
#include <cstdint>
#include "PartArray.h"
#include "CalculationParameter.h"

class MagnetisationLengthCore: public CalculationParameter
{

public:

    MagnetisationLengthCore(const std::string & parameterId,
        PartArray * prototype,
        const std::vector<uint64_t> & spins);

    virtual bool check(unsigned) const;
    virtual void printHeader(unsigned) const;
    virtual bool init(PartArray * sys);

    virtual void iterate(unsigned id);
    virtual void incrementTotal();
    virtual mpf_class getTotal(unsigned steps){ return this->mv / steps; }
    virtual mpf_class getTotal2(unsigned steps){ return this->mv2 / steps; }
    virtual mpf_class getTotal4(unsigned steps){ return this->mv4 / steps; }

    virtual MagnetisationLengthCore * copy() { return new MagnetisationLengthCore(*this); }

private:
    Vect method(unsigned spinId) const;

    double getFullTotal(Vect & val) const;

    std::vector<uint64_t> spins;

    Vect mOld;
    mpf_class mv;
    mpf_class mv2;
    mpf_class mv4;
};

#endif //MAGNETISATIONLENGTHCORE_H