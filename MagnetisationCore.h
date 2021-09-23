#ifndef MAGNETISATIONCORE_H
#define MAGNETISATIONCORE_H

#include <vector>
#include <string>
#include <map>
#include <gmpxx.h>
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
    virtual double getTotalDouble(unsigned steps){ mpf_class tmp = this->mv / steps; return tmp.get_d();}
    virtual double getTotal2Double(unsigned steps){ mpf_class tmp = this->mv2 / steps; return tmp.get_d();}

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
    double _sumModule;
};

#endif //MAGNETISATIONCORE_H