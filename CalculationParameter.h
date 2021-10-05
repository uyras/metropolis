#ifndef CALCULATIONPARAMETER_H
#define CALCULATIONPARAMETER_H

#include <string>
#include <gmpxx.h>
#include "PartArray.h"

class CalculationParameter
{
public:
    CalculationParameter(const std::string & parameterId, PartArray * prototype):
        _debug(false),_parameterId(parameterId),prototype(prototype) {};
    std::string parameterId() const {return this->_parameterId;}
    void setDebug(){this->_debug = true;}
    virtual bool check(unsigned) const = 0;
    virtual void printHeader(unsigned) const = 0;
    virtual bool init(PartArray * sys) {this->sys = sys; return true;};

    virtual void iterate(unsigned id) = 0;
    virtual void incrementTotal() = 0;
    virtual mpf_class getTotal(unsigned) = 0;
    virtual mpf_class getTotal2(unsigned) = 0;
    double getTotalDouble(unsigned steps) { return getTotal(steps).get_d(); };
    double getTotal2Double(unsigned steps) { return getTotal2(steps).get_d(); };

    virtual CalculationParameter * copy() = 0;

protected:
    bool _debug;
    PartArray * sys;
    const PartArray * prototype;

    void prototypeInit(PartArray* prototype){
        this->init(prototype);
        //this->sys = nullptr;
    }

private:
    std::string _parameterId;
};

#endif //CALCULATIONPARAMETER_H