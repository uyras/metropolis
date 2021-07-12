#ifndef CALCULATIONPARAMETER_H
#define CALCULATIONPARAMETER_H

#include <string>
#include "PartArray.h"

class CalculationParameter
{
public:
    CalculationParameter(const std::string & parameterId, PartArray * prototype):
        _inited(false),_debug(false),_parameterId(parameterId),prototype(prototype) {};
    std::string parameterId() const {return this->_parameterId;}
    void setDebug(){this->_debug = true;}
    virtual bool check(unsigned) const = 0;
    virtual void printHeader(unsigned) const = 0;
    virtual bool init(PartArray * sys) {return true;};
    virtual bool unInit() = 0;
    bool isInited() { return _inited; }

    virtual void iterate(unsigned id) = 0;
    virtual void incrementTotal() = 0;
    virtual double getTotalDouble(unsigned) = 0;
    virtual double getTotal2Double(unsigned) = 0;

    virtual CalculationParameter * copy() = 0;

protected:
    bool _inited;
    bool _debug;
    PartArray * sys;
    const PartArray * prototype;

    void prototypeInit(PartArray* prototype){
        this->init(prototype);
        this->sys = nullptr;
        this->_inited = false;
    }

private:
    std::string _parameterId;
};

#endif //CALCULATIONPARAMETER_H