#ifndef CALCULATIONPARAMETER_H
#define CALCULATIONPARAMETER_H

#include <string>
#include <gmpxx.h>
#include "PartArray.h"

class CalculationParameter
{
public:
    CalculationParameter(const std::string & parameterId, PartArray * prototype):
        _debug(false), _binder(false), _parameterId(parameterId), prototype(prototype) {};
    std::string parameterId() const {return this->_parameterId;}
    void setDebug(){this->_debug = true;}
    void setBinder(bool val){this->_binder = val;}
    virtual bool check(unsigned) const = 0;
    virtual void printHeader(unsigned) const = 0;
    virtual bool init(PartArray * sys) {this->sys = sys; return true;};

    virtual void iterate(unsigned id) = 0; // запускается при каждом успешном перевороте спина
    virtual void incrementTotal() = 0; // запускается после каждого шага Метрополиса
    virtual mpf_class getTotal(unsigned) = 0;
    virtual mpf_class getTotal2(unsigned) = 0;
    virtual mpf_class getTotal4(unsigned) = 0;
    double getTotalDouble(unsigned steps) { return getTotal(steps).get_d(); };
    double getTotal2Double(unsigned steps) { return getTotal2(steps).get_d(); };

    virtual CalculationParameter * copy() = 0;

    /**
     * @brief This function is running after completing calculations for the temperature,
     * and after total averages are printed on the screen, i.e. right before this parameter is deleted from memory.
     * You can save desired information to the file in this function
     * 
     * @param num sequential number of temperature in the list for which it calculates.
     */
    virtual void save(unsigned num){};

protected:
    bool _debug;
    bool _binder;
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