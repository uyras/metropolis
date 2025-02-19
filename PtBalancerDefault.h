#ifndef PTBLANCERDEFAULT_H
#define PTBLANCERDEFAULT_H

#include "PtBalancerInterface.h"
#include <string>

class PtBalancerDefault : public PtBalancerInterface
{
private:
    void setBaseTemperatures(std::vector<double>) {};

public:
    PtBalancerDefault(std::vector<double> baseTemperatures):PtBalancerInterface(0){
        PtBalancerInterface::setBaseTemperatures(baseTemperatures);
    };
    ~PtBalancerDefault(){};

    void parseConfig(inicpp::section configSection) {};
    unsigned size() { return this->base_temperatures.size(); }; //общее число температур
    unsigned sizeBase() { return this->size(); }; //число базовых температур
    unsigned size(unsigned baseNum) { return 1; }; //число реплик для базовой температуры baseNum
    
    /** получить температуру для базового номера baseNum и реплики replicaNum. 
     * В каждой реплике базовая температура считается нулевой.
     * Если replicaNum будет больше чем size(unsigned baseNum), то выдает ошибку
    */
    double at(unsigned baseNum, unsigned replicaNum) {
        if (replicaNum !=0)
            throw(std::invalid_argument("There is no replica "+std::to_string(replicaNum)+" for temperature "+std::to_string(baseNum)));
        else return this->base_temperatures.at(baseNum);
    };

    temp_t at(unsigned temperatureNum) {
        return {temperatureNum,0,this->base_temperatures.at(temperatureNum)};
    };
};


#endif // PTBLANCERDEFAULT_H