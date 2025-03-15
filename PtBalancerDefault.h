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

    string name(){ return "default"; }

    void parseConfig(inicpp::section configSection) {};
    size_t size() { return this->base_temperatures.size(); }; //общее число температур
    size_t sizeBase() { return this->size(); }; //число базовых температур
    size_t size(size_t baseNum) { return 1; }; //число реплик для базовой температуры baseNum
    
    /** получить температуру для базового номера baseNum и реплики replicaNum. 
     * В каждой реплике базовая температура считается нулевой.
     * Если replicaNum будет больше чем size(unsigned baseNum), то выдает ошибку
    */
    double at(size_t baseNum, size_t replicaNum) {
        if (replicaNum !=0)
            throw(std::string("There is no replica "+std::to_string(replicaNum)+" for temperature "+std::to_string(baseNum)));
        else return this->base_temperatures.at(baseNum);
    };

    temp_t at(size_t temperatureNum) {
        return {temperatureNum,0,this->base_temperatures.at(temperatureNum)};
    };

    size_t to(size_t baseNum, size_t replicaNum){ return baseNum; }
};


#endif // PTBLANCERDEFAULT_H