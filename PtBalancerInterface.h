#ifndef PTBLANCERINTERFACE_H
#define PTBLANCERINTERFACE_H

#include <inicpp/inicpp.h>
#include <vector>
#include "misc.h"

class PtBalancerInterface
{
private:
    unsigned each_step;

protected:
    std::vector<double> base_temperatures;

public:
    PtBalancerInterface(unsigned each_step): 
        each_step(each_step) {};
    virtual ~PtBalancerInterface() {};
    void setBaseTemperatures(std::vector<double> base_temperatures){ this->base_temperatures = base_temperatures;};

    unsigned get_each_step(){ return each_step; }
    auto cbeginBase() const { return base_temperatures.cbegin(); }
    auto cendBase() const { return base_temperatures.cend(); }

    /**
     * Внутренняя инициализация. Запускается после того, 
     * как установлены базовые температуры и прошел парсинг конфигурации
     */
    virtual void init() {};
    virtual void parseConfig(inicpp::section configSection) = 0;
    virtual void rebalance() {};
    virtual unsigned size() = 0; //общее число температур
    virtual unsigned sizeBase() = 0; //число базовых температур
    virtual unsigned size(unsigned baseNum) = 0; //число реплик для базовой температуры baseNum
    
    /** получить температуру для базового номера baseNum и реплики replicaNum. 
     * В каждой реплике базовая температура считается нулевой.
     * Если replicaNum будет больше чем size(unsigned baseNum), то выдает ошибку
    */
    virtual double at(unsigned baseNum, unsigned replicaNum) = 0;

    /** 
     * получить температуру из общего списка температур. В списке перемешаны и реплики и базовые температуры
    */
    virtual temp_t at(unsigned temperatureNum) = 0;

    /**
     * @brief Преобразовать базовый номер температуры и номер реплики в порядковый номер температуры
     * 
     * @param baseNum номер базовой температуры
     * @param replicaNum номер реплики в базовой температуре (реплика 0 - это и есть базовая)
     * @return unsigned порядковый номер температуры в общем списке температур
     */
    virtual unsigned to(unsigned baseNum, unsigned replicaNum) = 0;
};


#endif // PTBLANCERINTERFACE_H