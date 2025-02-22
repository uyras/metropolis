#ifndef PTBLANCERMANUAL_H
#define PTBLANCERMANUAL_H

#include <string>
#include <inicpp/inicpp.h>
#include <vector>
#include <algorithm>
#include "PtBalancerInterface.h"

using namespace std;
class PtBalancerManual : public PtBalancerInterface
{
private:
    vector<double> replicaTemperatures;
    vector<unsigned> replicaTemperatureCounts;
    vector<unsigned> replicaTemperaturePositions; //кумулятивная сумма от replicaTemperatureCounts

    struct parserParameterItem {
        bool checked;
        string name;
        vector<double> temps;
    };
    vector<parserParameterItem> parsedTemps;

public:
    PtBalancerManual(unsigned each_step):PtBalancerInterface(each_step){};
    ~PtBalancerManual(){};

    void init();
    void parseConfig(inicpp::section configSection);
    void rebalance();
    unsigned size(); //общее число температур
    unsigned sizeBase(); //число базовых температур
    unsigned size(unsigned baseNum); //число реплик для базовой температуры baseNum
    
    /** получить температуру для базового номера baseNum и реплики replicaNum. 
     * В каждой реплике базовая температура считается нулевой.
     * Если replicaNum будет больше чем size(unsigned baseNum), то выдает ошибку
    */
    double at(unsigned baseNum, unsigned replicaNum);
    temp_t at(unsigned temperatureNum);

    unsigned to(unsigned baseNum, unsigned replicaNum);
};


#endif // PTBLANCERMANUAL_H