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
    vector<size_t> replicaTemperatureCounts;
    vector<size_t> replicaTemperaturePositions; //кумулятивная сумма от replicaTemperatureCounts

    struct parserParameterItem {
        bool checked;
        string name;
        vector<double> temps;
    };
    vector<parserParameterItem> parsedTemps;

public:
    PtBalancerManual(unsigned each_step):PtBalancerInterface(each_step){};
    ~PtBalancerManual(){};

    string name(){ return "manual"; }

    void init();
    void parseConfig(inicpp::section configSection);
    void rebalance();
    void printHeader();
    size_t size(); //общее число температур
    size_t sizeBase(); //число базовых температур
    size_t size(size_t baseNum); //число реплик для базовой температуры baseNum
    
    /** получить температуру для базового номера baseNum и реплики replicaNum. 
     * В каждой реплике базовая температура считается нулевой.
     * Если replicaNum будет больше чем size(size_t baseNum), то выдает ошибку
    */
    double at(size_t baseNum, size_t replicaNum);
    temp_t at(size_t temperatureNum);

    size_t to(size_t baseNum, size_t replicaNum);
};


#endif // PTBLANCERMANUAL_H