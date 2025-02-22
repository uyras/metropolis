#include "PtBalancerManual.h"

void PtBalancerManual::init()
{
    replicaTemperatures.clear();
    replicaTemperatureCounts.resize(base_temperatures.size(),0);
    replicaTemperaturePositions.resize(base_temperatures.size(),0);
    for (int base_temp_num = 0; base_temp_num < base_temperatures.size(); base_temp_num++){
        std:string paramName = "t"+std::to_string(base_temp_num);

        //обязательно добавляем базовую температуру в список
        replicaTemperatures.push_back(base_temperatures[base_temp_num]);
        replicaTemperatureCounts[base_temp_num] += 1;

        for (auto parsedTemp : parsedTemps){
            if (parsedTemp.name == paramName){
                sort(parsedTemp.temps.begin(),parsedTemp.temps.end());
                if (parsedTemp.temps[0] <= base_temperatures[base_temp_num]){ // проверяем чтобы температура реплики была не ниже базовой
                    throw(std::invalid_argument("Temperature "+std::to_string(parsedTemp.temps[0])+" of parameter "+paramName+" should not be less than its base temperature "+std::to_string(base_temperatures[base_temp_num])));
                }
                replicaTemperatures.insert(replicaTemperatures.end(),parsedTemp.temps.begin(),parsedTemp.temps.end());
                replicaTemperatureCounts[base_temp_num] += parsedTemp.temps.size();
            }
        }
        if (base_temp_num>0){
            replicaTemperaturePositions[base_temp_num] = replicaTemperaturePositions[base_temp_num-1]+replicaTemperatureCounts[base_temp_num-1];
        }
        //if (_configSection.contains(paramName)) // TODO Добавить в этот список ещё и базовую температуру
    }
}

void PtBalancerManual::parseConfig(inicpp::section configSection)
{
    //save configSection for later
    for (auto configParam : configSection){
        parsedTemps.push_back({false, configParam.get_name(), configParam.get_list<double>()});
    }
}

void PtBalancerManual::rebalance()
{
}

unsigned PtBalancerManual::size()
{
    return this->replicaTemperatures.size();
}

unsigned PtBalancerManual::sizeBase()
{
    return base_temperatures.size();
}

unsigned PtBalancerManual::size(unsigned baseNum)
{
    return replicaTemperatureCounts[baseNum];
}

double PtBalancerManual::at(unsigned baseNum, unsigned replicaNum)
{
    return this->replicaTemperatures.at(this->to(baseNum, replicaNum));
}

temp_t PtBalancerManual::at(unsigned temperatureNum)
{
    unsigned i;
    for (i=1; i<replicaTemperaturePositions.size(); i++){
        if (temperatureNum<replicaTemperaturePositions[i]) break;
    }
    return {i-1,temperatureNum-replicaTemperaturePositions[i-1],this->replicaTemperatures.at(temperatureNum)};
}

unsigned PtBalancerManual::to(unsigned baseNum, unsigned replicaNum)
{
    if (replicaNum >= size(baseNum)){
        throw(std::invalid_argument("There is no replica "+std::to_string(replicaNum)+" for temperature "+std::to_string(baseNum)));
    } else  {
        return replicaTemperaturePositions[baseNum]+replicaNum;
    }
}
