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
                    throw(std::string("Temperature "+std::to_string(parsedTemp.temps[0])+" of parameter "+paramName+" should not be less than its base temperature "+std::to_string(base_temperatures[base_temp_num])));
                }
                replicaTemperatures.insert(replicaTemperatures.end(),parsedTemp.temps.begin(),parsedTemp.temps.end());
                replicaTemperatureCounts[base_temp_num] += parsedTemp.temps.size();
            }
        }
        if (base_temp_num>0){
            replicaTemperaturePositions[base_temp_num] = replicaTemperaturePositions[base_temp_num-1]+replicaTemperatureCounts[base_temp_num-1];
        }
    }
}

void PtBalancerManual::parseConfig(inicpp::section configSection)
{
    //save configSection for later
    for (auto configParam : configSection){
        parsedTemps.push_back({false, configParam.get_name(), configParam.get_list<double>()});
    }
}

void PtBalancerManual::printHeader()
{
    printf("# PT replica temperature list:\n");
    for (size_t r=0; r<this->sizeBase(); r++){
        printf("# rep. %4lu: ",r);
        size_t from = replicaTemperaturePositions[r];
        size_t cnt = replicaTemperatureCounts[r];
        size_t to = from + cnt;
        printf("%e", replicaTemperatures[from]);
        for (size_t bt=from+1; bt<to; bt++){
            printf(", %e", replicaTemperatures[bt]);
        }
        printf(" (%lu pcs. from %e to %e)\n", cnt,
        std::min_element(replicaTemperatures.cbegin()+from,replicaTemperatures.cbegin()+to).operator*(),
        std::max_element(replicaTemperatures.cbegin()+from,replicaTemperatures.cbegin()+to).operator*());
    }
}

size_t PtBalancerManual::size()
{
    return this->replicaTemperatures.size();
}

size_t PtBalancerManual::sizeBase()
{
    return base_temperatures.size();
}

size_t PtBalancerManual::size(size_t baseNum)
{
    return replicaTemperatureCounts[baseNum];
}

double PtBalancerManual::at(size_t baseNum, size_t replicaNum)
{
    return this->replicaTemperatures.at(this->to(baseNum, replicaNum));
}

temp_t PtBalancerManual::at(size_t temperatureNum)
{
    unsigned i;
    for (i=1; i<replicaTemperaturePositions.size(); i++){
        if (temperatureNum<replicaTemperaturePositions[i]) break;
    }
    return {i-1,temperatureNum-replicaTemperaturePositions[i-1],this->replicaTemperatures.at(temperatureNum)};
}

size_t PtBalancerManual::to(size_t baseNum, size_t replicaNum)
{
    if (replicaNum >= size(baseNum)){
        throw(std::string("There is no replica "+std::to_string(replicaNum)+" for temperature "+std::to_string(baseNum)));
    } else  {
        return replicaTemperaturePositions[baseNum]+replicaNum;
    }
}
