#ifndef CONFIGMANAGER_H
#define CONFIGMANAGER_H

#include "defines.h"

#include <inicpp/inicpp.h>
#include <vector>
#include <map>
#include <iostream>
#include <numeric>
#include <memory>
#include <omp.h>
#include "MagneticSystem.h"
#include "CommandLineParameters.h"
#include "CalculationParameter.h"
#include "CorrelationCore.h"
//#include "CorrelationPointCore.h"
#include "MagnetisationCore.h"
#include "MagnetisationLengthCore.h"
#include "misc.h"
#include "Hamiltonians.h"
#include "PtBalancerFactory.h"
#include "PtBalancerDefault.h"



class ConfigManager
{
private:
    std::string sysfile;
    unsigned heatup = 0;
    unsigned calculate = 0;
    int seed = 0;
    Vect field;
    bool restart = true;
    double restartThreshold = 1e-6;
    unsigned saveStates = 0;
    std::string saveStateFileBasename;
    std::string newGSFilename;
    std::vector<std::unique_ptr< CalculationParameter > > parameters;

    /**
     * @brief Если в конфигурации есть параметры, на основе которых должно создаваться несколько CalculationParameter,
     * то эта функция разбивает секцию конфига на несколько виртуальных секций, в каждой из которых присутствует
     * всего по одному параметру.
     * К имени нового параметра добавляется индекс "_N", где N - порядковый номер.
     * 
     * @param origin Оригинальная секция
     * @param params список параметров, по которым разбивать.
     * например, указаны 2 параметра foo, bar и для них указаны значения [1,2] и ['a','b']
     * тогда будет такой список секций: [1,'a'],[2,'a'],[1,'b'],[2,'b']
     * @return std::vector< inicpp::section > 
     */
    std::vector< config_section_t > explode(const inicpp::section& origin, const vector<string>& params);

    template<class C> size_t registerCalculationParameter(inicpp::config & conf, vector<string> explode_params);

public:
    shared_ptr<MagneticSystem> system;
    bool debug = false;
    int threadCount=0;
    double deltaEnergy;
    PtBalancerInterface * temperatures = nullptr; // указатель на класс балансира параллельного отжига

    ConfigManager(int argc, char *argv[]);
    bool check_config();

    void printHeader();
    void printColumnNames();
    void getParameters(std::vector< std::unique_ptr< CalculationParameter > > &);

    int getSeed() const { return this->seed; }
    inline size_t N() const { return this->system->N(); }
    unsigned getHeatup() const { return this->heatup; }
    unsigned getCalculate() const { return this->calculate; }
    std::string getSysfile() const { return this->sysfile; }
    const Vect & getField() const { return this->field; }
    bool isRestart() const {return this->restart; }
    double getRestartThreshold() const {return this->restartThreshold; }
    std::string getNewGSFilename() const {return this->newGSFilename; }
    inline unsigned getSaveStates() const { return this->saveStates; }
    std::string getSaveStateFileName(int temperature, int step) const { 
        return this->saveStateFileBasename+"_"+std::to_string(temperature)+"_"+std::to_string(step)+".mfsys"; 
    }
    bool pt_enabled() const { return this->temperatures->get_each_step(); } // включен ли параллельный отжиг (отключен когда each_step==0)
    
};

#endif //CONFIGMANAGER_H
