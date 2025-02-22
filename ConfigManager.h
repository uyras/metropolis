#ifndef CONFIGMANAGER_H
#define CONFIGMANAGER_H

#include "defines.h"

#include <inicpp/inicpp.h>
#include <vector>
#include <map>
#include <iostream>
#include <omp.h>
#include "PartArray.h"
#include "Vect.h"
#include "CommandLineParameters.h"
#include "CalculationParameter.h"
#include "CorrelationCore.h"
#include "CorrelationPointCore.h"
#include "MagnetisationCore.h"
#include "MagnetisationLengthCore.h"
#include "misc.h"
#include "PtBalancerFactory.h"
#include "PtBalancerDefault.h"

static const std::map<std::string, unsigned> methods = 
    {{"xor",0},{"energy",1},{"scalar",2}};


double hamiltonianDipolarPBC(Part* a, Part* b);
double hamiltonianDipolarCSV(Part* a, Part* b);
Vect radiusPBC(const Vect& a, const Vect& b);

class ConfigManager
{
public:
    ConfigManager(const CommandLineParameters & commandLineParameters, const inicpp::config & iniconfig);
    bool check_config();

    void printHeader();
    void getParameters(std::vector< std::unique_ptr< CalculationParameter > > &);

    const PartArray & getSystem() const {return this->system;}
    void saveSystem(std::string filename){ return this->system.save(filename); }
    void applyState(string s);


    int getSeed() const { return this->seed; }
    unsigned N() const { return this->system.size(); }
    unsigned getHeatup() const { return this->heatup; }
    unsigned getCalculate() const { return this->calculate; }
    std::string getSysfile() const { return this->sysfile; }
    const Vect & getField() const { return this->field; }
    bool isPBC() const { return this->pbc; }
    bool isCSV() const { return this->_csv;}
    bool isRestart() const {return this->restart; }
    double getRestartThreshold() const {return this->restartThreshold; }
    std::string getNewGSFilename() const {return this->newGSFilename; }
    inline unsigned getSaveStates() const { return this->saveStates; }
    std::string getSaveStateFileName(int temperature, int step) const{ return this->saveStateFileBasename+"_"+std::to_string(temperature)+"_"+std::to_string(step)+".mfsys"; }

    bool debug = false;
    int threadCount=0;
    static Vect size;
    static vector < vector < double > > energyTable;
    double deltaEnergy;

    static void setPBCEnergies(PartArray & sys);
    static void setCSVEnergies(PartArray & sys);

    bool pt_enabled(){ return this->temperatures->get_each_step(); } // включен ли параллельный отжиг (отключен когда each_step==0)
    PtBalancerInterface * temperatures = nullptr; // указатель на класс балансира параллельного отжига

private:
    ConfigManager(){};

    std::string sysfile;
    bool pbc = 0;
    bool _csv = 0;
    unsigned heatup = 0;
    unsigned calculate = 0;
    double range = 0;
    int seed = 0;
    Vect field;
    bool restart = true;
    double restartThreshold = 1e-6;
    unsigned saveStates = 0;
    std::string saveStateFileBasename;
    std::string newGSFilename;
    std::vector<std::unique_ptr< CalculationParameter > > parameters;
    PartArray system;

    static Vect strToVect(std::string val){
        Vect target;
        if (val=="x" || val=="X"){
            target = Vect(1,0,0);
        } else {
            if (val=="y" || val=="Y") 
                target = Vect(0,1,0); 
            else {
                if (val=="z" || val=="Z")
                    target = Vect(0,0,1);
                else
                {
                    auto pos = val.find('|');
                    if (pos!=std::string::npos){
                        double x = stod(val.substr(0,pos));
                        auto pos_second = val.find('|',pos+1);
                        if (pos_second!=std::string::npos){
                            double y = stod(val.substr(pos+1,pos_second));
                            double z = stod(val.substr(pos_second+1));
                            target = Vect(x,y,z);
                        } else {
                            double y = stod(val.substr(pos+1));
                            target = Vect(x,y,0);
                        }
                    } else {
                        throw(std::invalid_argument("Vector format is x|y or x|y|z"));
                    }
                }
            }
        }
        return target;
    }
};

#endif //CONFIGMANAGER_H