#include "PtBalancerFactory.h"

PtBalancerInterface * PtBalancerFactory(std::string type, unsigned each_step){
    if (type == "manual"){
        return new PtBalancerManual(each_step);
    }
    /*else if (type == "random"){ //добавляем для каждого алгоритма балансира свой if
        return new PtBalancerRandom();
    }*/
    else return nullptr;
}

PtBalancerInterface* PtParseConfig(const inicpp::config &iniconfig){
    PtBalancerInterface* balancer_instance = nullptr;
    inicpp::section sect = iniconfig["parallel_tempering"];
    int balancer_each_step;
    if (sect.contains("each_step")){
        balancer_each_step = sect["each_step"].get<inicpp::unsigned_ini_t>();
    } else {
        throw(std::invalid_argument("parameter 'each_step' in section 'parallel_tempering' of .ini file is mandatory"));
    }
    if (sect.contains("balancer")){
        auto balancer_name = sect["balancer"].get<inicpp::string_ini_t>();
        balancer_instance = PtBalancerFactory(balancer_name,balancer_each_step);
        if (!balancer_instance)
            throw(std::invalid_argument("Parallel tempering balancer called '"+balancer_name+"' not found"));
  
        string ini_section_name = "parallel_tempering.balancer."+balancer_name;
        balancer_instance->parseConfig(iniconfig[ini_section_name]); // парсим специфичную секцию для балансира
    } else {
        throw(std::invalid_argument("parameter 'balancer' in section 'parallel_tempering' of .ini file is mandatory"));
    }
    return balancer_instance;
}