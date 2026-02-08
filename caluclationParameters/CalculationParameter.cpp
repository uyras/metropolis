#include "CalculationParameter.h"
#include "ConfigManager.h"

CalculationParameter::CalculationParameter(const config_section_t &sect, const ConfigManager *conf):
sys(conf->system)
{
    this->_debug = (sect.data.contains("debug") && sect.data["debug"].get<inicpp::boolean_ini_t>()) ? true : false;

    const std::string parameterString = sect.name;
    std::size_t colpos=parameterString.find_first_of(':');
    if (colpos==std::string::npos) throw string("For the calculation parameter section should contain ':' symbol");
    this->_parameterId = parameterString.substr(colpos+1);
};
