#ifndef PTBALANCERFACTORY_H
#define PTBALANCERFACTORY_H

#include <string>
#include "PtBalancerInterface.h"
#include "ConfigManager.h"
#include <inicpp/inicpp.h>

class ConfigManager;

#include "PtBalancerManual.h"
//тут добавляем include для всех используемых балансеров, по аналогии с manual


/**
 * Factory function that creates and returns a pointer to a PtBalancerInterface
 * implementation based on the provided type.
 *
 * @param type A string representing the type of balancer to create. For example,
 * "manual" for a manual temperature balancer. Additional types can be added
 * with corresponding conditions.
 * 
 * @return A pointer to a PtBalancerInterface object of the specified type, or
 * nullptr if the type is not recognized.
 */
PtBalancerInterface * PtBalancerFactory(std::string type, unsigned each_step);

/**
 * Функция которая считывает параметры из ini-конфиг файла
 */
PtBalancerInterface * PtParseConfig(const inicpp::config &iniconfig);

#endif  // PTBALANCERFACTORY_H