#ifndef PARAMETERFACTORY_H
#define PARAMETERFACTORY_H

#include <unordered_map>
#include <memory>
#include <functional>
#include <inicpp/inicpp.h>
#include <vector>
#include <string>

class ConfigManager; //защита от cross-include
class CalculationParameter;

class ParameterFactory
{
public:
    // тип данных - функция, которая производит экземпляры класса
    using CreatorFunc = std::function<std::unique_ptr<CalculationParameter>(
        const config_section_t&, const ConfigManager*
    )>;

    static ParameterFactory& instance() {
        static ParameterFactory factory;
        return factory;
    }

    void registerParameter(const std::string& name, CreatorFunc creator, vector<string> explode_params) {
        creators_[name] = creator;
        explodes_[name] = explode_params;
    }

    std::unique_ptr<CalculationParameter> createParameter(
        const std::string& name, 
        const config_section_t& sect, 
        const ConfigManager* conf) {
        
        auto it = creators_.find(name);
        if (it != creators_.end()) {
            return it->second(sect, conf);  // Передаем параметры в конструктор
        }
        return nullptr;
    }

    const vector<string>& getExplodeList(const std::string& name) const {
        return explodes_.at(name);
    }

    // Получение списка доступных параметров
    std::vector<std::string> getAvailableParameters() const {
        std::vector<std::string> types;
        for (const auto& pair : creators_) {
            types.push_back(pair.first);
        }
        return types;
    }

private:
    ParameterFactory() = default;
    std::unordered_map<std::string, CreatorFunc> creators_;
    std::unordered_map<std::string, vector<string>> explodes_;
};



#endif //PARAMETERFACTORY_H