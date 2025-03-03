#ifndef CALCULATIONPARAMETER_H
#define CALCULATIONPARAMETER_H

#include <string>
#include <memory>
#include <gmpxx.h>
#include "misc.h"
#include "MagneticSystem.h"

/**
 * @brief Класс используется как шаблон для дополнительных вычисляемых параметров
 * Вычисляемый параметр - значение магнитной системы, которое собирается статистически во время МК-семплирования
 * Каждый вычисляемый параметр работает независимо для каждой температуры.
 * Ярким примером является намагниченность, определяемая как сумма векторов магнитных моментов 
 * вдоль определенного направления. На этом примере будет дальнейшее объяснение.
 * 
 * Параметр должен собирать среднее значение, и квадрат среднего значения (для удобного вычисления производных).
 * Так как расчет параметра бывает вычислительно-накладным, 
 * то добавлена возможность его итеративного вычисления.
 * 
 * При каждом успешном перевороте спина выполняется функция iterate(id) где id - номер перевернутого спина.
 * Шаг Метрополиса - N попыток перевернуть случайный спин. После каждого шага запускается
 * функция incrementTotal(). Она позволяет
 * 
 */
class CalculationParameter
{
public:

    /**
     * @brief Создает новый объект с вычисляемым параметром
     * 
     * @param parameterId id параметра. Обязателен для всех параметров.
     * @param prototype Прототип магнитной системы. Пока не понятно нафиг нужен
     */
    CalculationParameter(const std::string & parameterId, shared_ptr<MagneticSystem> prototype):
        _debug(false),_parameterId(parameterId),prototype(prototype) {};
    std::string parameterId() const {return this->_parameterId;}
    void setDebug(){this->_debug = true;}
    virtual bool check(unsigned) const = 0;
    virtual void printHeader(unsigned) const = 0;

    /**
     * @brief Инициирует внутренние значения класса для того чтобы работал iterate.
     * Запускается когда в системе изменилось более одного спина. Функция может быть запущена несколько раз
     * в любой момент работы программы. Например, при обмене репликами.
     * 
     * @param state конфигурация магнитной системы, для которой инициировать состояние
     */
    virtual void init(state_t state) = 0;

    virtual void iterate(unsigned id) = 0; // запускается при каждом успешном перевороте спина
    virtual void incrementTotal() = 0; // запускается после каждого шага Метрополиса
    virtual mpf_class getTotal(unsigned) = 0; //todo параметр функции - число шагов. Нужно избавиться от него, чтобы этот класс считал его сам
    virtual mpf_class getTotal2(unsigned) = 0;
    double getTotalDouble(unsigned steps) { return getTotal(steps).get_d(); };
    double getTotal2Double(unsigned steps) { return getTotal2(steps).get_d(); };

    virtual CalculationParameter * copy() = 0;

    /**
     * @brief This function is running after completing calculations for the temperature,
     * and after total averages are printed on the screen, i.e. right before this parameter is deleted from memory.
     * You can save desired information to the file in this function
     * 
     * @param num sequential number of temperature in the list for which it calculates.
     */
    virtual void save(unsigned num){};

protected:
    bool _debug;
    shared_ptr<MagneticSystem> prototype;

private:
    std::string _parameterId;
};

#endif //CALCULATIONPARAMETER_H