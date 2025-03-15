#ifndef PTBLANCERINTERFACE_H
#define PTBLANCERINTERFACE_H

#include <inicpp/inicpp.h>
#include <vector>
#include <string>
#include "misc.h"

class Worker; // хэдер с ним нельзя подключать из-за циклической зависимости 
// PtBalancerInterface->Worker->ConfigManager->PtBalancerDefault->PtBalancerInterface

class PtBalancerInterface
{
private:
    unsigned each_step;

protected:
    std::vector<double> base_temperatures;

public:
    PtBalancerInterface(unsigned each_step): 
        each_step(each_step) {};
    virtual ~PtBalancerInterface() {};
    void setBaseTemperatures(std::vector<double> base_temperatures){ this->base_temperatures = base_temperatures;};

    unsigned get_each_step() const { return each_step; }
    auto cbeginBase() const { return base_temperatures.cbegin(); } //todo подумать об удалении этого списка
    auto cendBase() const { return base_temperatures.cend(); }

    /**
     * @brief Имя класса балансера
     * 
     * @return string 
     */
    virtual string name() = 0;

    /**
     * @brief Функция отвечает за внутреннюю инициализацию алгоритма балансировки
     * Запускается после того, как установлены базовые температуры и прошел парсинг конфигурации
     * Но дизайн функции должен быть таким, чтобы повторный запуск init() повторно инициализировал класс
     */
    virtual void init() {};

    /**
     * @brief Тут выполняется парсинг раздела конфигурации, который относится к конкретному балансиру.
     * Эта функция вызывается после создания экземпляра класса, и до функции init().
     * 
     * @param configSection Секция конфигурации, относящаяся к данному балансиру.
     */
    virtual void parseConfig(inicpp::section configSection) = 0;

    /**
     * @brief Функция запускается сразу после обмена конфигурациями, только во время прогрева системы. 
     * В ней можно произвести ребалансировку температуры. Статистику по всем вычисляемым параметрам можно получить из
     * аргумента workers
     * 
     * @param workers массив объектов магнитной системы. Размер массива равен общему числу температур (параметр size)
     */
    virtual void rebalance(const vector <shared_ptr<Worker>> & workers) {};

    /**
     * @brief Выводит на экран (в поток cout) информационные данные о конфигурации балансира.
     * 
     * Запускается после чтения конфигурации, во время вывода общего хэдера на экран.
     * Все строки должны начинаться со знака #, чтобы не мешать основному выводу программы.
     * Вывод должен заканчиваться символом новой строки.
     */
    virtual void printHeader() {};

    /**
     * @brief Функция должна возвращать общее число температур.
     * На основе этого значения создается нужное количество worker-ов (потоков процессора).
     * 
     * @return size_t количество температур
     */
    virtual size_t size() = 0;

    /**
     * @brief Функция должна возвращать количество базовых температур.
     * Сами температуры по умолчанию определены в параметре this->base_temperatures.
     * 
     * @return size_t количество температур
     */
    virtual size_t sizeBase() = 0; //число базовых температур

    /**
     * @brief Число реплик для любой базовой температуры.
     * Сама базовая температура тоже считается репликой.
     * Поэтому, если для базовой температуры нет реплик, то должно возращаться значение 1
     * 
     * @param baseNum Номер базовой температуры. Нумерация с ноля.
     * @return size_t Число реплик для базовой температуры. 1 и более.
     */
    virtual size_t size(size_t baseNum) = 0; //число реплик для базовой температуры baseNum
    
    /** получить температуру для базового номера baseNum и реплики replicaNum. 
     * В каждой реплике базовая температура считается нулевой.
     * Если replicaNum будет больше чем size(size_t baseNum), то выдает exception
    */
    virtual double at(size_t baseNum, size_t replicaNum) = 0;

    /**
     * @brief получить все температуры (реплики и базовые) из общего списка.
     * Порядок нумерации и распределения не важен, главное чтобы
     * выполнение было согласовано с функцией to(i,j),
     * то есть выполнялось условие:
     * at( to(i,j) ) == at(i,j)
     * 
     * @param temperatureNum порядковый номер температуры в общем списке в диапазоне от 0 до size()
     * @return temp_t структура, в которой записана температура, ее базовый номер и номер реплики
     */
    virtual temp_t at(size_t temperatureNum) = 0;

    /**
     * @brief Преобразовать базовый номер температуры и номер реплики в порядковый номер температуры
     * 
     * @param baseNum номер базовой температуры
     * @param replicaNum номер реплики в базовой температуре (реплика 0 - это и есть базовая)
     * @return unsigned порядковый номер температуры в общем списке температур
     */
    virtual size_t to(size_t baseNum, size_t replicaNum) = 0;
};


#endif // PTBLANCERINTERFACE_H