#ifndef WORKER_H
#define WORKER_H

#include <gmpxx.h>
#include <chrono>
#include <random>
#include <memory>
#include <optional>
#include <utility>
#include <omp.h>
#include "misc.h"
#include "defines.h"
#include "ConfigManager.h"

using namespace std;

class Worker {
private:
    unsigned _num; //номер воркера
    unsigned stepsMadeDummy; // сколько шагов сделано с calculateStatistics=False // todo переделать на массив из двух элементов
    double stepsAcceptedDummy; // сколько успешных переворотов для прогревных шагов, считается в МК-шагах поэтому double
    unsigned stepsMadeCalculate; // сколько шагов сделано с calculateStatistics=True
    double stepsAcceptedCalculate; // сколько успешных переворотов для боевых шагов, считается в МК-шагах поэтому double
    int _seed;
    const ConfigManager* config;
    state_t state;
    default_random_engine generator;
    bool _previousCalculateStatistics;

    chrono::steady_clock::time_point _startTimePoint;
    void _startTimer();
    void _stopTimer();
    double fullRefreshEnergy();

public:
    Worker(unsigned num, int seed, const ConfigManager* _config);

    /**
     * @brief Основной процесс который выполняет вычисления. 
     * Потоко-безопасный, работает с локальной копией магнитной системы.
     * 
     * Если включен параллельный отжиг, то делает один "холостой" шаг вне зависимости от 
     * значения ключа calculateStatistics, то есть не учитывает этот шаг в статистике, в числе шагов,
     * после него потом не нужно делать полный пересчет всех конфигураций.
     * "Холостой" шаг выполняется первым при запуске функции work(), и только когда config->pt_enabled() == true,
     * Холостой шаг не делается при первом запуске work(), то есть когда stepsMadeDummy+stepsMadeCalculate>0.
     * 
     * @param steps Число МК-шагов которые должна сделать функция
     * @param calculateStatistics если true, то считает статистику шагов, но требует больше процессорного времени.
     * Значение false используется для прогревных шагов
     * @param lowestEnergy Задается минимальное значение энергии. 
     * Если найдется энергия ниже, то worker досрочно завершится при следующем полном обновлении энергии,
     * обновление энергии выполняется каждые FULL_REFRESH_EVERY шагов
     * @return optional< pair<double,string> > возвращает пустое значение если не нашлось энергии ниже, 
     * либо пару - минимальную энергию и соответствующую ей конфигурацию в виде строки
     */
    optional< pair<double,state_t> > work(unsigned steps, bool calculateStatistics, double lowestEnergy);

    void printout(temp_t temperature);

    mpf_class e;
    mpf_class e2;
    double eActual; // энергия текущей конфигурации системы. Обновляется итеративно при каждом перевороте. иногда полностью пересчитывается
    chrono::milliseconds duration;
    vector<std::unique_ptr<CalculationParameter>> calculationParameters;

    unsigned ptSuccessExchangesDummy; // количество успешных обменов с температурой ниже при прогреве
    unsigned ptSuccessExchangesCalculate; // количество успешных обменов с температурой ниже при основных шагах


    static bool exchange(shared_ptr<Worker> w1, shared_ptr<Worker> w2, double dBeta);
};

#endif //WORKER_H