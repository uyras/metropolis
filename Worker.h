#ifndef WORKER_H
#define WORKER_H

#include <gmpxx.h>
#include <chrono>
#include <random>
#include <memory>
#include <optional>
#include <utility>
#include <omp.h>
#include "PartArray.h"
#include "misc.h"
#include "defines.h"
#include "ConfigManager.h"

using namespace std;

class Worker {
private:
    unsigned _num; //номер воркера
    unsigned stepsMade;
    int _seed;
    const ConfigManager* config;
    PartArray sys;
    default_random_engine generator;

    chrono::steady_clock::time_point _startTimePoint;
    void _startTimer();
    void _stopTimer();
    double fullRefreshEnergy();

public:
    Worker(unsigned num, int seed, const ConfigManager* _config);
    Worker(const Worker&) = delete; //явно отключаем конструктор копирования, чтобы добавление этого класса в вектор работало корректно
    Worker(Worker&&) = default;

    /**
     * @brief Основной процесс который выполняет вычисления. 
     * Потоко-безопасный, работает с локальной копией магнитной системы
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
    optional< pair<double,string> > work(unsigned steps, bool calculateStatistics, double lowestEnergy);

    void printout(temp_t temperature);

    void printout_service();

    mpf_class e;
    mpf_class e2;
    chrono::milliseconds duration;
    vector<std::unique_ptr<CalculationParameter>> calculationParameters;


    static bool exchange(const Worker &w1, const Worker &w2);
};

#endif //WORKER_H