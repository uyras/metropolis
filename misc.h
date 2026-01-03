#ifndef MISC_H
#define MISC_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <inicpp/inicpp.h>

using namespace std;

//------------------структуры данных--------------------//

typedef vector<signed char> state_t;

/**
 * @brief структура-обертка для секции, так как в inicpp 
 * нельзя напрямую поменять имя секции. Оно задается только при создании новой секции, при копировании копируется старое
 */
struct config_section_t {
    string name;
    inicpp::section data;
};

/**
 * @brief Структура для сохранения статистики работы Монте-Карло.
 *  Эту структуру удобно возвращать из функции
 * 
 * initEnergy - начальная энергия вычислений
 * lowerEnergy - нижайшая энергия найденная при вычислениях
 * 
 */
struct monteCarloStatistics {
	double initEnergy; 
	double lowerEnergy;
	bool foundLowerEnergy;
	int temperatureOfLowerEnergy;
	state_t lowerEnergyState;
	int64_t time_proc_total;
};

/**
 * @brief Структура для описания температуры при МК-вычислениях
 * 
 * unsigned num_base - порядковый номер базовой температуры
 * unsigned num_replica - порядковый номер реплики
 * double t - значение температуры
 */
struct temp_t {
    size_t num_base;
    size_t num_replica;
    double t;
    
    std::string to_string() const {
        std::string text(30, '\0');
        std::snprintf(text.data(),30,"T%lu.%lu=%e",num_base,num_replica,t);
        return text;
    }
};

struct Vect{
    double x;
    double y;
    double z;
};

struct Part{
    Vect p;
    Vect m;
};


typedef double (*hamiltonian_t)(const Part&, const Part&);


//------------------полезные функции--------------------//

/**
 * @brief На вход получает две конфигурации, в которых записаны состояния спинов
 * и делает посимвольный xor для них
 * 
 * @param s1 
 * @param s2 
 * @return state_t
 */
state_t xorstate(const state_t &s1, const state_t &s2);

std::string stateToString(const state_t &state);


/**
 * @brief Считывает матрицы энергий из CSV файла
 * 
 * @param filename имя файла
 * @return vector < vector < double > > массив массивов - матрица энергий
 */
vector < vector < double > > readCSV(string filename);

Vect strToVect(std::string val);

/**
 * @brief Транслирует координаты периодических граничных условий и возвращает модифицированные координаты
 * вектора b
 * 
 * @param a основной вектор относительно которого делать преобразование координат
 * @param b вспомогательный вектор который двигаем в пространстве
 * @param size размер пространственного куба. Если одна из координат задана как 0, то не учитываем ПГУ в этом направлении
 * @return Vect Модифицированные координаты вектора b, транслированные с учетом ПГУ
 */
Vect translatePBC(const Vect &a, const Vect &b, const Vect size);

double distance(const Vect &a, const Vect &b);
double distance_2(const Vect &a, const Vect &b);

double scalar(const Vect &a, const Vect &b);

#endif