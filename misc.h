#ifndef MISC_H
#define MISC_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>

using namespace std;


//------------------структуры данных--------------------//

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
	string lowerEnergyState;
	int64_t time_proc_total;
};

/**
 * @brief Структура для описания температуры при МК-вычислениях
 * 
 * num_base - порядковый номер базовой температуры
 * num_replica - порядковый номер реплики
 */
struct temp_t {
    unsigned num_base;
    unsigned num_replica;
    double t;
    
    std::string to_string() const {
        std::string text(30, '\0');
        std::snprintf(text.data(),30,"T%u.%u=%e",num_base,num_replica,t);
        return text;
    }
};





//------------------полезные функции--------------------//

/**
 * @brief На вход получает две текстовые строки, в байтах которых записаны состояния спинов
 * и делает посимвольный xor для них
 * 
 * @param s1 
 * @param s2 
 * @return std::string 
 */
std::string xorstr(std::string s1,std::string s2);


/**
 * @brief Считывает матрицы энергий из CSV файла
 * 
 * @param filename имя файла
 * @return vector < vector < double > > массив массивов - матрица энергий
 */
vector < vector < double > > readCSV(string filename);


#endif