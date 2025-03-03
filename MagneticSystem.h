#ifndef MAGNETICSYSTEM_H
#define MAGNETICSYSTEM_H

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <optional>
#include <numeric>
#include "stringfunctions.h"
#include "misc.h"

using namespace std;

class MagneticSystem
{
private:
    string _fileContent;
    int _fileVersion;
    double _range;
    optional<Vect> _size;
    optional<hamiltonian_t> _hamiltonian;
    state_t _state; //текущее состояние системы, для которого расчитана матрица энергий
    optional<double> _E; //кэшированная энергия чтобы не пересчитывать постоянно


    void load_csv();
    void load_v1();
    void load_v2();
    void readFileToString(string fileName);
    void buildEnergyTable();


public:
    MagneticSystem(
        string fileName, 
        double range=0, 
        optional<Vect> size = nullopt, 
        optional<hamiltonian_t> hamiltonian = nullopt, 
        state_t state = {});
    ~MagneticSystem();

    void printHeader(Vect field);

    size_t size() const { return _state.size(); }

    double E() { return (_E) ? (*_E) : (_E = accumulate(eMatrix.begin(), eMatrix.end(), 0.0 ) / 2.0).value(); }

    int getFileVersion(){ return _fileVersion; }
    optional<Vect> getPBCSize(){ return this->_size; }
    double getRange() {return this->_range; }

    /**
     * @brief применить конфигурацию, то есть домножить все строки и столбцы матрицы энергии на соответствующие значения
     * 
     * по факту домножает значения спинов магнитной системы на заданную конфигурацию, и потом пересоздает матрицу энергий.
     * Дальнейшая работа с системой действует будто она в состоянии "111111..."
     * 
     * @param s массив значений, на которые домножает значения спинов
     */
    void applyState(const vector<signed char> s = {});

    /**
     * @brief Возвращает версию файла с магнитной системой
     * 
     * @param file путь до файла
     * @return int 
     * -1 - если формат не опознан
     * 0 - если это .csv файл (проверка только по расширению)
     * 1 - если это .mfsys версии 1 (https://github.com/uyras/partsEngine/wiki/Формат-файла-магнитной-системы#версия-1)
     * 1 - если это .mfsys версии 2 (https://github.com/uyras/partsEngine/wiki/Формат-файла-магнитной-системы#версия-2)
     */
    static int fileVersion(std::string file);

    void save(string filename, state_t state = {}) const;

    /**
     * для спина i энергии соседей будут в интервале от neighbours_from[i] до neighbours_to[i] 
     */
    vector< size_t > neighbours_from;
    vector< size_t > neighbours_count;
    vector< double > eMatrix;
    vector< size_t > neighbourNums;
    vector< Part >   parts;
};


#endif //MAGNETICSYSTEM_H