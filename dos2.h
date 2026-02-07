#ifndef DOS2_H
#define DOS2_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <unordered_map>

using namespace std;


template <typename T>
class Dos2
{
public:
    Dos2():min(0),max(0),intervals(0){}

    Dos2(T min, T max, unsigned intervals)
    {
        this->resize(min,max,intervals);
    }

    inline void resize(T min, T max, unsigned intervals){
        this->min = min;
        this->max= max;
        this->intervals = intervals;
        data.resize(intervals);
    }

    virtual void save(string file);
    virtual void load(string file);
    string toString();

    inline T& at(unsigned i){return data[i];}
    inline const T& at(unsigned i) const {return data[i];}
    inline T& operator[](T a){return data[num(a)];}
    inline const T& operator[](T a) const {return data[num(a)];}

    inline bool const operator==(const Dos2<T>& rhs) const {
        if (this->min != rhs.min) return false;
        if (this->max != rhs.max) return false;
        if (this->intervals != rhs.intervals) return false;
        if (this->data != rhs.data) return false;

        return true;
    }
    inline bool const operator!=(const Dos2<T>& rhs) const { return !(this->operator==(rhs)); }

    inline unsigned Intervals() const {return intervals;}
    inline T Min() const {return this->min;}
    inline T Max() const {return this->max;}


    /**
     * @brief val Возращает значение, с которого начинается i интервал
     * @param i Интервал, для которого надо получить соответствующее значение типа T
     * @return
     */
    inline T val(unsigned i) const{ return min + (max-min)*(i/(double)(intervals)); }

    /**
     * @brief num Возращает номер интервала для соответствующего значения
     * @param a
     */
    inline unsigned num(T a) const{
        if (a<min)
            return 0;
        if(a>=max)
            return intervals-1;
        unsigned r= (unsigned)floor(double(a-min)/(max-min)*(intervals));
        if (r>intervals-1)
            return unsigned(intervals-1);
        else
            return r;
    }

    /**
     * @brief clear Reset all values of DOS to 0
     */
    inline void clear(){
        for (auto i: data){
            i=0.;
        }
    }

protected:
    T min, max;
    unsigned intervals;
    vector< T > data;
};

#endif // DOS2_H
