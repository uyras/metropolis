#include "dos2.h"

template <typename T>
void Dos2<T>::save(string file)
{
    ofstream f(file);
    for (unsigned i=0; i<data.size()-1; i++){
            f<< val(i)<<"\t";
            f<< val(i+1)<<"\t";
            f<< data[i] <<endl;
    }


    f<< val(data.size()-1)<<"\t";
    f<< max <<"\t";
    f<< data[data.size()-1] <<endl;

    f.close();
}

template <typename T>
void Dos2<T>::load(string file)
{
    ifstream f(file);
    data.clear();
    double a=0, b=0, c=0;
    f>>a; f>>b; f>>c;
    data.push_back(c);
    this->min = a;
    while (f >> a >> b >> c){
        data.push_back(c);
    }
    this->max = b;
    this->intervals = data.size();
    f.close();
}

template <typename T>
string Dos2<T>::toString()
{
    stringstream f;
    for (unsigned i=0; i<data.size()-1; i++){
            f<< val(i)<<"\t";
            f<< val(i+1)<<"\t";
            f<< data[i] <<endl;
    }


    f<< val(data.size()-1)<<"\t";
    f<< max <<"\t";
    f<< data[data.size()-1] <<endl;
    return f.str();
}


template class Dos2<int>;
template class Dos2<double>;
template class Dos2<unsigned>;
