#include "CorrelationCore.h"

CorrelationCore::CorrelationCore(
    const config_section_t &sect,
    const ConfigManager *conf):
CalculationParameter(sect, conf),
cp(0,1024*8),
cp2(0,2048*8),
currentState(sys->N(),1)
{
    static const std::map<std::string, unsigned> methods_list = {{"xor",0},{"energy",1},{"scalar",2}};

    if (!sect.data.contains("method"))
        throw(std::string("Parameter 'method' is required in section "+sect.data.get_name()));
    const string method_name = sect.data["method"].get<inicpp::string_ini_t>();
    try {
        _methodVar = methods_list.at(method_name);
    } catch (std::out_of_range e){
        throw string("method '"+method_name+"' is unknown in section "+sect.data.get_name());
    }

    if (!sect.data.contains("maxrange")) throw(std::string("Parameter 'maxrange' is required in section "+sect.data.get_name()));
    _maxRange =  sect.data["maxrange"].get<inicpp::float_ini_t>();

    _minRange = (sect.data.contains("minrange")) ? sect.data["minrange"].get<inicpp::float_ini_t>() : 0;

    
    if (sect.data.contains("spins"))
        this->spins = sect.data["spins"].get_list<inicpp::unsigned_ini_t>();
    if (this->spins.size() == 0) {
        this->spins.resize(this->sys->N(),0);
        for (uint64_t i=0; i<this->sys->N(); ++i) this->spins[i]=i;
    }


    // далее определяем соседей
    correlationNeighbours.clear();
    correlationValues.clear();
    correlationPairsNum=0;
    
    correlationNeighbours.resize(sys->N());

    
    double minRange2 = _minRange*_minRange; // для ускорения расстояния считаются и сравниваются в квадрате
    double maxRange2 = _maxRange*_maxRange;

    for (size_t idA = 0; idA < sys->N(); idA++){
        for (size_t idB = 0; idB < sys->N(); idB++){
            if (idA==idB)
                continue;

            const auto& partA = sys->parts[idA];
            const auto& partB = sys->parts[idB];
            
            double space2 = distance_2(partA.p,partB.p);

            if (space2>=minRange2 && space2<=maxRange2){
                if (!areNeighbours(idA,idB)){
                    switch (this->_methodVar)
                    {
                        case 0: this->initMethod1(idA,idB); break;
                        case 1: this->initMethod2(idA,idB); break;
                        case 2: this->initMethod3(idA,idB); break;
                    }
                }
                if (!areNeighbours(idB,idA)){
                    switch (this->_methodVar)
                    {
                        case 0: this->initMethod1(idB,idA); break;
                        case 1: this->initMethod2(idB,idA); break;
                        case 2: this->initMethod3(idB,idA); break;
                    }
                }
            }
        }
    }

    this->correlationPairsNum/=2;

    if (isDebug()){
        for (size_t idA = 0; idA < sys->N(); idA++){
            if (!this->correlationNeighbours[idA].empty()){
                fprintf(stderr,"# neigh for %ld: ", idA);
                for (auto idB:this->correlationNeighbours[idA]){
                    fprintf(stderr,"%ld,", idB);
                }
                fprintf(stderr,"\n");
            }
        }
    }

}

bool CorrelationCore::check(unsigned N) const
{
    /*if (cctemp.correlationPairsNum==0){
            cerr<<"Check the correlation range #"<<i<<"."<<endl;
            cerr<<"Could not find any spin pairs within this distance!"<<endl;
            return 1;
        }*/
    return true;
}

void CorrelationCore::printHeader() const
{
    printf("# type: %s\n",CorrelationCore::name().c_str());
    printf("# id: %s\n",this->parameterId().c_str());
    
    int totalInteractingSpins = 0;
    for (int i=0; i<this->sys->N(); ++i){
        if (!this->correlationNeighbours[i].empty())
            ++totalInteractingSpins;
    }

    printf("# minimal interaction distance: %.2f\n",this->_minRange);
    printf("# maximal interaction distance: %.2f\n",this->_maxRange);
    printf("# average neighbours: %.2f\n",this->correlationPairsNum/double(totalInteractingSpins)*2);
    if (this->_methodVar==0)
        printf("# method: XOR\n");
    if (this->_methodVar==1)
        printf("# method: energy\n");
    if (this->_methodVar==2)
        printf("# method: scalar\n");
    printf("# initial value: %g\n",this->getFullTotal(state_t(sys->N(),1)) / double(sys->N()));
    
    printf("# spins: ");
    if (this->spins.size() == this->sys->N()){
        printf("All\n");
    } else {
        printf("%lu",this->spins[0]);
        for (int ss=1; ss<this->spins.size(); ++ss){
            printf(",%lu",this->spins[ss]);
        }
        printf("\n");
    }
    
}

void CorrelationCore::init(const state_t &state)
{
    this->currentState = state;
    this->cpOld = this->getFullTotal(state);
    
    return;
}

void CorrelationCore::iterate(size_t id){
    this->currentState[id] = - this->currentState[id];
    for (auto idB : this->correlationNeighbours[id]){
        this->cpOld += 2*this->method(id,idB);
    }


    if (this->isDebug()){
        long res = this->getFullTotal(this->currentState);
        if (res!=this->cpOld) 
            cerr<<"# (dbg CorrelationCore#"<<this->parameterId()<<") total value is different: iterative="<<this->cpOld<<", full="<<res<<endl;
    }
}

void CorrelationCore::incrementTotal(){
    this->cp += double(this->cpOld)/this->correlationPairsNum;
    this->cp2 += double(this->cpOld * this->cpOld)/(this->correlationPairsNum*this->correlationPairsNum);
}

long CorrelationCore::getFullTotal(const state_t &_state) const
{
    long res=0;
    for (size_t idA = 0; idA < sys->N(); idA++){
        for (auto idB : this->correlationNeighbours[idA]){
            // чтобы отключить занчение текущего состояния системы, домножаем на fxor
            res += this->method(idA, idB) * _state[idA] * _state[idB] * fxor(idA,idB);
        }
    }
    return res/2;
}

char CorrelationCore::method(const size_t idA, const size_t idB) const
{
    if (this->_methodVar==0) return fxor(idA,idB);
    if (this->_methodVar==1) return feSign(idA,idB);
    if (this->_methodVar==2) return fScalar(idA,idB);

    throw(string("Invalid _methodVar values found in CorrelationCore::method"));
    return 0;
}

// functions for correlation options
char CorrelationCore::fxor(const size_t idA, const size_t idB) const
{
    return currentState[idA] * currentState[idB]; // +1 если равны, -1 если не равны
}
char CorrelationCore::feSign(const size_t idA, const size_t idB) const
{
    return fxor(idA,idB) * 
        this->correlationValues.at(std::make_pair(idA,idB));
}
char CorrelationCore::fScalar(const size_t idA, const size_t idB) const
{
    return fxor(idA,idB) * 
        this->correlationValues.at(std::make_pair(idA,idB));
}

void CorrelationCore::initMethod1( size_t idA, size_t idB ){
    this->correlationNeighbours[idA].push_front(idB);
    ++correlationPairsNum;
}

void CorrelationCore::initMethod2( size_t idA, size_t idB ){
    double eTemp;
    
    this->correlationNeighbours[idA].push_front(idB);
    ++correlationPairsNum;

    //В матрицу надо помещать энергии только в неперевернутых состояниях
    eTemp = hamiltonian_dipolar(this->sys->parts[idA],this->sys->parts[idB]); 

    this->correlationValues[std::make_pair(idA,idB)] = (eTemp>0) ? 1 : -1;
}

void CorrelationCore::initMethod3( size_t idA, size_t idB ){
    double eTemp;
    
    this->correlationNeighbours[idA].push_front(idB);
    ++correlationPairsNum;

    eTemp = scalar(this->sys->parts[idA].m, this->sys->parts[idB].m); 

    this->correlationValues[std::make_pair(idA,idB)] = (eTemp>0) ? 1 : -1;
}