#include "MagnetisationLengthCore.h"

MagnetisationLengthCore::MagnetisationLengthCore(
        const config_section_t &sect,
        const ConfigManager *conf):
CalculationParameter(sect,conf),
mv(0,1024*8),
mv2(0,2048*8),
currentState(sys->N(),1),
enabledSpins(sys->N(),0)
{
    // получаем параметры из ini-конфига
    if (sect.data.contains("spins"))
        this->spins = sect.data["spins"].get_list<inicpp::unsigned_ini_t>();

    if (this->spins.size() == 0) {
        this->spins.resize(this->sys->N(),0);
        for (size_t i=0; i<this->sys->N(); ++i) 
            this->spins[i]=i;
    }

    for (auto spinId: spins){
        enabledSpins[spinId] = 1;
    }

}

bool MagnetisationLengthCore::check(unsigned N) const
{
    return true;
}

void MagnetisationLengthCore::printHeader() const
{
    printf("# type: %s\n", MagnetisationLengthCore::name().c_str());

    printf("# id: %s\n",this->parameterId().c_str());

    printf("# initial value %g\n",
        this->getFullTotal(state_t(sys->N(),1)).length() / double(this->spins.size()) );
    printf("# spins: ");
    if (this->spins.size()==this->sys->N()){
        printf("All\n");
    } else {
        printf("%lu",this->spins[0]);
        for (int ss=1; ss<this->spins.size(); ++ss){
            printf(",%lu",this->spins[ss]);
        }
        printf("\n");
    }

    printf("#\n");
    return;
}

void MagnetisationLengthCore::init(const state_t &state)
{
    this->currentState = state;
    this->mOld = getFullTotal(state);
    return;
}

void MagnetisationLengthCore::iterate(size_t id){
    this->currentState[id] = -this->currentState[id];
    this->mOld += this->method(id,this->currentState[id]) * 2;
    if (this->isDebug()){
        Vect res = this->getFullTotal(this->currentState);
        if ((res - this->mOld).length()>0.01) 
            cerr<<"# (dbg MagnetisationLengthCore#"<<this->parameterId()<<") total value is different: iterative="<<this->mOld<<", full="<<res<<endl;
    }
}

void MagnetisationLengthCore::incrementTotal(){
    double addVal = this->mOld.length()  / this->spins.size();
    
    this->mv += addVal;
    this->mv2 += addVal*addVal;
}

Vect MagnetisationLengthCore::getFullTotal(const state_t &_state) const
{
    Vect res = {0};
    for (auto spinId: spins){
        res += this->method(spinId, _state[spinId]);
    }
    return res;
}

Vect MagnetisationLengthCore::method(unsigned spinId, signed char spinState) const
{
        return sys->parts[spinId].m * spinState;
}