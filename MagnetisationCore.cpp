#include "MagnetisationCore.h"

MagnetisationCore::MagnetisationCore(
        const config_section_t &sect,
        const ConfigManager *conf):
CalculationParameter(sect,conf),
mv(0,1024*8),
mv2(0,2048*8),
currentState(sys->N(),1),
magnetisationValues(sys->N(),0),
_sumModule(false)
{
    // получаем параметры из ini-конфига
    if (sect.data.contains("spins"))
        this->spins = sect.data["spins"].get_list<inicpp::unsigned_ini_t>();

    if (this->spins.size() == 0) {
        this->spins.resize(this->sys->N(),0);
        for (size_t i=0; i<this->sys->N(); ++i) 
            this->spins[i]=i;
    }


    if (sect.data.contains("module"))
        this->_sumModule = sect.data["module"].get<inicpp::boolean_ini_t>();



    if (!sect.data.contains("axis"))
        throw(std::string("Parameter 'axis' is required in section "+sect.data.get_name()));
    std::string tmp_axis = sect.data["axis"].get<inicpp::string_ini_t>();
    this->vector = makeUnit(strToVect(tmp_axis));

    for (auto spinId: spins){
        magnetisationValues[spinId] = scalar(sys->parts[spinId].m, vector);        
    }

}

bool MagnetisationCore::check(unsigned N) const
{
    return true;
}

void MagnetisationCore::printHeader() const
{
    // get saturation magnetisation
    double saturation=0;
    for (auto s : this->spins){
        saturation += fabs(scalar(this->sys->parts[s].m, this->vector));
    }
    saturation /= this->spins.size();

    if (this->_sumModule)
        printf("# type: %s module\n", MagnetisationCore::name().c_str());
    else
        printf("# type: %s\n", MagnetisationCore::name().c_str());

    printf("# id: %s\n",this->parameterId().c_str());

    printf("# magnetisation vector: (%g|%g|%g); initial value %g\n",
        this->vector.x,
        this->vector.y,
        this->vector.z,
        this->getFullTotal(state_t(sys->N(),1)) / double(this->spins.size()) );
    printf("# saturation magnetisation / N: %g\n", saturation);
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

void MagnetisationCore::init(const state_t &state)
{
    this->currentState = state;
    this->mOld = getFullTotal(state);
    return;
}

void MagnetisationCore::iterate(size_t id){
    this->currentState[id] = -this->currentState[id];
    this->mOld += 2*this->method(id,this->currentState[id]);
    if (this->isDebug()){
        double res = this->getFullTotal(this->currentState);
        if (fabs(res-this->mOld)>0.01) 
            cerr<<"# (dbg MagnetisationCore#"<<this->parameterId()<<") total value is different: iterative="<<this->mOld<<", full="<<res<<endl;
    }
}

void MagnetisationCore::incrementTotal(){
    double addVal = double(this->mOld) / this->spins.size();
    if (this->_sumModule)
        this->mv += fabs(addVal);
    else
        this->mv += addVal;
    this->mv2 += addVal*addVal;
}

double MagnetisationCore::getFullTotal(const state_t &_state) const
{
    double res = 0;
    for (auto spinId: spins){
        res += this->method(spinId, _state[spinId]);
    }
    return res;
}

double MagnetisationCore::method(unsigned spinId, signed char spinState) const
{
        return magnetisationValues[spinId] * spinState;
}