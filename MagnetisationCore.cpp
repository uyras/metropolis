#include "MagnetisationCore.h"

MagnetisationCore::MagnetisationCore(const Vect & vector, const std::vector<unsigned> & spins):
vector(vector),
spins(spins),
mv(0,2048),
mv2(0,2048),
mOld(0),
_sys(nullptr)
{

}

void MagnetisationCore::init(PartArray & sys)
{
    this->_sys = &sys;
    magnetisationValues.resize(sys.size(),0);

    for (auto spinId: spins){
        Part* part = sys[spinId];
        magnetisationValues[spinId] = part->m.scalar(vector);
        if (!part->state)
            magnetisationValues[spinId] *= -1;
        
    }
}

long MagnetisationCore::getMOFull()
{
    double res = this->mOld;
    this->mOld = 0;
    for (auto spinId: spins){
        this->method(spinId);
    }
    std::swap(res,this->mOld);
    return res/2;
}

void MagnetisationCore::method(unsigned spinId)
{
    this->mOld += 2*magnetisationValues[spinId]*(this->_sys->getById(spinId)->state)?-1:+1;
}