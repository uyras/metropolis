#include "MagnetisationLengthCore.h"

MagnetisationLengthCore::MagnetisationLengthCore(const std::string & parameterId, 
    PartArray * prototype,
    const std::vector<uint64_t> & spins):
CalculationParameter(parameterId,prototype),
spins(spins),
mv(0,1024*8),
mv2(0,2048*8),
mOld(0,0,0)
{
    if (this->spins.size() == 0) {
        this->spins.resize(this->prototype->size(),0);
        for (uint64_t i=0; i<this->prototype->size(); ++i) this->spins[i]=i;
    }

    this->prototypeInit(prototype);
}

bool MagnetisationLengthCore::check(unsigned N) const
{
    return true;
}

void MagnetisationLengthCore::printHeader(unsigned num) const
{
    // get saturation magnetisation
    /*double saturation=0;
    for (auto s : this->spins){
        saturation += fabs(this->prototype->parts[s]->m.scalar(this->vector));
    }
    saturation /= this->spins.size();*/


    printf("##### calculation param #%d #####\n",num);

    printf("# type: magnetisation length\n");

    printf("# id: %s\n",this->parameterId().c_str());

    Vect tmp;
    printf("# initial value %f\n",
        this->getFullTotal(tmp)/this->spins.size());
    //printf("# saturation magnetisation / N: %f\n", saturation);
    printf("# spins: ");
    if (this->spins.size()==this->prototype->size()){
        printf("All\n");
    } else {
        printf("%d",this->spins[0]);
        for (int ss=1; ss<this->spins.size(); ++ss){
            printf(",%d",this->spins[ss]);
        }
        printf("\n");
    }

    printf("#\n");
    return;
}

bool MagnetisationLengthCore::init(PartArray * sys)
{  
    this->sys = sys;

    getFullTotal(this->mOld);

    return true;
}

void MagnetisationLengthCore::iterate(unsigned id){
    this->mOld += this->method(id)*2;
    if (_debug){
        Vect tmp;
        this->getFullTotal(tmp);
        if ((tmp - this->mOld).length()>0.001) 
            cerr<<"# (dbg MagnetisationLengthCore#"<<this->parameterId()<<") total vecto is different: iterative="<<this->mOld<<", full="<<tmp<<endl;
    }
}

void MagnetisationLengthCore::incrementTotal(){
    double addVal = this->mOld.length() / this->spins.size();
    this->mv += addVal;
    this->mv2 += addVal*addVal;
}

double MagnetisationLengthCore::getFullTotal(Vect & val) const
{
    val.setXYZ(0,0,0);
    for (auto spinId: spins){
        val += this->method(spinId);
    }
    return val.length();
}

Vect MagnetisationLengthCore::method(unsigned spinId) const
{
        return sys->parts[spinId]->m;
}