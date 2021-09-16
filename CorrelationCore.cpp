#include "CorrelationCore.h"

CorrelationCore::CorrelationCore(
    const std::string & parameterId,
    PartArray * prototype,
    double minRange, 
    double maxRange, 
    unsigned methodVar,
    const std::vector<uint64_t> & spins):
CalculationParameter(parameterId, prototype),
_minRange(minRange),
_maxRange(maxRange),
_methodVar(methodVar),
spins(spins),
cp(0,2048),
cp2(0,2048)
{
    minRange2 = _minRange*_minRange;
    maxRange2 = _maxRange*_maxRange;
    if (this->spins.size() == 0) {
        this->spins.resize(this->prototype->size(),0);
        for (uint64_t i=0; i<this->prototype->size(); ++i) this->spins[i]=i;
    }

    this->prototypeInit(prototype);
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

void CorrelationCore::printHeader(unsigned num) const
{
    printf("##### calculation param #%d #####\n",num);
    printf("# type: correlation\n");
    printf("# id: %s\n",this->parameterId().c_str());
    
    int totalSpins = this->spins.size();
    if (totalSpins==0) totalSpins=this->prototype->size();

    printf("# minimal interaction distance: %.2f\n",this->_minRange);
    printf("# maximal interaction distance: %.2f\n",this->_maxRange);
    printf("# average neighbours: %.2f\n",this->correlationPairsNum/double(totalSpins));
    if (this->_methodVar==0)
        printf("# method: XOR\n");
    if (this->_methodVar==1)
        printf("# method: energy\n");
    if (this->_methodVar==2)
        printf("# method: scalar\n");
    
    printf("# spins: ");
    if (this->spins.size() == this->prototype->size()){
        printf("All\n");
    } else {
        printf("%d",this->spins[0]);
        for (int ss=1; ss<this->spins.size(); ++ss){
            printf(",%d",this->spins[ss]);
        }
        printf("\n");
    }

    printf("#\n");
    
}

bool CorrelationCore::init(PartArray * sys)
{

    //reset
    correlationNeighbours.clear();
    correlationValues.clear();
    correlationPairsNum=0;
    
    this->sys = sys;
    correlationNeighbours.resize(sys->size());

    for (auto partA: sys->parts){
        for (auto partB: sys->parts){
            if (partA==partB)
                continue;
            
            double space2 = partA->pos.space_2(partB->pos);

            if (space2>=this->minRange2 && space2<=this->maxRange2){
                switch (this->_methodVar)
                {
                    case 0: this->initMethod1(partA,partB); break;
                    case 1: this->initMethod2(partA,partB); break;
                    case 2: this->initMethod3(partA,partB); break;
                }
            }
        }
    }
    this->correlationPairsNum/=2;

    if (_debug){
        for (auto partA : this->sys->parts){
            fprintf(stderr,"# neigh for %u: ",partA->Id());
            for (auto partB:this->correlationNeighbours[partA->Id()]){
                fprintf(stderr,"%u,",partB->Id());
            }
            fprintf(stderr,"\n");
        }
    }

    this->cpOld = this->getFullTotal(this->sys);

    return true;
}

void CorrelationCore::iterate(unsigned id){
    for (auto partB:this->correlationNeighbours[id]){
        this->cpOld += 2*this->method(this->sys->getById(id),partB);
    }


    if (_debug){
        long res = this->getFullTotal(this->sys);
        if (res!=this->cpOld) 
            cerr<<"# (dbg CorrelationCore) total value is different: iterative="<<this->cpOld<<", full="<<res<<endl;
    }
}

void CorrelationCore::incrementTotal(){
    this->cp += double(this->cpOld)/this->correlationPairsNum;
    this->cp2 += double(this->cpOld * this->cpOld)/(this->correlationPairsNum*this->correlationPairsNum);
}

long CorrelationCore::getFullTotal(const PartArray * _sys) const
{

    long res=0;
    for (auto partA : _sys->parts){
        for (auto partB:this->correlationNeighbours[partA->Id()]){
            res += this->method(partA, partB);
        }

    }
    return res/2;
}

char CorrelationCore::method(const Part* partA, const Part* partB) const
{
    if (this->_methodVar==0) return fxor(partA,partB);
    if (this->_methodVar==1) return feSign(partA,partB);
    if (this->_methodVar==2) return fScalar(partA,partB);

    throw(invalid_argument("Invalid _methodVar values found in CorrelationCore::method"));
    return 0;
}

// functions for correlation options
char CorrelationCore::fxor(const Part* partA,const Part* partB) const
{
    return (partA->state ^ partB->state)?-1:+1;
}
char CorrelationCore::feSign(const Part* partA, const Part* partB) const
{
    return fxor(partA,partB) * 
        this->correlationValues.at(std::make_pair(partA->Id(),partB->Id()));
}
char CorrelationCore::fScalar(const Part* partA, const Part* partB) const
{
    return fxor(partA,partB) * 
        this->correlationValues.at(std::make_pair(partA->Id(),partB->Id()));
}

void CorrelationCore::initMethod1(Part* partA, Part* partB){
    this->correlationNeighbours[partA->Id()].push_front(partB);
    ++correlationPairsNum;
}

void CorrelationCore::initMethod2( Part* partA, Part* partB){
    double eTemp;
    
    this->correlationNeighbours[partA->Id()].push_front(partB);
    ++correlationPairsNum;

    //В матрицу надо помещать энергии только в неперевернутых состояниях
    eTemp = hamiltonianDipolar(partA,partB)*-1; 
    if (partA->state!=partB->state)
        eTemp*=-1.;

    if (eTemp>0)
        this->correlationValues[std::make_pair(partA->Id(),partB->Id())] = 1;
    else
        this->correlationValues[std::make_pair(partA->Id(),partB->Id())] = -1;
}

void CorrelationCore::initMethod3( Part* partA, Part* partB){
    double eTemp;
    
    this->correlationNeighbours[partA->Id()].push_front(partB);
    ++correlationPairsNum;

    //В матрицу надо помещать энергии только в неперевернутых состояниях
    eTemp = partA->m.scalar(partB->m); 
    if (partA->state!=partB->state)
        eTemp*=-1.;


    if (eTemp>0)
        this->correlationValues[std::make_pair(partA->Id(),partB->Id())] = 1;
    else
        this->correlationValues[std::make_pair(partA->Id(),partB->Id())] = -1;
}