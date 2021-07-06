#include "CorrelationCore.h"

CorrelationCore::CorrelationCore(
    double minRange, 
    double maxRange, 
    unsigned methodVar,
    const std::vector<unsigned> & spins):
_minRange(minRange),
_maxRange(maxRange),
_methodVar(methodVar),
spins(spins),
cp(0,2048),
cp2(0,2048)
{
    minRange2 = _minRange*_minRange;
    maxRange2 = _maxRange*_maxRange;
}

void CorrelationCore::init(const PartArray & sys)
{
    correlationNeighbours.resize(sys.size());
    correlationPairsNum=0;

    for (auto partA: sys.parts){
        for (auto partB: sys.parts){
            if (partA==partB)
                continue;
            
            double space2 = partA->pos.space_2(partB->pos);

            if (space2>=this->minRange2 && space2<=this->maxRange2){
                switch (this->_methodVar)
                {
                    case 1: this->initMethod1(partA,partB); break;
                    case 2: this->initMethod2(partA,partB); break;
                    case 3: this->initMethod3(partA,partB); break;
                }
            }
        }
    }
    this->correlationPairsNum/=2;
}

long CorrelationCore::getCPFull(const PartArray & sys)
{
    long res = this->cpOld;
    this->cpOld = 0;
    for (auto partA:sys.parts){
        if (dbg) {
            fprintf(stderr,"# neigh for %u: ",partA->Id());
            for (auto partB:this->correlationNeighbours[partA->Id()]){
                fprintf(stderr,"%u,",partB->Id());
            }
            fprintf(stderr,"\n");
        }

        this->method(partA);
    }
    std::swap(res,this->cpOld);
    return res/4;
}

void CorrelationCore::method(const Part* partA)
{
    for (auto partB:this->correlationNeighbours[partA->Id()]){
        if (this->_methodVar==1) this->cpOld += 2*fxor(partA,partB);
        if (this->_methodVar==2) this->cpOld += 2*feSign(partA,partB);
        if (this->_methodVar==3) this->cpOld += 2*fScalar(partA,partB);
    }
}

// functions for correlation options
char CorrelationCore::fxor(const Part* partA,const Part* partB)
{
    return (partA->state ^ partB->state)?-1:+1;
}
char CorrelationCore::feSign(const Part* partA, const Part* partB)
{
    return fxor(partA,partB) * 
        this->correlationValues[std::make_pair(partA->Id(),partB->Id())];
}
char CorrelationCore::fScalar(const Part* partA, const Part* partB)
{
    return fxor(partA,partB) * 
        this->correlationValues[std::make_pair(partA->Id(),partB->Id())];
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