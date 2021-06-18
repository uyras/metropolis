#include "CorrelationCore.h"

CorrelationCore::CorrelationCore(double minRange, double maxRange, unsigned methodVar):
_minRange(minRange),_maxRange(maxRange),_methodVar(methodVar),
cp(0,2048),
cp2(0,2048)
{

}

void CorrelationCore::init(const PartArray & sys)
{
    correlationNeighbours.resize(sys.size());
    correlationPairsNum=0;

    double space2;
    const double minr = this->_minRange * this->_minRange;
    const double maxr = this->_maxRange * this->_maxRange;
    double eTemp;
    for (auto partA: sys.parts){
        for (auto partB: sys.parts){
            if (partA==partB)
                continue;
            
            space2 = partA->pos.space_2(partB->pos);
            if (space2>=minr && space2<=maxr){
                this->correlationNeighbours[partA->Id()].push_front(partB);
                ++correlationPairsNum;

                if (this->_methodVar==2 || this->_methodVar==3){
                    if (this->_methodVar==2){
                        //В матрицу надо помещать энергии только в неперевернутых состояниях
                        eTemp = hamiltonianDipolar(partA,partB)*-1; 
                        if (partA->state!=partB->state)
                            eTemp*=-1.;
                    }

                    if (this->_methodVar==3){
                        eTemp = partA->m.scalar(partB->m); 
                        if (partA->state!=partB->state)
                            eTemp*=-1.;
                    }


                    if (eTemp>0)
                        this->correlationValues[std::make_pair(partA->Id(),partB->Id())] = 1;
                    else
                        this->correlationValues[std::make_pair(partA->Id(),partB->Id())] = -1;
                }
            }
        }
    }
    this->correlationPairsNum/=2;
}

long CorrelationCore::getCPFull(const PartArray & sys)
{
    long res = 0;
    for (auto partA:sys.parts){
        if (dbg) fprintf(stderr,"# neigh for %u: ",partA->Id());
        for (auto partB:this->correlationNeighbours[partA->Id()]){
            res += this->method(partA,partB);
            if (dbg) fprintf(stderr,"%u,",partB->Id());
        }
        if (dbg) fprintf(stderr,"\n");
    }
    return res/2;
}

char CorrelationCore::method(const Part* partA,const Part* partB)
{
        if (this->_methodVar==1) return fxor(partA,partB);
        if (this->_methodVar==2) return feSign(partA,partB);
        if (this->_methodVar==3) return fScalar(partA,partB);
        return 0;
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