#include "CorrelationPointCore.h"

CorrelationPointCore::CorrelationPointCore(const std::vector<double> & X, 
        const std::vector<double> & Y, 
        double distance, double minRange, double maxRange):
X(X),Y(Y),
_distance(distance),
_minRange(minRange),
_maxRange(maxRange),
cp(0,2048),
cp2(0,2048)
{

}

void CorrelationPointCore::init(const PartArray & sys)
{
    correlationPointSpins.resize(X.size());
    correlationNeighbours.resize(sys.size());
    correlationPairsNum=0;
    spinsInPoint=0;
    double space2;

    //first find the spins around each point
    double dist2 = this->_distance * this->_distance;
    for (int i=0; i<X.size(); ++i){
        Vect point = Vect(X[i],Y[i],0);
        for (auto part: sys.parts){
            space2 = point.space_2(part->pos);
            if (space2<=dist2){
                this->correlationPointSpins[i].push_front(part);
                spinsInPoint += 1;
            }
        }
        if (std::distance(this->correlationPointSpins[i].begin(),
                          this->correlationPointSpins[i].end())==0){
            cerr << "# Corellation point "<<i<<" has no spins around. Check your config."<<endl;
        }
    }
    spinsInPoint /= this->pointCount();

    const double minr2 = this->_minRange * this->_minRange;
    const double maxr2 = this->_maxRange * this->_maxRange;
    double eTemp;
    //find neighbours in all corellation points
    for (auto cps: this->correlationPointSpins){
        for (auto partA: cps){
            for (auto partB: cps){
                if (partA==partB)
                    continue;
                
                space2 = partA->pos.space_2(partB->pos);
                if (space2>=minr2 && space2<=maxr2){
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
            }
        }
    }
    this->correlationPairsNum/=2;
}

long CorrelationPointCore::getCPFull(const PartArray & sys){
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

void CorrelationPointCore::method(const Part* partA)
{
    for (auto partB:this->correlationNeighbours[partA->Id()]){
        this->cpOld += 2*((partA->state ^ partB->state)?-1:+1) * this->correlationValues[std::make_pair(partA->Id(),partB->Id())];
    }
}