#include "CorrelationPointCore.h"

CorrelationPointCore::CorrelationPointCore(
        const std::string & parameterId,
        PartArray * prototype,
        const std::vector<double> & X, 
        const std::vector<double> & Y, 
        double distance, double minRange, double maxRange):
CalculationParameter(parameterId,prototype),
X(X),Y(Y),
_distance(distance),
_minRange(minRange),
_maxRange(maxRange),
cp(0,2048),
cp2(0,2048)
{
    this->prototypeInit(prototype);
}

bool CorrelationPointCore::check(unsigned N) const
{
    if (this->X.size() != this->Y.size()){
                cerr<<"Parameters -X, and -Y should \
                have the same count of values."<<endl;
                return false;
        }

    /*if (this->correlationPairsNum==0){
            cerr<<"Check the correlation point ranges."<<endl;
            cerr<<"Could not find any spin pairs within this distance!"<<endl;
            return false;
        }*/

    return true;
}

void CorrelationPointCore::printHeader(unsigned num) const
{

    printf("##### calculation param #%d #####\n",num);
    printf("# type: correlationpoint\n");
    printf("# id: %s\n",this->parameterId().c_str());

    printf("# distance aroun point: %.2f\n",this->_distance);
    printf("# minimal interaction distance: %.2f\n",this->_minRange);
    printf("# maximal interaction distance: %.2f\n",this->_maxRange);

    printf("# points: %d; (%.2f avg spins per point, %.2f avg. neighbours)\n",
            this->X.size(),
            this->spinsInPoint,
            this->correlationPairsNum/double(this->prototype->size()));
    printf("#    coordinates: format is <num:(x,y):spins>, <...>, ...\n");
    printf("# 0:(%f,%f):%d",
        this->X[0],
        this->Y[0],
        std::distance(this->correlationPointSpins[0].begin(),
                        this->correlationPointSpins[0].end()));
    for (int i=1; i<this->X.size(); ++i)
        printf(", %d:(%f,%f):%d",i,
            this->X[i],
            this->Y[i],
            std::distance(this->correlationPointSpins[i].begin(),
                        this->correlationPointSpins[i].end()));
    printf("\n");

    printf("#\n");
    return;
}

bool CorrelationPointCore::init(PartArray * sys)
{

    // clear variables
    correlationPointSpins.clear();
    correlationNeighbours.clear();
    correlationValues.clear();
    correlationPairsNum=0;
    spinsInPoint=0;

    
    this->sys = sys;

    correlationPointSpins.resize(X.size());
    correlationNeighbours.resize(sys->size());
    double space2;

    //first find the spins around each point
    double dist2 = this->_distance * this->_distance;
    for (int i=0; i<X.size(); ++i){
        Vect point = Vect(X[i],Y[i],0);
        for (auto part: sys->parts){
            space2 = point.space_2(part->pos);
            if (space2<=dist2){
                this->correlationPointSpins[i].push_front(part);
                spinsInPoint += 1;
            }
        }
        if (std::distance(this->correlationPointSpins[i].begin(),
                          this->correlationPointSpins[i].end())==0){
            throw(std::invalid_argument("# Corellation point "+std::to_string(i)+" has no spins around. Check your config."));
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

    if (_debug) {
        for (auto partA : sys->parts){
            fprintf(stderr,"# neigh for %u: ",partA->Id());
            for (auto partB : this->correlationNeighbours[partA->Id()]){
                fprintf(stderr,"%u,",partB->Id());
            }
            fprintf(stderr,"\n");
        }
    }

    this->cpOld = this->getFullTotal(this->sys);

    return true;
}

void CorrelationPointCore::iterate(unsigned id){
    for (auto partB:this->correlationNeighbours[id]){
        this->cpOld += 2*this->method(this->sys->getById(id),partB);
    }


    if (_debug){
        long res = this->getFullTotal(sys);
        if (res!=this->cpOld) 
            cerr<<"# (dbg CorrelationPointCore) total value is different: iterative="<<this->cpOld<<", full="<<res<<endl;
    }
}

void CorrelationPointCore::incrementTotal(){
    this->cp += double(this->cpOld)/this->pointCount();
    this->cp2 += double(this->cpOld * this->cpOld)/(this->pointCount() * this->pointCount());
}

long CorrelationPointCore::getFullTotal(const PartArray * _sys) const
{
    long res = 0;
    for (auto partA : _sys->parts){
            for (auto partB : this->correlationNeighbours[partA->Id()]){
                res += this->method(partA,partB);
            }
    }
    return res/2;
}

short CorrelationPointCore::method(const Part* partA, const Part* partB) const
{
    return ((partA->state ^ partB->state)?-1:+1) * this->correlationValues.at(std::make_pair(partA->Id(),partB->Id()));
}