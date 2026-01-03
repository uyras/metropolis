#include "CorrelationPointCore.h"

// @todo Сейчас в correlationNeighbours хранятся все соседи для каждого спина. 
// Но спин может принадлежать сразу нескольким кореляционным точкам, и тогда будет иметь разных соседей.

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
cp(0,1024*8),
cp2(0,2048*8),
cp4(0,3076*8),
_histogramEnabled(false),
_histogramFilename("")
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
    if (this->_histogramEnabled)
        printf("# histogram: enabled (slow)\n");
    else
        printf("# histogram: disabled\n");

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

    if (this->_histogramEnabled){
        correlationNeighboursByPoint.resize(X.size());
        for (auto &cps: this->correlationNeighboursByPoint){
            cps.resize(sys->size());
        }
    }

    //first find the spins around each point
    double dist2 = this->_distance * this->_distance;
    unsigned maxSpinsInPoint=0;
    for (int i=0; i<X.size(); ++i){
        Vect point = Vect(X[i],Y[i],0);
        unsigned spinsInCurrentPoint = 0;
        for (auto part: sys->parts){
            space2 = point.space_2(part->pos);
            if (space2<=dist2){
                this->correlationPointSpins[i].push_front(part);
                ++spinsInCurrentPoint;
                spinsInPoint += 1;
            }
        }
        if (spinsInCurrentPoint==0){
            throw(std::invalid_argument("# Corellation point "+std::to_string(i)+" has no spins around. Check your config."));
        }
        spinsInPoint += spinsInCurrentPoint;
        if (spinsInCurrentPoint > maxSpinsInPoint)
            maxSpinsInPoint = spinsInCurrentPoint;
    }
    spinsInPoint /= this->pointCount();

    const double minr2 = this->_minRange * this->_minRange;
    const double maxr2 = this->_maxRange * this->_maxRange;
    double eTemp;

    int i=0;
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
                    
                    if (this->_histogramEnabled){
                        this->correlationNeighboursByPoint[i][partA->Id()].push_front(partB);
                    }

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
        ++i;
    }
    this->correlationPairsNum/=2;

    if (_debug) {
        for (auto partA : sys->parts){
            fprintf(stderr,"# neigh for %zd: ",partA->Id());
            for (auto partB : this->correlationNeighbours[partA->Id()]){
                fprintf(stderr,"%zd,",partB->Id());
            }
            fprintf(stderr,"\n");
        }
    }

    this->cpOld = this->getFullTotal(this->sys);

    if (this->_histogramEnabled){
        dos.resize(-maxSpinsInPoint,maxSpinsInPoint,maxSpinsInPoint*2+1);
        dos.clear();
    }

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
    double addVal = double(this->cpOld)/this->pointCount();
    this->cp += addVal;
    this->cp2 += addVal*addVal;
    if (this->_binder)
        this->cp4 += addVal*addVal*addVal*addVal;

    if (this->_histogramEnabled){
        int i=0;
        for (auto cps: this->correlationPointSpins){
            int cpVal = 0;
            for (auto partA: cps){
                for (auto partB: this->correlationNeighboursByPoint[i][partA->Id()]){
                    cpVal += this->method(partA,partB);
                }
            }
            this->dos[cpVal]++;
            ++i;
        }
    }
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

void CorrelationPointCore::save(unsigned num){
    if (this->_histogramEnabled && !this->_histogramFilename.empty()){
        std::string fname = _histogramFilename;
        fname.replace(fname.find("$"),1,std::to_string(num));
        dos.save(fname);
    }
}