#include "CorrelationPointCore.h"

// @todo Сейчас в correlationNeighbours хранятся все соседи для каждого спина. 
// Но спин может принадлежать сразу нескольким кореляционным точкам, и тогда будет иметь разных соседей.

CorrelationPointCore::CorrelationPointCore(
    const config_section_t &sect, 
    const ConfigManager *conf):
CalculationParameter(sect, conf),
cp(0,1024*8),
cp2(0,2048*8),
currentState(sys->N(),1),
_histogramEnabled(false),
_histogramFilename(""),
_minRange(0)
{

    // парсим точки
    if (!sect.data.contains("points"))
        throw(std::string("Parameter 'points' is required in section "+sect.data.get_name()));
    for (auto & a : sect.data["points"].get_list<inicpp::string_ini_t>()){
        Vect t = strToVect(a);
        this->X.push_back(t.x);
        this->Y.push_back(t.y);
    }

    if (sect.data.contains("minrange")) 
        this->_minRange = sect.data["minrange"].get<inicpp::float_ini_t>();

    if (!sect.data.contains("maxrange"))
        throw(std::string("Parameter 'maxrange' is required in section "+sect.data.get_name()));
    this->_maxRange = sect.data["maxrange"].get<inicpp::float_ini_t>();

    if (!sect.data.contains("distance"))
        throw(std::string("Parameter 'distance' is required in section "+sect.data.get_name()));
    this->_distance = sect.data["distance"].get<inicpp::float_ini_t>();


    // if need to save the histogram to the file
    if (sect.data.contains("histogram") && sect.data["histogram"].get<inicpp::boolean_ini_t>()==true){
        //dollar is replaced by number in time of saving
        this->enableHistogram("histogram_"+this->parameterId()+"_$.txt"); 
    }

    this->findPointsArray();
}

void CorrelationPointCore::findPointsArray()
{
    // clear variables
    correlationPointSpins.clear();
    correlationNeighbours.clear();
    correlationValues.clear();
    correlationPairsNum=0;
    spinsInPoint=0;


    correlationPointSpins.resize(X.size());
    correlationNeighbours.resize(sys->N());
    double space2;

    if (this->_histogramEnabled){
        correlationNeighboursByPoint.resize(X.size());
        for (auto &cps: this->correlationNeighboursByPoint){
            cps.resize(sys->N());
        }
    }

    //first find the spins around each point
    double dist2 = this->_distance * this->_distance;
    unsigned maxSpinsInPoint=0;
    for (int i=0; i<X.size(); ++i){
        Vect point = {X[i], Y[i], 0};
        unsigned spinsInCurrentPoint = 0;
        for (size_t id=0; id < sys->N(); id++){
            space2 = distance_2(point, sys->parts[id].p);
            if (space2<=dist2){
                this->correlationPointSpins[i].push_front(id);
                ++spinsInCurrentPoint;
                spinsInPoint += 1;
            }
        }
        if (spinsInCurrentPoint==0){
            throw(std::string("# Corellation point "+std::to_string(i)+" has no spins around. Check your config."));
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
        for (auto idA: cps){
            for (auto idB: cps){
                if (idA==idB)
                    continue;
                
                space2 = distance_2(sys->parts[idA].p, sys->parts[idB].p);
                if (space2>=minr2 && space2<=maxr2){
                    this->correlationNeighbours[idA].push_front(idB);
                    ++correlationPairsNum;
                    
                    if (this->_histogramEnabled){
                        this->correlationNeighboursByPoint[i][idA].push_front(idB);
                    }

                    eTemp = hamiltonian_dipolar(sys->parts[idA], sys->parts[idB])*-1; 

                    if (eTemp>0)
                        this->correlationValues[std::make_pair(idA,idB)] = 1;
                    else
                        this->correlationValues[std::make_pair(idA,idB)] = -1;
                }
            }
        }
        ++i;
    }
    this->correlationPairsNum/=2;

    if (this->isDebug()) {
        for (size_t idA = 0; idA < sys->N(); idA++){
            fprintf(stderr,"# neigh for %zd: ",idA);
            for (auto idB : this->correlationNeighbours[idA]){
                fprintf(stderr,"%zd,",idB);
            }
            fprintf(stderr,"\n");
        }
    }

    if (this->_histogramEnabled){
        dos.resize(-maxSpinsInPoint,maxSpinsInPoint,maxSpinsInPoint*2+1);
        dos.clear();
    }

    return;
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

void CorrelationPointCore::printHeader() const
{

    printf("# type: %s\n",CorrelationPointCore::name().c_str());
    printf("# id: %s\n",this->parameterId().c_str());
    if (this->_histogramEnabled)
        printf("# histogram: enabled (slow)\n");
    else
        printf("# histogram: disabled\n");

    printf("# distance around point: %.2f\n",this->_distance);
    printf("# minimal interaction distance: %.2f\n",this->_minRange);
    printf("# maximal interaction distance: %.2f\n",this->_maxRange);

    printf("# points: %lu; (%.2f avg spins per point, %.2f avg. neighbours)\n",
            this->X.size(),
            this->spinsInPoint,
            this->correlationPairsNum/double(this->sys->N()));
    printf("#    coordinates: format is <num:(x,y):spins>, <...>, ...\n");
    printf("# 0:(%f,%f):%ld",
        this->X[0],
        this->Y[0],
        std::distance(this->correlationPointSpins[0].begin(),
                        this->correlationPointSpins[0].end()));
    for (int i=1; i<this->X.size(); ++i)
        printf(", %d:(%f,%f):%ld",i,
            this->X[i],
            this->Y[i],
            std::distance(this->correlationPointSpins[i].begin(),
                        this->correlationPointSpins[i].end()));
    printf("\n");

    printf("#\n");
    return;
}

void CorrelationPointCore::init(const state_t &state)
{
    
    this->currentState = state;
    this->cpOld = this->getFullTotal(this->currentState);
}

void CorrelationPointCore::iterate(size_t id){
    this->currentState[id] = - this->currentState[id];
    for (auto partB:this->correlationNeighbours[id]){
        this->cpOld += 2*this->method(id,partB);
    }


    if (this->isDebug()){
        long res = this->getFullTotal(this->currentState);
        if (res!=this->cpOld) 
            cerr<<"# (dbg CorrelationPointCore) total value is different: iterative="<<this->cpOld<<", full="<<res<<endl;
    }
}

void CorrelationPointCore::incrementTotal(){
    this->cp += double(this->cpOld)/this->pointCount();
    this->cp2 += double(this->cpOld * this->cpOld)/(this->pointCount() * this->pointCount());

    if (this->_histogramEnabled){
        int i=0;
        for (auto cps: this->correlationPointSpins){
            int cpVal = 0;
            for (auto idA: cps){
                for (auto idB: this->correlationNeighboursByPoint[i][idA]){
                    cpVal += this->method(idA,idB);
                }
            }
            this->dos[cpVal]++;
            ++i;
        }
    }
}

long CorrelationPointCore::getFullTotal(const state_t &_state) const
{
    long res = 0;
    for (size_t idA = 0; idA < sys->N(); idA++){
        for (auto idB : this->correlationNeighbours[idA]){
            // чтобы отключить занчение текущего состояния системы,
            // которое учитывает method(idA, idB), домножаем на currentState обоих спинов
            res += this->method(idA, idB) * currentState[idA] * currentState[idB] * _state[idA] * _state[idB];
        }
    }
    return res/2;
}

short CorrelationPointCore::method(const size_t idA, const size_t idB) const
{
    return currentState[idA] * currentState[idB] * this->correlationValues.at(std::make_pair(idA, idB));
}

void CorrelationPointCore::save(unsigned num){
    if (this->_histogramEnabled && !this->_histogramFilename.empty()){
        std::string fname = _histogramFilename;
        fname.replace(fname.find("$"),1,std::to_string(num));
        dos.save(fname);
    }
}