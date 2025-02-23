#include "Worker.h"

void Worker::_startTimer()
{
    _startTimePoint = chrono::steady_clock::now();
}

void Worker::_stopTimer()
{
    duration += chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - _startTimePoint);
}

double Worker::fullRefreshEnergy()
{   
    double eOld = sys.E();
    // add external field
    for (auto p : sys.parts)
    {
        eOld -= p->m.scalar(config->getField());
    }
    return eOld;
}

Worker::Worker(unsigned num, int seed, const ConfigManager* _config) : 
    _num(num),
    stepsMade(0),
    _seed(seed),
    config(_config),
    _previousCalculateStatistics(false),
    e(0, 1024 * 8),
    e2(0, 2048 * 8),
    duration(0)
{
    generator.seed(seed);

    this->sys = PartArray(config->getSystem());

    if (config->isCSV()){
        ConfigManager::setCSVEnergies(sys);
    } else {
        if (config->isPBC())
        {
            ConfigManager::setPBCEnergies(sys);
        }
    }

    // print neighbours and energies
	/*cout<<sys.E()<<endl;
	for (unsigned i=0; i<sys.size(); i++){
		cout<<i<<": ";
		unsigned j=0;
		for (auto p: sys.neighbours[i]){
			cout<<p->Id()<<"("<<sys.eAt(i,j)<<"), ";
			++j;
		}
		cout<<endl;
	}*/
}

optional< pair<double,string> > Worker::work(unsigned steps, bool calculateStatistics, double lowestEnergy)
{
    _startTimer();

    temp_t temperature = config->temperatures->at(_num);

    uniform_int_distribution<int> intDistr(0, sys.size() - 1); // including right edge
    uniform_real_distribution<double> doubleDistr(0, 1);	   // right edge is not included

    //если вычисление статистики было выключено и включилось, инициировать 
    if (calculateStatistics && !_previousCalculateStatistics){ 
        for (auto &cp : calculationParameters){
            cp->init(&sys); // attach the system and calculate the init value
        }
    }

    bool isFoundLowest = false; // нашлась или нет более низкая энергия
    string lowestState;

    bool swapRes;
    unsigned swapNum;
    const Vect field = config->getField();
    double eOld;

    double dE, p, randNum;

    bool acceptSweep;

    for (unsigned step = 0; step < steps; ++step)
    {
        if ((this->stepsMade + step) == 0 || (this->stepsMade + step) % FULL_REFRESH_EVERY == 0)
        {
            eOld = this->fullRefreshEnergy();
            if (config->isRestart() && isFoundLowest){
                break;
            }
        }

        for (unsigned sstep = 0; sstep < sys.size(); ++sstep)
        {

            dE = 0;
            swapNum = intDistr(generator);
            Part *partA = sys.getById(swapNum);

            { // get dE
                unsigned j = 0;

                if (sys.interactionRange() != 0.0)
                {
                    for (Part *neigh : sys.neighbours[swapNum])
                    {
                        if (neigh->state == partA->state) // assume it is rotated, inverse state in mind
                            dE -= 2. * sys.eAt(swapNum, j);
                        else
                            dE += 2. * sys.eAt(swapNum, j);
                        ++j;
                    }
                }
                else
                {
                    for (Part *neigh : sys.parts)
                    {
                        if (partA != neigh)
                        {
                            if (neigh->state == partA->state)
                                dE -= 2. * sys.eAt(swapNum, j);
                            else
                                dE += 2. * sys.eAt(swapNum, j);
                            ++j;
                        }
                    }
                }

                dE += 2 * partA->m.scalar(field);
            }

            acceptSweep = false;
            if (dE < 0 || temperature.t == 0)
            {
                acceptSweep = true;
            }
            else
            {
                p = exp(-dE / temperature.t);
                randNum = doubleDistr(generator);
                if (randNum <= p)
                {
                    acceptSweep = true;
                }
            }

            if (acceptSweep)
            {
                sys.parts[swapNum]->rotate(false);
                eOld += dE;

                for (auto &cp : calculationParameters)
                {
                    cp->iterate(partA->Id());
                }

                if (config->debug)
                {
                    // recalc energy
                    double eTmp = sys.E();

                    // add external field
                    for (auto pt : sys.parts)
                    {
                        eTmp -= pt->m.scalar(field);
                    }

                    if (fabs(eTmp - eOld) > 0.00001)
                    {
                        cerr << "# (dbg main#" << calculateStatistics << ") energy is different. iterative: " << eOld << "; actual: " << eTmp << endl;
                    }
                }

                
                if (config->isRestart() && (eOld - lowestEnergy) < -config->deltaEnergy) // if found lower energy
                {
                    isFoundLowest = true;
                    lowestEnergy = eOld;
                    lowestState = sys.state.toString();
                }
            }
        }

        // update thermodynamic averages (porosyenok ;)
        if (calculateStatistics)
        {
            e += eOld;
            e2 += eOld * eOld;
            for (auto &cp : calculationParameters)
            {
                cp->incrementTotal();
            }

            if (config->getSaveStates()>0 && (stepsMade+step) % config->getSaveStates() == 0){
                sys.save( config->getSaveStateFileName(_num,step) );
            }
        }
    }

    stepsMade += steps;
    _previousCalculateStatistics = calculateStatistics;
    _stopTimer();

    if (isFoundLowest) return make_pair(lowestEnergy,lowestState); else return nullopt;
}

void Worker::printout(temp_t temperature)
{
    mpf_class ee = e / config->getCalculate();
    mpf_class ee2 = e2 / config->getCalculate();

    mpf_class cT = (ee2 - (ee * ee)) / (temperature.t * temperature.t * sys.size());

    {
        gmp_printf("%e %.30Fe %.30Fe %.30Fe %d %d",
                temperature.t, cT.get_mpf_t(), ee.get_mpf_t(), ee2.get_mpf_t(),
                omp_get_thread_num(), _seed);
        for (auto &cp : calculationParameters)
        {
            gmp_printf(" %.30Fe %.30Fe",
                    cp->getTotal(config->getCalculate()).get_mpf_t(),
                    cp->getTotal2(config->getCalculate()).get_mpf_t());
        }
        auto rtime = duration.count();
        printf(" %f", rtime / 1000.);
        printf("\n");
        fflush(stdout);
        for (auto &cp : calculationParameters)
        {
            cp->save(_num);
        }
    }
}

void Worker::printout_service()
{
    printf("#%d, time=%fs, %s, E=%e, final state: %s\n",
        _num,
        duration.count() / 1000.,
        config->temperatures->at(_num).to_string().c_str(),
        this->fullRefreshEnergy(),
        sys.state.toString().c_str());
}

bool Worker::exchange(const shared_ptr<Worker> w1, const shared_ptr<Worker> w2)
{
    return false;
}
