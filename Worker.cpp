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
    double eOld = 0;
    for (size_t i=0; i<config->system->size(); i++){
        size_t nf = config->system->neighbours_from[i];
        size_t nc = config->system->neighbours_count[i];
        for (size_t k = nf; k < nf + nc; k++){
            eOld +=  config->system->eMatrix[k] * state[i] * state[config->system->neighbourNums[k]];
        }

        eOld -= scalar(config->system->parts[i].m, config->getField()) * state[i];
    }

    return eOld;
}

Worker::Worker(unsigned num, int seed, const ConfigManager* _config) : 
    _num(num),
    stepsMadeDummy(0),
    stepsAcceptedDummy(0),
    stepsMadeCalculate(0),
    stepsAcceptedCalculate(0),
    _seed(seed),
    config(_config),
    _previousCalculateStatistics(false),
    e(0, 1024 * 8),
    e2(0, 2048 * 8),
    duration(0)
{
    generator.seed(seed);
    state.clear();
    state.resize(config->N(),1);

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

optional< pair<double,state_t> > Worker::work(unsigned steps, bool calculateStatistics, double lowestEnergy)
{
    _startTimer();

    temp_t temperature = config->temperatures->at(_num);

    uniform_int_distribution<int> intDistr(0, config->N() - 1); // including right edge
    uniform_real_distribution<double> doubleDistr(0, 1);	   // right edge is not included

    //если вычисление статистики было выключено и включилось, инициировать 
    if (calculateStatistics && !_previousCalculateStatistics){ 
        for (auto &cp : calculationParameters){
            cp->init(state); // attach the system and calculate the init value
        }
    }

    bool isFoundLowest = false; // нашлась или нет более низкая энергия
    state_t lowestState;

    bool swapRes;
    unsigned swapNum;
    const Vect field = config->getField();

    double dE, p, randNum;

    bool acceptSweep;
    const unsigned stepsMade = (calculateStatistics) ? this->stepsMadeCalculate : this->stepsMadeDummy;

    for (unsigned step = 0; step < steps; ++step)
    {
        if ((stepsMade + step) == 0 || (stepsMade + step) % FULL_REFRESH_EVERY == 0)
        {
            eActual = this->fullRefreshEnergy();
            if (config->isRestart() && isFoundLowest){
                break;
            }
        }

        for (unsigned sstep = 0; sstep < config->N(); ++sstep)
        {

            dE = 0;
            swapNum = intDistr(generator);

            { // get dE
                for (
                    size_t neigh = config->system->neighbours_from[swapNum]; 
                    neigh < config->system->neighbours_from[swapNum] + config->system->neighbours_count[swapNum]; 
                    neigh++)
                {
                    if (state[config->system->neighbourNums[neigh]] == state[swapNum]) // assume it is rotated, inverse state in mind
                        dE -= 2. * config->system->eMatrix[neigh];
                    else
                        dE += 2. * config->system->eMatrix[neigh];
                }

                dE += 2 * scalar(config->system->parts[swapNum].m, field);
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
                if (calculateStatistics) this->stepsAcceptedCalculate++; else this->stepsAcceptedDummy++;
                state[swapNum] *= -1;
                eActual += dE;

                for (auto &cp : calculationParameters)
                {
                    cp->iterate(swapNum);
                }

                if (config->debug)
                {
                    // recalc energy
                    double eTmp = fullRefreshEnergy();

                    if (fabs(eTmp - eActual) > 0.00001)
                    {
                        cerr << "# (dbg main#" << calculateStatistics << ") energy is different. iterative: " << eActual << "; actual: " << eTmp << endl;
                    }
                }

                
                if (config->isRestart() && (eActual - lowestEnergy) < -config->deltaEnergy) // if found lower energy
                {
                    isFoundLowest = true;
                    lowestEnergy = eActual;
                    lowestState = state;
                }
            }
        }

        // update thermodynamic averages (porosyenok ;)
        if (calculateStatistics)
        {
            e += eActual;
            e2 += eActual * eActual;
            for (auto &cp : calculationParameters)
            {
                cp->incrementTotal();
            }

            if (config->getSaveStates()>0 && (stepsMade+step) % config->getSaveStates() == 0){
                config->system->save( config->getSaveStateFileName(_num,step), state );
            }
        }
    }

    if (calculateStatistics) this->stepsMadeCalculate += steps; else this->stepsMadeDummy += steps;
    _previousCalculateStatistics = calculateStatistics;
    _stopTimer();

    if (isFoundLowest) return make_pair(lowestEnergy,lowestState); else return nullopt;
}

void Worker::printout(temp_t temperature)
{
    mpf_class ee = e / config->getCalculate();
    mpf_class ee2 = e2 / config->getCalculate();

    mpf_class cT = (ee2 - (ee * ee)) / (temperature.t * temperature.t * config->N());

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
        stateToString(state).c_str());
}

bool Worker::exchange(const shared_ptr<Worker> w1, const shared_ptr<Worker> w2, double dT)
{
    uniform_real_distribution<double> doubleDistr(0, 1);	   // right edge is not included
    double dE = w2->eActual - w1->eActual;
    double p = exp(dE/dT);
    double randNum = doubleDistr(w1->generator);
    if (randNum <= p)
    {
        //todo тут добавить обмен конфигурациями между репликами. Но для этого сперва надо переделать всю магнитную систему
        return true;
    } else return false;
}
