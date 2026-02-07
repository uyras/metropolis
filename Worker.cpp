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
    for (size_t i=0; i<config->N(); i++){
        size_t nf = config->system->neighbours_from[i];
        size_t nc = config->system->neighbours_count[i];
        for (size_t k = nf; k < nf + nc; k++){
            eOld +=  config->system->eMatrix[k] * state[i] * state[config->system->neighbourNums[k]] / 2;
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
    ptSuccessExchangesDummy(0),
    ptSuccessExchangesCalculate(0),
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

    //создаем алиас для объекта магнитной системы
    const MagneticSystem &sys = *(config->system);

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
    size_t neighbours_from; //порядковый номер в массиве соседей, начало
    size_t neighbours_to; //порядковый номер в массиве соседей, конец
    const Vect &field = config->getField();

    double dE, p, randNum;

    bool acceptSweep;
    unsigned acceptedCount;
    const unsigned stepsMade = (calculateStatistics) ? this->stepsMadeCalculate : this->stepsMadeDummy;

    /**
     * Счетчик, сколько шагов нужно сделать в начале
     */
    unsigned ptTrialStepsMade = (config->pt_enabled()) ? 1 : 0;

    for (int step = 0; step < steps; ++step)
    {
        acceptedCount = 0;
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
            neighbours_from = sys.neighbours_from[swapNum];
            neighbours_to = sys.neighbours_from[swapNum] + sys.neighbours_count[swapNum];

            { // get dE
                for (size_t neigh = neighbours_from; neigh < neighbours_to; neigh++)
                {
                    if (state[sys.neighbourNums[neigh]] == state[swapNum]) // assume it is rotated, inverse state in mind
                        dE -= 2. * sys.eMatrix[neigh];
                    else
                        dE += 2. * sys.eMatrix[neigh];
                }

                dE += 2 * scalar(sys.parts[swapNum].m, field) * state[swapNum];
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
                acceptedCount++;
                state[swapNum] *= -1;
                eActual += dE;

                if (calculateStatistics) {
                    for (auto &cp : calculationParameters)
                    {
                        cp->iterate(swapNum);
                    }
                }

                if (config->debug)
                {
                    // recalc energy
                    double eTmp = fullRefreshEnergy();

                    if (fabs(eTmp - eActual) > 0.00001)
                    {
                        cerr << "# (dbg main#" << _num << ") energy is different. iterative: " << eActual << "; actual: " << eTmp << endl;
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

        // если нужно сделать прогревной шаг И это первый МК-шаг в текущей функции И это НЕ первый запуск функции worker
        if (ptTrialStepsMade && !step && (this->stepsMadeCalculate+this->stepsMadeDummy)) {
            step--; //делаем число шагов -1, чтобы он обнулился в следующем цикле
            ptTrialStepsMade--;
            //cerr<<"# make one trial step worker="<<this->_num<<endl;
        } else {
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
                    sys.save( config->getSaveStateFileName(_num,step), state );
                }
            }

            if (calculateStatistics) this->stepsAcceptedCalculate += double(acceptedCount) / config->N(); 
            else this->stepsAcceptedDummy += double(acceptedCount) / config->N();
        }
    }

    if (calculateStatistics) this->stepsMadeCalculate += steps; else this->stepsMadeDummy += steps;
    _previousCalculateStatistics = calculateStatistics;
    _stopTimer();

    if (isFoundLowest) return make_pair(lowestEnergy,lowestState); else return nullopt;
}

void Worker::printout(temp_t temperature)
{

    //print lines
    mpf_class ee = e / config->getCalculate();
    mpf_class ee2 = e2 / config->getCalculate();

    mpf_class cT = (ee2 - (ee * ee)) / (temperature.t * temperature.t * config->N());

    gmp_printf("%e %.10Fe %.10Fe %.10Fe %d %d",
            temperature.t, cT.get_mpf_t(), ee.get_mpf_t(), ee2.get_mpf_t(),
            omp_get_thread_num(), _seed);

    for (auto &cp : calculationParameters)
    {
        gmp_printf(" %.10Fe %.10Fe",
                cp->getTotal(config->getCalculate()).get_mpf_t(),
                cp->getTotal2(config->getCalculate()).get_mpf_t());
    }

    printf(" %f %f",
        this->stepsAcceptedDummy/this->stepsMadeDummy, 
        this->stepsAcceptedCalculate/this->stepsMadeCalculate
    );

    if (config->pt_enabled()){
        double exchange_count_d = ceil(config->getHeatup()/config->temperatures->get_each_step());
        double exchange_count_c = ceil(config->getCalculate()/config->temperatures->get_each_step());
        printf(" %f %f",
            this->ptSuccessExchangesDummy/exchange_count_d, 
            this->ptSuccessExchangesCalculate/exchange_count_c
        );
    }

    printf(" %lu %lu %f %e",
        temperature.num_base,
        temperature.num_replica,
        duration.count() / 1000.,
        this->fullRefreshEnergy()
    );
    printf("\n");
    fflush(stdout);
    for (auto &cp : calculationParameters)
    {
        cp->save(_num);
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

bool Worker::exchange(shared_ptr<Worker> w1, shared_ptr<Worker> w2, double dBeta)
{
    uniform_real_distribution<double> doubleDistr(0, 1);	   // right edge is not included
    double dE = w2->eActual - w1->eActual;
    double p = exp(dE*dBeta);
    double randNum = doubleDistr(w1->generator);
    if (randNum <= p)
    {
        swap(w1->state,w2->state);
        swap(w1->eActual,w2->eActual);
        //todo тут добавить обмен между доп параметрами, так как у основных систем резко сменилось состояние
        return true;
    } else return false;
}
