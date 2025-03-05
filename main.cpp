#include "defines.h"

#include <string_view>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <bitset>
#include <gmpxx.h>
#include <chrono>
#include <omp.h>
#include <memory>
#include <argumentum/argparse.h>
#include "CommandLineParameters.h"
#include "ConfigManager.h"
#include "CalculationParameter.h"
#include <inicpp/inicpp.h>
#include "misc.h"
#include "Worker.h"
#include "MagneticSystem.h"

monteCarloStatistics montecarlo(ConfigManager *config){
	unsigned temperatureCount = config->temperatures->size();

	monteCarloStatistics statData;
	statData.foundLowerEnergy = false;
	statData.time_proc_total = 0;

	{ // block to get initial energy
		const Vect field = config->getField();
		statData.initEnergy = config->system->E();
		for (auto &p : config->system->parts)
		{
			statData.initEnergy -= scalar(p.m, field);
		}
		statData.lowerEnergy = statData.initEnergy;
		config->deltaEnergy = fabs(statData.initEnergy * config->getRestartThreshold()); //todo переделать так чтобы эта дельта расчитывалась всего один раз
	}



	vector< optional< pair<double,state_t> > > workerResults(config->temperatures->size());
	vector <shared_ptr<Worker>> workers(config->temperatures->size());
	for (int tt = 0; tt < config->temperatures->size(); ++tt){
		workers[tt] = make_shared<Worker>(tt,config->getSeed() + tt, config);
		config->getParameters(workers[tt]->calculationParameters);
	}

	//новый кусок кода
	for (unsigned phase = 0; phase <= 1; ++phase)
	{
		unsigned calculateSteps;
		if (phase == 0)
			calculateSteps = config->getHeatup();
		else
			calculateSteps = config->getCalculate();

		unsigned eachsteps = config->temperatures->get_each_step();
		if (eachsteps==0) eachsteps = calculateSteps; //случай когда параллельный отжиг отключен
		
		unsigned stepResidual = calculateSteps % eachsteps;
		
		for (unsigned step = 0; step < calculateSteps + stepResidual; step += eachsteps){
			#pragma omp parallel for
			for (int tt = 0; tt < config->temperatures->size(); ++tt)
			{
				workerResults[tt] = workers[tt]->work((step<calculateSteps)?eachsteps:stepResidual, (bool)phase, statData.lowerEnergy);
			}

			// тут проверить нашлась ли во всех потоках энергия ниже начальной или нет
			for (int tt = 0; tt < config->temperatures->size(); ++tt)
			{
				if (workerResults[tt]){
					statData.foundLowerEnergy = true;
					if (workerResults[tt]->first < statData.lowerEnergy){
						statData.lowerEnergy = workerResults[tt]->first;
						statData.lowerEnergyState = workerResults[tt]->second;
						statData.temperatureOfLowerEnergy = tt;
					}
				}
			}
			if (statData.foundLowerEnergy){
                break; //break up 
            }

			for (int tt = 0; tt < config->temperatures->sizeBase(); ++tt)
			{
				if (config->temperatures->size(tt)>1){
					for (int r=0; r<config->temperatures->size(tt)-1;++r){
						double dt = config->temperatures->at(tt,r+1) - config->temperatures->at(tt,r);
						Worker::exchange(workers[config->temperatures->to(tt,r)],workers[config->temperatures->to(tt,r+1)], dt);
					}
				}
			}
		}

		
		if (statData.foundLowerEnergy){
			break; //break up 
		}
	}

	if (!statData.foundLowerEnergy) {
		for (int tt = 0; tt < config->temperatures->size(); ++tt)
		{
			workers[tt]->printout(config->temperatures->at(tt));
			statData.time_proc_total += workers[tt]->duration.count();
		}

		// print out the states and times of running
		printf("###########  end of calculations #############\n");
		printf("#\n");
		printf("###########     final notes:     #############\n");
		for (int tt = 0; tt < config->temperatures->size(); ++tt)
		{
			workers[tt]->printout_service();
		}
		printf("#\n");
	}
	

	return statData;
} 

int main(int argc, char *argv[])
{
	auto time_start = std::chrono::steady_clock::now();

	ConfigManager *config;
	try {
		config = new ConfigManager(argc,argv);
	} catch (const std::string& s){
		cerr<<s<<endl;
		return 0;
	}

	config->printHeader();

	bool programRestarted = false;
	monteCarloStatistics statData;
	state_t finalState = state_t(config->N(),1); // начальное состояние - все единицы. Остальные состояния собираем относительно него
	do {
		statData = montecarlo(config); // запуск самих вычислений
		if (statData.foundLowerEnergy){
			printf("# -- restart MC: found lower energy %g < %g, at %s new state: %s\n",
			   statData.lowerEnergy,
			   statData.initEnergy,
			   config->temperatures->at(statData.temperatureOfLowerEnergy).to_string().c_str(),
			   stateToString(statData.lowerEnergyState).c_str());
			programRestarted = true;
			finalState = xorstate(finalState,statData.lowerEnergyState);
			config->system->applyState(finalState);
		}
	} while(statData.foundLowerEnergy);

	

	auto time_end = std::chrono::steady_clock::now();

	
	int64_t time_total = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
	double speedup = double(statData.time_proc_total) / time_total;
	printf("# total time: %fs, speedup: %f%%, efficiency: %f%%\n", time_total / 1000., speedup * 100, speedup / config->threadCount * 100);

	if (programRestarted){
		printf("\n##### Warning! The program was restarted because it found the lower energy.\n");
		printf("##### But the console output before this moment can not be wiped!\n");
		printf("##### Remove all the result lines before the last line starting with:\n");
		printf("# -- restart MC:\n");
		//printf("# Command to delete this lines:\n");
		//printf("#    perl -p0e -i 's/(# 1:T[^\\n]+\\n).+# -- restart MC: found[^\\n]+\\n/$1/s'   <filename>\n#\n");

		printf("# configuration of the lowest energy: %s\n",stateToString(finalState).c_str());
		if (!config->getNewGSFilename().empty()){
			config->system->save(config->getNewGSFilename());
			printf("# system with found lowest energy is saved to file %s\n",config->getNewGSFilename().c_str());
		}
	}

	delete config;

	return 0;
}