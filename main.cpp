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
#include <argumentum/argparse.h>
#include "PartArray.h"
#include "Part.h"
#include "CorrelationCore.h"
#include "CorrelationPointCore.h"
#include "MagnetisationCore.h"
#include "MagnetisationLengthCore.h"
#include "CommandLineParameters.h"
#include "ConfigManager.h"
#include "CalculationParameter.h"
#include <inicpp/inicpp.h>
#include "misc.h"
#include "Worker.h"

ConfigManager* readParameters(int argc, char *argv[]){

	// get file name
	bool parse_failed = false;

	auto parser = argumentum::argument_parser{};
	parser.config().program("metropolis").description("Program for calculating heat capacity \
            and magnetisation of spin system with dipole-dipole hamiltonian v." +
													  std::string(METROPOLIS_VERSION));
	auto commandLineParameters = std::make_shared<CommandLineParameters>();
	parser.params().add_parameters(commandLineParameters);

	auto parseResult = parser.parse_args(argc, argv, 1);

	if (!parseResult)
	{
		if (commandLineParameters && commandLineParameters->showExample)
		{
			std::cout << endl;
			std::cout << "##########################################" << endl;
			std::cout << "######## contents of example.ini: ########" << endl;
			std::cout << "##########################################" << endl;
			std::cout << endl;
			std::cout << example_string << endl;
		}
		return nullptr;
	}

	inicpp::config iniconfig;
	if (!commandLineParameters->inifilename.empty())
	{
		iniconfig = inicpp::parser::load_file(commandLineParameters->inifilename);
	}

	ConfigManager *config = new ConfigManager(*(commandLineParameters.get()), iniconfig);

	bool configError = config->check_config();
	if (!configError)
	{
		cerr << "Program stopped with error" << endl;
		return nullptr;
	} else {
		config->printHeader();
	}

	return config;
}

monteCarloStatistics montecarlo(ConfigManager *config){
	unsigned temperatureCount = config->temperatures->size();

	monteCarloStatistics statData;
	statData.foundLowerEnergy = false;

	{ // block to get initial energy
		const Vect field = config->getField();
		PartArray sys(config->getSystem());
		if (config->isCSV()){
			ConfigManager::setCSVEnergies(sys);
		} else {
			if (config->isPBC())
			{
				ConfigManager::setPBCEnergies(sys);
			}
		}
		statData.initEnergy = sys.E();
		for (auto p : sys.parts)
		{
			statData.initEnergy -= p->m.scalar(field);
		}
		statData.lowerEnergy = statData.initEnergy;
		config->deltaEnergy = fabs(statData.initEnergy * config->getRestartThreshold()); //todo переделать так чтобы эта дельта расчитывалась всего один раз
	}



	vector< optional< pair<double,string> > > workerResults(config->temperatures->size());
	vector <Worker> workers;
	for (int tt = 0; tt < config->temperatures->size(); ++tt){
		workers.emplace_back(tt,config->getSeed() + tt, config);
		config->getParameters(workers[tt].calculationParameters);
	}

	//новый кусок кода
	for (unsigned phase = 0; phase <= 1; ++phase)
	{
		unsigned eachsteps = config->temperatures->get_each_step();
		unsigned calculateSteps;
		if (phase == 0)
			calculateSteps = config->getHeatup();
		else
			calculateSteps = config->getCalculate();
		
		unsigned stepResidual = calculateSteps % eachsteps;
		
		for (unsigned step = 0; step < calculateSteps + stepResidual; step += eachsteps){
			#pragma omp parallel for
			for (int tt = 0; tt < config->temperatures->size(); ++tt)
			{
				workerResults[tt] = workers[tt].work((step<calculateSteps)?eachsteps:stepResidual, (bool)phase, statData.lowerEnergy);
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

			for (int tt = 0; tt < config->temperatures->size()-1; ++tt)
			{
				Worker::exchange(workers[tt],workers[tt+1]);
			}
		}

		
		if (statData.foundLowerEnergy){
			break; //break up 
		}
	}

	if (!statData.foundLowerEnergy) {
		for (int tt = 0; tt < config->temperatures->size(); ++tt)
		{
			workers[tt].printout(config->temperatures->at(tt));
			statData.time_proc_total += workers[tt].duration.count();
		}
	}
	

	return statData;
} 

int main(int argc, char *argv[])
{
	auto time_start = std::chrono::steady_clock::now();

	ConfigManager *config = readParameters(argc,argv);
	if (!config){
		return 0;
	}

	bool programRestarted = false;
	monteCarloStatistics statData;
	std::string finalState = config->getSystem().state.toString();
	do {
		statData = montecarlo(config); // запуск самих вычислений
		if (statData.foundLowerEnergy){
			config->applyState(statData.lowerEnergyState);
			printf("# -- restart MC: found lower energy %g < %g, at %s new state: %s\n",
			   statData.lowerEnergy,
			   statData.initEnergy,
			   config->temperatures->at(statData.temperatureOfLowerEnergy).to_string().c_str(),
			   statData.lowerEnergyState.c_str());
			programRestarted = true;
			finalState = xorstr(finalState,statData.lowerEnergyState);
		}
	} while(statData.foundLowerEnergy);

	

	auto time_end = std::chrono::steady_clock::now();

	// print out the states and times of running
	printf("###########  end of calculations #############\n");
	printf("#\n");
	printf("###########     final notes:     #############\n");

	printf("#\n");
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

		printf("# configuration of the lowest energy: %s\n",finalState.c_str());
		if (!config->getNewGSFilename().empty()){
			config->saveSystem(config->getNewGSFilename());
			printf("# system with found lowest energy is saved to file %s\n",config->getNewGSFilename().c_str());
		}
	}

	return 0;
}