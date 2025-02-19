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

std::optional<ConfigManager> readParameters(int argc, char *argv[]){

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
		return {};
	}

	inicpp::config iniconfig;
	if (!commandLineParameters->inifilename.empty())
	{
		iniconfig = inicpp::parser::load_file(commandLineParameters->inifilename);
	}

	ConfigManager config = ConfigManager::init(*(commandLineParameters.get()), iniconfig);

	bool configError = config.check_config();
	if (!configError)
	{
		cerr << "Program stopped with error" << endl;
		return {};
	} else {
		config.printHeader();
	}

	return config;
}

monteCarloStatistics montecarlo(ConfigManager &config){
	unsigned temperatureCount = config.temperatures->size();

	monteCarloStatistics statData;
	statData.foundLowerEnergy = false;
	statData.finalStates.resize(temperatureCount);
	statData.finalEnergies.resize(temperatureCount);
	statData.temperature_times_start.resize(temperatureCount);
	statData.temperature_times_end.resize(temperatureCount);

	{ // block to get initial energy
		const Vect field = config.getField();
		PartArray sys(config.getSystem());
		if (config.isCSV()){
			ConfigManager::setCSVEnergies(sys);
		} else {
			if (config.isPBC())
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
		statData.deltaEnergy = fabs(statData.initEnergy * config.getRestartThreshold());
	}

#pragma omp parallel
	{
#pragma omp for
		for (int tt = 0; tt < config.temperatures->size(); ++tt)
		{
			{
				statData.temperature_times_start[tt] = std::chrono::steady_clock::now();

				std::vector<std::unique_ptr<CalculationParameter>> calculationParameters;
				config.getParameters(calculationParameters);

				const temp_t t = config.temperatures->at(tt);
				const unsigned trseed = config.getSeed() + tt;
				default_random_engine generator;
				generator.seed(trseed);
				uniform_int_distribution<int> intDistr(0, config.N() - 1); // including right edge
				uniform_real_distribution<double> doubleDistr(0, 1);	   // right edge is not included

				
				mpf_class e(0, 1024 * 8);
				mpf_class e2(0, 2048 * 8);

				/////////// duplicate the system
				PartArray sys(config.getSystem());
				if (config.isCSV()){
					ConfigManager::setCSVEnergies(sys);
				} else {
					if (config.isPBC())
					{
						ConfigManager::setPBCEnergies(sys);
					}
				}

				// print neighbours and energies
				/*sys.E();
				for (unsigned i=0; i<sys.size(); i++){
					cout<<i<<": ";
					unsigned j=0;
					for (auto p: sys.neighbours[i]){
						cout<<p->Id()<<"("<<sys.eAt(i,j)<<"), ";
						++j;
					}
					cout<<endl;
				}*/

				const unsigned N = sys.size();

				bool swapRes;
				unsigned swapNum;
				const Vect field = config.getField();
				double eOld;

				double dE, p, randNum;

				double mxOld;
				double myOld;

				bool acceptSweep;

				// phase=0 is the heatup, phase=1 is calculate
				for (unsigned phase = 0; phase <= 1; ++phase)
				{

					// full recalculte energy
					eOld = sys.E();
					// add external field
					for (auto p : sys.parts)
					{
						eOld -= p->m.scalar(field);
					}

					if (phase == 1)
					{
						for (auto &cp : calculationParameters)
						{
							cp->init(&sys); // attach the system and calculate the init value
						}
					}

					unsigned calculateSteps;
					if (phase == 0)
						calculateSteps = config.getHeatup();
					else
						calculateSteps = config.getCalculate();

					for (unsigned step = 0; step < calculateSteps; ++step)
					{
						// full recalculte energy every to avoid FP error collection
						if (step != 0 && step % FULL_REFRESH_EVERY == 0)
						{
							eOld = sys.E();
							// add external field
							for (auto p : sys.parts)
							{
								eOld -= p->m.scalar(field);
							}

							if (statData.foundLowerEnergy){
								//cancel the calculations
								phase = 1; //force go to the phase
								break; //break up the main for loop
							}
						}

						for (unsigned sstep = 0; sstep < N; ++sstep)
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
							if (dE < 0 || t.t == 0)
							{
								acceptSweep = true;
							}
							else
							{
								p = exp(-dE / t.t);
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

								if (phase == 1)
								{
									for (auto &cp : calculationParameters)
									{
										cp->iterate(partA->Id());
									}
								}

								if (config.debug)
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
										cerr << "# (dbg main#" << phase << ") energy is different. iterative: " << eOld << "; actual: " << eTmp << endl;
									}
								}

								
								if (config.isRestart() && (eOld - statData.lowerEnergy) < -statData.deltaEnergy) // if found lower energy
								{
#pragma omp critical
									{
										statData.foundLowerEnergy = 1;
										statData.lowerEnergy = eOld;
										statData.lowerEnergyState = sys.state.toString();
										statData.temperatureOfLowerEnergy = tt;
									}
								}
							}
						}

						// update thermodynamic averages (porosyenok ;)
						if (phase == 1)
						{
							e += eOld;
							e2 += eOld * eOld;
							for (auto &cp : calculationParameters)
							{
								cp->incrementTotal();
							}

							if (config.getSaveStates()>0 && step % config.getSaveStates() == 0){
								sys.save( config.getSaveStateFileName(tt,step) );
							}
						}
					}
				}

				if (!statData.foundLowerEnergy) {
					e /= config.getCalculate();
					e2 /= config.getCalculate();

					mpf_class cT = (e2 - (e * e)) / (t.t * t.t * N);

					statData.finalStates[tt] = sys.state.toString();
					statData.finalEnergies[tt] = eOld;
					statData.temperature_times_end[tt] = std::chrono::steady_clock::now();

	#pragma omp critical
					{
						gmp_printf("%e %.30Fe %.30Fe %.30Fe %d %d",
								t.t, cT.get_mpf_t(), e.get_mpf_t(), e2.get_mpf_t(),
								omp_get_thread_num(), trseed);
						for (auto &cp : calculationParameters)
						{
							gmp_printf(" %.30Fe %.30Fe",
									cp->getTotal(config.getCalculate()).get_mpf_t(),
									cp->getTotal2(config.getCalculate()).get_mpf_t());
						}
						auto rtime = std::chrono::duration_cast<std::chrono::milliseconds>(statData.temperature_times_end[tt] - statData.temperature_times_start[tt]).count();
						printf(" %f", rtime / 1000.);
						printf("\n");
						fflush(stdout);
						for (auto &cp : calculationParameters)
						{
							cp->save(tt);
						}
					}
				}
			}
		}
	}

	return statData;
} 

int main(int argc, char *argv[])
{
	auto time_start = std::chrono::steady_clock::now();

	auto config = readParameters(argc,argv);
	if (!config){
		return 0;
	}

	bool programRestarted = false;
	monteCarloStatistics statData;
	std::string finalState = config->getSystem().state.toString();
	do {
		statData = montecarlo(*config); // запуск самих вычислений
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
	int64_t time_proc_total = 0;
	printf("###########  end of calculations #############\n");
	printf("#\n");
	printf("###########     final notes:     #############\n");
	for (int tt = 0; tt < config->temperatures->size(); ++tt)
	{
		auto rtime = std::chrono::duration_cast<std::chrono::milliseconds>(statData.temperature_times_end[tt] - statData.temperature_times_start[tt]).count();
		printf("#%d, time=%fs, %s, E=%e, final state: %s\n",
			   tt,
			   rtime / 1000.,
			   config->temperatures->at(tt).to_string().c_str(),
			   statData.finalEnergies[tt],
			   statData.finalStates[tt].c_str());
		time_proc_total += rtime;
	}

	printf("#\n");
	int64_t time_total = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
	double speedup = double(time_proc_total) / time_total;
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