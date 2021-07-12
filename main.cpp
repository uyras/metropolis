#include <string_view>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <bitset>
#include <gmpxx.h>
#include <omp.h>
#include <argumentum/argparse.h>
#include "PartArray.h"
#include "Part.h"
#include "CorrelationCore.h"
#include "CorrelationPointCore.h"
#include "MagnetisationCore.h"
#include "CommandLineParameters.h"
#include "ConfigManager.h"
#include "CalculationParameter.h"
#include <inicpp/inicpp.h>

#define VERSION "0.0.3"

bool dbg=false;

int main(int argc, char* argv[])
{
    //get file name
    bool parse_failed=false;

    auto parser = argumentum::argument_parser{};
    parser.config().program("metropolis").description("Program for calculating heat capacity \
            and magnetisation of spin system with dipole-dipole hamiltonian");
    auto commandLineParameters = std::make_shared<CommandLineParameters>();
    parser.params().add_parameters( commandLineParameters );

    auto parseResult = parser.parse_args( argc, argv, 1 );
    if ( !parseResult )
      return 1;
    
    
    inicpp::config iniconfig;
    if (!commandLineParameters->inifilename.empty()){
        iniconfig = inicpp::parser::load_file(commandLineParameters->inifilename);
    }

    ConfigManager config = ConfigManager::init(*(commandLineParameters.get()),iniconfig);

    bool configError = config.check_config();
    if (!configError){
        cerr<<"Program stopped with error"<<endl;
        return 1;
    }

    config.printHeader();



    #pragma omp parallel
    {        
        #pragma omp for
        for (int tt=0; tt<config.temperatures.size();  ++tt){

            std::vector< std::unique_ptr< CalculationParameter > > calculationParameters;
            config.getParameters(calculationParameters);

            const double t = config.temperatures[tt];
            const unsigned trseed = config.getSeed()+tt;
            default_random_engine generator;
            generator.seed(trseed);
            uniform_int_distribution<int> intDistr(0, config.N()-1); // including right edge
            uniform_real_distribution<double> doubleDistr(0,1); // right edge is not included

            /////////// import the system
            mpf_class e(0,1024);
            mpf_class e2(0,2048);

            PartArray sys(config.getSystem());

            bool swapRes;
            unsigned swapNum;
            double eOld = sys.E();
            double dE,p,randNum;


            double mxOld;
            double myOld;

            for (unsigned step=0; step<config.getHeatup(); ++step){
                for (unsigned sstep=0; sstep<sys.size(); ++sstep){
                    swapNum = intDistr(generator);
                    sys.parts[swapNum]->rotate(true);
                    dE = eOld - sys.E();
                    p = exp(dE/t);
                    randNum = doubleDistr(generator);
                    if (dE>0 || (t>0 && randNum <= p)){
                        eOld = sys.E();
                    } else {
                        sys.parts[swapNum]->rotate(true);
                    }
                }
            }


            eOld = sys.E();

            for (auto &cp: calculationParameters){
                cp->init(&sys); //attach the system and calculate the init value
            }

            for (unsigned step=0; step<config.getCalculate(); ++step){
                for (unsigned sstep=0; sstep<sys.size(); ++sstep){
                    swapNum = intDistr(generator);
                    sys.parts[swapNum]->rotate(true);
                    dE = eOld - sys.E();
                    p = exp(dE/t);
                    randNum = doubleDistr(generator);
                    if (dE>0 || (t>0 && randNum <= p)){
                        Part* partA = sys.parts[swapNum];
                        eOld = sys.E();

                        for (auto &cp: calculationParameters){
                            cp->iterate(partA->Id());
                        }
                    } else {
                        sys.parts[swapNum]->rotate();
                    }
                    e += eOld;
                    e2 += eOld*eOld;
                    for (auto &cp: calculationParameters){
                            cp->incrementTotal();
                    }
                }
            }

            unsigned totsteps = config.getCalculate()*sys.size();

            e /= totsteps;
            e2 /= totsteps;

            mpf_class cT = (e2 - (e * e))/(t * t * sys.size());

            #pragma omp critical
            {
                printf("%e %e %e %e %d %d",
                    t, cT.get_d(), e.get_d(), e2.get_d(), 
                    omp_get_thread_num(), trseed);
                for (auto &cp: calculationParameters){
                    printf(" %e %e",cp->getTotalDouble(totsteps),cp->getTotal2Double(totsteps));
                }
                printf("\n");
                fflush(stdout);
            }
        }
        

    }
    

}