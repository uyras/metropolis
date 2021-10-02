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

int main(int argc, char* argv[])
{
    //get file name
    bool parse_failed=false;

    auto parser = argumentum::argument_parser{};
    parser.config().program("metropolis").description("Program for calculating heat capacity \
            and magnetisation of spin system with dipole-dipole hamiltonian v."+std::string(METROPOLIS_VERSION));
    auto commandLineParameters = std::make_shared<CommandLineParameters>();
    parser.params().add_parameters( commandLineParameters );

    auto parseResult = parser.parse_args( argc, argv, 1 );

    if ( !parseResult ){
        if (commandLineParameters && commandLineParameters->showExample){
            std::cout<<endl;
            std::cout<<"##########################################"<<endl;
            std::cout<<"######## contents of example.ini: ########"<<endl;
            std::cout<<"##########################################"<<endl;
            std::cout<<endl;
            std::cout<<example_string<<endl;
        }
      return 1;
    }
    
    
    
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

            const unsigned N = sys.size();

            bool swapRes;
            unsigned swapNum;
            const Vect field = config.getField();
            double eOld;

            double dE,p,randNum;

            double mxOld;
            double myOld;

            bool acceptSweep;

            //phase=0 is the heatup, phase=1 is calculate
            for (unsigned phase=0; phase<=1; ++phase){

                // full recalculte energy
                eOld = sys.E();
                // add external field
                for (auto p: sys.parts){
                    eOld -= p->m.scalar(field);
                }

                if (phase==1){
                    sys.save("tmp.mfsys");
                    for (auto &cp: calculationParameters){
                        cp->init(&sys); //attach the system and calculate the init value
                    }
                }

                unsigned calculateSteps;
                if (phase == 0) 
                    calculateSteps = config.getHeatup();
                else
                    calculateSteps = config.getCalculate();


                for (unsigned step=0; step<calculateSteps; ++step){

                    // full recalculte energy every to avoid FP error collection
                    if (step != 0 && step % FULL_REFRESH_EVERY == 0){
                        eOld = sys.E();
                        // add external field
                        for (auto p: sys.parts){
                            eOld -= p->m.scalar(field);
                        }
                    }

                    for (unsigned sstep=0; sstep<N; ++sstep){
                        dE = 0;
                        swapNum = intDistr(generator);
                        Part *partA = sys.getById(swapNum);

                        { // get dE
                            unsigned j=0;

                            if (sys.interactionRange()!=0.0){
                                for (Part* neigh : sys.neighbours[swapNum]){
                                    if (neigh->state==partA->state) //assume it is rotated, inverse state in mind
                                        dE -=  2. * sys.eAt(swapNum,j);
                                    else
                                        dE += 2. * sys.eAt(swapNum,j);
                                    ++j;
                                }
                            } else {
                                for (Part* neigh : sys.parts){
                                    if (partA!=neigh){
                                        if (neigh->state==partA->state)
                                            dE -=  2. * sys.eAt(swapNum,j);
                                        else
                                            dE += 2. * sys.eAt(swapNum,j);
                                        ++j;
                                    }
                                }
                            }

                            dE += 2 * partA->m.scalar(field);
                        }

                        acceptSweep = false;
                        if (dE<0 || t == 0){
                            acceptSweep = true;
                        } else {
                            p = exp(-dE/t);
                            randNum = doubleDistr(generator);
                            if (randNum <= p){
                                acceptSweep = true;
                            }
                        }

                        if (acceptSweep){
                            sys.parts[swapNum]->rotate();
                            eOld += dE;

                            if (phase==1){
                                for (auto &cp: calculationParameters){
                                    cp->iterate(partA->Id());
                                }
                            }

                            if (config.debug){
                                //recalc energy
                                double eTmp = sys.E();

                                // add external field
                                for (auto pt: sys.parts){
                                    eTmp -= pt->m.scalar(field);
                                }

                                if (fabs( eTmp-eOld ) > 0.00001 ){
                                    cerr<<"# (dbg main#"<<phase<<") energy is different. iterative: "<<eOld<<"; actual: "<<eTmp<<endl;
                                }
                            }
                        }

                        if (phase==1){
                            e += eOld;
                            e2 += eOld*eOld;
                            for (auto &cp: calculationParameters){
                                    cp->incrementTotal();
                            }
                        }
                    }
                }
            }

            unsigned totsteps = config.getCalculate()*N;

            e /= totsteps;
            e2 /= totsteps;

            mpf_class cT = (e2 - (e * e))/(t * t * N);

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