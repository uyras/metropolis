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

int main(int argc, char* argv[])
{
    auto time_start = std::chrono::steady_clock::now();

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


    vector<string> finalStates(config.temperatures.size());
    vector<double> finalEnergies(finalStates.size());
    vector< std::chrono::time_point<std::chrono::steady_clock> > temperature_times_start ( finalStates.size() );
    vector< std::chrono::time_point<std::chrono::steady_clock> > temperature_times_end ( finalStates.size() );



    #pragma omp parallel
    {        
        #pragma omp for
        for (int tt=0; tt<config.temperatures.size();  ++tt){

            temperature_times_start[tt] = std::chrono::steady_clock::now();

            std::vector< std::unique_ptr< CalculationParameter > > calculationParameters;
            config.getParameters(calculationParameters);

            const double t = config.temperatures[tt];
            const unsigned trseed = config.getSeed()+tt;
            default_random_engine generator;
            generator.seed(trseed);
            uniform_int_distribution<int> intDistr(0, config.N()-1); // including right edge
            uniform_real_distribution<double> doubleDistr(0,1); // right edge is not included

            /////////// import the system
            mpf_class e(0,1024*8);
            mpf_class e2(0,2048*8);

            PartArray sys(config.getSystem());
            if (config.isPBC()){
                ConfigManager::setPBCEnergies(sys);
            }



            // write neighbours and energies
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
                    }

                    // update thermodynamic averages (porosyenok ;) 
                    if (phase==1){
                        e += eOld;
                        e2 += eOld*eOld;
                        for (auto &cp: calculationParameters){
                                cp->incrementTotal();
                        }
                    }

                }
            }

            e /= config.getCalculate();
            e2 /= config.getCalculate();

            mpf_class cT = (e2 - (e * e))/(t * t * N);


            finalStates[tt] = sys.state.toString();
            finalEnergies[tt] = eOld;
            temperature_times_end[tt] = std::chrono::steady_clock::now();

            #pragma omp critical
            {
                gmp_printf("%e %.30Fe %.30Fe %.30Fe %d %d",
                    t, cT.get_mpf_t(), e.get_mpf_t(), e2.get_mpf_t(), 
                    omp_get_thread_num(), trseed);
                for (auto &cp: calculationParameters){
                    gmp_printf(" %.30Fe %.30Fe",
                        cp->getTotal (config.getCalculate()).get_mpf_t(),
                        cp->getTotal2(config.getCalculate()).get_mpf_t());
                }
                auto rtime = std::chrono::duration_cast<std::chrono::milliseconds>(temperature_times_end[tt]-temperature_times_start[tt]).count();
                printf(" %f",rtime/1000.);
                printf("\n");
                fflush(stdout);
            }
        }
        

    }

    auto time_end = std::chrono::steady_clock::now();
    
    //print out the states and times of running
    int64_t time_proc_total=0;
    printf("###########  end of calculations #############\n");
    printf("#\n");
    printf("###########     final notes:     #############\n");
    for (int tt=0; tt<config.temperatures.size(); ++tt){
        auto rtime = std::chrono::duration_cast<std::chrono::milliseconds>(temperature_times_end[tt]-temperature_times_start[tt]).count();
        printf("#%d, time=%fs, T=%e, E=%e, final state: %s\n",
        tt,
        rtime/1000.,
        config.temperatures[tt],
        finalEnergies[tt],
        finalStates[tt].c_str());
        time_proc_total += rtime;
    }

    printf("#\n");
    int64_t time_total = std::chrono::duration_cast<std::chrono::milliseconds>(time_end-time_start).count();
    double speedup = double(time_proc_total)/time_total;
    printf("# total time: %fs, speedup: %f%, efficiency: %f%\n",time_total/1000., speedup*100, speedup/config.threadCount*100 );

}