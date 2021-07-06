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
#include "CorrelationParameters.h"
#include "CorrelationPointParameters.h"
#include "MagnetisationParameters.h"

#define VERSION "0.0.3"

std::string filename;
std::vector<double> temperatures; // system temperature
unsigned hSteps;        // heatup steps
unsigned cSteps;        // calculate steps
unsigned rseed;         // random seed value
double iRange;          // interaction range
bool dbg=false;
unsigned N;

using namespace std;
using namespace argumentum;

int main(int argc, char* argv[])
{
    bool parse_failed=false;

    auto parser = argumentum::argument_parser{};
    auto params = parser.params();

    parser.config().program("metropolis")
        .description("Program for calculating heat capacity \
        and magnetisation of spin system with dipole-dipole hamiltonian");
    params.add_parameter(temperatures,"t","temperature").minargs(1)
        .help("Temperature of the system, in units D. \
            May be the space-separated list of temperatures, \
            in this case it will be created several threads for each temperature.\
            Another variant is to save temperatures to the text file, \
            one temperature per line and use this as: metropolis @param.txt");
    params.add_parameter(filename,"-f","--filename").nargs(1).required().metavar("FILE.mfsys")
        .help("Path to text file with structure of the system. \
            Format is the mfsys file.");
    params.add_parameter(hSteps,"-p","--heatup").nargs(1).required().metavar("STEPS")
        .help("Number of steps needed to bring the configuration to \
            the stable state. The single step here consists of several MC sweeps (trials)");
    params.add_parameter(cSteps,"-c","--calculate").nargs(1).required().metavar("STEPS")
        .help("Number of steps performed to Metropolis algorithm.\
            The single step here consists of several MC sweeps (trials)");
    params.add_parameter(iRange,"-r","--range").absent(0.).nargs(1).metavar("RANGE")
        .help("Maximal distance between particles which interaction counts in total energy. \
            Partices located further are considered non-interacting.\
            Default value is 0 which means all-to-all interaction.");
    params.add_parameter(rseed,"-s","--seed").absent(0.).nargs(1).metavar("SEED")
        .help("Random seed number. In case of parallel execution\
        it will be SEED+TEMPNUM, where TEMPNUM is the sequential number \
        of temperature in list. Default is 0.");

    params.add_command<CorrelationParameters>("correlation").help("Gather also correlation parameter. \
        This option has additional parameters. For information print: metropolis correlation --help.");
    params.add_command<CorrelationPointParameters>("corrpoint")
    .help("Gather also correlation parameter with spins around set of points. \
        This option has additional parameters. For information print: metropolis corrpoint --help.");
    params.add_command<MagnetisationParameters>("magnetisation")
    .help("Gather also the magnetisation information along selected axis. \
        This option has additional parameters. For information print: metropolis magnetisation --help.");
    params.add_command<MagnetisationParameters>("magnetisation2")
    .help("The same as magnetisation, ndeeded for gathering several parameters in one run.");
    params.add_command<MagnetisationParameters>("magnetisation3")
    .help("The same as magnetisation, ndeeded for gathering several parameters in one run.");
    params.add_command<MagnetisationParameters>("magnetisation4")
    .help("The same as magnetisation, ndeeded for gathering several parameters in one run.");
    params.add_command<MagnetisationParameters>("magnetisation5")
    .help("The same as magnetisation, ndeeded for gathering several parameters in one run.");


    auto res = parser.parse_args( argc, argv, 1 );

    if ( !res )
      return 1;

    auto correlationOptions = static_pointer_cast<CorrelationParameters>( res.findCommand("correlation") );
    auto correlationPointOptions = static_pointer_cast<CorrelationPointParameters>( res.findCommand("corrpoint") );
    std::vector< std::shared_ptr< MagnetisationParameters> > magnetisationOptions;
    for (std::string commandName: {"magnetisation", 
                                   "magnetisation2", 
                                   "magnetisation3", 
                                   "magnetisation4", 
                                   "magnetisation5"}){
        auto tmp = static_pointer_cast<MagnetisationParameters>( res.findCommand(commandName) );
        if (tmp)
            magnetisationOptions.push_back(tmp);
    }
    
    ifstream f(filename);
    if (!f.is_open()) {
        cout<<"file read error: "<<filename<<endl;
        return 1;
    }
    f.close();
    
    if (cSteps<1){
        cerr<<"error! --calculate should be greather than 0!"<<endl;
        return 1;
    }


    PartArray bigsys;
    bigsys.load(filename.c_str());
    bigsys.state.hardReset();
    bigsys.setInteractionRange(iRange);
    N=bigsys.size();
    if (dbg){
        cerr<<"# file '"<<filename<<"': imported "<<bigsys.size()<<" parts"<<endl;
        cerr<<"# Energy: "<<bigsys.E()<<" range: "<<iRange<<endl;
    }



    if (correlationOptions){
        if (!correlationOptions->check(N)){
                return 1;
        }
    }

    if (correlationPointOptions){
        if (correlationPointOptions->X.size() != correlationPointOptions->Y.size()){
                cerr<<"Parameters -X, and -Y should \
                have the same count of values."<<endl;
                return 1;
        }
    }

    for (auto mo: magnetisationOptions){
        if (!mo->check(N)){
            return 1;
        }
    }

    //obtain the avaliable number of threads
    int threadCount = 0;
    #pragma omp parallel
    {
        #pragma omp single
        {
            threadCount = omp_get_num_threads();
        }
    }

    printf("# Metropolis algorithm for calculating heating capacity v%s\n",VERSION);
    printf("#  filename: %s\n",filename.c_str());
    printf("#     rseed: %d+<tempNum>\n",rseed);
    printf("#    System: %d spins, %f interaction range\n",N,iRange);
    printf("#        MC: %u heatup, %u compute steps\n",hSteps,cSteps);
    printf("#   threads: %d\n",threadCount);
    printf("#    temps.: %d pcs. from %e to %e\n",
        temperatures.size(),
        std::min_element(temperatures.begin(),temperatures.end()).operator*(),
        std::max_element(temperatures.begin(),temperatures.end()).operator*());
    printf("# temp list: %e",temperatures[0]);
    for (int tt=1; tt<temperatures.size(); ++tt)
        printf(",%e",temperatures[tt]);
    printf("\n");


    if (correlationOptions){
        printf("#   correl.: %d cores\n",correlationOptions->minRange.size());
    
        for (int i=0; i<correlationOptions->minRange.size(); ++i){
            CorrelationCore cctemp(
                correlationOptions->minRange[i], 
                correlationOptions->maxRange[i],
                correlationOptions->methodVar[i],
                correlationOptions->spins[i]);
            cctemp.init(bigsys);
            if (cctemp.correlationPairsNum==0){
                cerr<<"Check the correlation range #"<<i<<"."<<endl;
                cerr<<"Could not find any spin pairs within this distance!"<<endl;
                return 1;
            }
            int totalSpins = correlationOptions->spins[i].size();
            if (totalSpins==0) totalSpins=N;
            printf("# core #%03d: %.2f min range, %.2f max range, %.2f avg. neigh.",i,
                correlationOptions->minRange[i], correlationOptions->maxRange[i],
                cctemp.correlationPairsNum/double(totalSpins));
            if (correlationOptions->methodVar[i]==1)
                printf(" method: XOR");
            if (correlationOptions->methodVar[i]==2)
                printf(" method: E/|E|");
            if (correlationOptions->methodVar[i]==3)
                printf(" method: scalar");
            
            printf(" spins: ");
            if (correlationOptions->spins[i].size()==0){
                printf("All\n");
            } else {
                printf("%d",correlationOptions->spins[i][0]);
                for (int ss=1; ss<correlationOptions->spins[i].size(); ++ss){
                    printf(",%d",correlationOptions->spins[i][ss]);
                }
                printf("\n");
            }
        }
    }

    if (correlationPointOptions){
        CorrelationPointCore cpcTemp(
            correlationPointOptions->X,
            correlationPointOptions->Y,
            correlationPointOptions->distance,
            correlationPointOptions->minRange,
            correlationPointOptions->maxRange);
        
        cpcTemp.init(bigsys);
        if (cpcTemp.correlationPairsNum==0){
            cerr<<"Check the correlation point ranges."<<endl;
            cerr<<"Could not find any spin pairs within this distance!"<<endl;
            return 1;
        }
        printf("# corrPoint: %.2f distance, %.2f min range, %.2f max range\n",
                correlationPointOptions->distance,
                correlationPointOptions->minRange, 
                correlationPointOptions->maxRange);
        printf("#            %d points, %.2f avg spins per point, %.2f avg. neigh\n",
                correlationPointOptions->X.size(),
                cpcTemp.spinsInPoint,
                cpcTemp.correlationPairsNum/double(N));
        printf("#    coord.: format is <num:(x,y):spins>, <...>, ...\n# 0:(%f,%f):%d",
            correlationPointOptions->X[0],
            correlationPointOptions->Y[0],
            std::distance(cpcTemp.correlationPointSpins[0].begin(),
                          cpcTemp.correlationPointSpins[0].end()));
        for (int i=1; i<correlationPointOptions->X.size(); ++i)
            printf(", %d:(%f,%f):%d",i,
                correlationPointOptions->X[i],
                correlationPointOptions->Y[i],
                std::distance(cpcTemp.correlationPointSpins[i].begin(),
                          cpcTemp.correlationPointSpins[i].end()));
        printf("\n");

    }

    for (auto mo: magnetisationOptions){
        MagnetisationCore mcTemp(mo->axis,mo->_spins);
        mcTemp.init(bigsys);
        printf("#     magn.: v(%f;%f), start val %f,",mo->axis.x,mo->axis.y,mcTemp.getMOFull()/mo->_spins.size());
        printf(" spins: ");
            if (mo->_spins.size()==N){
                printf("All\n");
            } else {
                printf("%d",mo->_spins[0]);
                for (int ss=1; ss<mo->_spins.size(); ++ss){
                    printf(",%d",mo->_spins[ss]);
                }
                printf("\n");
            }
    }

    printf("# T C(T)/N <E> <E^2> <mx> <mx^2> <my> <my^2> threadId seed");
    if (correlationOptions)
        for (int i=0; i<correlationOptions->minRange.size(); ++i) 
            printf(" <cc#%03d> <cc^2#%03d>",i,i);
    if (correlationPointOptions)
        printf(" <cp> <cc^2>");
    for (int i=0; i<magnetisationOptions.size(); ++i)
        printf(" <mo#%03d> <mo#%03d^2>",i,i);
    printf("\n");
    fflush(stdout);


return 0;
    #pragma omp parallel
    {        
        #pragma omp for
        for (int tt=0; tt<temperatures.size();  ++tt){

            std::vector<CorrelationCore> correlationCores;
            std::shared_ptr<CorrelationPointCore> correlationPointCore;
            std::vector<MagnetisationCore> magnetisationCores;

            const double t = temperatures[tt];
            const unsigned trseed = rseed+tt;
            default_random_engine generator;
            generator.seed(trseed);
            uniform_int_distribution<int> intDistr(0,N-1); // including right edge
            uniform_real_distribution<double> doubleDistr(0,1); // right edge is not included

            /////////// import the system
            mpf_class mx(0,1024);
            mpf_class my(0,1024);
            mpf_class e(0,1024);
            mpf_class mx2(0,2048);
            mpf_class my2(0,2048);
            mpf_class e2(0,2048);

            PartArray sys = bigsys;

            ///////// define neighbours for correlation parameter
            if (correlationOptions){
                for (int i=0; i<correlationOptions->maxRange.size(); ++i){
                    correlationCores.push_back(CorrelationCore(correlationOptions->minRange[i],
                                                                correlationOptions->maxRange[i],
                                                                correlationOptions->methodVar[i],
                                                                correlationOptions->spins[i]));
                    correlationCores[i].init(sys);
                    //correlationCores[i].dbg = dbg;
                }
            }

            ///////// setup correlation point parameter
            if (correlationPointOptions){
                correlationPointCore = make_shared<CorrelationPointCore>(
                    correlationPointOptions->X,
                    correlationPointOptions->Y,
                    correlationPointOptions->distance,
                    correlationPointOptions->minRange,
                    correlationPointOptions->maxRange);
                
                correlationPointCore->init(sys);
                //correlationPointCore->dbg = dbg;
            }


            ///////// setup magnetisation parameters
            for (int i=0; i<magnetisationOptions.size(); ++i){
                magnetisationCores.push_back(
                    MagnetisationCore(magnetisationOptions[i]->axis,(magnetisationOptions[i]->_spins)));
                magnetisationCores[i].init(sys);
            }

            bool swapRes;
            unsigned swapNum;
            double eOld = sys.E();
            double dE,p,randNum;


            double mxOld;
            double myOld;

            for (unsigned step=0; step<hSteps; ++step){
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
            mxOld = sys.M().x;
            myOld = sys.M().y;

            for (auto &cc: correlationCores){
                cc.cpOld = cc.getCPFull(sys);
                if (dbg){
                    cout<<"cpOld: "<<cc.cpOld<<endl;
                }
            }

            if (correlationPointCore){
                correlationPointCore->cpOld = correlationPointCore->getCPFull(sys);
            }


            for (auto &mc: magnetisationCores){
                mc.mOld = mc.getMOFull();
            }

            for (unsigned step=0; step<cSteps; ++step){
                for (unsigned sstep=0; sstep<sys.size(); ++sstep){
                    swapNum = intDistr(generator);
                    sys.parts[swapNum]->rotate(true);
                    dE = eOld - sys.E();
                    p = exp(dE/t);
                    randNum = doubleDistr(generator);
                    if (dE>0 || (t>0 && randNum <= p)){
                        Part* partA = sys.parts[swapNum];
                        eOld = sys.E();
                        mxOld += 2 * (partA->m.x);
                        myOld += 2 * (partA->m.y);

                        for (auto &cc: correlationCores){
                            cc.method(partA);

                            if (dbg){
                                long tmp = cc.getCPFull(sys);
                                if (tmp!=cc.cpOld){
                                    sys.save("tmp.mfsys");
                                    cout<<"CC err: "<<cc.cpOld<<" | "
                                    <<tmp<<endl;
                                }
                            }
                        }

                        if (correlationPointCore){
                            correlationPointCore->method(partA);

                            if (dbg){
                            
                                long tmp = correlationPointCore->getCPFull(sys);
                                if (tmp!=correlationPointCore->cpOld){
                                    sys.save("tmp.mfsys");
                                    cout<<"CP err: "<<correlationPointCore->cpOld<<" | "
                                    <<tmp<<endl;
                                }
                            }
                        }

                        for (auto &mc: magnetisationCores){
                            mc.method(partA->Id());

                            if (dbg){
                                long tmp = mc.getMOFull();
                                if (tmp!=mc.mOld){
                                    sys.save("tmp.mfsys");
                                    cout<<"MO err: "<<mc.mOld<<" | "
                                    <<tmp<<endl;
                                }
                            }
                        }
                    } else {
                        sys.parts[swapNum]->rotate();
                    }
                    e += eOld;
                    mx += mxOld;
                    my += myOld;
                    e2 += eOld*eOld;
                    mx2 += mxOld*mxOld;
                    my2 += myOld*myOld;
                    for (auto &cc: correlationCores){
                        cc.cp += double(cc.cpOld) / cc.correlationPairsNum;
                        cc.cp2 += (double(cc.cpOld*cc.cpOld) / (cc.correlationPairsNum * cc.correlationPairsNum));
                    }

                    if (correlationPointCore){
                        double addVal = correlationPointCore->cpOld / correlationPointCore->pointCount();
                        correlationPointCore->cp += addVal;
                        correlationPointCore->cp2 += addVal*addVal;
                    }

                    for (auto &mc: magnetisationCores){
                        double addVal = double(mc.mOld) / mc.spins.size();
                        mc.mv += addVal;
                        mc.mv2 += addVal*addVal;
                    }
                }
            }

            mx /= cSteps*N;
            mx2 /= cSteps*N;
            my /= cSteps*N;
            my2 /= cSteps*N;
            e /= cSteps*N;
            e2 /= cSteps*N;
            for (auto &cc: correlationCores){
                cc.cp /= cSteps*N;
                cc.cp2 /= cSteps*N;
            }

            if (correlationPointCore){
                correlationPointCore->cp /= cSteps*N;
                correlationPointCore->cp2 /= cSteps*N;
            }

            for (auto &mc: magnetisationCores){
                mc.mv /= cSteps*N;
                mc.mv2 /= cSteps*N;
            }

            mpf_class cT = (e2 - (e * e))/(t * t * sys.size());

            #pragma omp critical
            {
                printf("%e %e %e %e %e %e %e %e %d %d",
                    t, cT.get_d(), e.get_d(), e2.get_d(), 
                    mx.get_d(), mx2.get_d(), my.get_d(), my2.get_d(),
                    omp_get_thread_num(), trseed);
                for (auto &cc: correlationCores) 
                    printf(" %e %e",cc.cp.get_d(),cc.cp2.get_d());
                if (correlationPointCore){
                    printf(" %e %e",
                    correlationPointCore->cp.get_d(),
                    correlationPointCore->cp2.get_d());
                }
                for (auto &mc: magnetisationCores){
                    printf(" %e %e",mc.mv.get_d(),mc.mv2.get_d());
                }
                printf("\n");
                fflush(stdout);
            }
        }
        

    }
    

}