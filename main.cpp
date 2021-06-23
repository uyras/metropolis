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
#include "PartArray.h"
#include "Part.h"
#include <argumentum/argparse.h>
#include "CorrelationCore.h"

#define VERSION "0.0.3"

std::string filename;
std::vector<double> temperatures; // system temperature
unsigned hSteps;        // heatup steps
unsigned cSteps;        // calculate steps
unsigned rseed;         // random seed value
double iRange;          // interaction range
bool dbg=false;
unsigned N;

class correlationParameters: public argumentum::CommandOptions
{
    public:
    std::vector<double> minRange;
    std::vector<double> maxRange;
    std::vector<unsigned> methodVar;

    correlationParameters():CommandOptions("correlation"){}
    correlationParameters(std::string_view name):
        CommandOptions(name)
        {}


    

    protected:
    void add_parameters(argumentum::ParameterConfig& params) override
    {
        params.add_parameter(minRange,"","--minRange").minargs(1).required().metavar("RANGE")
            .help("minimal distance where to calculate the corelation parameter");
        params.add_parameter(maxRange,"","--maxRange").minargs(1).required().metavar("RANGE")
            .help("maximal distance where to calculate the corelation parameter");
        params.add_parameter(methodVar,"","--method")
            .choices({"xor","eSign","scalar"})
            .metavar("xor|eSign|scalar")
            .minargs(1)
            .required()
            .help("set the method of how to calculate the correlation. \
            Possible variants are \"xor\" or \"eSign\" or \"scalar\":\n\
            xor - the initial direction of M for each spin is multiplied by Si={-1;+1}.\
            xor is the production of Si*Sj. At the start of the program all values are +1.\n\
            eSign - here we use Si=H/|H|, where H is hamiltonian. \
            This function is the production of Si*Sj.\n\
            scalar - just a scalar production of two unit vectors.")
            .action([&](auto& target, const std::string& val){
                if (val=="xor") target.push_back(1);
                if (val=="eSign") target.push_back(2);
                if (val=="scalar") target.push_back(3);
            });
    }
};

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

    params.add_command<correlationParameters>("correlation").help("Gather also correlation parameter. \
        This option has additional parameters. For information print: metropolis correlation --help.");

    auto res = parser.parse_args( argc, argv, 1 );

    if ( !res )
      return 1;

    auto correlationOptions = static_pointer_cast<correlationParameters>( res.findCommand("correlation") );
    
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
        if (correlationOptions->maxRange.size() != correlationOptions->minRange.size() ||
            correlationOptions->maxRange.size() != correlationOptions->methodVar.size()){
                cerr<<"Parameters --minRange, --maxRange and --method should \
                have the same values count."<<endl;
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
                correlationOptions->methodVar[i]);
            cctemp.init(bigsys);
            if (cctemp.correlationPairsNum==0){
                cerr<<"Check the correlation range #"<<i<<". \
                    Could not find any spin pairs within this distance!"<<endl;
                return 1;
            }
            printf("# core #%03d: %.2f min range, %.2f max range, %.2f avg. neigh",i,
                correlationOptions->minRange[i], correlationOptions->maxRange[i],
                cctemp.correlationPairsNum/double(N));
            if (correlationOptions->methodVar[i]==1)
                printf(" method: XOR\n");
            if (correlationOptions->methodVar[i]==2)
                printf(" method: E/|E|\n");
            if (correlationOptions->methodVar[i]==3)
                printf(" method: scalar\n");
        }
    }

    printf("# T C(T)/N <E> <E^2> <mx> <mx^2> <my> <my^2> threadId seed");
    if (correlationOptions)
        for (int i=0; i<correlationOptions->minRange.size(); ++i) 
            printf(" <cp#%03d> <cp^2#%03d>",i,i);
    printf("\n");
    fflush(stdout);



    #pragma omp parallel
    {        
        #pragma omp for
        for (int tt=0; tt<temperatures.size();  ++tt){

            std::vector<CorrelationCore> correlationCores;

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
                                                                correlationOptions->methodVar[i]));
                    correlationCores[i].init(sys);
                    correlationCores[i].dbg = dbg;
                }
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

            for (unsigned step=0; step<cSteps; ++step){
                for (unsigned sstep=0; sstep<sys.size(); ++sstep){
                    swapNum = intDistr(generator);
                    sys.parts[swapNum]->rotate(true);
                    dE = eOld - sys.E();
                    p = exp(dE/t);
                    randNum = doubleDistr(generator);
                    if (dE>0 || (t>0 && randNum <= p)){
                        eOld = sys.E();
                        mxOld += 2 * (sys.parts[swapNum]->m.x);
                        myOld += 2 * (sys.parts[swapNum]->m.y);

                        Part* partA = sys.parts[swapNum];
                        for (auto &cc: correlationCores){
                            for (auto partB:cc.correlationNeighbours[swapNum]){
                                cc.cpOld += (cc.method(partA,partB));
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
            mpf_class cT = (e2 - (e * e))/(t * t * sys.size());

            #pragma omp critical
            {
                printf("%e %e %e %e %e %e %e %e %d %d",
                    t, cT.get_d(), e.get_d(), e2.get_d(), 
                    mx.get_d(), mx2.get_d(), my.get_d(), my2.get_d(),
                    omp_get_thread_num(), trseed);
                for (auto &cc: correlationCores) 
                    printf(" %e %e",cc.cp.get_d(),cc.cp2.get_d());
                printf("\n");
                fflush(stdout);
            }
        }
        

    }
    

}