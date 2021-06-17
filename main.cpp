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
#include "argumentum/include/argumentum/argparse.h"

#define VERSION "0.0.2"

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
    double minRange;
    double maxRange;
    unsigned methodVar;

    correlationParameters():CommandOptions("correlation"){}
    correlationParameters(std::string_view name):
        CommandOptions(name)
        {}


    std::vector< std::forward_list < Part* > > correlationNeighbours;
    std::map< std::pair<unsigned, unsigned>, short > correlationValues;
    unsigned correlationPairsNum;
    /**
     * @brief Сaches the neighbours and energies for further calculations
     * 
     * @param sys system
     */
    void init(const PartArray & sys){
        correlationNeighbours.resize(sys.size());
        correlationPairsNum=0;

        double space2;
        const double minr = this->minRange * this->minRange;
        const double maxr = this->maxRange * this->maxRange;
        double eTemp;
        for (auto partA: sys.parts){
            for (auto partB: sys.parts){
                if (partA==partB)
                    continue;
                
                space2 = partA->pos.space_2(partB->pos);
                if (space2>=minr && space2<=maxr){
                    correlationNeighbours[partA->Id()].push_front(partB);
                    ++correlationPairsNum;

                    if (this->methodVar==2 || this->methodVar==3){
                        if (this->methodVar==2){
                            //В матрицу надо помещать энергии только в неперевернутых состояниях
                            eTemp = hamiltonianDipolar(partA,partB)*-1; 
                            if (partA->state!=partB->state)
                                eTemp*=-1.;
                        }

                        if (this->methodVar==3){
                            eTemp = partA->m.scalar(partB->m); 
                            if (partA->state!=partB->state)
                                eTemp*=-1.;
                        }


                        if (eTemp>0)
                            correlationValues[std::make_pair(partA->Id(),partB->Id())] = 1;
                        else
                            correlationValues[std::make_pair(partA->Id(),partB->Id())] = -1;
                    }
                }
            }
        }
        correlationPairsNum/=2;
    }


    char method(const Part* partA,const Part* partB){
        if (this->methodVar==1) return fxor(partA,partB);
        if (this->methodVar==2) return feSign(partA,partB);
        if (this->methodVar==3) return fScalar(partA,partB);
        return 0;
    }
    // functions for correlation options
    char fxor(const Part* partA,const Part* partB){
        return (partA->state ^ partB->state)?-1:+1;
    }
    char feSign(const Part* partA, const Part* partB){
        return fxor(partA,partB) * 
            this->correlationValues[std::make_pair(partA->Id(),partB->Id())];
    }
    char fScalar(const Part* partA, const Part* partB){
        return fxor(partA,partB) * 
            this->correlationValues[std::make_pair(partA->Id(),partB->Id())];
    }

    protected:
    void add_parameters(argumentum::ParameterConfig& params) override
    {
        params.add_parameter(minRange,"","--minRange").nargs(1).required().metavar("RANGE")
            .help("minimal distance where to calculate the corelation parameter");
        params.add_parameter(maxRange,"","--maxRange").nargs(1).required().metavar("RANGE")
            .help("maximal distance where to calculate the corelation parameter");
        params.add_parameter(methodVar,"","--method")
            .choices({"xor","eSign","scalar"})
            .metavar("xor|eSign|scalar")
            .nargs(1)
            .help("set the method of how to calculate the correlation. \
            Possible variants are \"xor\" or \"eSign\" or \"scalar\":\n\
            xor - the initial direction of M for each spin is multiplied by Si={-1;+1}.\
            xor is the production of Si*Sj. At the start of the program all values are +1.\n\
            eSign - here we use Si=H/|H|, where H is hamiltonian. \
            This function is the production of Si*Sj.\n\
            scalar - just a scalar production of two unit vectors\
            Default is eSign.")
            .absent(2)
            .action([&](auto& target, const std::string& val){
                if (val=="xor") target = 1;
                if (val=="eSign") target = 2;
                if (val=="scalar") target = 3;
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


    PartArray bigsys;
    bigsys.setInteractionRange(iRange);
    bigsys.load(filename.c_str());
    N=bigsys.size();
    if (dbg){
        cerr<<"# file '"<<filename<<"': imported "<<bigsys.size()<<" parts"<<endl;
        cerr<<"# Energy: "<<bigsys.E()<<" range: "<<iRange<<endl;
    }



    ///////// define neighbours for correlation parameter
    if (correlationOptions){
        correlationOptions->init(bigsys);
        if (correlationOptions->correlationPairsNum==0){
            cerr<<"Check the correlation range. \
                Could not find any spin pairs within this distance!"<<endl;
            return 1;
        }
    }




    #pragma omp parallel
    {
        #pragma omp single
        {
            printf("# Metropolis algorithm for calculating heating capacity v%s\n",VERSION);
            printf("#  filename: %s\n",filename.c_str());
            printf("#    System: %d spins, %f interaction range\n",N,iRange);
            printf("#        MC: %u heatup, %u compute steps\n",hSteps,cSteps);
            printf("#   threads: %d\n",omp_get_num_threads());

            if (correlationOptions){
                printf("#   correl.: %f min range, %f max range, %f avg. neighbours\n",
                    correlationOptions->minRange, correlationOptions->maxRange, 
                    double(correlationOptions->correlationPairsNum)/N*2);
                if (correlationOptions->methodVar==1)
                    printf("#    method: XOR\n");
                if (correlationOptions->methodVar==2)
                    printf("#    method: E/|E|\n");
                if (correlationOptions->methodVar==3)
                    printf("#    method: scalar\n");
            }

            printf("# T C(T)/N <E> <E^2> <mx> <mx^2> <my> <my^2> threadId seed");
            if (correlationOptions) printf(" <cp> <cp^2>");
            printf("\n");
            fflush(stdout);
        }

        #pragma omp for
        for (int tt=0; tt<temperatures.size();  ++tt){
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
            mpf_class cp(0,2048);
            mpf_class cp2(0,2048);

            PartArray sys = bigsys;

            //sys.state.reset();

            bool swapRes;
            unsigned swapNum;
            double eOld = sys.E();
            double dE,p,randNum;


            double mxOld;
            double myOld;
            long cpOld;

            for (unsigned step=0; step<hSteps; ++step){
                for (unsigned sstep=0; sstep<sys.size(); ++sstep){
                    swapNum = intDistr(generator);
                    sys.parts[swapNum]->rotate(true);
                    dE = eOld - sys.E();
                    p = exp(dE/t);
                    randNum = doubleDistr(generator);
                    if (dE>0 || (t>0 && randNum <= p)){
                        //do nothing
                    } else {
                        sys.parts[swapNum]->rotate(true);
                    }
                }
            }


            mxOld = sys.M().x;
            myOld = sys.M().y;
            cpOld=0;
            if (correlationOptions){
                for (auto partA:sys.parts){
                    if (dbg) fprintf(stderr,"# neigh for %u: ",partA->Id());
                    for (auto partB:correlationOptions->correlationNeighbours[partA->Id()]){
                        cpOld += correlationOptions->method(partA,partB);
                        if (dbg) fprintf(stderr,"%u,",partB->Id());
                    }
                    if (dbg) fprintf(stderr,"\n");
                }
                cpOld /= 2;
            }
            if (dbg){
                cout<<"cpOld: "<<cpOld<<endl;
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
                        if (correlationOptions){
                            Part* partA = sys.parts[swapNum];
                            for (auto partB:correlationOptions->correlationNeighbours[swapNum]){
                                cpOld += (correlationOptions->method(partA,partB));
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
                    if (correlationOptions){
                        cp += double(cpOld) / correlationOptions->correlationPairsNum;
                        cp2 += (double(cpOld*cpOld) / (correlationOptions->correlationPairsNum*correlationOptions->correlationPairsNum));
                    }
                }
            }

            mx /= cSteps*N;
            mx2 /= cSteps*N;
            my /= cSteps*N;
            my2 /= cSteps*N;
            e /= cSteps*N;
            e2 /= cSteps*N;
            cp /= cSteps*N;
            cp2 /= cSteps*N;
            mpf_class cT = (e2 - (e * e))/(t * t * sys.size());

            #pragma omp critical
            {
                printf("%e %e %e %e %e %e %e %e %d %d",
                    t, cT.get_d(), e.get_d(), e2.get_d(), 
                    mx.get_d(), mx2.get_d(), my.get_d(), my2.get_d(),
                    omp_get_thread_num(), trseed);
                if (correlationOptions) printf(" %e %e",cp.get_d(),cp2.get_d());
                printf("\n");
                fflush(stdout);
            }
        }

    }

}