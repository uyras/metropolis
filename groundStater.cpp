#include <iostream>
#include <string>
#include <random>
#include <chrono>
#include "PartArray.h"
#include "Part.h"
#include <argumentum/argparse.h>

using namespace std;
using namespace argumentum;

int main(int argc, char* argv[])
{

    int time_secs; //one minute to run
    int heatupSteps = 10000;
    double t = 0.000001;
    auto parser = argumentum::argument_parser{};
    auto params = parser.params();

    std::string filename,newFilename;
    bool rewrite=false;
    bool shuffle = false;

    parser.config().program("groundStater")
        .description("Program to find the deeper groundState of magnetic system");
    params.add_parameter(filename,"-f","--filename").nargs(1).required().metavar("FILE.mfsys")
        .help("Path to text file with structure of the system. \
            Format is the mfsys file.");
    params.add_parameter(newFilename,"-n","--newfilename").absent("").nargs(1).metavar("FILE.mfsys")
        .help("New file where to save the result. By default it is the old file with '_gs' suffix.");
    params.add_parameter(rewrite,"-r","--rewrite").nargs(0)
        .help("Rewrite old file");
    params.add_parameter(time_secs,"-t","--time").nargs(1).absent(60)
        .help("time to run the code in seconds. Default is 60s - one minute");
    params.add_parameter(shuffle,"-s","--shuffle").nargs(0)
        .help("shuffle system states at start");

    auto res = parser.parse_args( argc, argv, 1 );

    if ( !res )
      return 1;

    if (rewrite){
        newFilename = filename;
    } else {
        if (newFilename==""){
            newFilename = filename;
            std::string::size_type pos = newFilename.rfind('.');
            newFilename.insert(pos,"_gs");
        }
    }

    //cout<<filename<<endl<<newFilename<<endl; //return 0;

    PartArray sys;
    sys.load(filename);

    const int N = sys.size();
    const double eInit = sys.E();

    if (shuffle){
        sys.state.randomize(N);
        cout<<"system shuffled"<<endl;
    }

    bool swapRes;
    unsigned swapNum;
    double eOld = sys.E();
    double dE,p,randNum;

    default_random_engine generator;
    generator.seed(time(NULL));
    uniform_int_distribution<int> intDistr(0, N-1); // including right edge
    uniform_real_distribution<double> doubleDistr(0,1); // right edge is not included

    if (N==0){
        throw(std::string("System size is 0 or file not found"));
    }

    cout<<"Run monte-carlo for "<<time_secs<<" seconds"<<endl;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;
    do {
        for (unsigned sstep=0; sstep<N; ++sstep){
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
        end_time = std::chrono::high_resolution_clock::now();
    } while(std::chrono::duration_cast<std::chrono::seconds>(end_time-start_time).count() < time_secs);

    printf("old E: %f, new E: %f\n", eInit, sys.E());
    if (sys.E()<eInit){
        cout<<"new state is: "<<sys.state.toString()<<endl;
        sys.state.hardReset();
        sys.save(newFilename);
        printf("Saved to %s\n",newFilename.c_str());
    } else {
        printf("Could not find deeper GS. Nothing to save. Exiting.\n");
    }

    return 0;
}