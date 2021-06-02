#include <iostream>
#include <sstream>
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
#include "cxxopts.hpp"
#include "argumentum/include/argumentum/argparse.h"

#define VERSION "0.0.1"

std::string filename;
std::vector<double> temperatures; // system temperature
unsigned hSteps;        // heatup steps
unsigned cSteps;        // calculate steps
double iRange;          // interaction range
bool dbg=true;
unsigned N;


int main(int argc, char* argv[])
{

    bool parse_failed=false;

    /** Parse Options */
    cxxopts::Options options("metropolis","Program for calculating heat capacity and magnetisation of spin system with dipole-dipole hamiltonian");
    options.add_options()
        ("f,filename", 
            "Path to text file with structure of the system. \
            Format is the mfsys file.", 
            cxxopts::value<std::string>())
        ("t,temperature", "Temperature of the system, in units D. \
            May be the comma-separated list of temperatures, \
            in this case it will be created several threads for each temperature",
            cxxopts::value<std::vector<double>>())
        ("s,heatup","Number of steps needed to bring the configuration to \
            the stable state. The single step here consists of several MC sweeps (trials)",
            cxxopts::value<unsigned>())
        ("c,calculate",
            "Number of steps performed to Metropolis algorithm.",
            cxxopts::value<unsigned>())
        ("r,range",
            "Maximal distance between particles which interaction counts in total energy. \
            Partices located further are considered non-interacting.\
            Default value is 0 which means all-to-all interaction.",
            cxxopts::value<double>()->default_value("0"))
        ("h,help", "Print usage");
    try{
        auto result = options.parse(argc, argv);
    
        if (result.count("help"))
        {
            parse_failed=true;
        }
    
        filename = result["filename"].as<std::string>();
        temperatures = result["temperature"].as<std::vector<double>>();
        hSteps = result["heatup"].as<unsigned>();
        cSteps = result["calculate"].as<unsigned>();
        iRange = result["range"].as<double>();
    } catch(std::domain_error e){
        std::cout<<"Error parsing command-line arguments.\nParameters f,t,s,c are required."<<std::endl;
        parse_failed=true;
    } catch(cxxopts::OptionException e){
        std::cout<<"Error! "<<e.what()<<std::endl;
        parse_failed=true;
    }

    if (parse_failed){
        std::cout << options.help() << std::endl;
        exit(0);
    }

    ifstream f(filename);
    if (!f.is_open()) {
        cout<<"file read error: "<<filename<<endl;
        exit(0);
    }
    f.close();

    {
        PartArray sys;
        sys.load(filename.c_str());
        N=sys.size();
        if (dbg){
            cerr<<"# file '"<<filename<<"': imported "<<sys.size()<<" parts"<<endl;
            cerr<<"# Energy: "<<sys.E()<<" range: "<<iRange<<endl;
        }
    }
    #pragma omp parallel
    {
        #pragma omp single
        {
            printf("# Metropolis algorithm for calculating heating capacity v%s\n",VERSION);
            printf("# filename: %s\n",filename.c_str());
            printf("#   System: %d spins, %f interaction range\n",N,iRange);
            printf("#       MC: %u heatup, %u compute steps\n",hSteps,cSteps);
            printf("#  threads: %d\n",omp_get_num_threads());
            printf("# T C(T)/N <E> <E^2> <mx> <mx^2> <my> <my^2> threadId\n");
        }
    }

    #pragma omp parallel for
    for (int tt=0; tt<temperatures.size();  ++tt){
        const double t = temperatures[tt]; 
        /////////// import the system
        mpf_class mx(0,1024);
        mpf_class my(0,1024);
        mpf_class e(0,1024);
        mpf_class mx2(0,2048);
        mpf_class my2(0,2048);
        mpf_class e2(0,2048);

        PartArray sys;
        sys.setInteractionRange(iRange);
        sys.load(filename.c_str());

        //sys.state.reset();

        bool swapRes;
        unsigned swapNum;
        double eOld = sys.E();
        double mxOld = sys.M().x;
        double myOld = sys.M().y;
        double dE,p,randNum;

        for (unsigned step=0; step<hSteps; ++step){
            for (unsigned sstep=0; sstep<sys.size(); ++sstep){
                swapNum = sys.state.randomize();
                dE = eOld - sys.E();
                p = exp(dE/t);
                randNum = (double)config::Instance()->rand()/(double)config::Instance()->rand_max;
                if (dE>0 || (t>0 && randNum <= p)){
                    //do nothing
                } else {
                    sys.parts[swapNum]->rotate();
                }
            }
        }

        for (unsigned step=0; step<cSteps; ++step){
            for (unsigned sstep=0; sstep<sys.size(); ++sstep){
                swapNum = sys.state.randomize();
                dE = eOld - sys.E();
                p = exp(dE/t);
                randNum = (double)config::Instance()->rand()/(double)config::Instance()->rand_max;
                if (dE>0 || (t>0 && randNum <= p)){
                    eOld = sys.E();
                    mxOld += 2 * (sys.parts[swapNum]->m.x);
                    myOld += 2 * (sys.parts[swapNum]->m.y);
                } else {
                    sys.parts[swapNum]->rotate();
                }
                e += eOld;
                mx += mxOld;
                my += myOld;
                e2 += eOld*eOld;
                mx2 += mxOld*mxOld;
                my2 += myOld*myOld;
            }
        }

        mx /= cSteps*N;
        mx2 /= cSteps*N;
        my /= cSteps*N;
        my2 /= cSteps*N;
        e /= cSteps*N;
        e2 /= cSteps*N;
        mpf_class cT = (e2 - (e * e))/(t * t * sys.size());

        #pragma omp critical
        {
            printf("%f %e %e %e %e %e %e %e %d\n",
                t, cT.get_d(), e.get_d(), e2.get_d(), 
                mx.get_d(), mx2.get_d(), my.get_d(), my2.get_d(),
                omp_get_thread_num());
            fflush(stdout);
        }
    }

}