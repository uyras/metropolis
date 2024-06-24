#ifndef COMMANDLINEPARAMETERS_H
#define COMMANDLINEPARAMETERS_H

#include <argumentum/argparse.h>

class CommandLineParameters : public argumentum::Options
{
public:
    std::string inifilename;
    std::string sysfilename;          //mfsys
    std::vector<double> temperatures; // system temperature
    int hSteps;                  // heatup steps
    int cSteps;                  // calculate steps
    int rseed;                   // random seed value
    double iRange;                    // interaction range
    bool showExample;
    int saveStates;

protected:
    void add_parameters(argumentum::ParameterConfig &params) override
    {

        params.add_parameter(inifilename,"i","inifile").nargs(1).absent("").metavar("FILE.ini")
            .help("Path to the configuration file. \
                It is not required if you set the basic parameters via command line.");
        params.add_parameter(sysfilename,"-f","--filename").nargs(1).metavar("FILE.mfsys")
            .help("Path to text file with structure of the system. \
                Format is the mfsys file.");
        params.add_parameter(hSteps,"-p","--heatup").nargs(1).absent(-1).metavar("STEPS")
            .help("Number of steps needed to bring the configuration to \
                the stable state. The single step here consists of several MC sweeps (trials)");
        params.add_parameter(cSteps,"-c","--calculate").nargs(1).absent(-1).metavar("STEPS")
            .help("Number of steps performed to Metropolis algorithm.\
                The single step here consists of several MC sweeps (trials)");
        params.add_parameter(iRange,"-r","--range").absent(NAN).nargs(1).metavar("RANGE")
            .help("Maximal distance between particles which interaction counts in total energy. \
                Partices located further are considered non-interacting.\
                Default value is 0 which means all-to-all interaction.");
        params.add_parameter(rseed,"-s","--seed").absent(-1).nargs(1).metavar("SEED")
            .help("Random seed number. In case of parallel execution\
            it will be SEED+TEMPNUM, where TEMPNUM is the sequential number \
            of temperature in list. Default is 0.");
        params.add_parameter(temperatures,"-t","--temperature").minargs(1)
            .help("Temperature of the system, in units D. \
                May be the space-separated list of temperatures, \
                in this case it will be created several threads for each temperature.\
                Another variant is to save temperatures to the text file, \
                one temperature per line and use this as: metropolis @param.txt");
        params.add_parameter(saveStates,"","--save").maxargs(1).absent(0).metavar("N")
            .help("Enable saving system configurations every N-th MC step. \
                Dont work on heatup phase. Saves only on calculate. \
                Default value is 0 means do not save the data.\
                File name is \"<input filename>_<temperature_number>_<step>.mfsys\"");
        params.add_parameter(showExample,"-e","--example")
            .help("Print out the example of ini-file and exit.");
        params.add_default_help_option();
    }
};

#endif //COMMANDLINEPARAMETERS_H