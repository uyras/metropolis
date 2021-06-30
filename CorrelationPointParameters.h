#ifndef CORELLATIONPOINTPARAMETERS_H
#define CORELLATIONPOINTPARAMETERS_H

#include <argumentum/argparse.h>

class CorrelationPointParameters: public argumentum::CommandOptions
{
    public:
    std::vector<double> X;
    std::vector<double> Y;
    double minRange;
    double maxRange;
    double distance;

    CorrelationPointParameters():CommandOptions("corrpoint"){}
    CorrelationPointParameters(std::string_view name):
        CommandOptions(name)
        {}


    

    protected:
    void add_parameters(argumentum::ParameterConfig& params) override
    {
        params.add_parameter(X,"-x").minargs(1).required().metavar("POS")
            .help("set of X points where to find spins for this parameters. \
            Should me the same count of points as in -y.");
        params.add_parameter(Y,"-y").minargs(1).required().metavar("POS")
            .help("set of Y points where to find spins for this parameters. \
            Should me the same count of points as in -x.");
        params.add_parameter(distance,"-d","--distance").nargs(1).required().metavar("POS")
            .help("the distance around points where to find spins linked to point");
        params.add_parameter(minRange,"","--minRange").nargs(1).absent(0.).metavar("RANGE")
            .help("minimal distance where to calculate the corelation parameter, default is 0.");
        params.add_parameter(maxRange,"","--maxRange").nargs(1).required().metavar("RANGE")
            .help("maximal distance where to calculate the corelation parameter");
        /*params.add_parameter(methodVar,"","--method")
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
            });*/
    }
};

#endif //CORELLATIONPOINTPARAMETERS_H