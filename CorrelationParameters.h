#ifndef CORELLATIONPARAMETERS_H
#define CORELLATIONPARAMETERS_H

#include <argumentum/argparse.h>

class CorrelationParameters: public argumentum::CommandOptions
{
    public:
    std::vector<double> minRange;
    std::vector<double> maxRange;
    std::vector<unsigned> methodVar;

    CorrelationParameters():CommandOptions("correlation"){}
    CorrelationParameters(std::string_view name):
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

#endif //CORELLATIONPARAMETERS_H