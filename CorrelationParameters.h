#ifndef CORELLATIONPARAMETERS_H
#define CORELLATIONPARAMETERS_H

#include <argumentum/argparse.h>

class CorrelationParameters: public argumentum::CommandOptions
{
    public:
    std::vector<double> minRange;
    std::vector<double> maxRange;
    std::vector<unsigned> methodVar;
    std::vector< std::vector <unsigned> > spins;

    CorrelationParameters():CommandOptions("correlation"){}
    CorrelationParameters(std::string_view name):
        CommandOptions(name)
        {}

    bool check(const unsigned systemSize){
        bool res=true;

        //check equality of minrange and maxrange
        if (this->minRange.size() != this->maxRange.size()){
            cerr<<"Parameters --minRange and --maxRange should have the same count of values."<<endl;
            res=false;
        }

        //normalize the methodVar
        if (this->methodVar.size()==1 && this->maxRange.size()>1){
            methodVar.resize(this->maxRange.size(),this->methodVar.back());
        }
        //check equality of methodVar and maxrange
        if (this->methodVar.size() != this->maxRange.size()){
            cerr<<"Parameter --method should have one value or the same count of values as --maxRange."<<endl;
            res=false;
        }

        //normalize spin numbers
        this->spins.push_back({});
        if (this->spins.size()==1 && this->maxRange.size()>1){
            spins.resize(this->maxRange.size(),this->spins.back());
        }
        //check equality of spin count and maxrange
        if (this->spins.size() != this->maxRange.size()){
            cerr<<"Parameter --spins should have one value or the same count of values as --maxRange."<<endl;
            res=false;
        } else {
            //check that all spins are in range
            for (int i=0;i<this->spins.size(); ++i){
                for (auto ss : this->spins[i]){
                    if (ss >= systemSize){
                        fprintf(stderr,"There is no spin with ID #%u in system of %u spins. ",ss,systemSize);
                        fprintf(stderr,"You setted this in --spins parameter for correlation #%d\n",i);
                        res=false;
                    }
                }
            }
        }

        return res;

    }

    void dump(){
        cout<<"minRange: "; for (auto i : this->minRange) cout<<i<<","; cout<<endl;
        cout<<"maxRange: "; for (auto i : this->maxRange) cout<<i<<","; cout<<endl;
        cout<<"methodVar: "; for (auto i : this->methodVar) cout<<i<<","; cout<<endl;
        cout<<"spins: "; 
        for (auto i : this->spins){
            cout<<"(";
            for (auto j:i) cout<<j<<",";
            cout<<"), ";
        }
        cout<<endl;
    }


    

    protected:
    void add_parameters(argumentum::ParameterConfig& params) override
    {
        params.add_parameter(minRange,"","--minRange").minargs(1).required().metavar("RANGE")
            .help("minimal distance where to calculate the corelation parameter");
        params.add_parameter(maxRange,"","--maxRange").minargs(1).required().metavar("RANGE")
            .help("maximal distance where to calculate the corelation parameter");
        /*params.add_parameter(spins,"","--spins").minargs(1).metavar("ID")
            .help("Set the list of spin IDs for which the selected parameter wil be calculated.\
            IDs should be divided by comma. If you need all spins, write the dot symbol.\
            The amount of this arguments should be the same as --maxRange, \
            or specify just one value and it will be applied to all correlation cores.\
            Default value is \".\" (all spins).")
            .absent(std::vector< std::vector < unsigned> > ({{}}))
            .action([&](auto& target, const std::string& val){
                target.push_back({});
                if (val!="."){
                    istringstream f(val);
                    string tmp;
                    while (getline(f,tmp,',')){
                        target.back().push_back(stoi(tmp));
                    }
                }
            });*/
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
                if (val=="m") target.push_back(4);
            });
    }
};

#endif //CORELLATIONPARAMETERS_H