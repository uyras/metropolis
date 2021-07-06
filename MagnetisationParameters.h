#ifndef MAGNETISATIONPARAMETERS_H
#define MAGNETISATIONPARAMETERS_H

#include <vector>
#include <argumentum/argparse.h>

class MagnetisationParameters: public argumentum::CommandOptions
{
    public:
    std::vector<unsigned> _spins;
    Vect axis;

    MagnetisationParameters():CommandOptions("magnetisation"){}
    MagnetisationParameters(std::string_view name):
        CommandOptions(name)
        {}

    bool check(const unsigned systemSize){
        bool res=true;


        for (auto ss : this->_spins){
            if (ss >= systemSize){
                fprintf(stderr,"There is no spin with ID #%u in system of %u spins. ",ss,systemSize);
                res=false;
            }
        }


        if (this->_spins.size()==0){
            for (int i=0; i<systemSize; ++i){
                this->_spins.push_back(i);
            }
        }



        this->axis.setUnitary();

        return res;

    }

    void dump(){
        cout<<"axis: "<<this->axis<<endl;
        cout<<"spins: ("; 
        for (auto i : this->_spins){
            cout<<i<<",";
        }
        cout<<")";
        cout<<endl;
    }


    

    protected:
    void add_parameters(argumentum::ParameterConfig& params) override
    {
        params.add_parameter(_spins,"","--spins").nargs(1)
            .metavar("ID1,ID2,...")
            .help("Set the list of spin IDs for which the magnetisation will be calculated.\
            IDs should be divided by comma. If you need all spins, skip this parameter.\
            By deafult all spins are included.")
            .action([&](auto& target, const std::string& val){
                istringstream f(val);
                string tmp;
                while (getline(f,tmp,',')){
                    target.push_back(stoi(tmp));
                }
            });
        params.add_parameter(axis,"","--axis")
            .nargs(1)
            .metavar("X,Y")
            .absent(Vect(1,0,0))
            .help("Set the normal vector along which the magnetisation will be calculated. \
            The format is X,Y coordinates. \
            For example: 1,0 will measure magnetisation along X axis; \
            1,1 will measure diagonal magnetisation. \
            This vector will be automatically normalised.\
            Aliases X and Y instead of coordinates are also avaliable.\
            Default is X vector.")
            .action([&](auto& target, const std::string& val){
                if (val=="x" || val=="X") 
                    target = Vect(1,0,0); 
                else
                    if (val=="y" || val=="Y") 
                        target = Vect(0,1,0); 
                    else {
                        cerr<<val<<endl;
                        auto pos = val.find(',');
                        if (pos!=std::string::npos){
                            double x = stod(val.substr(0,pos));
                            double y = stod(val.substr(pos+1));
                            target = Vect(x,y,0);
                        } else {
                            throw(std::invalid_argument("format is x,y"));
                        }
                    }
                
            });
    }
};

#endif //MAGNETISATIONPARAMETERS_H