#include "ConfigManager.h"

Vect ConfigManager::size;

bool ConfigManager::check_config()
{

    //check main parameters
    {
        if (this->calculate<1){
            cerr<<"error! --calculate should be greather than 0!"<<endl;
            return false;
        }

    }

    //check parameters of the core
    for (auto & co : parameters){
        if (!co->check(this->system.size())) return false;
    }

    return true;
}

ConfigManager ConfigManager::init(
    const CommandLineParameters & commandLineParameters, 
    const inicpp::config & iniconfig)
{
    ConfigManager tmp;

    if (iniconfig.contains("main")){
        inicpp::section sect = iniconfig["main"];
        if (sect.contains("file")) tmp.sysfile = sect["file"].get<inicpp::string_ini_t>();
        if (sect.contains("heatup")) tmp.heatup = sect["heatup"].get<inicpp::unsigned_ini_t>();
        if (sect.contains("calculate")) tmp.calculate = sect["calculate"].get<inicpp::unsigned_ini_t>();
        if (sect.contains("range")) tmp.range = sect["range"].get<inicpp::float_ini_t>();
        if (sect.contains("seed")) tmp.seed = sect["seed"].get<inicpp::unsigned_ini_t>();
        if (sect.contains("temperature")) tmp.temperatures = sect["temperature"].get_list<inicpp::float_ini_t>();
        if (sect.contains("boundaries") && sect["boundaries"].get<inicpp::string_ini_t>()=="periodic") 
            tmp.pbc = true;
        if (sect.contains("size")) {
            ConfigManager::size =
             ConfigManager::strToVect(sect["size"].get<inicpp::string_ini_t>());
        } else if (tmp.pbc) {
            throw std::invalid_argument("You have to set the \"size\" parameter when using boundaries=periodic");
        }
        
        if (sect.contains("field")) 
            tmp.field = ConfigManager::strToVect(sect["field"].get<inicpp::string_ini_t>());
        else
            tmp.field.setXYZ(0,0,0);

        if (sect.contains("debug")) tmp.debug = sect["debug"].get<inicpp::boolean_ini_t>();

        if (sect.contains("restart")) tmp.restart = sect["restart"].get<inicpp::boolean_ini_t>();
        if (sect.contains("restartThreshold")) tmp.restartThreshold = sect["restartThreshold"].get<inicpp::float_ini_t>();
        if (sect.contains("restartthreshold")) tmp.restartThreshold = sect["restartthreshold"].get<inicpp::float_ini_t>();
        if (sect.contains("saveGS")) tmp.newGSFilename = sect["saveGS"].get<inicpp::string_ini_t>();
        if (sect.contains("savegs")) tmp.newGSFilename = sect["savegs"].get<inicpp::string_ini_t>();
    }
    
    if (!commandLineParameters.sysfilename.empty())
        tmp.sysfile = commandLineParameters.sysfilename;
    if (commandLineParameters.hSteps != -1)
        tmp.heatup = commandLineParameters.hSteps;
    if (commandLineParameters.cSteps != -1)
        tmp.calculate = commandLineParameters.cSteps;
    if (!isnan(commandLineParameters.iRange))
        tmp.range = commandLineParameters.iRange;
    if (commandLineParameters.rseed!=-1)
        tmp.seed = commandLineParameters.rseed;
    if (commandLineParameters.temperatures.size()>0)
        tmp.temperatures = commandLineParameters.temperatures;

    
    tmp.system.load(tmp.sysfile);
    tmp.system.state.hardReset();
    tmp.system.setInteractionRange(tmp.range);
    if (tmp.isPBC()){
        ConfigManager::setPBCEnergies(tmp.system);
    }

    for (auto & sect: iniconfig){
        const std::string parameterString = sect.get_name();

        if (parameterString == "main") continue;

        std::size_t colpos=parameterString.find_first_of(':');
        if (colpos==std::string::npos)
            throw(std::invalid_argument("Section name in ini file should contane colon symbol"));

        const std::string parameterName = parameterString.substr(0,colpos);
        const std::string parameterId = parameterString.substr(colpos+1);

        bool setDebug = false;
        if (sect.contains("debug") && sect["debug"].get<inicpp::boolean_ini_t>()==true)
            setDebug = true;


        if (parameterName == "correlation") {
            if (!sect.contains("method"))
                throw(std::invalid_argument("Parameter " + parameterString + " should have method field."));
            const auto opt = sect["method"];

            std::vector<uint64_t> spins;
            if (sect.contains("spins"))
                spins = sect["spins"].get_list<inicpp::unsigned_ini_t>();

            std::string t_parameterId(parameterId);
            bool is_list = opt.is_list();

            int i=0;
            for (auto & a : opt.get_list<inicpp::string_ini_t>()){
                if (is_list) t_parameterId = parameterId+"_"+to_string(i);
                std::unique_ptr<CorrelationCore> core = 
                    make_unique<CorrelationCore>(t_parameterId,
                        &tmp.system,
                        sect["minrange"].get<inicpp::float_ini_t>(),
                        sect["maxrange"].get<inicpp::float_ini_t>(),
                        methods.at(a),
                        spins);
                
                if (setDebug) core->setDebug();

                tmp.parameters.push_back(std::move(core));
                ++i;
            }


        } else if (parameterName == "correlationpoint") {
                std::vector<double> X,Y;

                for (auto & a : sect["points"].get_list<inicpp::string_ini_t>()){
                    Vect t = ConfigManager::strToVect(a);
                    X.push_back(t.x);
                    Y.push_back(t.y);
                }

                double minRange = 0;
                if (sect.contains("minrange")) 
                    minRange = sect["minrange"].get<inicpp::float_ini_t>();

                std::unique_ptr<CorrelationPointCore> core = 
                    make_unique<CorrelationPointCore>(parameterId,
                        &tmp.system, X, Y,
                        sect["distance"].get<inicpp::float_ini_t>(),
                        minRange,
                        sect["maxrange"].get<inicpp::float_ini_t>());
                //core._methodVar = m;
                if (setDebug) core->setDebug();

                // if need to save the histogram to the file
                if (sect.contains("histogram") && sect["histogram"].get<inicpp::boolean_ini_t>()==true){
                    core->enableHistogram("histogram_"+parameterId+"_$.txt"); //replace dollar by number in time of saving
                }

                tmp.parameters.push_back(std::move(core));

        } else if (parameterName == "magnetisation") {
            int i=0;

            std::vector<uint64_t> spins;
            if (sect.contains("spins"))
                spins = sect["spins"].get_list<inicpp::unsigned_ini_t>();

            bool setModule = false;;
            if (sect.contains("module"))
                setModule = sect["module"].get<inicpp::boolean_ini_t>();


            std::string t_parameterId(parameterId);
            bool is_list = sect["axis"].is_list();

            for (auto & a : sect["axis"].get_list<inicpp::string_ini_t>()){
                if (is_list) t_parameterId = parameterId+"_"+to_string(i);

                std::unique_ptr<MagnetisationCore> core = 
                    make_unique<MagnetisationCore>(t_parameterId,
                        &tmp.system,
                        ConfigManager::strToVect(a),
                        spins);

                if (setDebug) core->setDebug();
                if (setModule) core->setModule(true);

                tmp.parameters.push_back(std::move(core));
                ++i;
            }
        } else if (parameterName == "magnetisationlength") {

            std::vector<uint64_t> spins;
            if (sect.contains("spins"))
                spins = sect["spins"].get_list<inicpp::unsigned_ini_t>();

            std::unique_ptr<MagnetisationLengthCore> core = 
                make_unique<MagnetisationLengthCore>(parameterId,
                    &tmp.system,
                    spins);

            if (setDebug) core->setDebug();

            tmp.parameters.push_back(std::move(core));

        } else {
            throw(std::invalid_argument("Parameter " + parameterName + " is unknown"));
        }
    }

    return tmp;
}

void ConfigManager::printHeader()
{
    //obtain the avaliable number of threads
    #pragma omp parallel
    {
        #pragma omp single
        {
            threadCount = omp_get_num_threads();
        }
    }

    double avgNeighb = 0;
    for (int i=0; i<this->system.size(); ++i)
        avgNeighb += this->system.neighbourSize(i);
    avgNeighb /= this->system.size();

    double e = system.E();
    for (auto p : this->system.parts)
    {
        e -= p->m.scalar(field);
    }

    printf("# Metropolis algorithm for calculating heating capacity v%s\n",METROPOLIS_VERSION);

    printf("#   sysfile: %s\n",this->sysfile.c_str());
    printf("#    system: %d spins, %f interaction range, %f avg. neighbours\n",this->system.size(),this->range, avgNeighb);
    printf("#   physics: energy: %g, ext.filed: (%g,%g,%g), hamiltonian: dipole, space: 2D\n",e,this->field.x,this->field.y,this->field.z);
    printf("#    bounds: ");
    if (this->isPBC()){
        printf("periodic, system size: (%g,%g,%g)\n",ConfigManager::size.x,ConfigManager::size.y,ConfigManager::size.z);
    } else {
        printf("open\n");
    }
    printf("#        MC: %u heatup, %u compute steps\n",this->heatup,this->calculate);
    if (this->isRestart())
        printf("#   restart: enabled, delta E threshold: %g*energy=%g\n",
            this->getRestartThreshold(),
            fabs(this->getRestartThreshold()*e));
    else
        printf("#   restart: disabled\n");
    printf("#   threads: %d\n",threadCount);
    printf("#     rseed: %d+<temperature number>\n",this->seed);
    printf("#    temps.: %zd pcs. from %e to %e\n",
        temperatures.size(),
        std::min_element(temperatures.begin(),temperatures.end()).operator*(),
        std::max_element(temperatures.begin(),temperatures.end()).operator*());
    printf("# temp list: %e",temperatures[0]);
    for (int tt=1; tt<temperatures.size(); ++tt)
        printf(",%e",temperatures[tt]);
    printf("\n");
    printf("#   params.: %zd\n",this->parameters.size());
    printf("#\n");

    unsigned i=0;
    for (auto & co : parameters){
        co->printHeader(i); ++i;
    }

    printf("# legend (column names):\n");
    printf("# 1:T 2:C(T)/N 3:<E> 4:<E^2> 5:threadId 6:seed");
    i=7;
    for (auto & co : parameters){
        printf(" %d:<%s> %d:<%s^2>",i,co->parameterId().c_str(),i+1,co->parameterId().c_str());
        i+=2;
    }
    printf(" %d:time,s",i);
    printf("\n");
    fflush(stdout);
}

void ConfigManager::applyState(string s)
{ 
    this->system.state.fromString(s);
    this->system.state.hardReset();
    this->system.setInteractionRange(this->range);
    if (this->isPBC()){
        ConfigManager::setPBCEnergies(this->system);
    } 
}

void ConfigManager::getParameters(std::vector< std::unique_ptr< CalculationParameter > > & calculationParameters)
{
    for (auto & co : parameters){
        calculationParameters.push_back(std::unique_ptr<CalculationParameter>(co->copy()));
    }
    return;
}

void ConfigManager::setPBCEnergies(PartArray & sys)
{
    // first update all neighbours
    sys.neighbours.clear();

    //определяем соседей частицы
    if (sys.interactionRange() != 0.){ //только если не все со всеми
        sys.neighbours.resize(sys.size());
        Part *part, *temp;
        for (unsigned i=0; i<sys.size(); i++){
            sys.neighbours[i].clear();
            part = sys[i];
            vector<Part*>::iterator iter = sys.parts.begin();
            while(iter!=sys.parts.end()){
                temp = *iter;
                if (temp != part && radiusPBC(part->pos,temp->pos).length() < sys.interactionRange()){
                    sys.neighbours[i].push_front(temp);
                }
                iter++;
            }
        }
    }
    sys.changeSystem();

    //then set the hamiltonian
    sys.setHamiltonian(hamiltonianDipolarPBC);
}



double hamiltonianDipolarPBC(Part *a, Part *b)
{
    Vect rij = radiusPBC(b->pos,a->pos);
    double r2, r, r5,E;
    r2 = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
    r = sqrt(r2); //трудное место, заменить бы
    r5 = r2 * r2 * r; //радиус в пятой
    
    E = //энергия считается векторным методом, так как она не нужна для каждой оси
            (( (a->m.x * b->m.x + a->m.y * b->m.y + a->m.z * b->m.z) * r2)
                -
                (3 * (b->m.x * rij.x + b->m.y * rij.y + b->m.z * rij.z) * (a->m.x * rij.x + a->m.y * rij.y + a->m.z * rij.z)  )) / r5;
    return E;
}

Vect radiusPBC(const Vect& a, const Vect& b){
    Vect dist = a - b;
    if (fabs(dist.x) > ConfigManager::size.x/2){
        if (a.x < b.x){
            dist.x += ConfigManager::size.x;
        } else {
            dist.x -= ConfigManager::size.x;
        }
    }

    if (fabs(dist.y) > ConfigManager::size.y/2){
        if (a.y < b.y){
            dist.y += ConfigManager::size.y;
        } else {
            dist.y -= ConfigManager::size.y;
        }
    }

    if (fabs(dist.z) > ConfigManager::size.z/2){
        if (a.z < b.z){
            dist.z += ConfigManager::size.z;
        } else {
            dist.z -= ConfigManager::size.z;
        }
    }
    
    return dist;
}