#include "ConfigManager.h"

ConfigManager::ConfigManager(int argc, char *argv[])
{
    //--------------- старая функция readParameters
    // get file name
	bool parse_failed = false;

	auto parser = argumentum::argument_parser{};
	parser.config().program("metropolis").description("Program for calculating heat capacity \
            and magnetisation of spin system with dipole-dipole hamiltonian v." +
													  std::string(METROPOLIS_VERSION));
	auto commandLineParameters = std::make_shared<CommandLineParameters>();
	parser.params().add_parameters(commandLineParameters);

	auto parseResult = parser.parse_args(argc, argv, 1);

	if (!parseResult)
	{
		if (commandLineParameters && commandLineParameters->showExample)
		{
			std::cout << endl;
			std::cout << "##########################################" << endl;
			std::cout << "######## contents of example.ini: ########" << endl;
			std::cout << "##########################################" << endl;
			std::cout << endl;
			std::cout << example_string << endl;
		}
		throw string("# There is nothing to calculate. Stopping.");
	}

	inicpp::config iniconfig;
	if (!commandLineParameters->inifilename.empty())
	{
		iniconfig = inicpp::parser::load_file(commandLineParameters->inifilename);
	}
    //--------------- старая функция readParameters, конец

    vector<double> base_temperatures;

    Vect size = {0,0,0};
    double range = 0;
    hamiltonian_t ham = hamiltonian_selector("dipolar");

    if (iniconfig.contains("main")){
        inicpp::section sect = iniconfig["main"];
        if (sect.contains("boundaries")){
            cerr<<"# Warning! The option 'main/boundaries' is no longer in use. "<<endl;
            cerr<<"#          Set non-zero value of 'main/size' instead to enable PBC long desired axis."<<endl;
            cerr<<"#          Otherwise the OBC will be used."<<endl;
        }
        if (sect.contains("file")) this->sysfile = sect["file"].get<inicpp::string_ini_t>();
        if (sect.contains("heatup")) this->heatup = sect["heatup"].get<inicpp::unsigned_ini_t>();
        if (sect.contains("calculate")) this->calculate = sect["calculate"].get<inicpp::unsigned_ini_t>();
        if (sect.contains("range")) range = sect["range"].get<inicpp::float_ini_t>();
        if (sect.contains("seed")) this->seed = sect["seed"].get<inicpp::unsigned_ini_t>();
        if (sect.contains("temperature")) base_temperatures = sect["temperature"].get_list<inicpp::float_ini_t>();
        if (sect.contains("size")) size = strToVect(sect["size"].get<inicpp::string_ini_t>());
        
        if (sect.contains("field")) 
            this->field = strToVect(sect["field"].get<inicpp::string_ini_t>());
        else
            this->field = {0,0,0};

        if (sect.contains("debug")) this->debug = sect["debug"].get<inicpp::boolean_ini_t>();

        if (sect.contains("restart")) this->restart = sect["restart"].get<inicpp::boolean_ini_t>();
        if (sect.contains("restartThreshold")) this->restartThreshold = sect["restartThreshold"].get<inicpp::float_ini_t>();
        if (sect.contains("restartthreshold")) this->restartThreshold = sect["restartthreshold"].get<inicpp::float_ini_t>();
        if (sect.contains("saveGS")) this->newGSFilename = sect["saveGS"].get<inicpp::string_ini_t>();
        if (sect.contains("savegs")) this->newGSFilename = sect["savegs"].get<inicpp::string_ini_t>();
    }
    
    if (!commandLineParameters->sysfilename.empty())
        this->sysfile = commandLineParameters->sysfilename;
    if (commandLineParameters->hSteps != -1)
        this->heatup = commandLineParameters->hSteps;
    if (commandLineParameters->cSteps != -1)
        this->calculate = commandLineParameters->cSteps;
    if (!isnan(commandLineParameters->iRange))
        range = commandLineParameters->iRange;
    if (commandLineParameters->rseed!=-1)
        this->seed = commandLineParameters->rseed;
    if (commandLineParameters->temperatures.size()>0)
        base_temperatures = commandLineParameters->temperatures;

    if (iniconfig.contains("parallel_tempering")){
        this->temperatures = PtParseConfig(iniconfig);
        this->temperatures->setBaseTemperatures(base_temperatures);
    } else {
        this->temperatures = new PtBalancerDefault(base_temperatures);
    }
    this->temperatures->init();

    this->system = make_shared<MagneticSystem>(this->sysfile,range,size,ham);
    
    this->saveStates = commandLineParameters->saveStates;
    this->saveStateFileBasename = this->sysfile.substr(0, this->sysfile.find_last_of("."));
/*
    for (auto & sect: iniconfig){
        const std::string parameterString = sect.get_name();

        // все вычисляемые параметры должны иметь в имени секции символ двоеточие
        std::size_t colpos=parameterString.find_first_of(':');
        if (colpos==std::string::npos) continue;

        const std::string parameterName = parameterString.substr(0,colpos);
        const std::string parameterId = parameterString.substr(colpos+1);

        bool setDebug = false;
        if (sect.contains("debug") && sect["debug"].get<inicpp::boolean_ini_t>()==true)
            setDebug = true;


        if (parameterName == "correlation") {
            if (!sect.contains("method"))
                throw(std::string("Parameter " + parameterString + " should have method field."));
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
                        &this->system,
                        sect["minrange"].get<inicpp::float_ini_t>(),
                        sect["maxrange"].get<inicpp::float_ini_t>(),
                        methods.at(a),
                        spins);
                
                if (setDebug) core->setDebug();

                this->parameters.push_back(std::move(core));
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
                        &this->system, X, Y,
                        sect["distance"].get<inicpp::float_ini_t>(),
                        minRange,
                        sect["maxrange"].get<inicpp::float_ini_t>());
                //core._methodVar = m;
                if (setDebug) core->setDebug();

                // if need to save the histogram to the file
                if (sect.contains("histogram") && sect["histogram"].get<inicpp::boolean_ini_t>()==true){
                    core->enableHistogram("histogram_"+parameterId+"_$.txt"); //replace dollar by number in time of saving
                }

                this->parameters.push_back(std::move(core));

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
                        &this->system,
                        ConfigManager::strToVect(a),
                        spins);

                if (setDebug) core->setDebug();
                if (setModule) core->setModule(true);

                this->parameters.push_back(std::move(core));
                ++i;
            }
        } else if (parameterName == "magnetisationlength") {

            std::vector<uint64_t> spins;
            if (sect.contains("spins"))
                spins = sect["spins"].get_list<inicpp::unsigned_ini_t>();

            std::unique_ptr<MagnetisationLengthCore> core = 
                make_unique<MagnetisationLengthCore>(parameterId,
                    &this->system,
                    spins);

            if (setDebug) core->setDebug();

            this->parameters.push_back(std::move(core));

        } else {
            throw(std::string("Parameter " + parameterName + " is unknown"));
        }
    }
*/

    //obtain the avaliable number of threads
    #pragma omp parallel
    {
        #pragma omp single
        {
            threadCount = omp_get_num_threads();
        }
    }

    if (!this->check_config()){
        throw string("Configuration error");
    }

    return;
}

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
        if (!co->check(this->system->size())) return false;
    }

    return true;
}

void ConfigManager::printHeader()
{
    double e = this->system->E();
    for (auto p : this->system->parts)
    {
        e -= scalar(p.m, this->field);
    }

    printf("# Metropolis algorithm for calculating heating capacity v%s\n",METROPOLIS_VERSION);
    printf("#   sysfile: %s\n",this->sysfile.c_str());

    this->system->printHeader(this->field);

    printf("#        MC: %u heatup, %u compute steps\n",this->heatup,this->calculate);
    if (this->isRestart())
        printf("#   restart: enabled, delta E threshold: %g*energy=%g\n",
            this->getRestartThreshold(),
            fabs(this->getRestartThreshold()*e));
    else
        printf("#   restart: disabled\n");
    printf("#   threads: %d\n",threadCount);
    printf("#     rseed: %d+<temperature number>\n",this->seed);
    printf("#    temps.: %u pcs. from %e to %e\n",
        temperatures->sizeBase(),
        std::min_element(temperatures->cbeginBase(),temperatures->cendBase()).operator*(),
        std::max_element(temperatures->cbeginBase(),temperatures->cendBase()).operator*());
    printf("# temp list: %e",temperatures->at(0,0));
    for (int tt=1; tt<temperatures->sizeBase(); ++tt)
        printf(",%e",temperatures->at(tt,0));
    printf("\n");
    printf("#   params.: %zd\n",this->parameters.size());
    printf("#\n");

    unsigned i=0;
    for (auto & co : parameters){
        co->printHeader(i); ++i;
    }

    
}

void ConfigManager::printColumnNames()
{
    // print header
    printf("# legend (column names):\n");
    printf("# 1:T 2:C(T)/N 3:<E> 4:<E^2> 5:threadId 6:seed");
    unsigned i=7;
    for (auto & co : parameters){
        printf(" %d:<%s> %d:<%s^2>",i,co->parameterId().c_str(),i+1,co->parameterId().c_str());
        i+=2;
    }
    printf(" %d:acc_rate_heatup %d:acc_rate_calculate",i,i+1); i+=2; // внутри printf i++ работает некорректно
    if (this->pt_enabled()){
        printf(" %d:pt_acc_rate_heatup %d:pt_acc_rate_calculate",i,i+1); i+=2;
    }
    printf(" %d:id(T) %d:id(replica)",i,i+1); i+=2;
    printf(" %d:time,s %d:E_final",i,i+1);
    printf("\n");
    fflush(stdout);
    // end print header
}

void ConfigManager::getParameters(std::vector< std::unique_ptr< CalculationParameter > > & calculationParameters)
{
    for (auto & co : parameters){
        calculationParameters.push_back(std::unique_ptr<CalculationParameter>(co->copy()));
    }
    return;
}