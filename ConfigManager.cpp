#include "ConfigManager.h"

std::vector<config_section_t> ConfigManager::explode(const inicpp::section& origin, const vector<string>& params)
{
    std::vector<config_section_t> list_to_explode; //список секций которые нужно расширить
    list_to_explode.push_back({origin.get_name(),origin});

    std::vector<config_section_t> list_to_explode_copy;
    for (const auto & param : params){
        if (origin.contains(param) && origin[param].is_list()){
            // чтобы не изменять список по которому итерируем, копируем его в list_to_explode_copy
            list_to_explode_copy = list_to_explode;
            list_to_explode.clear();
            for (auto val : origin[param].get_list<inicpp::string_ini_t>()){
                // копируем список list_to_explode где каждый элемент list_to_explode[i][param] больше не массив
                for (auto origin_list_item : list_to_explode_copy){
                    origin_list_item.data[param].set<inicpp::string_ini_t>(val);
                    origin_list_item.name = origin.get_name() + "_" + to_string(list_to_explode.size());
                    list_to_explode.push_back(origin_list_item);
                }
            }
        }
    }
    // если ни одного параметра не раскрыли, то возвращаем оригинальный параметр в виде списка
    return list_to_explode; 
}

void ConfigManager::parseCalculationParameters(inicpp::config & conf)
{
    for (auto & sect: conf){
        const std::string parameterString = sect.get_name();

        // все вычисляемые параметры должны иметь в имени секции символ двоеточие
        std::size_t colpos=parameterString.find_first_of(':');
        if (colpos==std::string::npos) continue;

        const std::string parameterName = parameterString.substr(0,colpos);
        const std::string parameterId = parameterString.substr(colpos+1);

        const vector<string> & explode_params = ParameterFactory::instance().getExplodeList(parameterName);
        std::vector< config_section_t > sections = explode(sect,explode_params);
        for (auto &s : sections){
            this->parameters.push_back(ParameterFactory::instance().createParameter(parameterName, s, this));
            this->parameters.back()->setBinder(this->isBinder());
        }
    }
    return;
}

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
        if (sect.contains("binder") && 
            sect["binder"].get<inicpp::boolean_ini_t>()) this->_binder = 1;
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
    if (commandLineParameters->binder)
        this->_binder = 1;

    if (iniconfig.contains("parallel_tempering")){
        this->temperatures = PtParseConfig(iniconfig);
        this->temperatures->setBaseTemperatures(base_temperatures);
    } else {
        this->temperatures = new PtBalancerDefault(base_temperatures);
    }
    this->temperatures->init();

    this->system = make_shared<MagneticSystem>(this->sysfile,range,size,ham);
    
    this->saveStates = commandLineParameters->saveStates;
    this->saveShort = commandLineParameters->saveShort;
    this->saveStateFileBasename = this->sysfile.substr(0, this->sysfile.find_last_of("."));


    parseCalculationParameters(iniconfig);


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
        if (!co->check(this->system->N())) return false;
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
    printf("#    temps.: %lu pcs. from %e to %e\n",
        temperatures->sizeBase(),
        std::min_element(temperatures->cbeginBase(),temperatures->cendBase()).operator*(),
        std::max_element(temperatures->cbeginBase(),temperatures->cendBase()).operator*());
    printf("# temp list: %e",temperatures->at(0,0));
    for (int tt=1; tt<temperatures->sizeBase(); ++tt)
        printf(",%e",temperatures->at(tt,0));
    printf("\n");

    if (this->pt_enabled()){
        printf("# par.temp.: enabled, balancer: %s, each_step: %u\n",this->temperatures->name().c_str(),this->temperatures->get_each_step());
        this->temperatures->printHeader();
    } else {
        printf("# par.temp.: disabled\n");
    }

    printf("#   params.: %zd\n",this->parameters.size());
    printf("#\n");

    unsigned i=0;
    for (auto & co : parameters){
        printf("##### calculation param #%3d #####\n",i); ++i;
        co->printHeader(); 
        printf("##### ---------------------- #####\n#\n");
    }

    
}

void ConfigManager::printColumnNames()
{
    // print header
    printf("# legend (column names):\n");
    printf("# 1:T 2:C(T)/N 3:<E> 4:<E^2>");
    unsigned i=5;
    if(this->isBinder()){
        printf(" %d:<E^4>",i);
        i+=1;
    }
    printf(" %d:threadId %d:seed",i,i+1);
    i+=2;

    for (auto & co : parameters){
        printf(" %d:<%s> %d:<%s^2>",i,co->parameterId().c_str(),i+1,co->parameterId().c_str());
        i+=2;
        if(this->isBinder()){
            printf(" %d:<%s^4>",i,co->parameterId().c_str());
            i+=1;
        }
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