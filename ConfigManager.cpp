#include "ConfigManager.h"

Vect ConfigManager::size;

bool ConfigManager::check_config()
{
    //check main parameters
    {
        if (this->calculate<1){
            cerr<<"error! --calculate should be greather than 0!"<<endl;
            return 0;
        }

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

    printf("# Metropolis algorithm for calculating heating capacity v%s\n",METROPOLIS_VERSION);

    printf("#   sysfile: %s\n",this->sysfile.c_str());
    printf("#    system: %d spins, %f interaction range, %f avg. neighbours\n",this->system.size(),this->range, avgNeighb);
    printf("#   physics: ext.filed: (%g,%g), hamiltonian: dipole, space: 2D\n",this->field.x,this->field.y);
    printf("#    bounds: ");
    if (this->isPBC()){
        printf("periodic, system size: (%g,%g)\n",ConfigManager::size.x,ConfigManager::size.y);
    } else {
        printf("open\n");
    }
    printf("#        MC: %u heatup, %u compute steps\n",this->heatup,this->calculate);
    printf("#   threads: %d\n",threadCount);
    printf("#     rseed: %d+<temperature number>\n",this->seed);
    printf("#    temps.: %d pcs. from %e to %e\n",
        temperatures.size(),
        std::min_element(temperatures.begin(),temperatures.end()).operator*(),
        std::max_element(temperatures.begin(),temperatures.end()).operator*());
    printf("# temp list: %e",temperatures[0]);
    for (int tt=1; tt<temperatures.size(); ++tt)
        printf(",%e",temperatures[tt]);
    printf("\n");
    printf("# legend (column names):\n");
    printf("# 1:T 2:C(T)/N 3:<E> 4:<E^2> 5:threadId 6:seed");
    printf(" %d:time,s",7);
    printf("\n");
    fflush(stdout);
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
    r2 = rij.x * rij.x + rij.y * rij.y;
    r = sqrt(r2); //трудное место, заменить бы
    r5 = r2 * r2 * r; //радиус в пятой
    
    E = //энергия считается векторным методом, так как она не нужна для каждой оси
            (( (a->m.x * b->m.x + a->m.y * b->m.y) * r2)
                -
                (3 * (b->m.x * rij.x + b->m.y * rij.y) * (a->m.x * rij.x + a->m.y * rij.y)  )) / r5;
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
    
    return dist;
}