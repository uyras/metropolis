#include "MagneticSystem.h"

MagneticSystem::MagneticSystem(
    string fileName, 
    double range, 
    optional<Vect> size, 
    optional<hamiltonian_t> hamiltonian, 
    state_t state) :
    _range(range),
    _size(size),
    _hamiltonian(hamiltonian)
{
    _fileVersion = fileVersion(fileName); 
    _fileContent.clear();
    readFileToString(fileName); //кэшируем содержимое файла во внутреннюю строку
    applyState(state);
}

MagneticSystem::~MagneticSystem()
{
}

void MagneticSystem::printHeader(Vect field)
{
    double avgNeighb = std::accumulate(this->neighbours_count.begin(), this->neighbours_count.end(), 0.0 ) / this->size();

    double e = this->E();
    for (auto p : this->parts)
    {
        e -= scalar(p.m, field);
    }

    printf("#    system: %lu spins, ", this->size());
    if (_fileVersion!=0) printf("%f interaction range, ", this->getRange());
    printf("%f avg. neighbours\n", avgNeighb);
    printf("#   physics: energy: %g, ext.filed: (%g,%g,%g), ",e,field.x,field.y,field.z);
    if (_fileVersion==0)
        printf("hamiltonian: csv, ");
    else
        printf("hamiltonian: dipole, ");
    printf("space: 2D\n");
    if (_fileVersion!=0) {
        printf("#    bounds: ");
        if (_size->x == 0 && _size->y == 0 && _size->z == 0){
            printf("open\n");
            
        } else {
            printf("periodic, system size: (%g,%g,%g)\n",_size->x,_size->y,_size->z);
        }
    }
}

void MagneticSystem::load_csv()
{
    char delimiter = ';';

    stringstream file(_fileContent);

    vector < vector < double > > result; //todo переписать считывание csv сразу в итоговый массив

    int linenum = 0;
    int linecount = -1;
    do {
        string line;
        getline(file,line);
        trim(line);
        if (line.length()<2 || line[0]=='#') continue;

        if (linecount==-1){ //read count of columns from the first line
            linecount = count(line.begin(), line.end(), delimiter)+1;
            result.resize(linecount);
            for (int i = 0; i < linecount; ++i)
                result[i].resize(linecount);
        }
        
        int colnum = 0;
        size_t pos = 0;
        std::string sval;
        double dval;
        do {
            pos = line.find(delimiter);
            sval = (pos != std::string::npos) ? line.substr(0, pos) : line;
            if (sval.length()>0)
                dval = stod(sval);
            else
                dval = 0;
            line.erase(0, pos + 1);
            result[linenum][colnum] = dval;
            colnum++;
            if (colnum>linecount) throw(string("Too much columns in CSV file "));
        } while (pos != std::string::npos);
        linenum++;
        if (linenum>linecount) throw(string("Too much lines in CSV file"));
    } while (!file.eof());

    //проверяем конфигурацию системы
    if (_state.size()==0) _state.resize(result.size(),1);
    if (_state.size() != result.size()){
        throw(string("Your configuration with size ")+to_string(_state.size())+" can not be applied to system with "+to_string(result.size())+" spins");
    }

    // делаем на основе считанных данных матрицу энергий
    this->parts.clear(); //удаляем все частицы
    this->eMatrix.clear();
    this->neighbourNums.clear();
    this->neighbours_count.clear();
    this->neighbours_from.clear();
    
    this->eMatrix.reserve(result.size() * result.size());
    this->neighbourNums.reserve(result.size() * result.size());
    this->neighbours_count.reserve(result.size());
    this->neighbours_from.reserve(result.size());

    size_t i=0;
    for (auto rr : result){
        {
            // создание частиц нужно для того чтобы осталась возможность применять поле во время вычислений.
            double si = (double)_state[i];
            Part pt = {p:{(double)i,0,0},m:{si,si,si}}; 
            this->parts.emplace_back(pt);
        }
        this->neighbours_from.push_back(eMatrix.size());
        size_t nCount = 0;
        size_t j = 0;
        for (auto r: rr){
            if (r !=0){
                this->eMatrix.push_back(r * _state[i] * _state[j]);
                this->neighbourNums.push_back(j);
                nCount++;
            }
            j++;
        }
        this->neighbours_count.push_back(nCount);
        i++;
    }

    this->eMatrix.shrink_to_fit();
    this->neighbourNums.shrink_to_fit();
    this->neighbours_count.shrink_to_fit();
    this->neighbours_from.shrink_to_fit();
}

void MagneticSystem::load_v1()
{
    stringstream f(_fileContent);

    this->parts.clear(); //удаляем все частицы

    //сначала сохраняем xyz
    double dummy;
    f >> dummy;
    f >> dummy;
    f >> dummy;

    int i=0;

    //пропускаем строку с заголовками
    char c[256];
    f.getline(c,256,'\n');
    f.getline(c,256,'\n');

    //затем читаем все магнитные моменты системы и положения точек
    double radius = 0;
    string shape;
    while (!f.eof()) {
        Part temp;
        if (!(f >> temp.p.x).good()) break; //если не получилось считать - значит конец файла
        f >> temp.p.y;
        f >> temp.p.z;
        f >> temp.m.x;
        f >> temp.m.y;
        f >> temp.m.z;
        f >> dummy; //w
        f >> dummy; //h
        //f >> temp.sector; для MPI реализации, @todo потом перегрузить
        f >> dummy; //r

        f >> shape;

        parts.push_back(temp);
        i++;
    }
}

void MagneticSystem::load_v2()
{
    this->parts.clear(); //удаляем все частицы

    double dummy;

    stringstream f(_fileContent);
  
    f.seekg(0);
    string section = "[parts]";
    std::string str;
    while (!f.eof() && str != section){
        std::getline(f,str);
        rtrim(str);
    }
    if (str!=section)
        throw(string("section [parts] not found in file")); //todo сделать нормальные классы для исключений


    while (!f.eof()){
        getline(f,str);
        trim(str);
        if (str.empty()) continue;
        if (str[0]=='[' && str[str.length()-1]==']') break;
        stringstream helper(str);

        Part temp;
        helper >> dummy; //id

        helper >> temp.p.x;
        helper >> temp.p.y;
        helper >> temp.p.z;
        helper >> temp.m.x;
        helper >> temp.m.y;
        helper >> temp.m.z;
        helper >> dummy; //state

        parts.push_back(temp);
    }
}

void MagneticSystem::readFileToString(string fileName)
{
    ifstream inFile;
    inFile.open(fileName);
    if (!inFile.good()) throw(string("Error reading file ") + fileName);
    std::stringstream strStream;
    strStream << inFile.rdbuf(); //read the file
    _fileContent = strStream.str(); //str holds the content of the file
    inFile.close();
}

void MagneticSystem::buildEnergyTable()
{

    // делаем на основе считанных данных матрицу энергий
    this->eMatrix.clear(); //удаляем все частицы
    this->neighbourNums.clear();
    this->neighbours_count.clear();
    this->neighbours_from.clear();
    
    this->eMatrix.reserve(parts.size() * parts.size());
    this->neighbourNums.reserve(parts.size() * parts.size());
    this->neighbours_count.reserve(parts.size());
    this->neighbours_from.reserve(parts.size());

    Part tmp;
    for (size_t i=0; i<parts.size(); i++){
        this->neighbours_from.push_back(this->eMatrix.size());
        size_t nCount = 0;
        for (size_t j=0; j<parts.size(); j++){
            tmp.p = (_size) ? translatePBC(parts[i].p, parts[j].p, *_size) : parts[j].p;
            if (i!=j){
                if (_range==0 || distance(parts[i].p,tmp.p) < _range){ // если range==0 то считаем все со всеми
                    tmp.m = parts[j].m;
                    double e = (*_hamiltonian)(parts[i],tmp);
                    this->eMatrix.push_back( e );
                    this->neighbourNums.push_back(j);
                    nCount++;
                }
            }
        }
        this->neighbours_count.push_back(nCount);
    }    
    
    this->eMatrix.shrink_to_fit();
    this->neighbourNums.shrink_to_fit();
    this->neighbours_count.shrink_to_fit();
    this->neighbours_from.shrink_to_fit();
    
}

void MagneticSystem::save(string filename, state_t state) const
{
    if (state.size()==0){
        state.resize(this->size());
    }

    ofstream f;
    f.open(filename, ios_base::out|ios_base::trunc);
    if (f.fail())
        throw(string("saveHelper: file "+filename+" is unwritable or not found"));
    
    f<<"[header]"<<endl;
    f<<"version=2"<<endl;
    f<<"dimensions=3"<<endl;
    f<<"type=standart"<<endl;
    f<<"size="+std::to_string(this->size())<<endl;
    f<<"state="+string(this->size(),'0')<<endl;
    f<<"interactionrange="+std::to_string(this->_range)<<endl;
    f<<"sizescale=1"<<endl;
    f<<"magnetizationscale=1"<<endl;
    f<<"[parts]"<<endl;

    for (size_t i=0; i<this->size(); i++) {
        f << i << "\t";
        f << parts[i].p.x << "\t";
        f << parts[i].p.y << "\t";
        f << parts[i].p.z << "\t";
        f << parts[i].m.x * state[i] << "\t";
        f << parts[i].m.y * state[i] << "\t";
        f << parts[i].m.z * state[i] << "\t";
        f << "0";
        f<<endl;
    }

    f.close();
}

void MagneticSystem::applyState(const vector<signed char> s)
{
    _E = nullopt;
    this->_state = s;
    if (this->_fileVersion == 0){
        // в случае с csv таблица энергий строится на лету
        // на систему не влияют граничные условия, гамильтониан и радиус взаимодействия, 
        // какими бы они не были. Влияет только конфигурация спинов
        if (this->_hamiltonian)
            throw(string("Your file is CSV, hamiltonian parameter is not applicable"));
        if (this->_range)
            throw(string("Your file is CSV, range parameter is not applicable"));
        if (this->_size)
            throw(string("Your file is CSV, PBC is not applicable"));
        this->load_csv(); 
    } else if (this->_fileVersion == 1 || this->_fileVersion == 2){
        if (this->_fileVersion == 1) this->load_v1();
        if (this->_fileVersion == 2) this->load_v2();

        //проверяем конфигурацию системы
        if (_state.size()==0) {
            _state.resize(parts.size(),1);
        } else {
            if (_state.size() != parts.size()){
                throw(string("Your configuration with size ")+to_string(_state.size())+" can not be applied to system with "+to_string(parts.size())+" spins");
            }
            // переворачиваем спины как задано в конфигурации
            for (size_t i=0; i<parts.size(); i++){
                parts[i].m.x *= _state[i];
                parts[i].m.y *= _state[i];
                parts[i].m.z *= _state[i];
            }
        }

        buildEnergyTable();
    }
}

int MagneticSystem::fileVersion(std::string file)
{
    if (ends_with(file,".csv")){ //if filename ends with .csv
        return 0;
    } else if (ends_with(file,".mfsys")){
        std::ifstream f(file);
        if (f.good()) {
            std::string s;
            std::getline(f,s);
            rtrim(s);
            if (s=="[header]"){
                f.close();
                return 2;
            } else {
                std::getline(f,s); //read 2 line
                std::getline(f,s); //read 3 line
                std::getline(f,s); //read 4 line
                if (s=="x\ty\tz\tMx\tMy\tMz\tr"){
                    f.close();
                    return 1;
                } else {
                    f.close();
                    return -1;
                }
            }
            f.close();
        } else {
            std::string("file "+file+" not found"); //todo сделать нормальные классы для исключений
        }
    } else {
        throw(std::string("Workg input file extention. Only mfsys and csv files are supported!"));
    }
    return 0;
}
