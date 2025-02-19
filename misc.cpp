#include "misc.h"

#include <algorithm> 
#include <cctype>
#include <locale>

/*================ inner functions ===================*/
// trim from start (in place)
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
inline void trim(std::string &s) {
    rtrim(s);
    ltrim(s);
}

std::string xorstr(std::string s1, std::string s2)
{
	std::string s(s1);
	for (int i=0; i<s1.size(); i++){
		s[i] = (s1[i]==s2[i])?'0':'1';
	}
	return s;
}

/*================ outer functions ===================*/
vector < vector < double > > readCSV(string filename){
    char delimiter = ';';

    ifstream file(filename);
    if (!file.is_open()) throw(string("Error reading csv file ") + filename);

    vector < vector < double > > result;

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
            if (colnum>linecount) throw(string("Too much columns in file ") + filename);
        } while (pos != std::string::npos);
        linenum++;
        if (linenum>linecount) throw(string("Too much lines in file ") + filename);
    } while (!file.eof());

    return result;
}