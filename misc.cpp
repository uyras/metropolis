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

state_t xorstate(const state_t &s1, const state_t &s2)
{
    state_t s(s1);
	for (int i=0; i<s1.size(); i++){
		s[i] = s1[i] * s2[i];
	}
	return s;
}

std::string stateToString(const state_t &state)
{
    std::vector<Part*>::const_iterator iter;
    std::string str="";
    for (auto s : state){
        if (s>0){
            str+="0";
        } else {
            str+="1";
        }
    }
    return str;
}

Vect strToVect(std::string val){
    Vect target;
    if (val=="x" || val=="X"){
        target = {1,0,0};
    } else {
        if (val=="y" || val=="Y") 
            target = {0,1,0}; 
        else {
            if (val=="z" || val=="Z")
                target = {0,0,1};
            else
            {
                auto pos = val.find('|');
                if (pos!=std::string::npos){
                    double x = stod(val.substr(0,pos));
                    auto pos_second = val.find('|',pos+1);
                    if (pos_second!=std::string::npos){
                        double y = stod(val.substr(pos+1,pos_second));
                        double z = stod(val.substr(pos_second+1));
                        target = {x,y,z};
                    } else {
                        double y = stod(val.substr(pos+1));
                        target = {x,y,0};
                    }
                } else {
                    throw(std::invalid_argument("Vector format is x|y or x|y|z"));
                }
            }
        }
    }
    return target;
}

Vect translatePBC(const Vect &a, const Vect &b, const Vect size)
{
    Vect res = b;
    if (size.x && fabs(a.x - b.x) > size.x/2){
        if (a.x < b.x){
            res.x -= size.x;
        } else {
            res.x += size.x;
        }
    }

    
    if (size.y && fabs(a.y - b.y) > size.y/2){
        if (a.y < b.y){
            res.y -= size.y;
        } else {
            res.y += size.y;
        }
    }

    if (size.z && fabs(a.z - b.z) > size.z/2){
        if (a.z < b.z){
            res.z -= size.z;
        } else {
            res.z += size.z;
        }
    }
    return res;
}

double distance(const Vect &a, const Vect &b)
{
    Vect d = {a.x-b.x,a.y-b.y,a.z-b.z};
    return sqrt(
        d.x * d.x +
        d.y * d.y +
        d.z * d.z
        );
}

double scalar(const Vect &a, const Vect &b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}
