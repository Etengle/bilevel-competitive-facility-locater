#ifndef UTIL_H
#define UTIL_H


#include <bits/stdc++.h>
#include <getopt.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <regex>

#ifdef __linux__
#define bold_on  "\e[1m"
#define bold_off "\e[0m"
#else
#define bold_on  ""
#define bold_off ""
#endif

#define DEBUG_BREAK_POINT \
if (_gp->debugLevel() > 0) std::cerr << "DEBUG BREAK POINT HERE!\n";

#define error_msg(str) \
std::cout << "Error: " << str << std::endl;

#define warning_msg(str) \
std::cout << "Warning: " << str << std::endl;

#define info_msg(str) \
std::cout << "Info: " << str << std::endl;

#define warning_if_not(cond,str) \
if (!(cond)) { warning_msg(str); }

#define warning_if(cond,str) \
if ((cond)) { warning_msg(str); }

#define error_if_not(cond,str) \
if (!(cond)) { error_msg(str); assert(cond); }

#define error_if(cond,str) \
if ((cond)) { error_msg(str); assert(!(cond)); }

#define error_assert(str) \
error_if_not(false,str);


using PartSolution = boost::dynamic_bitset<>;
using FullSolution = std::pair<PartSolution, PartSolution>;

// for cout enum
template<typename T>
std::ostream& operator<<(typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& out, const T& e)
{
    return out << static_cast<typename std::underlying_type<T>::type>(e);
}

// for underlying type enum
template <typename E>
constexpr typename std::underlying_type<E>::type V(E e) {
    return static_cast<typename std::underlying_type<E>::type>(e);
}

// hash enum cless
struct EnumClassHash
{
    template <typename T>
    std::size_t operator()(T t) const
    {
        return static_cast<std::size_t>(t);
    }
};

enum DAY_TYPE {DAY_UNINIT = -1, WEEKDAY, WEEKEND, DAY_SIZE};
const std::map<DAY_TYPE, std::string> _dtName = {
    {DAY_UNINIT, "UNINIT"},
    {WEEKDAY, "WEEKDAY"}, {WEEKEND, "WEEKEND"},
    {DAY_SIZE, "END"}
};
const std::vector<DAY_TYPE> _vecDt = {WEEKDAY, WEEKDAY, WEEKDAY, WEEKDAY, WEEKDAY, WEEKEND, WEEKEND};
const std::map<std::string, DAY_TYPE> _nameDt = {
    {"UNINIT", DAY_UNINIT},
    {"WEEKDAY", WEEKDAY}, {"WEEKEND", WEEKEND},
    {"END", DAY_SIZE}
};
enum MARKET_TYPE {MARKET_UNINIT = -1, SUPER, HYPER, MARKET_SIZE};
const std::map<MARKET_TYPE, std::string> _mtName = {
    {MARKET_UNINIT, "UNINIT"},
    {SUPER, "SUPER"}, {HYPER, "HYPER"},
    {MARKET_SIZE, "END"}
};
const std::map<std::string, MARKET_TYPE> _nameMt = {
    {"UNINIT", MARKET_UNINIT},
    {"SUPER", SUPER}, {"HYPER", HYPER},
    {"END", MARKET_SIZE}
};
// using FAMILY_TYPE = int;
enum FAMILY_TYPE {FAMILY_UNINIT = -1, POOR, MEDIUM, RICH, FAMILY_SIZE};
const std::map<FAMILY_TYPE, std::string> _ftName = {
    {FAMILY_UNINIT, "UNINIT"},
    {POOR, "POOR"}, {MEDIUM, "MEDIUM"}, {RICH, "RICH"},
    {FAMILY_SIZE, "END"}
};
const std::map<std::string, FAMILY_TYPE> _nameFt = {
    {"UNINIT", FAMILY_UNINIT},
    {"POOR", POOR}, {"MEDIUM", MEDIUM}, {"RICH", RICH},
    {"END", FAMILY_SIZE}
};

class Layer;
class Family;

inline void printBanner(const std::string& str, std::ostream& out = std::cout, const std::string& prefix = ""){
    out << prefix << "===== " << str << " =====" << std::endl;
}

inline void printBanner(const std::string& str, const int& n = 5, const char& c = '=', const std::string& prefix = "", std::ostream& out = std::cout){
    out << prefix;
    for (int i = 0; i < n; ++i)
        out << c;
    out << " " << str << " ";
    for (int i = 0; i < n; ++i)
        out << c;
    out << std::endl;
}

inline bool myGetline(std::istream& fin, std::stringstream& ss){
    std::string buf;
    if (std::getline(fin, buf)) {
        ss.str("");
        ss.clear();
        ss.str(buf);
        // std::cout << "ss.str(buf) = " << ss.str() << std::endl;
        return true;
    }
    return false;
}

template <typename T>
inline bool myGet(std::stringstream& ss, T& d) {
    // std::cout << "in.ss.str() = " << ss.str() << std::endl;
    if (ss >> d) {
        // std::cout << "in.d = " << d << std::endl;
        return true;
    }
    ss.str("");
    ss.clear();
    return false;
}
template <typename T>
inline bool myGet(std::stringstream& ss, T& x, T& y) {
    char _;
    return ss >> _ >> x >> _ >> y >> _;
}

template <class T>
class R {
    public:
        R() : _rd(new std::random_device()), _X(std::mt19937((*_rd)())) {}
        ~R() { /* delete _rd; */ }
        T Ui(const T& s, const T& e) {
            std::uniform_int_distribution<T> U(s, e);
            return U(_X);
        }
        T Ui(const T& e) { return Ui(1, e); }
        T Ur(const T& s, const T& e) {
            std::uniform_real_distribution<T> U(s, e);
            return U(_X);
        }
        T Ur(const T& e) { return Ur(std::numeric_limits<T>::epsilon(), e); }
        T Nr(const T& mu, const T& sigma)  {
            std::normal_distribution<T> N(mu, sigma);
            return N(_X);
        }
        T Nr(const T& sigma) { return Nr(0, sigma); }
        T Bi(const T& n, const T& p) {
            std::uniform_real_distribution<T> B(n, p);
            return B(_X);
        }
        void saveSeed(std::ostream& out) { out << _X; }
        void loadSeed(std::istream& in) { in >> _X; }
        std::mt19937& X() { return _X; }
    private:
        std::random_device *_rd = nullptr;
        std::mt19937 _X;
};


class GlobalParameters {

    friend std::ostream& operator << (std::ostream& out, const GlobalParameters& gp);

public:

    void show(std::ostream& out = std::cout) const { out << *this << std::endl; }
    void setCP(const MARKET_TYPE& mt, const double& c) {
        warning_if_not(c > 0, "commodity price should be a positive number instead of " + std::to_string(c));
        error_if_not(mt >= SUPER && mt <= HYPER, "market type should be SUPER or HYPER " + _mtName.at(mt));
        _cp[mt] = c;
    }
    void setPA(const MARKET_TYPE& mt, const double& p) {
        warning_if_not(p >= 0, "parking attraction should be a non-negative number instead of " + std::to_string(p));
        error_if_not(mt >= SUPER && mt <= HYPER, "market type should be SUPER or HYPER " + _mtName.at(mt));
        _pa[mt] = p;
    }
    void setBP(const FAMILY_TYPE& ft, const double& b) {
        warning_if_not(b > 0, "buying power should be a positive number instead of " + std::to_string(b));
        error_if_not(ft >= POOR && ft <= RICH, "family type should be POOR or MEDIUM or RICH " + _ftName.at(ft));
        _bp[ft] = b;
    }
    void setAT(const MARKET_TYPE& mt, const FAMILY_TYPE& ft, const DAY_TYPE& dt, const double& a) {
        warning_if_not(a > 0, "attraction should be a positive number instead of " + std::to_string(a));
        error_if_not(mt >= SUPER && mt <= HYPER, "market type should be SUPER or HYPER " + _mtName.at(mt));
        error_if_not(ft >= POOR && ft <= RICH, "family type should be POOR or MEDIUM or RICH " + _ftName.at(ft));
        error_if_not(dt >= WEEKDAY && dt <= WEEKEND, "day type should be WEEKDAY or WEEKEND " + _dtName.at(dt));
        _at[mt][ft][dt] = a;
    }
    void setTheta(const int& i, const double& t){
        warning_if_not(t > 0, "theta should be a positive number instead of " + std::to_string(t));
        error_if_not(i >= 1 && i <= 4, "index should be 1, 2, 3, or 4" + std::to_string(i));
        _theta[i] = t;
    }
    void setNumFamily(const int& numFamily) {
        warning_if_not(numFamily > 0, "#family should be a positive integer instead of " + std::to_string(numFamily));
        _numFamily = numFamily;
    }
    void setMapVecFamily(const FAMILY_TYPE& ft, Family* fam) {
        error_if_not(ft >= POOR && ft <= RICH, "family type should be POOR or MEDIUM or RICH " + _ftName.at(ft));
        _mapVecFamily[ft].push_back(fam);
    }
    void setCache(bool c) { _is_cache = c; }
    void saveSeeds() {
        std::ofstream fi("seed_i.txt", std::ios_base::out);
        std::ofstream fr("seed_r.txt", std::ios_base::out);
        if (fi.is_open()) {
            _ri.saveSeed(fi);
            fi.close();
        }
        if (fr.is_open()) {
            _rr.saveSeed(fr);
            fr.close();
        }
    }
    void saveCSeeds() {
        std::ofstream fic("seed_i_" + currentDateTime() + ".txt", std::ios_base::out);
        std::ofstream frc("seed_r_" + currentDateTime() + ".txt", std::ios_base::out);
        if (fic.is_open()) {
            _ri.saveSeed(fic);
            fic.close();
        }
        if (frc.is_open()) {
            _rr.saveSeed(frc);
            frc.close();
        }
    }
    void loadSeeds() {
        std::ifstream fi("seed_i.txt", std::ios_base::in);
        std::ifstream fr("seed_r.txt", std::ios_base::in);
        if (fi.is_open()) {
            _ri.loadSeed(fi);
            fi.close();
        }
        if (fr.is_open()) {
            _rr.loadSeed(fr);
            fr.close();
        }
    }

    const std::string currentDateTime() const {
        time_t     now = time(0);
        struct tm  tstruct;
        char       buf[80];
        tstruct = *localtime(&now);
        // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
        // for more information about date/time format
        strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

        return buf;
    }

    inline const std::array<double, MARKET_SIZE>& cp() const { return _cp; }
    inline const std::array<double, MARKET_SIZE>& pa() const { return _pa; }
    inline const std::array<double, FAMILY_SIZE>& bp() const { return _bp; }
    inline const std::array<std::array<
                 std::array<double, MARKET_SIZE>,
                 FAMILY_SIZE>, DAY_SIZE>& at() const { return _at; }
    inline const std::array<double, 5>& theta() const { return _theta; }
    inline const std::vector<Family*>& vecFamily() const { return _vecFamily; }
    inline const std::map<FAMILY_TYPE, std::vector<Family*>>& mapVecFamily() const { return _mapVecFamily; }
    inline const int& numFamily() const { return _numFamily; }
    inline bool isVerbose() const { return _is_verbose; }
    inline bool isQt() const { return _is_qt; }
    inline bool isCache() const { return _is_cache; }
    inline bool isInputFile() const { return _is_input_file; }
    inline bool isOutputFile() const { return _is_output_file; }
    inline bool isSaveSeeds() const { return _is_save_seeds; }
    inline bool isLoadSeeds() const { return _is_load_seeds; }
    inline bool isTimed() const { return _is_timed; }
    inline bool isCircle() const { return _is_circle; }
    inline int debugLevel() const { return _debug_level; }
    inline int randomTimes() const { return _random_times; }
    inline const Layer* upperLayer() const { return _upperLayer; }
    inline const Layer* lowerLayer() const { return _lowerLayer; }
    inline const int& bitDiff() const { return _bit_diff; }
    inline const double& maxbp() const { return _maxbp; }
    inline bool drawOnly() const { return _drawOnly; }
    inline bool isIgnoreInitial() const { return _is_ignore_initial; }
    inline bool isNeighborCache() const { return _is_neighbor_cache; }
    inline void startTiming() { _g_start = clock(); return; }
    inline void showElapsedTime() { std::cout << "Elapsed Time = " << ((double) (clock() -_g_start) / CLOCKS_PER_SEC) << " sec." << std::endl; }


    std::array<double, MARKET_SIZE> _cp = {{}}; // commodity price
    std::array<double, MARKET_SIZE> _pa = {{}}; // parking attraction
    std::array<double, FAMILY_SIZE> _bp = {{}}; // buying power
    std::array<double, 5> _theta = {{}}; // theta [1]~[4], theta [0] don't care
    std::array<std::array<std::array<double, MARKET_SIZE>, FAMILY_SIZE>, DAY_SIZE> _at{{}}; // attraction
    std::vector<Family*> _vecFamily = std::vector<Family*>();
    // std::array<int, FAMILY_SIZE> _numFamily{{}}; // number of family each type
    std::map<FAMILY_TYPE, std::vector<Family*>> _mapVecFamily = std::map<FAMILY_TYPE, std::vector<Family*>>();

    int _numFamily = 0;
    bool _is_verbose = false, _is_input_file = false, _is_output_file = false;
    bool _is_exact = false, _is_qt = false, _is_cache = false, _is_timed = false;
    bool _is_save_seeds = false, _is_load_seeds = false, _is_circle = false;
    int _random_times = 0;
    int _debug_level = 0;
    int _bit_diff = 4;
    std::string _input_file_name = "parameter.txt";
    std::string _output_file_name = "parameter_out.txt";
    Layer *_upperLayer = nullptr, *_lowerLayer = nullptr;
    double _maxbp = 0.0;
    bool _drawOnly = false;
    std::string _initial_sol_file_name = "";
    bool _is_ignore_initial = false;
    bool _is_neighbor_cache = false;
    clock_t _g_start = 0;

    R<int> _ri;
    R<double> _rr;
};

class FuncLog {
    public:
        FuncLog(GlobalParameters* gp, const std::string& msg = "", const std::string& prefix = "",
                bool enforce = false) :
            _gp(gp), _msg(msg), _prefix(prefix), _isEnforced(enforce) {}
        ~FuncLog() {
            if (_isEnforced || (_gp && _gp->isTimed())){
                printBanner("Function " + _msg + " end", 5, '-', _prefix);
                printBanner("Elapsed Time " +
                            std::to_string((double) (clock() - _start) / CLOCKS_PER_SEC) + " sec.",
                            5, '-', _prefix);
                std::cout << std::endl;
            }
        }
        void startLog() {
            _start = clock();
            if (_isEnforced || (_gp && _gp->isTimed()))
                printBanner("Function " + _msg + " start", 5, '-', _prefix);
        }

        void endLog() {
            if (_isEnforced || (_gp && _gp->isTimed())){
                printBanner("Function " + _msg + " end", 5, '-', _prefix);
                printBanner("Elapsed Time " +
                            std::to_string((double) (clock() - _start) / CLOCKS_PER_SEC) + " sec.",
                            5, '-', _prefix);
                std::cout << std::endl;
            }
        }
    private:
        GlobalParameters* _gp = nullptr;
        std::string _msg = "", _prefix = "";
        bool _isEnforced = false;
        clock_t _start;
};


inline std::ostream& operator << (std::ostream& out, const PartSolution& bs) {
    out << std::noboolalpha;
    for (int i = 0, sz = bs.size(); i < sz; ++i)
        out << bs[i];
    return out;
}

inline std::ostream& operator << (std::ostream& out, const FullSolution& fs) {
    out << std::noboolalpha;
    out << "S(Y = " << fs.first << ", X = " << fs.second << ")";
    return out;
}

inline int dist(const PartSolution& a, const PartSolution& b) {
    return (a ^ b).count();
}

namespace std {

  template <>
  struct hash<PartSolution>
  {
    std::size_t operator()(const PartSolution& p) const
    {
      return std::hash<unsigned long>()(p.to_ulong());
    }
  };

  template <>
  struct hash<FullSolution>
  {
    std::size_t operator()(const FullSolution& f) const
    {
      return (std::hash<PartSolution>()(f.first) << 1) ^ (std::hash<PartSolution>()(f.second) >> 1);
    }
  };

  template <>
  struct hash<DAY_TYPE>
  {
    std::size_t operator()(const DAY_TYPE& d) const
    {
      return std::hash<int>()(d);
    }
  };

}

#endif // UTIL_H
