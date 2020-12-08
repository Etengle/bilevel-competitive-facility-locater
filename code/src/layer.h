#ifndef LAYER_H
#define LAYER_H
#include "facility.h"
#include "family.h"

class Layer{
	
	friend std::ostream& operator << (std::ostream&, const Layer&);

    class Tabu {
        friend class Layer;
    public:
        inline const int& size(void) const { return _size; }
        inline const int& gmax(void) const { return _gMax; }
        inline const int& lmax(void) const { return _lMax; }
        inline const int& cmax(void) const { return _cMax; }
        inline const double& itmax(void) const { return _itMax; }
    private:
        int _size = 0;
        int _gMax = 0, _lMax = 0, _cMax = 0;
        double _itMax = 0.0;

    };

	public:

        Layer(GlobalParameters* gp, const std::string & name = "Unknown Layer",
              bool isHyper = false, bool isUpper = false)
            : _gp(gp), _name(name), _isHyper(isHyper), _isUpper(isUpper), _mt((_isHyper) ? HYPER : SUPER), _tabu(new Tabu()) {}
        /*  _numFacility(0), _numExistedFacility(-1),
            _budget(0), _gMax(0), _kMax(0), _lMax(0),
            _vecFacility(std::vector<Facility*>()),
            _vecExistedFacility(std::vector<Facility*>()),
            _vecFamily(std::vector<Family*>()) */

        ~Layer() {}

        void show(std::ostream& out = std::cout) { out << *this << std::endl; }
        void randomGenerateOrAdjust();
        void setName(const std::string& name) { _name = name; }
        void setBudget(const double& budget) {
            warning_if_not(budget > 0, _name + " budget should be a positive number instead of " + std::to_string(budget));
            _budget = budget;
		}
        void setNumFacility(const int& numFacility) {
            warning_if_not(numFacility > 0, _name + " #candidateFacility should be a positive integer instead of " + std::to_string(numFacility));
            _numFacility = numFacility;
		}
        void setNumExistedFacility(const int& numExistedFacility) {
            warning_if_not(numExistedFacility >= 0, _name + " #existedFacility should be a non-negative integer instead of " + std::to_string(numExistedFacility));
            _numExistedFacility = numExistedFacility;
        }
        void setGMax(const int& gMax) {
            warning_if_not(gMax > 0, _name + " gmax should be a positive integer instead of " + std::to_string(gMax));
            _gMax = gMax;
        }
        void setKMax(const int& kMax) {
            warning_if_not(kMax > 0, _name + " kmax should be a positive integer instead of " + std::to_string(kMax));
            _kMax = kMax;
        }
        void setLMax(const int& lMax) {
            warning_if_not(lMax > 0, _name + " kmax should be a positive integer instead of " + std::to_string(lMax));
            _lMax = lMax;
        }
        void setItMax(const double& itMax) {
            warning_if_not(itMax > 0, _name + " itmax should be a positive number instead of " + std::to_string(itMax));
            _itMax = itMax;
        }
        void setVecFacility(const std::vector<Facility*>& vecFacility) {
            warning_if_not((int)vecFacility.size() == _numFacility, _name +
                           " vecFacility.size() = " + std::to_string(vecFacility.size()) +
                           " != _numFacility = " + std::to_string(_numFacility));
            _vecFacility = vecFacility;
        }
        void setVecExistedFacility(const std::vector<Facility*>& vecExistedFacility) {
            warning_if_not((int)vecExistedFacility.size() == _numExistedFacility, _name +
                           " vecFacility.size() = " + std::to_string(vecExistedFacility.size()) +
                           " != _numExistedFacility = " + std::to_string(_numExistedFacility));
            _vecExistedFacility = vecExistedFacility;
        }
        void setVecFamily(const std::vector<Family*>& vecFamily) {
            _vecFamily = vecFamily;
        }
        void setTabuSize(const int& size) {
            warning_if_not(size > 0, _name + " tabu size should be a positive integer instead of " + std::to_string(size));
            _tabu->_size = size;
        }
        void setTabuGMax(const int& gMax) {
            warning_if_not(gMax > 0, _name + " tabu gmax should be a positive integer instead of " + std::to_string(gMax));
            _tabu->_gMax = gMax;
        }
        void setTabuLMax(const int& lMax) {
            warning_if_not(lMax > 0, _name + " tabu lmax should be a positive integer instead of " + std::to_string(lMax));
            _tabu->_lMax = lMax;
        }
        void setTabuCMax(const int& cMax) {
            warning_if_not(cMax > 0, _name + " tabu cmax should be a positive integer instead of " + std::to_string(cMax));
            _tabu->_cMax = cMax;
        }
        void setTabuItMax(const double& itMax) {
            warning_if_not(itMax > 0, _name + " itmax should be a positive number instead of " + std::to_string(itMax));
            _tabu->_itMax = itMax;
        }
        void setVecNSS(const std::vector<double>& vecNSS) {
            warning_if_not((int)vecNSS.size() == _kMax+1, _name +
                           " vecNSS.size()-1 = " + std::to_string((int)(vecNSS.size())-1) +
                           " != _kMax = " + std::to_string(_kMax));
            _vecNSS = vecNSS;
            for (int k = 1; k <= _kMax; ++k)
                warning_if_not(_vecNSS[k] > _vecNSS[k-1], _name + " k = " + std::to_string(k) +
                        ", NSS[k] = " + std::to_string(_vecNSS[k]) + " <= NSS[k-1] = " + std::to_string(_vecNSS[k-1]));
        }

        void allocateSolutions();
        void initSolutions();
        void random(PartSolution& sol);
        // inline void randomTry(const int&);
        bool isBudgetEnough(const PartSolution&, double&);
        bool isBudgetEnough(const PartSolution&);

        FullSolution tabuSearch(Layer*, const FullSolution&);
        FullSolution changeLocation(const FullSolution&, const int& = 1);

        inline double utility(Facility*, Family*, const DAY_TYPE&, bool isDebug = false);
        inline double estUtility(Facility*);
        double revenue(Layer*, const FullSolution&, bool = false, bool = false);
        inline double totalExistedUtilityPerFamily(Family*, const DAY_TYPE&, bool = false);
        int findNeighbors(const PartSolution&,
                          std::unordered_set<PartSolution>&,
                          const int& = 1, bool isInitial = true);
        inline void recursive(const PartSolution& curSol,
                              const PartSolution& flipedSol,
                              PartSolution& newSol,
                              std::unordered_set<PartSolution>& setSols,
                              int s0, const int& e0, int n0, const int& m0,
                              int s1, const int& e1, int n1, const int& m1);

        inline const Layer& getValue(void) const { return *this; }
        inline const std::string& name(void) const { return _name; }
        inline const double& budget(void) const { return _budget; }
        inline const int& numFacility(void) const { return _numFacility; }
        inline const int& numExistedFacility(void) const { return _numExistedFacility; }
        inline const int& gmax(void) const { return _gMax; }
        inline const int& kmax(void) const { return _kMax; }
        inline const int& lmax(void) const { return _lMax; }
        inline const double& itmax(void) const { return _itMax; }
        inline const std::vector<Facility*>& vecFacility() const { return _vecFacility; }
        inline const std::vector<Facility*>& vecExistedFacility() const { return _vecExistedFacility; }
        inline const std::vector<Family*>& vecFamily() const { return _vecFamily; }
        inline const MARKET_TYPE& mt() const { return _mt; }
        inline PartSolution& cur_sol(){ return _cur_sol; }
        inline const Tabu* tabu() const { return _tabu; }
        inline const double & tabuRuntime() const { return _tabuRuntime; }
        inline const double & maxTabuRuntime() const { return _maxTabuRuntime; }
        inline const double & minTabuRuntime() const { return _minTabuRuntime; }
        inline const int & tabuIter() const { return _tabuIter; }
        inline const double& nss(const int& k) const {
            error_if_not(k >= 0 && k <= _kMax, "access k out of range [0, kmax]");
            warning_if(k == 0, "nss[0] is always 0");
            return _vecNSS[k];
        }
        inline const std::vector<double>& vecNSS() const { return _vecNSS; }
        const GlobalParameters* gp() { return _gp; }
        void show(std::ostream& out = std::cout) const { out << *this << std::endl; }
        void showSolutions(std::ostream& out = std::cout) const {
            printBanner(_name + " Solution Info", out, "\t");
            out << "\tcurrentSolution = " << _cur_sol << std::endl
                << "\tpreviousSolution = " << _prev_sol << std::endl
                << "\tbestSolution = " << _best_sol << std::endl;
        }
        inline std::unordered_map<PartSolution, PartSolution>& tabuBestMap() { return _tabuBestMap; }
        inline std::unordered_map<FullSolution, double>& revenueCache() { return _revenueCache; }
        void setBestSol(const PartSolution& par) {
            error_if_not(par.size() == _best_sol.size(), "best solution size differs!");
            _best_sol = par;
        }
        inline const PartSolution& best_sol() const { return _best_sol; }
        inline bool isHyper() const { return _isHyper; }

	private:
        GlobalParameters* _gp = nullptr;
        std::string _name = "";
        int _numFacility = 0, _numExistedFacility = -1;
        double _budget = 0.0;
        int _gMax = 0, _kMax = 0, _lMax = 0;
        bool _isHyper = false, _isUpper = false;
        MARKET_TYPE _mt = MARKET_UNINIT;
        double _itMax = 0.0;
        std::vector<double> _vecNSS = std::vector<double>();

        std::vector<Facility*> _vecFacility = std::vector<Facility*>();
        std::vector<Facility*> _vecExistedFacility = std::vector<Facility*>();
        std::vector<Family*> _vecFamily = std::vector<Family*>(); // shared for each layer
        PartSolution _cur_sol, _prev_sol, _best_sol;


        std::unordered_map<Facility*, std::unordered_map<Family*, std::unordered_map<DAY_TYPE, double>>> _utilCache;
        std::unordered_map<Family*, std::unordered_map<DAY_TYPE, double>> _utilPerFamCache;
        std::unordered_map<PartSolution, bool> _budgetCache;
        std::unordered_map<FullSolution, double> _revenueCache;

        Tabu* _tabu = nullptr;
        std::queue<PartSolution> _tabuList;
        std::unordered_map<PartSolution, bool> _tabuMap;


        std::unordered_map<PartSolution, PartSolution> _tabuBestMap;
        std::unordered_map<PartSolution, std::unordered_set<PartSolution>*> _neighborSetMap;

        int _tabuIter = 0;
        double _tabuRuntime = 0.0, _maxTabuRuntime = DBL_MIN, _minTabuRuntime = DBL_MAX;
        bool _isDebug = false;
        // std::map<PartSolution, double> _tabuBestRev;
        // std::valarray<bool> _cur_sol, _prev_sol, _best_sol;
};

#endif // __LAYER_H__
