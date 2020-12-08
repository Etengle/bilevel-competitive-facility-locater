#include "layer.h"
#define W1 0.12
#define W2 0.2

std::ostream& operator << (std::ostream& out, const Layer& layer)
{

    out << "\tname = \"" << layer._name << "\"\n"
        << "\tmarketType = " << _mtName.at(layer._mt) << std::endl
        << "\tisUpper = " << layer._isUpper << std::endl
        << "\tisHyper = " << layer._isHyper << std::endl
        << "\tbudget = " << layer._budget << std::endl
        << "\tGmax = " << layer._gMax << std::endl
        << "\tKmax = " << layer._kMax << std::endl
        << "\tLmax = " << layer._lMax << std::endl
        << "\tItmax = " << layer._itMax << std::endl;

    out << "\tNSS = ";
    for (const auto& nss : layer._vecNSS)
        out << nss << " ";
    out << std::endl;

    printBanner("Facility Info", out, "\t");
    out << "\t#candidateFacility = " << layer._numFacility << std::endl;
    out << "\tvecFacility(" << layer._vecFacility.size() << ") = \n\t{" << std::endl;
    int i = 0;
    for (const Facility* fac : layer._vecFacility) {
        if (i > 3 && layer._gp->debugLevel() < 2) {
            out << "\t\t...\n";
            break;
        }
        out << "\t\t" << layer._name << "::CandidateFacility[" << i++ << "]: " << *fac << " " << std::endl;
    }
    out << "\t}" << std::endl;

    out << "\t#existedFacility = " << layer._numExistedFacility << std::endl;
    out << "\tvecExistedFacility(" << layer._vecExistedFacility.size() << ") = \n\t{" << std::endl;
    i = 0;
    for (const Facility* fac : layer._vecExistedFacility) {
        if (i > 3 && layer._gp->debugLevel() < 2) {
            out << "\t\t...\n";
            break;
        }
        out << "\t\t" << layer._name << "::ExistedFacility[" << i++ << "]: " << *fac << " " << std::endl;
    }
    out << "\t}" << std::endl;

    printBanner("Tabu Info", out, "\t");
    out << "\tsize = " << layer._tabu->size() << std::endl
        << "\tGmax = " << layer._tabu->gmax() << std::endl
        << "\tCmax = " << layer._tabu->cmax() << std::endl
        << "\tLmax = " << layer._tabu->lmax() << std::endl
        << "\tItmax = " << layer._tabu->itmax() << std::endl;

    layer.showSolutions(out);

    return out;
}


void Layer::recursive(const PartSolution& curSol,
                      const PartSolution& flippedSol,
                      PartSolution& newSol,
                      std::unordered_set<PartSolution>& setSols,
                      int s0, const int& e0, int n0, const int& m0,
                      int s1, const int& e1, int n1, const int& m1) {
    int j = 0;
    while (s0 < e0 && n0 < m0) {
        if (_gp->debugLevel() >= 2)
            std::cout << "n0 = " << n0 << ", s0 = " << s0 << ", flipped sol = "<< flippedSol << std::endl;
        size_t i = (s0 == 0 && j++ == 0 && n0 == 0) ? flippedSol.find_first() : flippedSol.find_next(s0);
        if (_gp->debugLevel() >= 2)
            std::cout << "find 0 i = " << i << ", newSol = " << newSol << std::endl;
        if (i == PartSolution::npos) return;
        newSol.flip(i);
        if (_gp->debugLevel() >= 2)
            std::cout << "find 0 i = " << i << ", changed newSol = " << newSol << std::endl;
        s0 = i;
        recursive(curSol, flippedSol, newSol, setSols,
                  s0, e0, n0+1, m0, s1, e1, n1, m1);
        newSol.flip(i);
    }
    int k = 0;
    while (s1 < e1 && n1 < m1) {
        if (_gp->debugLevel() >= 2)
            std::cout << "n1 = " << n1 << ", s1 = " << s1 << ", cur sol = " << curSol << std::endl;
        size_t i = (s1 == 0 && k++ == 0 && n1 == 0) ? curSol.find_first() : curSol.find_next(s1);
        if (_gp->debugLevel() >= 2)
            std::cout << "find 1 i = " << i << ", newSol = " << newSol << std::endl;
        if (i == PartSolution::npos) return;
        newSol.flip(i);
        if (_gp->debugLevel() >= 2)
            std::cout << "changed new sol = " << newSol << std::endl;
        s1 = i;
        recursive(curSol, flippedSol, newSol, setSols,
                  s0, e0, n0, m0, s1, e1, n1+1, m1);
        newSol.flip(i);
    }
    setSols.insert(newSol);
}

int Layer::findNeighbors(const PartSolution& curSol,
                         std::unordered_set<PartSolution>& setSols,
                         const int& k, bool isInitial) {
    FuncLog f(_gp, "find " + std::string((isInitial) ? "Initial" : (_isUpper) ? "Upper Tabu" : "Lower Tabu") + " Neighbors", "\t\t");
    if (1 || isInitial) {
        f.startLog();
    }

    if (_gp->isNeighborCache() && !isInitial){
        auto it = _neighborSetMap.find(curSol);
        if (it != _neighborSetMap.end() && it->second && (*it->second).size() > 0) {
            setSols = (*it->second);
            // if (_gp->isVerbose())
                std::cout << "Found cached "
                          << ((_isUpper) ? "Upper" : "Lower")
                          << " #tabu neighbors (+itself) = " << setSols.size() << std::endl;
            return (int) setSols.size();
        }
    }

    setSols.clear();
    std::unordered_set<PartSolution> setAllSols;
    setAllSols = setSols;
    setAllSols.insert(curSol);
    PartSolution flipedSol = curSol, newSol = curSol;
    flipedSol.flip();

    if (isInitial) {
        warning_if_not(k == 1, "find Initial Neighbors 1 != k = " + std::to_string(k));
        int numNSS = std::ceil(_vecNSS[k]*curSol.size());
        if (numNSS <= 0) {
            warning_msg("#change bit = 0, auto increase to 1");
            numNSS = 1;
        }
        if (_gp->isVerbose())
            std::cout << "Mode is "<< (_gp->isCircle() ? "circle" : "ring") << ".\n"
                      << "It will enumerate the changes containing " << (_gp->isCircle() ? "0 ~ " : "")
                      << numNSS << " (= " << _vecNSS[k]*100 << "% * "
                      << curSol.size() << ") 1's to 0's and also 0's to 1's\n";

        if (_gp->isCircle()) {
            for (int i = 0; i <= numNSS; ++i)
                for (int j = 0; j <= numNSS; ++j)
                    if (!(i == 0 && j == 0))
                        recursive(curSol, flipedSol, newSol, setAllSols,
                                  0, curSol.size(), 0, i,
                                  0, curSol.size(), 0, j);
        } else if (numNSS != 0)
            recursive(curSol, flipedSol, newSol, setAllSols,
                      0, curSol.size(), 0, numNSS,
                      0, curSol.size(), 0, numNSS);

        if (_gp->isVerbose())
            std::cout << "Found #VNS neighbors (+itself) = " << setAllSols.size() << "\n";

        while (numNSS < (int) curSol.size() && (setAllSols.size() <= 1 || setSols.size() <= 1)) {

            if (setAllSols.size() <= 1 || setSols.size() == 1) {
                setSols.clear();
                setAllSols.clear();
                setAllSols.insert(curSol);

                if (!_gp->isCircle()) {
                    warning_msg("None neighbor is found, auto change to circle search.");
                    _gp->_is_circle = true;
                }
                else if (numNSS < (int) curSol.size()) {
                    warning_msg("None neighbor is found, auto increase circle radius by 1.");
                    numNSS++;
                }

                for (int i = 0; i <= numNSS; ++i)
                    for (int j = 0; j <= numNSS; ++j)
                        if (!(i == 0 && j == 0))
                            recursive(curSol, flipedSol, newSol, setAllSols,
                                      0, curSol.size(), 0, i,
                                      0, curSol.size(), 0, j);
                if (_gp->isVerbose())
                    std::cout << "New found #VNS neighbors (+itself) = " << setAllSols.size() << "\n";
            }

            if (_gp->isVerbose()) {
                int j = 0;
                for (auto& i : setAllSols) {
                    if (j >= 3 && _gp->debugLevel() < 2) {
                        std::cout << "...\n";
                        break;
                    }
                    std::cout << "[" << ++j << "] = " << i << std::endl;
                }
                std::cout << "(Verbose) Mode is "<< (_gp->isCircle() ? "circle" : "ring") << ".\n"
                          << "(Verbose) It enumerated the changes containing " << (_gp->isCircle() ? "0 ~ " : "")
                          << numNSS << " (= " << _vecNSS[k]*100 << "%*"
                          << curSol.size() << ") 1's to 0's and also 0's to 1's.\n";
                std::cout << "(Verbose) Initial solution = " << curSol << "\n";
                std::cout << "(Verbose) Found total #VNS neighbors (+itself) = " << setAllSols.size() << "\n\n";
            }
            if (_gp->isVerbose())
                std::cout << "Excluding each solution that is not affordable...\n";
            for (auto& b : setAllSols){
                if (b.any() && isBudgetEnough(b))
                    setSols.insert(b);
            }
            if (_gp->isVerbose())
                std::cout << "Finalized #VNS neighbors (+itself) = " << setSols.size() << "\n\n";

            if (_gp->isVerbose()) {
                int j = 0;
                for (auto& i : setSols) {
                    if (j >= 3 && _gp->debugLevel() < 2) {
                        std::cout << "...\n";
                        break;
                    }
                    std::cout << "[" << ++j << "] = " << i << std::endl;
                }
                std::cout << "(Verbose) Mode is "<< (_gp->isCircle() ? "circle" : "ring") << ".\n"
                          << "(Verbose) It enumerated the changes containing " << (_gp->isCircle() ? "0 ~ " : "")
                          << numNSS << " (= " << _vecNSS[k]*100 << "%*"
                          << curSol.size()
                          << ") 1's to 0's and also 0's to 1's.\n";
                std::cout << "(Verbose) Initial solution = " << curSol << "\n";
                std::cout << "(Verbose) Originally found total #VNS neighbors (+itself) = " << setAllSols.size() << "\n";
                std::cout << "(Verbose) Finalized #VNS neighbors = " << setSols.size() << "\n\n";
            }
        }

        if (setSols.size() == 1)
            warning_msg("After several tries, still none neighbor is found. Please check your budget setting!\n");


    }
    else { // tabu
        std::unordered_set<PartSolution> addSols,swapSols,reduceSols;
        std::vector<PartSolution> vecAddSols, vecSwapSols, vecReduceSols;

        // AddLocation
        recursive(curSol, flipedSol, newSol, addSols,
                  0, curSol.size(), 0, 1,
                  0, curSol.size(), 0, 0);
        // SwapLocation
        recursive(curSol, flipedSol, newSol, swapSols,
                  0, curSol.size(), 0, 1,
                  0, curSol.size(), 0, 1);
        // ReduceLocation
        recursive(curSol, flipedSol, newSol, reduceSols,
                  0, curSol.size(), 0, 1,
                  0, curSol.size(), 0, 2);

        std::string ul = std::string(((_isUpper) ? "U" : "L"));
        int as = addSols.size(), ss = swapSols.size(), rs = reduceSols.size();
        int totalsize = as + ss + rs;
        if (_tabu->cmax() > totalsize) _tabu->_cMax = totalsize;
        double ratio = _tabu->cmax()/(double)totalsize;
        std::cout << "Found " + ul + "#tabu neighbors (+itself) = " << totalsize + setAllSols.size() <<
                     ", (#add = " << as << ", #swap = " << ss << ", #reduce = " << rs <<
                     ", cmax = " << _tabu->cmax() << ", ratio = " << ratio << ")" << std::endl;

        std::copy(addSols.begin(), addSols.end(), std::back_inserter(vecAddSols));
        std::copy(swapSols.begin(), swapSols.end(), std::back_inserter(vecSwapSols));
        std::copy(reduceSols.begin(), reduceSols.end(), std::back_inserter(vecReduceSols));

        std::random_shuffle(vecAddSols.begin(), vecAddSols.end());
        std::random_shuffle(vecSwapSols.begin(), vecSwapSols.end());
        std::random_shuffle(vecReduceSols.begin(), vecReduceSols.end());

        std::copy(vecAddSols.begin(), vecAddSols.begin() + std::ceil(ratio*as), std::inserter(setAllSols, setAllSols.end()));
        std::copy(vecSwapSols.begin(), vecSwapSols.begin() + std::ceil(ratio*ss), std::inserter(setAllSols, setAllSols.end()));
        std::copy(vecReduceSols.begin(), vecReduceSols.begin() + std::ceil(ratio*rs), std::inserter(setAllSols, setAllSols.end()));

        std::cout << "Random trimmed " + ul + "#tabu neighbors (+itself) = " << setAllSols.size() <<
                     ", (#add = " << std::ceil(ratio*as) << ", #swap = " << std::ceil(ratio*ss) << ", #reduce = " <<
                     std::ceil(ratio*rs) << ")" << std::endl;


        if (_gp->isVerbose()) {
            int j = 0;
            for (auto& i : setAllSols) {
                if (j >= 3 && _gp->debugLevel() < 2) {
                    std::cout << "...\n";
                    break;
                }
                std::cout << "[" << ++j << "] = " << i << std::endl;
            }
            std::cout << "(Verbose) Initial solution = " << curSol << "\n";
            std::cout << "(Verbose) Found total " + std::string(((_isUpper) ? "U" : "L")) + " #tabu neighbors (+itself) = " << setAllSols.size() << "\n\n";
        }
        if (_gp->isVerbose())
            std::cout << "Excluding each solution that is not affordable...\n";
        for (auto& b : setAllSols)
            if (b.any() && isBudgetEnough(b)) {
                auto it = _tabuMap.find(b);
                if (it == _tabuMap.end() || it->second == false)
                   setSols.insert(b);
            }
        setSols.insert(curSol);
        if (_gp->isVerbose())
            std::cout << "Finalized " + std::string(((_isUpper) ? "U" : "L")) + " #tabu neighbors (+itself) = " << setSols.size() << "\n\n";

        if (_gp->isVerbose()) {
            int j = 0;
            for (auto& i : setSols) {
                if (j >= 3 && _gp->debugLevel() < 2) {
                    std::cout << "...\n";
                    break;
                }
                std::cout << "[" << ++j << "] = " << i << std::endl;
            }

            std::cout << "(Verbose) Initial solution = " << curSol << "\n";
            std::cout << "(Verbose) Originally found total " + std::string(((_isUpper) ? "U" : "L")) + " #tabu neighbors (+itself) = " << setAllSols.size() << "\n";
            std::cout << "(Verbose) Finalized " + std::string(((_isUpper) ? "U" : "L")) + " #tabu neighbors = " << setSols.size() << "\n\n";
        }
        if (!_neighborSetMap.count(curSol))
            _neighborSetMap[curSol] = new std::unordered_set<PartSolution>(setSols);
    }

    return (int) setSols.size();
}

FullSolution Layer::tabuSearch(Layer* aLayer, const FullSolution& curFsol) {

    FuncLog f(_gp, std::string(((_isUpper) ? "U" : "L")) + "TabuSearch", "\t\t");
    if (1 || _isUpper) {
        f.startLog();
    }

    error_if_not(aLayer, "aLayer is a nullptr");
    clock_t start = clock();
    std::unordered_set<FullSolution> setBestSols;
    std::unordered_set<PartSolution> setSols;
    setBestSols.insert(curFsol);
    PartSolution max_my;
    _tabuList = std::queue<PartSolution>();
    int g = 0, l = 1, h = 0;
    bool isFirst = true;
    FullSolution iterBestFsol, lastIterBestFsol, firstIterFSol;
    std::cout << std::string(((_isUpper) ? "U" : "L")) << "TabuSearch" << ", initFSol = " << curFsol  << std::endl;
    while (g < _tabu->gmax()) {
        _tabuIter++;
        clock_t istart = clock();

        if (isFirst) {
            firstIterFSol = lastIterBestFsol = iterBestFsol = curFsol;
            isFirst = false;
        }

        const PartSolution& x = iterBestFsol.second;
        const PartSolution& y = iterBestFsol.first;
        const PartSolution& my = (_isUpper) ? y : x;
        const PartSolution& his = (_isUpper) ? x : y;
        PartSolution max_his = his;

        std::cout << "one loop start" << std::endl;
        findNeighbors(my, setSols, 0, false);

        if (_gp->isVerbose()) {
            int i = 0;
            for (auto& sol : setSols) {
                if (i >= 3 && _gp->debugLevel() < 2) {
                    std::cout << "...\n";
                    break;
                }
                FullSolution fs = (_isUpper) ? std::make_pair(sol, his) : std::make_pair(his, sol);
                std::cout << "sol[" << ++i << "] = " << fs << ", ulRev = " <<
                             revenue(aLayer, fs) << ", llRev = " << aLayer->revenue(this, fs) << std::endl;
            }
        }

        error_if_not(setSols.size() > 0, "the set of solutions is empty!");
        std::cout << std::string(((_isUpper) ? "Upper" : "Lower")) << " tabu iter l = " << l << ", found #neighbor = " << setSols.size() << std::endl;
        int useMap = 0, calcu = 0;
        // Upper tabu
        if (_isUpper) {
            std::unordered_set<FullSolution> setNewSols;
            for (auto& sol : setSols) {
                std::cout << "neighbor sol = " << sol << ", curFSol = " << curFsol  << std::endl;
                FullSolution fs = std::make_pair(sol, his);
                bool isMapFound = false;
                FullSolution mapSol;
                FullSolution newSol;
                if (_gp->isCache()) {
                    auto it = _tabuBestMap.find(sol);
                    if (it != _tabuBestMap.end()) {
                        isMapFound = true;
                        mapSol = std::make_pair(sol, it->second);
                        if (dist(it->second, his) < _gp->bitDiff()) {
                            useMap++;
                            newSol = std::make_pair(sol, it->second);
                        }
                        else {
                            calcu++;
                            newSol = aLayer->tabuSearch(this, fs);
                        }
                        /*
                         *
                        if (sol == my) {
                            FullSolution mapSol = std::make_pair(sol, it->second);
                            FullSolution hisSol = std::make_pair(my, his);
                            warning_msg("x should be the same for an old y");
                            std::cout << "y = " << sol << ", x in map = " << it->second
                                      << ", ulRev = " << revenue(aLayer, mapSol)
                                      << ", llRev = " << aLayer->revenue(this, mapSol)
                                      << std::endl;

                            std::cout << "y = " << sol << ", x last time = " << his
                                      << ", ulRev = " << revenue(aLayer, hisSol)
                                      << ", llRev = " << aLayer->revenue(this, hisSol)
                                      << std::endl;
                        } */
                        if (it->second != his && _gp->isVerbose()) {
                            FullSolution mapSol = std::make_pair(sol, it->second);
                            FullSolution hisSol = std::make_pair(sol, his);
                            warning_msg("x should be the same for a new y");
                            std::cout << "y = " << sol << ", x in map = " << it->second
                                      << ", ulRev = " << revenue(aLayer, mapSol)
                                      << ", llRev = " << aLayer->revenue(this, mapSol)
                                      << std::endl;

                            std::cout << "y = " << sol << ", x last time = " << his
                                      << ", ulRev = " << revenue(aLayer, hisSol)
                                      << ", llRev = " << aLayer->revenue(this, hisSol)
                                      << std::endl;
                        }
                    }
                    else {
                        calcu++;
                        newSol = aLayer->tabuSearch(this, fs);
                    }
                }
                else {
                    calcu++;
                    newSol = aLayer->tabuSearch(this, fs);
                }
                error_if_not(fs.first == newSol.first, "y should be the same!");
                if (!isMapFound || (isMapFound && newSol != mapSol &&
                                    aLayer->revenue(this, newSol) >
                                    aLayer->revenue(this, mapSol)))
                    _tabuBestMap[newSol.first] = newSol.second;
                setNewSols.insert(newSol);
            }
            error_if_not(setNewSols.size() > 0, "the set of solutions is empty!");
            FullSolution max_fsol = *std::max_element(setNewSols.begin(), setNewSols.end(),
                             [&](const FullSolution& fa, const FullSolution& fb) -> bool {
                return revenue(aLayer, fa) < revenue(aLayer, fb);
            });
            max_my = max_fsol.first;
            max_his = max_fsol.second;
        } else
            max_my = *std::max_element(setSols.begin(), setSols.end(),
                                       [&](const PartSolution& a, const PartSolution& b) -> bool {
                FullSolution fa = (_isUpper) ? std::make_pair(a, his) : std::make_pair(his, a);
                FullSolution fb = (_isUpper) ? std::make_pair(b, his) : std::make_pair(his, b);
                return revenue(aLayer, fa) < revenue(aLayer, fb);
            });
        if (_isUpper)
        std::cout << "total #neighbor = " << setSols.size() << ", #cache = " << useMap << ", #tabuCal = " << calcu
                  << std::endl;
        warning_if_not(max_my.count() >= my.count(),
                       "best #selected locations = " + std::to_string(max_my.count()) +
                       " < original #selected locations = " + std::to_string(my.count()));

        FullSolution fsol_max = (_isUpper) ? std::make_pair(max_my, max_his) : std::make_pair(max_his, max_my);
        FullSolution fsol_my = (_isUpper) ? std::make_pair(my, his) : std::make_pair(his, my);

        if (_gp->isVerbose() && !_isUpper && revenue(aLayer, fsol_max) < revenue(aLayer, fsol_my)) {
            FullSolution mYLX = (_isUpper) ? std::make_pair(max_my, his) : std::make_pair(his, max_my);
            FullSolution mXLY = (_isUpper) ? std::make_pair(my, max_his) : std::make_pair(max_his, my);

            std::cout << "last iter sol (last x last y) = " << iterBestFsol << ", rev = " << revenue(aLayer, iterBestFsol)
                      << ", aRev = " << aLayer->revenue(this, iterBestFsol) << std::endl;
            std::cout << "last iter fsol_my             = " << fsol_my << ", rev = " << revenue(aLayer, fsol_my)
                      << ", aRev = " << aLayer->revenue(this, fsol_my) << std::endl;
            std::cout << "max y last x                  = " << mYLX << ", rev = " << revenue(aLayer, mYLX)
                      << ", aRev = " << aLayer->revenue(this, mYLX) << std::endl;
            std::cout << "max x last y                  = " << mXLY << ", rev = " << revenue(aLayer, mXLY)
                      << ", aRev = " << aLayer->revenue(this, mXLY) << std::endl;
            std::cout << "max x max y                   = " << fsol_max << ", rev = " << revenue(aLayer, fsol_max)
                      << ", aRev = " << aLayer->revenue(this, fsol_max) << std::endl;
        }

        if (!_isUpper)
            error_if(revenue(aLayer, fsol_max) < revenue(aLayer, fsol_my), "max revenue = " + std::to_string(revenue(aLayer, fsol_max)) +
                 " < current revenue = " + std::to_string(revenue(aLayer, fsol_my)));
        if (_gp->isNeighborCache() &&
                ((!_isUpper && max_my == my) ||
                 (_isUpper && fsol_max == fsol_my) ||
                 revenue(aLayer, fsol_max) == revenue(aLayer, fsol_my))) {
            warning_msg("best solution is equal to last iter solution! leave the loop!");
            std::cout << "last iter = " << iterBestFsol << ", rev = " << revenue(aLayer, iterBestFsol)
                      << ", aRev = " << aLayer->revenue(this, iterBestFsol) << std::endl;
            std::cout << "this iter = " << fsol_max << ", rev = " << revenue(aLayer, fsol_max)
                      << ", aRev = " << aLayer->revenue(this, fsol_max) << std::endl;
            break;
        }

        if (revenue(aLayer, fsol_max) > revenue(aLayer, iterBestFsol))
            iterBestFsol = fsol_max;

        if ((int )_tabuList.size() == _tabu->size()) {
            _tabuMap[_tabuList.front()] = false;
            _tabuList.pop();
        }
        _tabuList.push(max_my);
        _tabuMap[max_my] = true;
        FullSolution tmpFsol = iterBestFsol;
        if (l > 1) {
            if (revenue(aLayer, iterBestFsol) < revenue(aLayer, lastIterBestFsol)*1.01) {
                std::cout << "improve less than 1% => g = g+1" << std::endl;
                g++;
                h++;

            }
            else {
                std::cout << "improve more than 1% => g = 0 and h = 0" << std::endl;
                h = 0;
                g = 0;
            }
            if (!_isUpper && h >= _tabu->lmax()) {
                tmpFsol = changeLocation(iterBestFsol, _kMax);
                h = 0;
            }
        }

        // if (_gp->isVerbose()) {
            std::cout << std::string(((_isUpper) ? "Upper" : "Lower")) << " tabu now " << "g = " << g
                      << ", l = " << l << ", h = " << h << std::endl;
            std::cout << "lastIterFSol = " << lastIterBestFsol<< ", myRev = "
                  << revenue(aLayer, lastIterBestFsol) << ", hisRev = " << aLayer->revenue(this, lastIterBestFsol)
                  << std::endl;
            std::cout << "iterBestSol = " << iterBestFsol << ", myRev = "
                  << revenue(aLayer, iterBestFsol) << ", hisRev = " << aLayer->revenue(this, iterBestFsol)
                  << std::endl;
            if (!_isUpper)
                std::cout << "tmpFsol = " << tmpFsol << ", myRev = "
                      << revenue(aLayer, tmpFsol) << ", hisRev = " << aLayer->revenue(this, tmpFsol)
                      << std::endl;
            _gp->showElapsedTime();
             std::cout << "one loop end" << std::endl << std::endl;
        // }
        setBestSols.insert(iterBestFsol);
        lastIterBestFsol = iterBestFsol;
        if (l == 1) firstIterFSol = iterBestFsol;

        if (l > 1 && !_isUpper) iterBestFsol = tmpFsol;

        double curUsedTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        double iterTime = (double) (clock() - istart) / CLOCKS_PER_SEC;

        l++;

        _maxTabuRuntime = std::max(iterTime, _maxTabuRuntime);
        _minTabuRuntime = std::min(iterTime, _minTabuRuntime);

        if (curUsedTime > _tabu->itmax()) {
            warning_msg("curUsedTime = " + std::to_string(curUsedTime) + " > max = " + std::to_string(_tabu->itmax()) + " sec., terminated.");
            break;
        }
    }

    if (_gp->isVerbose())
        std::cout << std::string(((_isUpper) ? "Upper" : "Lower")) << " tabu used #iterations = " << l << std::endl;

    FullSolution max_fsol = *std::max_element(setBestSols.begin(), setBestSols.end(),
      [&](const FullSolution& fa, const FullSolution& fb) -> bool {
      return revenue(aLayer, fa) < revenue(aLayer, fb);
    });

    double tabuTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    _tabuRuntime += tabuTime;

    std::cout << "================================" << std::endl;
    std::cout << std::string(((_isUpper) ? "Upper" : "Lower")) << " tabu "
              << "final total g = " << g << ", total l = " << l << ", total time = " << tabuTime
              << " sec." << std::endl;
    std::cout << "input Sol = " << curFsol << ", myRev = "
          << revenue(aLayer, curFsol) << ", hisRev = " << aLayer->revenue(this, curFsol)
          << std::endl;
    std::cout << "find maxSol = " << max_fsol << ", myRev = "
          << revenue(aLayer, max_fsol, true, _isDebug) << ", hisRev = " << aLayer->revenue(this, max_fsol, true, _isDebug)
          << std::endl;
    std::cout << std::string(((_isUpper) ? "Upper" : "Lower")) << " all loops end" << std::endl;
    _gp->showElapsedTime();
    std::cout << "================================" << std::endl << std::endl;
    _isDebug = false;

    return max_fsol;
}

FullSolution Layer::changeLocation(const FullSolution& curFsol, const int& k) {
    FuncLog f(_gp, "changeLocation", "\t\t");
    f.startLog();


    FullSolution fsol = curFsol;
    PartSolution& sol = (!_isUpper) ? fsol.second : fsol.first;

    if (sol.count() == 0) return fsol;

    int numNSS = std::ceil(_vecNSS[k]*sol.count());

    // std::ceil(_vecNSS[k]*sol.size());
    if (numNSS <= 0) {
        warning_msg("#change bit = 0, auto increase to 1");
        numNSS = 1;
    }
    // std::cout << "enter changeLocation: #change 1 = " << numNSS << std::endl;
    error_if_not(isBudgetEnough(sol), "budget is not enough before change location");
    int iter = 0;
    while (fsol == curFsol || !isBudgetEnough(sol)) {
        std::vector<size_t> vec1Loc, vec0Loc;
        PartSolution flipSol = sol;
        flipSol.flip();

        // std::cout << "sol = " << sol << std::endl;
        size_t i = sol.find_first();
        while (i != PartSolution::npos){
            vec1Loc.push_back(i);
            i = sol.find_next(i);
        }

        /*
        for (auto & i : vec1Loc)
            std::cout << i << " ";
        std::cout << std::endl;
        */
        std::random_shuffle(vec1Loc.begin(), vec1Loc.end());
        /*
        for (auto & i : vec1Loc)
            std::cout << i << " ";
        std::cout << std::endl;
        */
        // std::cout << "flipped sol = " << flipSol << std::endl;
        size_t j = flipSol.find_first();
        while (j != PartSolution::npos){
            vec0Loc.push_back(j);
            j = flipSol.find_next(j);
        }

        /*
        for (auto & i : vec0Loc)
            std::cout << i << " ";
        std::cout << std::endl;
        */

        std::random_shuffle(vec0Loc.begin(), vec0Loc.end());
        /*
        for (auto & i : vec0Loc)
            std::cout << i << " ";
        std::cout << std::endl;
        */
        warning_if_not((int)vec0Loc.size() >= numNSS, "#0's = " + std::to_string(vec0Loc.size()) +
                       " is not enough for #nss = " + std::to_string(numNSS));
        warning_if_not((int)vec1Loc.size() >= numNSS, "#1's = " + std::to_string(vec1Loc.size()) +
                       " is not enough for #nss = " + std::to_string(numNSS));
        for (int i = 0, sz = std::min((int)vec0Loc.size(), numNSS); i < sz; i++){
            int j = vec0Loc[i];
            // std::cout << sol << " " << j << " " << sol.test(j) << std::endl;
            error_if_not(!sol.test(j), "position " + std::to_string(j) + " is not 0");
            sol.flip(j);
        }
        for (int i = 0, sz = std::min((int)vec1Loc.size(), numNSS); i < sz; i++){
            int j = vec1Loc[i];
            // std::cout << sol << " " << j << " " << sol.test(j) << std::endl;
            error_if_not(sol.test(j), "position " + std::to_string(j) + " is not 1");
            sol.flip(j);
        }
        // std::cout << "original fsol = " << curFsol << ", changed = " << fsol << std::endl;
        iter++;
        if (iter > 1000) {
            warning_msg("Try change location more than 500 failures!");
            fsol = curFsol;
            break;
        }
    }

    std::cout << "#Tries of change location = " << iter << std::endl;

    return fsol;
}

double Layer::estUtility(Facility* fac) {
    error_if_not(fac, "facility is a nullptr");
    double total = 0;
    // for (auto& fam : _vecFamily)
    //    total += fam->bp()*(W1*utility(fac, fam, WEEKDAY) + W2*utility(fac, fam, WEEKEND));
    for (auto& fam : _vecFamily)
        total += fam->dist(fac);
    // fam->bp()*(W1*utility(fac, fam, WEEKDAY) + W2*utility(fac, fam, WEEKEND));
    return total;
}

void Layer::allocateSolutions() {
    warning_if_not(_numFacility > 0, _name + " #candidateFacility should be a positive integer instead of " + std::to_string(_numFacility));
    _cur_sol = _prev_sol = _best_sol = PartSolution(_numFacility, 0);
}

void Layer::initSolutions() {

    std::vector<Facility*> vecSortedFacility = _vecFacility;

    for (auto& fac : vecSortedFacility)
        fac->setEstUtil(estUtility(fac));

    std::sort(vecSortedFacility.begin(), vecSortedFacility.end(),
    [](const Facility* a, const Facility* b) {
        // sort by cp-value
        // return a->estUtil()/a->cost() > b->estUtil()/b->cost();

        // sort by est utility
        return a->estUtil() < b->estUtil();
    });

    if (_gp->debugLevel() >= 2)
        for (auto& fac : vecSortedFacility)
           fac->show();

    double remain = _budget;
    for (auto& fac : vecSortedFacility) {
        if (fac->cost() > remain) continue;
        else {
            remain -= fac->cost();
            _cur_sol[fac->id()] = 1;
        }
    }
    _prev_sol = _best_sol = _cur_sol;
}

void Layer::random(PartSolution& sol) {
    for (int i = 0; i < _numFacility; ++i)
        sol[i] = _gp->_ri.Ui(0,1);
    return;
}

/*
void Layer::randomTry(const int& n) {
    int t = 0;
    while (++t <= n) {
        random(_cur_sol);
        double r = 0.0;
        bool isBE = isBudgetEnough(_cur_sol, r);
        if (_gp->debugLevel() > 0)
            std::cout << "Try " << t << ": cur_sol = " << _cur_sol << ", isBudgetEnough = "
                      << std::boolalpha << isBE << ", remain = " << r << std::endl;
        if (isBE)
            std::cout << "Revenue = " << revenue << std::endl;

    }
}
*/
bool Layer::isBudgetEnough(const PartSolution& sol, double& remains) {

    if (_gp->isCache()) {
        auto it = _budgetCache.find(sol);
        if (it != _budgetCache.end())
            return it->second;
    }

    double total = _budget;
    if (true){
        size_t i = sol.find_first();
        while(i != PartSolution::npos){
            total -= _vecFacility[i]->cost();
            i = sol.find_next(i);
        }
    }
    else {
        for (int i = 0; i < _numFacility ; ++i) {
            // if (_gp->debugLevel() > 0)
            //    std::cout << "i = " << i << ", v = " << _vecFacility[i]->cost() << std::endl;
            if (sol.test(i))
                total -= _vecFacility[i]->cost();
        }
    }
    remains = total;

    bool isBudgetEnough = (total >= 0);
    _budgetCache[sol] = isBudgetEnough;
    // if (_gp->debugLevel() > 0) std::cout << "remains = " << remains << std::endl;
    return isBudgetEnough;
}

bool Layer::isBudgetEnough(const PartSolution& sol) {
    double tmp = 0;
    return isBudgetEnough(sol, tmp);
}

void Layer::randomGenerateOrAdjust() {
    if (_budget <= 0) {
        _budget = _gp->_rr.Ur(500.0, 3000.0);
        if (_gp->isVerbose()) info_msg(_name + " randomed budget = " + std::to_string(_budget));
    }
    // paper: _gp->_rr.Ur(100.0, 125000.0);
    if (_gMax <= 0) {
        _gMax = (_isUpper) ? _gp->_ri.Ui(5) : 1;
        if (_gp->isVerbose()) info_msg(_name + " randomed gmax = " + std::to_string(_gMax));
    }
    if (_kMax <= 0) {
        _kMax = _gp->_ri.Ui(5);
        if (_gp->isVerbose()) info_msg(_name + " randomed kmax = " + std::to_string(_kMax));
    }
    if ((int) _vecNSS.size() != _kMax+1) {
        _vecNSS.resize(_kMax+1);
        if (_vecNSS[0] != 0) {
            _vecNSS[0] = 0;
            if (_gp->isVerbose()) info_msg(_name + " set _vecNSS[0] = 0");
        }
        double max_step = 1.0/_vecNSS.size();
        for (int k = 1; k <= _kMax; ++k)
            if (_vecNSS[k] <= _vecNSS[k-1] || _vecNSS[k] >= 1.0) {
                _vecNSS[k] = _vecNSS[k-1] + _gp->_rr.Ur(0.2*std::fmin(max_step, 1.0-_vecNSS[k-1])); ;
                if (_gp->isVerbose()) info_msg(_name + " randomed NSS[" + std::to_string(k) + "] = " +
                            std::to_string(_vecNSS[k]));
            }
    }
    if (_lMax <= 0) {
        _lMax = _gp->_ri.Ui(5);
        if (_gp->isVerbose()) info_msg(_name + " randomed lmax = " + std::to_string(_lMax));
    }
    if (_itMax <= 0) {
        _itMax = _gp->_rr.Ur(20, 30);
        if (_gp->isVerbose()) info_msg(_name + " randomed itmax = " + std::to_string(_itMax));
    }
    if (_numFacility <= 0) {
        // _numFacility = (_isHyper) ? _gp->_ri.Ui(5, 15) : _gp->_ri.Ui(15, 35);
        _numFacility = (_isHyper) ? _gp->_ri.Ui(5) : _gp->_ri.Ui(5);
        if (_gp->isVerbose()) info_msg(_name + " randomed #candidate facility = " + std::to_string(_numFacility));
    }
    if (_numExistedFacility < 0) {
        _numExistedFacility = _gp->_ri.Ui(0, 5);
        if (_gp->isVerbose()) info_msg(_name + " randomed #existed facility = " + std::to_string(_numExistedFacility));
    }
    int id = 0;
    if ((int) _vecFacility.size() != _numFacility) {
        _vecFacility.resize(_numFacility);
        for (auto& fac : _vecFacility)
            if (!fac)
                fac = new Facility(_gp->_rr.Nr(50), _gp->_rr.Nr(50), _gp->_rr.Ur(100.0, 500.0), _gp->_cp[_mt], _gp->_pa[_mt], id++);
            else if (fac->cost() <= 0) {
                fac->setCost(_gp->_rr.Ur(100.0, 500.0));
                if (_gp->isVerbose()) info_msg(_name + " randomed facility[" + std::to_string(fac->id()) + "].cost() = " +
                            std::to_string(fac->cost()));
            }
    }
    if ((int) _vecExistedFacility.size() != _numExistedFacility) {
        _vecExistedFacility.resize(_numExistedFacility);
        for (auto& fac : _vecExistedFacility)
            if (!fac)
                fac = new Facility(_gp->_rr.Nr(50), _gp->_rr.Nr(50), -1, _gp->_cp[_mt], _gp->_pa[_mt]);

    }

    if (_cur_sol.size() == 0 || _best_sol.size() == 0 || _prev_sol.size() == 0)
        allocateSolutions();

    if (!_tabu) _tabu = new Tabu();
    if (_tabu->_gMax <= 0) {
        _tabu->_gMax = _gp->_ri.Ui(5);
        if (_gp->isVerbose()) info_msg(_name + " randomed tabu gmax = " + std::to_string(_tabu->_gMax));
    }
    if (_tabu->_lMax <= 0) {
        _tabu->_lMax = _gp->_ri.Ui(50);
        if (_gp->isVerbose()) info_msg(_name + " randomed tabu lmax = " + std::to_string(_tabu->_lMax));
    }
    if (_tabu->_cMax <= 0) {
        _tabu->_cMax = _gp->_ri.Ui(100);
        if (_gp->isVerbose()) info_msg(_name + " randomed tabu cmax = " + std::to_string(_tabu->_cMax));
    }
    if (_tabu->_size <= 0) {
        _tabu->_size = _gp->_ri.Ui(20);
        if (_gp->isVerbose()) info_msg(_name + " randomed tabu list size = " + std::to_string(_tabu->_size));
    }
    if (_tabu->_itMax <= 0) {
        _tabu->_itMax = (!_isUpper) ? 5 : _gp->_rr.Ur(5,20);
        if (_gp->isVerbose()) info_msg(_name + " randomed tabu itmax = " + std::to_string(_tabu->_itMax));
    }

}

double Layer::utility(Facility* fac, Family* fam, const DAY_TYPE& dt, bool isDebug) {
    error_if_not(fac, "facility is a nullptr");
    error_if_not(fam, "family is a nullptr");
    error_if_not(_gp, "global parameter is a nullptr");

    if (_gp->isCache()) {
        auto it = _utilCache.find(fac);
        if (it != _utilCache.end()) {
            auto it2 = it->second.find(fam);
            if (it2 != it->second.end()) {
                auto it3 = it2->second.find(dt);
                if (it3 != it2->second.end()) {
                    // if (_gp->debugLevel() > 0)
                    //    std::cout << "found cache untility: " << it3->second << std::endl;
                    return it3->second;
                }
            }
        }
    }
    const double &i1 = std::pow(_gp->at()[mt()][fam->type()][dt], _gp->theta()[1]),
            &i2 = std::pow(1+_gp->pa()[mt()], _gp->theta()[2]),
            &i3 = std::exp(_gp->theta()[3]*_gp->cp()[mt()]),
            &i4 = std::exp(_gp->theta()[4]*fac->dist(fam));
    const double &util = (i1*i2)/(i3*i4);
    const double &exputil = std::exp(util);


    // double util = std::pow(_gp->at()[mt()][fam->type()][dt], _gp->theta()[1])/
    //        std::exp(_gp->theta()[3]*_gp->cp()[mt()])*
    //        std::pow(1+_gp->pa()[mt()], _gp->theta()[2])/
    //        std::exp(_gp->theta()[4]*fac->dist(fam));
    if (isDebug) {
        std::cout << "------------------------------" << std::endl;
        fac->show();
        fam->show();
        std::cout << "Day " << _dtName.at(dt) << ", isUpper = " << _isUpper << std::endl;
        std::cout << "util = " << util << ", i1 = " << i1
                << ", i2 = " << i2 << ", i3 = " << i3
                << ", i4 = " << i4 << ", i3*i4 = " << i3*i4
                << ", exp(util) = " << std::exp(util) << std::endl;
        std::cout << "------------------------------" << std::endl;
    }
    _utilCache[fac][fam][dt] = exputil;
    return exputil;
}

double Layer::revenue(Layer* aLayer, const FullSolution& curFsol, bool isTest, bool isDebug) {
    error_if_not(aLayer, "aLayer is a nullptr");

    const PartSolution& x = curFsol.second;
    const PartSolution& y = curFsol.first;

    const PartSolution& my = (_isUpper) ? y : x;
    const PartSolution& his = (_isUpper) ? x : y;

    if (isTest){
        bool isBE = true;
        double r;
        if (!isBudgetEnough(my, r)) {
            std::string str;
            to_string(my, str);
            // if (_gp->isVerbose())
                warning_msg(_name + " budget is not enough, cur_sol = " + str + ", budget = " +
                        std::to_string(_budget) + ", remain = " + std::to_string(r));
            isBE = false;
        }
        if (_gp->debugLevel() >= 3) std::cout << "r = " << r << std::endl;
        if (!aLayer->isBudgetEnough(his, r)) {
            std::string str;
            to_string(his, str);
            // if (_gp->isVerbose())
                warning_msg(aLayer->name() + " budget is not enough, cur_sol = " + str + ", budget = " +
                        std::to_string(aLayer->budget()) + ", remain = " + std::to_string(r));
            isBE = false;
        }
        if (_gp->debugLevel() >= 3) std::cout << "r2 = " << r << std::endl;
        if (!isBE) return -1;
    }

    double rev = 0.0, rcrev = 0.0, aRev = 0.0;

    if (_gp->isCache()) {
        auto it = _revenueCache.find(curFsol);
        if (it != _revenueCache.end()) {
            if (isDebug)
                rcrev = it->second;
            else
                return it->second;
        }
    }
    if (isDebug)
        std::cout << "Enter " << ((_isUpper) ? "U" : "L") << " revenue!\n";

    int k = 0;
    for (auto& fam : _vecFamily) {
        int t = 0;
        double all_week_rev = 0.0, all_week_arev = 0.0, alz1 = 0.0, alz2 = 0.0, alz3 = 0.0;
        for (auto& dt : _vecDt) {
            double teu = totalExistedUtilityPerFamily(fam, dt, isDebug),
                   aTeu = aLayer->totalExistedUtilityPerFamily(fam, dt, isDebug);
            double tcu = 0.0, aTcu = 0.0;
            // for (int i = 0; i < _numFacility ; ++i)
            //    if (my.test(i))
            size_t i = my.find_first();
            while (i != PartSolution::npos) {
                double subutil = utility(_vecFacility[i], fam, dt, isDebug);

                tcu += subutil;
                if (isDebug && k < 10) {
                    std::cout << "sub util = " << subutil << std::endl;
                }
                i = my.find_next(i);
            }
            // for (int i = 0; i < aLayer->numFacility() ; ++i)
            //    if (aLayer->cur_sol().test(i))
            size_t j = his.find_first();
            while (j != PartSolution::npos) {
                double subutil = aLayer->utility(aLayer->vecFacility()[j], fam, dt, isDebug);
                aTcu += subutil;
                if (isDebug && k < 10) {
                    std::cout << "sub util = " << subutil << std::endl;
                }
                j = his.find_next(j);
            }
            const double& wt = (dt == WEEKDAY) ? W1 : W2;
            const double& nume = (tcu+teu), &anume = (aTcu+aTeu), &deno = (tcu+teu+aTcu+aTeu);
            const double& subrev = wt*nume/deno;
            const double& subarev = wt*anume/deno;
            const double& subz1 = wt*tcu/deno, &subz2 = wt*anume/deno, &subz3 = wt*teu/deno;


            all_week_rev += subrev;
            all_week_arev += subarev;
            alz1 += subz1;
            alz2 += subz2;
            alz3 += subz3;
            // if (_gp->debugLevel() > 0)
            //    std::cout << "bp = " << fam->bp() << ", rev = " << rev << std::endl;
            if (isDebug && k < 10) {
                // this->show();
                std::cout << "t = " << t << ", dt = " << _dtName.at(dt)
                          << ", fam[" << k << "].bp = " << fam->bp()
                          << ", tcu = " << tcu << ", teu = " << teu << ", aTcu = " << aTcu
                          << ", aTeu = " << aTeu << ", nume (tcu+teu) = " << nume
                          << ", anume (aTcu+aTeu) = " << anume
                          << ", deno (tcu+teu+aTcu+aTeu) = " << deno
                          << ", frac (nume/deno) = " << nume/deno
                          << ", afrac (anume/deno) = " << anume/deno
                          << ", wt = " << wt
                          << ", subrev (bp*wt*frac) = " << subrev
                          << ", subarev (bp*wt*afrac) = " << subarev
                          << std::endl;
                std::cout << "sub z1 = " << subz1
                          << ", sub z2 = " << subz2
                          << ", sub z3 = " << subz3
                          << ", sub z1+z3 = " << subz1+subz3
                          << std::endl << std::endl;
            }
            t++;
        }
        rev += fam->bp()*all_week_rev;
        aRev += fam->bp()*all_week_arev;
        if (isDebug && k < 10) {
            // this->show();
            std::cout << "+++++++++\n";
            std::cout << "fam[" << k << "].bp = " << fam->bp()
                      << ", 7 days z1 = " << alz1 << ", z2 = " << alz2 << ", z3 = " << alz3
                      << ", z1+13 = " << alz1 + alz3
                      << ", 7 days rev = " << all_week_rev
                      << ", 7 days arev = " << all_week_arev
                      << ", rev = " << rev << "arev = " << aRev << std::endl;
            std::cout << "+++++++++\n";
        }
        if (isDebug) isDebug = false;
        k++;
    }
    /*
    for (auto& fam : _vecFamily) {
        double teu = totalExistedUtilityPerFamily(fam, isDebug),
               aTeu = aLayer->totalExistedUtilityPerFamily(fam, isDebug);
        double tcu = 0.0, aTcu = 0.0;
        // for (int i = 0; i < _numFacility ; ++i)
        //    if (my.test(i))
        size_t i = my.find_first();
        while (i != PartSolution::npos) {
            double weekday = utility(_vecFacility[i], fam, WEEKDAY, isDebug),
                   weekend = utility(_vecFacility[i], fam, WEEKEND, isDebug);

            tcu += W1*weekday+W2*weekend;
            if (isDebug && k < 10) {
                std::cout << "my weekday = " << weekday << ", W1 = " << W1 << ", weekend = " << weekend << ", W2 = " << W2 << std::endl;
                std::cout << "my subtotal = W1*weekday+W2*weekend = " << W1*weekday+W2*weekend << std::endl;
            }
            i = my.find_next(i);
        }
        // for (int i = 0; i < aLayer->numFacility() ; ++i)
        //    if (aLayer->cur_sol().test(i))
        size_t j = his.find_first();
        while (j != PartSolution::npos) {
            double weekday = aLayer->utility(aLayer->vecFacility()[j], fam, WEEKDAY, isDebug),
                   weekend = aLayer->utility(aLayer->vecFacility()[j], fam, WEEKEND, isDebug);
            aTcu += W1*weekday+W2*weekend;
            if (isDebug && k < 10) {
                std::cout << "competitor weekday = " << weekday << ", W1 = " << W1 << ", weekend = " << weekend << ", W2 = " << W2 << std::endl;
                std::cout << "competitor subtotal = W1*weekday+W2*weekend = " << W1*weekday+W2*weekend << std::endl;
            }
            j = his.find_next(j);
        }
        rev += fam->bp()*(tcu+teu)/(tcu+teu+aTcu+aTeu);
        aRev += fam->bp()*(aTcu+aTeu)/(tcu+teu+aTcu+aTeu);

        if (isDebug && k < 10) {
            // this->show();
            std::cout << "fam[" << k << "] = " << fam->bp() << ", tcu = " << tcu
                      << ", teu = " << teu << ", aTcu = " << aTcu
                      << ", aTeu = " << aTeu << ", nume = " << (tcu+teu) << ", anume = " << (aTcu+aTeu)
                      << ", deno = " << (tcu+teu+aTcu+aTeu) << ", frac = " << (tcu+teu)/(tcu+teu+aTcu+aTeu)
                      << ", afrac = " << (aTcu+aTeu)/(tcu+teu+aTcu+aTeu)
                      << ", bp = " << fam->bp() << ", total = "
                      << fam->bp()*(tcu+teu)/(tcu+teu+aTcu+aTeu)
                      << ", atotal = " << fam->bp()*(aTcu+aTeu)/(tcu+teu+aTcu+aTeu)
                      << std::endl;
            std::cout << "z1 = " << (tcu)/(tcu+teu+aTcu+aTeu)
                      << ", z2 = " << (aTcu+aTeu)/(tcu+teu+aTcu+aTeu)
                      << ", z3 = " << (teu)/(tcu+teu+aTcu+aTeu) << std::endl << std::endl;
        }


        // if (_gp->debugLevel() > 0)
        //    std::cout << "bp = " << fam->bp() << ", rev = " << rev << std::endl;
        k++;
    }
    */
    // if (_gp->debugLevel() > 0)
    //    std::cout << ", rev = " << rev << std::endl;

    if (isDebug)
        std::cout << "\nsol = " << curFsol << ", rev = " << rev << ", aRev = " << aRev << ", rcrev = " << rcrev
                  << std::endl;
    error_if(rcrev > 0.0 && rev != rcrev, "rev recorded is not equal to calculated");

    _revenueCache[curFsol] = rev;
    aLayer->revenueCache()[curFsol] = aRev;

    return rev;
}

double Layer::totalExistedUtilityPerFamily(Family* fam, const DAY_TYPE& dt, bool isDebug) {
    error_if_not(fam, "family is a nullptr");
    // std::cout << "is cache = " << _gp->isCache() << std::endl;
    if (_gp->isCache()) {
        // std::cout << "is cache " << std::endl;
        auto it = _utilPerFamCache.find(fam);
        if (it != _utilPerFamCache.end()) {
            auto it2 = it->second.find(dt);
            if (it2 != it->second.end()) {
                // if (_gp->debugLevel() > 0)
                //    std::cout << "found cache existed untility: " << it->second << std::endl;
                return it2->second;
            }
        }
    }

    double total = 0.0;
    for (auto& efac : _vecExistedFacility) {
        double subutil = utility(efac, fam, dt, isDebug);
        if (isDebug) {
            std::cout << "exist subutil = " << subutil << std::endl;
        }
        total += subutil;
    }
    /*
    for (auto& efac : _vecExistedFacility) {
        double weekday = utility(efac, fam, WEEKDAY, isDebug),
               weekend = utility(efac, fam, WEEKEND, isDebug);

        if (isDebug) {
            std::cout << "my exist weekday = " << weekday << ", W1 = " << W1 << ", exist weekend = " << weekend << ", W2 = " << W2 << std::endl;
            std::cout << "my exist subtotal = W1*weekday+W2*weekend = " << W1*weekday+W2*weekend << std::endl;
        }
        total += W1*weekday+W2*weekend;
    }
    // total += W1*utility(efac, fam, WEEKDAY, isDebug) + W2*utility(efac, fam, WEEKEND, isDebug);
    */
    _utilPerFamCache[fam][dt] = total;
    return total;
}
