//#include "../include/gurobi_c++.h"
#include "util.h"
#include "location.h"
#include "family.h"
#include "facility.h"
#include "layer.h"
#include "parser.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>



#ifdef USE_QT
#include "mainwindow.h"
#include <QApplication>
#endif

void testPointInsidePolygon() {
    FuncLog f(nullptr, "testPointInsidePolygon", "\t", true);
    f.startLog();

    typedef boost::geometry::model::d2::point_xy<double> point_type;
    typedef boost::geometry::model::polygon<point_type> polygon_type;

    polygon_type poly;
    boost::geometry::read_wkt(
       "POLYGON((2 1.3,2.4 1.7,2.8 1.8,3.4 1.2,3.7 1.6,3.4 2,4.1 3,5.3 2.6,5.4 1.2,4.9 0.8,2.9 0.7,2 1.3)"
           "(4.0 2.0, 4.2 1.4, 4.8 1.9, 4.4 2.2, 4.0 2.0))", poly);

    point_type p(4, 1);

    std::cout << "within: " << (boost::geometry::within(p, poly) ? "yes" : "no") << std::endl;
}

void setInitialSolution(const GlobalParameters *gp, FullSolution& curFsol){
    std::ifstream inf(gp->_initial_sol_file_name, std::ios::in);
    if (gp->_initial_sol_file_name.empty() || !inf.is_open()){
        warning_msg("Cannot open initial solution file \"" +
                    gp->_initial_sol_file_name + "\", use original greedy solution.");
        return;
    }
    std::string y, x;
    inf >> y >> x;
    std::reverse(y.begin(),y.end());
    std::reverse(x.begin(),x.end());
    curFsol = {PartSolution(y), PartSolution(x)};
    return;
}

void randomTry(Layer* layer, Layer* aLayer, const int& n) {
    FuncLog f(nullptr, "randomTry", "\t\t", true);
    f.startLog();

    int t = 0;
    int tn = n/10;
    int rejectTimes = 0, totalTimes = 0;
    while (t <= n) {
        totalTimes++;
        layer->random(layer->cur_sol());
        aLayer->random(aLayer->cur_sol());

        FullSolution fsol = std::make_pair(layer->cur_sol(), aLayer->cur_sol());

        double r = 0.0, ar = 0.0;
        bool isBE = layer->isBudgetEnough(layer->cur_sol(), r),
             isABE = aLayer->isBudgetEnough(aLayer->cur_sol(), ar);
        if (layer->gp()->debugLevel() >= 2)
            std::cout << "Try " << t << ": \ncur_sol = " << layer->cur_sol() << ", enough = "
                      << std::boolalpha << isBE << ", remain = " << r
                      << "\na.cur_sol = " << aLayer->cur_sol() << ", enough = "
                      << std::boolalpha << isABE << ", remain = " << ar <<  std::endl;
        if (isBE && isABE) {
            double r = layer->revenue(aLayer, fsol, 1), ar = aLayer->revenue(layer, fsol, 1);
            if (layer->gp()->debugLevel() >= 2) {
                std::cout << "Revenue = " << r << std::endl;
                std::cout << "a.Revenue = " << ar << std::endl;
            }
            if (t % tn == 0) {
                std::cout << "Randomed " << t << "/" << n << " solutions" << std::endl;
                std::cout << "\tRejection rate = " << rejectTimes << "/" << totalTimes
                          << " (" << (double) rejectTimes/totalTimes*100 << "%)"  << std::endl;
            }
            t++;
        } else rejectTimes++;
    }
}

void parse(int argc, char *argv[], Parser& p, GlobalParameters& gp) {
    FuncLog f(&gp, "parse", "\t", true);
    f.startLog();

    p.parseCommand(argc, argv);
    if (gp.isSaveSeeds()) gp.saveSeeds();
    if (gp.isLoadSeeds()) gp.loadSeeds();
    gp.saveCSeeds();

    if (gp.isInputFile()) p.readInputParameterFile();
    else p.initLayers();
    error_if_not(gp.upperLayer() && gp.lowerLayer(), "Layers are invalid");
    if (!gp.isInputFile()) {
      p.randomGenerateOrAdjustGlobalParameters();
      gp._upperLayer->randomGenerateOrAdjust();
    }
    gp._upperLayer->allocateSolutions();
    if (!gp.isInputFile()) {
      gp._lowerLayer->randomGenerateOrAdjust();
    }
    gp._lowerLayer->allocateSolutions();

    if (gp.isOutputFile()) p.writeOutputParameterFile();

    // if (gp.isVerbose())
    gp.show();
}

FullSolution completeAlgorithm(GlobalParameters& gp){
    FuncLog f(&gp, "CompleteAlgorithm", "\t");
    f.startLog();

    FuncLog fi(&gp, "initial part", "\t\t");
    fi.startLog();
    clock_t i_start = clock();

    int g = 0, l = 1, k = 1;
    Layer *ul = gp._upperLayer, *ll = gp._lowerLayer;
    ul->initSolutions();
    ll->initSolutions();

    if (gp.isVerbose())
        gp.show();

    FullSolution initFsol = std::make_pair(ul->cur_sol(), ll->cur_sol());
    if (gp.isVerbose()){
        ul->showSolutions();
        ll->showSolutions();
    }

    std::unordered_set<PartSolution> setSols;
    std::unordered_set<FullSolution> setFullSols;

    std::cout << "Initial upper layer solution = " << ul->cur_sol() << std::endl;
    std::cout << "Initial lower layer solution = " << ll->cur_sol() << std::endl;

    if (!gp.isIgnoreInitial()) ul->findNeighbors(ul->cur_sol(), setSols, k, true /* Initial */);
    int llg = ll->tabu()->gmax();

    // ll->setTabuGMax(1);

    std::cout << "Initial VNS Y's #neighbor = " << setSols.size() << std::endl;

    if (!gp.isIgnoreInitial())
        for (auto& sol : setSols) {
            double remain = 0;
            if (!ul->isBudgetEnough(sol,remain)) {
                std::string str;
                to_string(sol, str);
                warning_msg("ignore a sol " + str + " is not enough budget, remain = "
                            + std::to_string(remain));
                continue;
            }
            FullSolution tabuFsol;
            FullSolution fsol = std::make_pair(sol, ll->cur_sol());
            if (gp.isCache()) {
                auto it = ul->tabuBestMap().find(sol);
                if (it != ul->tabuBestMap().end())
                    tabuFsol = std::make_pair(sol, it->second);
                else
                    tabuFsol = ll->tabuSearch(ul, fsol);
            } else
                tabuFsol = ll->tabuSearch(ul, fsol);

            error_if_not(fsol.first == tabuFsol.first, "y should be the same!");
            error_if(ul->tabuBestMap().count(tabuFsol.first), "a y is already in map");
            ul->tabuBestMap()[tabuFsol.first] = tabuFsol.second;
            setFullSols.insert(tabuFsol);
        }
        if (gp.isVerbose()) {
            std::cout << "M = {\n";
            int j = 0;
            for (auto& sol : setFullSols) {
                if (j >= 10 && gp.debugLevel() < 2) {
                    std::cout << "...\n";
                    break;
                }
                std::cout << "\t" << sol << std::endl;
            }
            std::cout << "}\n";
        }

    FullSolution curFsol = initFsol;

    if (!gp.isIgnoreInitial())
        curFsol = *std::max_element(setFullSols.begin(), setFullSols.end(),
          [&](const FullSolution& fa, const FullSolution& fb) -> bool {
          return ul->revenue(ll, fa) < ul->revenue(ll, fb);
        });



    std::cout << "initSol                = " << initFsol
              << ", initULRev = " << ul->revenue(ll, initFsol)
              << ", initLLRev = " << ll->revenue(ul, initFsol)
              << std::endl;

    PartSolution initY = initFsol.first;
    PartSolution bestX = initFsol.second;
    bool isXFound = false;
    for (auto & sol : setFullSols)
        if (sol.first == initY) {
             bestX = sol.second;
             isXFound = true;
             break;
        }
    if (!gp.isIgnoreInitial())
        warning_if_not(isXFound, "best x not found for initial y!");
    FullSolution iYbX = std::make_pair(initY, bestX);

    if (!gp.isIgnoreInitial())
    std::cout << "initSol Y with Best X  = " << iYbX
              << ", initULRev = " << ul->revenue(ll, iYbX)
              << ", initLLRev = " << ll->revenue(ul, iYbX)
              << std::endl;

    PartSolution initX = initFsol.second;
    PartSolution bestY = initFsol.first;
    PartSolution worstY = initFsol.first;
    // bool isYFound = false;
    double minRev = DBL_MAX, maxRev = DBL_MIN;
    for (auto & sol : setFullSols)
        if (sol.second == initX) {
             const double &yRev = ul->revenue(ll, sol);
             if (yRev > maxRev) {
                 maxRev = yRev;
                 bestY = sol.first;
             }
             if (yRev < minRev) {
                 minRev = yRev;
                 worstY = sol.first;
             }
        }
    // warning_if_not(isYFound, "best y not found for initial x!");
    FullSolution iXbY = std::make_pair(bestY, initX);
    FullSolution iXwY = std::make_pair(worstY, initX);

     if (!gp.isIgnoreInitial())
    std::cout << "initSol X with Best Y  = " << iXbY
              << ", initULRev = " << ul->revenue(ll, iXbY)
              << ", initLLRev = " << ll->revenue(ul, iXbY)
              << std::endl;
      if (!gp.isIgnoreInitial())
    std::cout << "initSol X with Worst Y = " << iXwY
              << ", initULRev = " << ul->revenue(ll, iXwY)
              << ", initLLRev = " << ll->revenue(ul, iXwY)
              << std::endl;


    if (gp.isIgnoreInitial()) setInitialSolution(ul->gp(), curFsol);
    FullSolution initialBestSol = curFsol;

    std::cout << "initialVNSBestSol      = " << curFsol
              << ", bestULRev = " << ul->revenue(ll, curFsol)
              << ", bestLLRev = " << ll->revenue(ul, curFsol)
              << std::endl << std::endl;

    fi.endLog();
    double i_time = (double) (clock() - i_start) / CLOCKS_PER_SEC;
    clock_t u_start = clock();

    bool isFirst = true;
    FullSolution iterBestFsol, lastIterBestFsol;
    std::unordered_set<FullSolution> setBestSols;
    FullSolution trashSol = curFsol;
    ll->setTabuGMax(llg);
    while (g < ul->gmax()) {
        iterBestFsol = curFsol;
        k = 1;
        while (k <= ul->kmax()) {
            if (!isFirst) {
                trashSol = ul->changeLocation(curFsol, k);
                std::cout << "random jump solution from cur = " << curFsol << " to new = " << trashSol << std::endl;
            }
            isFirst = false;
            FullSolution maxFsol = ul->tabuSearch(ll, trashSol);
            // std::cout << "maxFsol = " << maxFsol << ", curFsol = " << curFsol << std::endl;
            if (ul->revenue(ll, maxFsol) > ul->revenue(ll, curFsol)) {
                curFsol = maxFsol;
                k = 1;
            } else k++;

            if (ul->revenue(ll, curFsol) > ul->revenue(ll, iterBestFsol)) {
                iterBestFsol = curFsol;
            }
            std::cout << "main k = " << k << ", l = " << l << std::endl;
        }
        if (l > 1) {
            if (ul->revenue(ll, iterBestFsol) < ul->revenue(ll, lastIterBestFsol)*1.01)
                g++;
            else
                g = 0;
        }
        l++;

        if (gp.isVerbose())
            std::cout << "main " << "g = " << g << ", l = " << l << std::endl;
        std::cout << "Complete Algorithm iteration l = " << l << ", iterBestFsol = " << iterBestFsol << std::endl;
        std::cout << "ulRev = " << ul->revenue(ll, iterBestFsol) << ", llRev = "
                  << ll->revenue(ul, iterBestFsol) << std::endl;
        std::cout << "main upper loop once k = " << k << ", g = " << g << ", l = " << l << std::endl;
        gp.showElapsedTime();
        setBestSols.insert(iterBestFsol);
        lastIterBestFsol = iterBestFsol;
    }

    error_if_not(setBestSols.size() > 0, "the set of solutions is empty!");
    FullSolution max_fsol = *std::max_element(setBestSols.begin(), setBestSols.end(),
        [&](const FullSolution& fa, const FullSolution& fb) -> bool {
        return ul->revenue(ll, fa) < ul->revenue(ll, fb);
    });

    ul->setBestSol(max_fsol.first);
    ll->setBestSol(max_fsol.second);
    std::cout << "ul best = " << ul->best_sol() << std::endl;
    std::cout << "ll best = " << ll->best_sol() << std::endl;

    gp.show();

    std::cout << "completeAlgorithm used #iterations = " << l << std::endl;
    std::cout << std::endl;
    std::cout << "Upper Tabu used #iterations        = " << ul->tabuIter() << std::endl;
    std::cout << "Upper Tabu used runtime            = " << ul->tabuRuntime() << std::endl;
    std::cout << "Upper Tabu avg runtime per iter    = " << ul->tabuRuntime()/ul->tabuIter() << std::endl;
    std::cout << "Upper Tabu max runtime per iter    = " << ul->maxTabuRuntime() << std::endl;
    std::cout << "Upper Tabu min runtime per iter    = " << ul->minTabuRuntime() << std::endl;
    std::cout << std::endl;
    std::cout << "Lower Tabu used #iterations        = " << ll->tabuIter() << std::endl;
    std::cout << "Lower Tabu used runtime            = " << ll->tabuRuntime() << std::endl;
    std::cout << "Lower Tabu avg runtime per iter    = " << ll->tabuRuntime()/ll->tabuIter() << std::endl;
    std::cout << "Lower Tabu max runtime per iter    = " << ll->maxTabuRuntime()  << std::endl;
    std::cout << "Lower Tabu min runtime per iter    = " << ll->minTabuRuntime()  << std::endl;


    std::cout << "Initial solution       (init x init y)  = " << initFsol <<
                 ", ulRev = " << ul->revenue(ll, initFsol, 1) <<
                 ", llRev = " << ll->revenue(ul, initFsol, 1) << std::endl;
     if (!gp.isIgnoreInitial())
    std::cout << "initSol Y with Best X  (best x init y)  = " << iYbX
              << ", ulRev = " << ul->revenue(ll, iYbX)
              << ", llRev = " << ll->revenue(ul, iYbX)
              << std::endl;
     if (!gp.isIgnoreInitial())
     std::cout << "initSol X with Best Y  (init x best y)  = " << iXbY
              << ", ulRev = " << ul->revenue(ll, iXbY)
              << ", llRev = " << ll->revenue(ul, iXbY)
              << std::endl;
      if (!gp.isIgnoreInitial())
     std::cout << "initSol X with Worst Y (init x worst y) = " << iXwY
              << ", ulRev = " << ul->revenue(ll, iXwY)
              << ", llRev = " << ll->revenue(ul, iXwY)
              << std::endl;
if (gp.isIgnoreInitial()) std::cout << "Input initial sol = \n";
    std::cout << "Initial VNS best sol   (best x best y)  = " << initialBestSol <<
                 ", ulRev = " << ul->revenue(ll, initialBestSol, 1) <<
                 ", llRev = " << ll->revenue(ul, initialBestSol, 1) << std::endl;
    std::cout << "Global best sol        (tabu x tabu y)  = " << max_fsol <<
                 ", ulRev = " << ul->revenue(ll, max_fsol, 1) <<
                 ", llRev = " << ll->revenue(ul, max_fsol, 1) << std::endl;

    double u_time = (double) (clock() - u_start) / CLOCKS_PER_SEC;
    double t_time = (double) (clock() - i_start) / CLOCKS_PER_SEC;

    std::cout << "Initial VNS time               = " << i_time << " sec." << std::endl;
    std::cout << "Total tabu time                = " << u_time << " sec." << std::endl;
    std::cout << "Total complete algo. time      = " << t_time << " sec." << std::endl;

    return max_fsol;
}

void testRandom(GlobalParameters& gp){
    FuncLog f(&gp, "testRandom", "\t");
    f.startLog();
    if (gp.randomTimes() > 0) {
        gp._upperLayer->initSolutions();
        gp._lowerLayer->initSolutions();
        if (gp.isVerbose()){
            gp._upperLayer->showSolutions();
            gp._lowerLayer->showSolutions();
        }
        gp.setCache(true);
        randomTry(gp._upperLayer, gp._lowerLayer, gp.randomTimes());
        gp.setCache(false);
        randomTry(gp._upperLayer, gp._lowerLayer, gp.randomTimes());
    }
}

int main(int argc, char *argv[])
{
    FuncLog f(nullptr, "main", "", true);
    f.startLog();
    GlobalParameters gp;
    gp.startTiming();
    Parser p(&gp);

    parse(argc, argv, p, gp);

#ifdef USE_QT
    if (p.isQT() && gp.drawOnly()) {
        if (gp.isIgnoreInitial()) {
            FullSolution bestSol;
            setInitialSolution(&gp, bestSol);
            gp._upperLayer->setBestSol(bestSol.first);
            gp._lowerLayer->setBestSol(bestSol.second);
        }
        QApplication a(argc, argv);
        MainWindow w(gp._upperLayer, gp._lowerLayer, &gp, gp._vecFamily);
        w.show();
        return a.exec();
    }
#endif

    // test runtime
    testRandom(gp);

    // complete algorithm
    FullSolution allBest = completeAlgorithm(gp);


    // testPointInsidePolygon();

#ifdef USE_QT
    if (p.isQT()) {
        QApplication a(argc, argv);
        MainWindow w(gp._upperLayer, gp._lowerLayer, &gp, gp._vecFamily);
        w.show();
        return a.exec();
    }
#endif

}

std::ostream& operator << (std::ostream& out, const GlobalParameters& gp) {
    printBanner("Program Parameter Info", out);
    out << "_is_verbose = " << std::boolalpha << gp._is_verbose << std::endl
        << "_is_input_file = " << gp._is_input_file << std::endl
        << "_is_output_file = " << gp._is_output_file << std::endl
        << "_is_exact = " << gp._is_exact << std::endl
        << "_is_cache = " << gp.isCache() << std::endl
        << "_is_circle = " << gp.isCircle() << std::endl
        << "_is_timed = " << gp.isTimed() << std::endl
        << "_is_save_seed = " << gp.isSaveSeeds() << std::endl
        << "_is_load_seed = " << gp.isLoadSeeds() << std::endl
        << "_is_qt = " << gp._is_qt << std::endl
        << "_debug_level = " << std::noboolalpha << gp._debug_level << std::endl
        << "_bit_diff = " << gp.bitDiff() << std::endl
        << "_random_times = " << gp.randomTimes() << std::endl
        << "_input_file_name = " << gp._input_file_name << std::endl
        << "_output_file_name = " << gp._output_file_name << std::endl;

    printBanner("Statistic Parameter Info", out);
    int i = 0, j = 0, k = 0;
    out << "Commodity Price: " << std::endl;
    for (auto& c : gp._cp)
        out << "\t" << _mtName.at(static_cast<MARKET_TYPE>(i++)) << " = " << c << std::endl;
    i = 0;
    out << "Parking Attraction: " << std::endl;
    for (auto& p : gp._pa)
        out << "\t" << _mtName.at(static_cast<MARKET_TYPE>(i++)) << " = " << p << std::endl;
    i = 0;
    out << "Buying Power: " << std::endl;
    for (auto& b : gp._bp)
        out << "\t" << _ftName.at(static_cast<FAMILY_TYPE>(i++)) << " = " << b << std::endl;

    i = 0;
    out << "Attraction: " << std::endl;
    for (auto& a1 : gp._at) {
        j = 0;
        for (auto& a2 : a1) {
            k = 0;
            for (auto& a3 : a2) {
                out << "\t" << _mtName.at(static_cast<MARKET_TYPE>(i)) <<
                       " " << _ftName.at(static_cast<FAMILY_TYPE>(j)) <<
                       " " << _dtName.at(static_cast<DAY_TYPE>(k++)) << " = " << a3 << std::endl;
            }
            ++j;
        }
        ++i;
    }

    for (int i = 1; i < 5; ++i)
        out << "Theta[" << i << "] = " << gp._theta[i] << std::endl;

    printBanner("Family Info", out);
    out << "#Family: " << std::endl;
    for (auto& el : gp.mapVecFamily())
       out << "\t" << _ftName.at(el.first) << " = " << el.second.size() << std::endl;
    i = 0;
    out << "_vecFamily(" << gp._vecFamily.size() << ") = \n{" << std::endl;
    for (auto& fam : gp._vecFamily) {
        if (i > 3 && gp.debugLevel() < 2) {
            out << "\t...\n";
            break;
        }
        out << "\tFamily[" << i++ << "] = " << *fam << std::endl;
    }
    out << "}" << std::endl;

    printBanner("Upper Layer Info", out);
    out << "_upperLayer = \n{\n" << *gp._upperLayer << "}" << std::endl;
    printBanner("Lower Layer Info", out);
    out << "_lowerLayer = \n{\n" << *gp._lowerLayer << "}" << std::endl;

    return out;
}
