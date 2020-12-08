#include "parser.h"

void Parser::readInputParameterFile(const std::string& str)
{
    if (!str.empty()) _gp->_input_file_name = str;
    std::ifstream fin(_gp->_input_file_name, std::ios::in);
    error_if_not(fin.is_open(), "Cannot open the input parameter file \"" + _gp->_input_file_name + "\", terminated");

    std::stringstream ss;
    bool isHyper = false, isUpper = false;
    double num = 0.0;
    int n = 0;
    double x = 0.0, y = 0.0;
    std::string specifier, title, mt, ft, dt, name, type;
    Layer** cLayer = nullptr;
    std::vector<FAMILY_TYPE> vecFt;
    std::vector<double> vecCost[2], vecNSS[2], vecFamBP;
    std::vector<Location> vecFamLoc, vecFacLoc[2], vecEFacLoc[2];
    while (myGetline(fin, ss)) {
        if (_gp->debugLevel() > 0)
            std::cout << "1st ss.str(): " <<  ss.str() << std::endl;
        if (myGet(ss, specifier)) {
            if (_gp->debugLevel() > 0)
                std::cout << "1st specifier: " <<  specifier << std::endl;
            if (boost::iequals(specifier, "begin") && myGet(ss, title)) {
                if (boost::iequals(title, "global")) {
                    while (myGetline(fin, ss)) {
                        if (_gp->debugLevel() > 0)
                            std::cout << "2nd ss.str(): " <<  ss.str() << std::endl;
                        if (myGet(ss, specifier)) {
                            if (_gp->debugLevel() > 0)
                                std::cout << "2nd specifier: " <<  specifier << std::endl;
                            if (boost::iequals(specifier, "end")) break;
                            else if (boost::iequals(specifier, "commodityPrices")) {
                                int i = 0;
                                while (i++ < MARKET_SIZE && myGetline(fin, ss) && myGet(ss, mt) && myGet(ss, num))
                                    _gp->setCP(_nameMt.at(mt), num);
                            }
                            else if (boost::iequals(specifier, "parkingAttractions")) {
                                int i = 0;
                                while (i++ < MARKET_SIZE && myGetline(fin, ss) && myGet(ss, mt) && myGet(ss, num))
                                    _gp->setPA(_nameMt.at(mt), num);
                            }
                            else if (boost::iequals(specifier, "buyingPowers")) {
                                int i = 0;
                                while (i++ < FAMILY_SIZE && myGetline(fin, ss) && myGet(ss, ft) && myGet(ss, num)) {
                                    _gp->_maxbp = std::max(_gp->_maxbp, num);
                                    _gp->setBP(_nameFt.at(ft), num);
                                }
                            }
                            else if (boost::iequals(specifier, "attractions")) {
                                int i = 0;
                                while (i++ < FAMILY_SIZE*MARKET_SIZE*DAY_SIZE &&
                                        myGetline(fin, ss) && myGet(ss, mt) && myGet(ss, ft) &&
                                        myGet(ss, dt) && myGet(ss, num))
                                    _gp->setAT(_nameMt.at(mt), _nameFt.at(ft), _nameDt.at(dt), num);
                            }
                            else if (boost::iequals(specifier, "theta")) {
                                int i = 0;
                                while (++i < 5 && myGet(ss, num))
                                    _gp->setTheta(i, num);
                            }
                            else if (boost::iequals(specifier, "#family") && myGet(ss, n))
                                _gp->setNumFamily(n);
                            else if (boost::iequals(specifier, "vecfamilyTypes")) {
                                int i = 0;
                                while (i++ < _gp->numFamily() && myGet(ss, title)) {
                                    if (_gp->debugLevel() > 0) std::cout << title << std::endl;
                                    if (boost::iequals(title.substr(0, 1), "p"))
                                        vecFt.push_back(POOR);
                                    else if (boost::iequals(title.substr(0, 1), "m"))
                                        vecFt.push_back(MEDIUM);
                                    else if (boost::iequals(title.substr(0, 1), "r"))
                                        vecFt.push_back(RICH);
                                }
                            }
                            else if (boost::iequals(specifier, "vecfamilyBuyingPowers")) {
                                int i = 0;
                                while (i++ < _gp->numFamily() && myGet(ss, num)) {
                                    if (_gp->debugLevel() > 0)
                                        std::cout << num << std::endl;
                                    _gp->_maxbp = std::max(_gp->_maxbp, num);
                                    vecFamBP.push_back(num);
                                }
                            }
                            else if (boost::iequals(specifier, "vecfamilyLocations")) {
                                int i = 0;
                                while (i++ < _gp->numFamily() && myGet(ss, x, y)) {
                                    if (_gp->debugLevel() > 0)  std::cout << x << " " << y << std::endl;
                                    vecFamLoc.push_back(Location{x, y});
                                }
                            }
                        }
                    }
                }
                else if (boost::iequals(title, "layer")) {
                    while (myGetline(fin, ss))
                        if (myGet(ss, specifier)) {
                            if (_gp->debugLevel() > 0)
                                std::cout << "2nd specifier: " <<  specifier << std::endl;
                            if (boost::iequals(specifier, "end")) break;
                            else if (boost::iequals(specifier, "layerType") && myGet(ss, type)) {
                                if (boost::iequals(type.substr(0, 1), "u")) {
                                    isUpper = true;
                                    cLayer = &(_gp->_upperLayer);
                                }
                                else if (boost::iequals(type.substr(0, 1), "l")) {
                                    isUpper = false;
                                    cLayer = &(_gp->_lowerLayer);
                                }
                            }
                            else if (boost::iequals(specifier, "name") && myGet(ss, name)) {
                                if (*cLayer) (*cLayer)->setName(name);
                            }
                            else if (boost::iequals(specifier, "marketType") && myGet(ss, mt)) {
                                if (boost::iequals(mt.substr(0, 1), "s")) isHyper = false;
                                else if (boost::iequals(mt.substr(0, 1), "h")) isHyper = true;
                                (*cLayer) = new Layer(_gp, name, isHyper, isUpper);
                            }
                            else if (boost::iequals(specifier, "gmax") && myGet(ss, n)) {
                                if (*cLayer) (*cLayer)->setGMax(n);
                            }
                            else if (boost::iequals(specifier, "kmax") && myGet(ss, n)) {
                                if (*cLayer) (*cLayer)->setKMax(n);
                            }
                            else if (boost::iequals(specifier, "nss")) {
                                int i = 0;
                                vecNSS[isUpper].push_back(0.0);
                                while (i++ < (*cLayer)->kmax() && myGet(ss, num))
                                    vecNSS[isUpper].push_back(num);
                                if (*cLayer) (*cLayer)->setVecNSS(vecNSS[isUpper]);
                            }
                            else if (boost::iequals(specifier, "lmax") && myGet(ss, n)) {
                                if (*cLayer) (*cLayer)->setLMax(n);
                            }
                            else if (boost::iequals(specifier, "itmax") && myGet(ss, num)) {
                                if (*cLayer) (*cLayer)->setItMax(num);
                            }
                            else if (boost::iequals(specifier, "budget") && myGet(ss, num)) {
                                if (*cLayer) (*cLayer)->setBudget(num);
                            }
                            else if (boost::iequals(specifier, "#candidatefacility") && myGet(ss, n)) {
                                if (*cLayer) (*cLayer)->setNumFacility(n);
                            }
                            else if (boost::iequals(specifier, "vecCandidateLocations")) {
                                int i = 0;
                                while (i++ < (*cLayer)->numFacility() && myGet(ss, x, y)) {
                                    if (_gp->debugLevel() > 0)  std::cout << x << " " << y << std::endl;
                                    vecFacLoc[isUpper].push_back(Location{x, y});
                                }
                            }
                            else if (boost::iequals(specifier, "vecCanidateCosts")) {
                                int i = 0;
                                while (i++ < (*cLayer)->numFacility() && myGet(ss, num))
                                    vecCost[isUpper].push_back(num);
                            }
                            else if (boost::iequals(specifier, "#existedfacility") && myGet(ss, n))
                                (*cLayer)->setNumExistedFacility(n);
                            else if (boost::iequals(specifier, "vecExistedLocations")) {
                                int i = 0;
                                while (i++ < (*cLayer)->numExistedFacility() && myGet(ss, x, y)) {
                                    if (_gp->debugLevel() > 0)  std::cout << x << " " << y << std::endl;
                                    vecEFacLoc[isUpper].push_back(Location{x, y});
                                }
                            }
                            else if (boost::iequals(specifier, "begin") && myGet(ss, title)) {
                                if (boost::iequals(title, "tabu")) {
                                    while (myGetline(fin, ss)) {
                                        if (myGet(ss, specifier)) {
                                            if (boost::iequals(specifier, "end")) break;
                                            else if (boost::iequals(specifier, "size") && myGet(ss, n))
                                                (*cLayer)->setTabuSize(n);
                                            else if (boost::iequals(specifier, "gmax") && myGet(ss, n))
                                                (*cLayer)->setTabuGMax(n);
                                            else if (boost::iequals(specifier, "lmax") && myGet(ss, n))
                                                (*cLayer)->setTabuLMax(n);
                                            else if (boost::iequals(specifier, "cmax") && myGet(ss, n))
                                                (*cLayer)->setTabuCMax(n);
                                            else if (boost::iequals(specifier, "itmax") && myGet(ss, num))
                                                (*cLayer)->setTabuItMax(num);
                                        }
                                    }
                                }
                            }
                        }
                }
            }
        }
    }

    for (int i = 0; i < _gp->numFamily(); ++i) {
        // if (_gp->debugLevel() > 0)  std::cout << vecFamLoc[i].x() << " " << vecFamLoc[i].y() << std::endl;
        const double& bp = (vecFamBP.size() > 0) ? vecFamBP[i] : _gp->bp()[vecFt[i]];
        Family* fam = new Family(vecFamLoc[i].x(), vecFamLoc[i].y(), vecFt[i], bp);
        _gp->_vecFamily.push_back(fam);
        _gp->setMapVecFamily(vecFt[i], fam);
    }

    for (int j = 0; j < 4; ++j) {
        Layer* cLayer = (j >= 2) ? _gp->_upperLayer : _gp->_lowerLayer;
        if (j % 2) cLayer->setVecFamily(_gp->vecFamily());
        int n = (j % 2) ? cLayer->numFacility() : cLayer->numExistedFacility();
        std::vector<Location>& vecLoc = (j % 2) ? vecFacLoc[j >= 2] : vecEFacLoc[j >= 2];
        const std::vector<double>& vecRCost = (j % 2) ? vecCost[j >= 2] : std::vector<double>(n, 0.0);
        std::vector<Facility*> vecFac;
        for (int i = 0; i < n; ++i) {
            Facility* fac = new Facility(
                vecLoc[i].x(), vecLoc[i].y(), vecRCost[i],
                _gp->cp()[cLayer->mt()], _gp->pa()[cLayer->mt()], i);
            vecFac.push_back(fac);
        }
        (j % 2) ? cLayer->setVecFacility(vecFac) : cLayer->setVecExistedFacility(vecFac);
    }

}

void Parser::writeOutputParameterFile(const std::string& str)
{
    if (!str.empty()) _gp->_output_file_name = str;
    std::ofstream fout(_gp->_output_file_name, std::ios::out);
    warning_if_not(fout.is_open(), "Cannot open the output parameter file \"" + _gp->_output_file_name + "\", skipped");
    if (!fout.is_open()) return;
    // global
    fout << "BEGIN GLOBAL\n";
    fout << "\tcommodityPrices\n";
    int i = 0;
    for (auto& n : _gp->cp())
        fout << "\t\t" << _mtName.at(static_cast<MARKET_TYPE>(i++)) << "\t" << n << std::endl;
    fout << "\tparkingAttractions\n";
    i = 0;
    for (auto& n : _gp->pa())
        fout << "\t\t" << _mtName.at(static_cast<MARKET_TYPE>(i++)) << "\t" << n << std::endl;
    fout << "\tbuyingPowers\n";
    i = 0;
    for (auto& n : _gp->bp())
        fout << "\t\t" << _ftName.at(static_cast<FAMILY_TYPE>(i++)) << "\t" << n << std::endl;
    fout << "\tattractions\n";
    i = 0;
    for (auto& a1 : _gp->at()) {
        int j = 0;
        for (auto& a2 : a1) {
            int k = 0;
            for (auto& a3 : a2) {
                fout << "\t\t" << _mtName.at(static_cast<MARKET_TYPE>(i)) <<
                     "\t" << _ftName.at(static_cast<FAMILY_TYPE>(j)) <<
                     "\t" << _dtName.at(static_cast<DAY_TYPE>(k++)) << "\t" << a3 << std::endl;
            }
            ++j;
        }
        ++i;
    }
    fout << "\ttheta\t";
    for (int i = 1; i < 5; ++i)
        fout << _gp->theta()[i] << " ";
    fout << std::endl;
    fout << "\t#family\t" << _gp->numFamily() << std::endl;
    fout << "\tvecfamilyTypes\t";
    for (auto& f : _gp->vecFamily())
        fout << _ftName.at(f->type()) << " ";
    fout << std::endl;
    fout << "\tvecfamilyBuyingPowers\t";
    for (auto& f : _gp->vecFamily())
        fout << f->bp() << " ";
    fout << std::endl;
    fout << "\tvecfamilyLocations\t";
    for (auto& f : _gp->vecFamily())
        fout << "(" << f->x() << ", " << f->y() << ") ";
    fout << std::endl;
    fout << "END\n";
    // layers
    for (int i = 0; i < 2; ++i) {
        fout << std::endl;
        const Layer* cLayer = (i) ? _gp->lowerLayer() : _gp->upperLayer();
        fout << "BEGIN LAYER\n";
        fout << "\tlayerType\t" << ((i) ? "lower" : "upper") << std::endl;
        fout << "\tname\t" << cLayer->name() << std::endl;
        fout << "\tmarketType\t" << _mtName.at(cLayer->mt()) << std::endl;
        fout << "\tgMax\t" << cLayer->gmax() << std::endl;
        fout << "\tkMax\t" << cLayer->kmax() << std::endl;
        fout << "\tnss\t";
        int k = 0;
        for (auto& nss : cLayer->vecNSS())
            if (k++ > 0) {
                fout << nss << " ";
            }
        fout << std::endl;
        fout << "\tlMax\t" << cLayer->lmax() << std::endl;
        fout << "\titMax\t" << cLayer->itmax() << std::endl;
        fout << "\tbudget\t" << cLayer->budget() << std::endl;
        fout << "\t#candidatefacility\t" << cLayer->numFacility() << std::endl;
        fout << "\tvecCandidateLocations\t";
        for (auto& f : cLayer->vecFacility())
            fout << "(" << f->x() << ", " << f->y() << ") ";
        fout << std::endl;
        fout << "\tvecCanidateCosts\t";
        for (auto& f : cLayer->vecFacility())
            fout << f->cost() << " ";
        fout << std::endl;
        fout << "\t#existedfacility\t" << cLayer->numExistedFacility() << std::endl;
        fout << "\tvecExistedLocations\t";
        for (auto& f : cLayer->vecExistedFacility())
            fout << "(" << f->x() << ", " << f->y() << ") ";
        fout << std::endl;
        fout << "\tBEGIN TABU\n";
        fout << "\t\tsize\t" << cLayer->tabu()->size() << std::endl;
        fout << "\t\tgMax\t" << cLayer->tabu()->gmax() << std::endl;
        fout << "\t\tcMax\t" << cLayer->tabu()->cmax() << std::endl;
        fout << "\t\tlMax\t" << cLayer->tabu()->lmax() << std::endl;
        fout << "\t\titMax\t" << cLayer->tabu()->itmax() << std::endl;
        fout << "\tEND\n";
        fout << "END\n";
    }
}

// Print cmd usage
void Parser::printUsage(const std::string& str)
{
    std::cerr << bold_on << "\nNAME\n\t" << str.substr(2) << bold_off
              << 	" - a bi-level competitive facility locater" << std::endl;

    std::cerr << bold_on << "SYNOPSIS\n\t" << bold_off << str
              // << " [-b <str>] [-d <num>] [-g <str>] [-i <name>] [-k <str>] [-o <name>] [-s <str>] [-ehqv]" << std::endl;
                 << " [-I <name>] [-b <num>][-d <num>] [-i <name>] [-o <name>] [-r <num>] [-CDchlnqstv]" << std::endl;

    std::cerr << bold_on << "DESCRIPTION\n" << bold_off
              << bold_on <<	"\t-C,         --Circle\n\t\t" << bold_off
              << 			"init by circle (default is ring, circle = off)"<< std::endl
              << bold_on <<	"\t-D,   --Draw\n\t\t" << bold_off
              << 			"draw only (default = off)"<< std::endl
              << bold_on <<	"\t-I <name>,  --Initial_sol=<file name>\n\t\t" << bold_off
              << 			"set initial solution file name, and ignore initial VNS (default <name> = \"\")"<< std::endl
              << bold_on <<	"\t-b <num>,   --bit_diff=<number>\n\t\t" << bold_off
              << 			"set bit difference (default <num> = 4)"<< std::endl
              << bold_on <<	"\t-c,         --cache\n\t\t" << bold_off
              << 			"cache utility (default off)"<< std::endl
              << bold_on <<	"\t-d <num>,   --debug_level=<number>\n\t\t" << bold_off
              << 			"set debug level (default <num> = 0)"<< std::endl
              // << bold_on <<	"\t-e,         --exact\n\t\t" << bold_off
              // << 			"exact inputs (default off)"<< std::endl
              // << bold_on <<	"\t-g <str>,   --gmax=<string>\n\t\t" << bold_off
              // << 			"set gmax via string (default <str> = \"\")"<< std::endl
              << bold_on <<	"\t-h,         --help\n\t\t" << bold_off
              << 			"print the command usage (default = off)"<< std::endl
              << bold_on <<	"\t-i <name>,  --input_file=<file name>\n\t\t" << bold_off
              << 			"set input parameter file name (default <name> = \"parameter.txt\")"<< std::endl
              << bold_on <<	"\t-l <str>,   --load_seeds\n\t\t" << bold_off
              << 			"load saved random seeds (default = off)"<< std::endl
              << bold_on <<	"\t-o <name>,  --output_file=<file name>\n\t\t" << bold_off
              << 			"set output parameter file name (default <name> = \"parameter_out.txt\")"<< std::endl
              << bold_on <<	"\t-n,         --neighbor\n\t\t" << bold_off
              << 			"neighbor cache (default = off)"<< std::endl
              << bold_on <<	"\t-q,         --qt\n\t\t" << bold_off
              << 			"draw results using qt (default = off)"<< std::endl
              << bold_on <<	"\t-s,         --save_seeds\n\t\t" << bold_off
              << 			"save random seeds (default = off)"<< std::endl
              << bold_on <<	"\t-r <num>,   --random_times=<number>\n\t\t" << bold_off
              << 			"set random times (default <num> = 0)"<< std::endl
              << bold_on <<	"\t-t,         --time\n\t\t" << bold_off
              << 			"time each function (default = off)"<< std::endl
              // << bold_on <<	"\t-s <str>,   --size=<string>\n\t\t" << bold_off
              // << 			"set size via string (default <str> = \"\")"<< std::endl
              << bold_on <<	"\t-v,         --verbose\n\t\t" << bold_off
              << 			"print the verbose usage (default = off)"<< std::endl;

    std::cerr << bold_on << "EXAMPLES" << std::endl << bold_off
              << "\t" << str << "" << std::endl
              << "\t" << str << " --i=parameter.txt" << std::endl
              // << "\t" << str << " -b \"100\" -s \"10 13\" -k 3 -g 6"  << std::endl
              // << "\t" << str << " -e -b \"100\" -s \"10 13\" -k 3 -g 6"  << std::endl
              << "\t" << str << " -o" << std::endl;

    std::cerr << bold_on << "AUTHOR" << std::endl << bold_off
              << "\tCode\n\t\tYu-Hsiang Cheng (鄭宇翔) <slightencheng@gmail.com>" << std::endl
              << "\tAlgorithm\n\t\tKeng-Hua Chuang (莊耿樺) <az41304az@gmail.com>" << std::endl;

    std::cerr << bold_on << "COPYRIGHT" << std::endl << bold_off
              << "\tCopyright © 2019 NCTU" << std::endl << std::endl;
}

// Parse command
void Parser::parseCommand(int argc, char** argv)
{
    char opt = 0;
    int index;
    static const struct option longopts[] = {
        {"Circle",          no_argument,        0, 'C'},
        {"Draw",            no_argument,        0, 'D'},
        {"Initial_sol", 	optional_argument,  0, 'I'},
    //    {"budget",          required_argument,  0, 'b'},
        {"bit_diff",        required_argument,  0, 'b'},
        {"cache",           no_argument,        0, 'c'},
        {"debug_level", 	optional_argument,  0, 'd'},
    //    {"exact",          no_argument,        0, 'e'},
    //    {"gmax",            required_argument,  0, 'g'},
        {"help",       		no_argument,        0, 'h'},
        {"input_file", 		optional_argument,  0, 'i'},
        {"load_seeds",      no_argument,        0, 'l'},
    //    {"kmax",            required_argument,  0, 'k'},
        {"neighbor",        no_argument,        0, 'n'},
        {"output_file",     optional_argument,  0, 'o'},

        {"qt",              no_argument,        0, 'q'},
        {"random_times",    required_argument,  0, 'r'},
        {"save_seeds",      no_argument,        0, 's'},
        {"time",            no_argument,        0, 't'},
    //    {"size",            required_argument,  0, 's'},
        {"verbose",    		no_argument,        0, 'v'},
        {0, 0, 0, 0}
    };
    // double ub = 0, lb = 0;
    // int ug = 0, lg = 0;
    // int uk = 0, lk = 0;
    // while ((opt = getopt_long(argc, argv, "b:d::eg:hi::k:o::qs:v", longopts, &index)) != -1) {
    while ((opt = getopt_long(argc, argv, "CDI::b:cd::hi::lno::qr:stv", longopts, &index)) != -1) {
        switch (opt) {
            // case 'b': {
            //    std::stringstream ss(optarg);
            //    ss >> ub >> lb;
            //    _gp->_upperLayer->setBudget(ub > 0 || _gp->_is_exact ? ub : lb);
            //    _gp->_lowerLayer->setBudget(lb > 0 || _gp->_is_exact ? lb : ub);
            //    break;
            // }
            case 'C':
                _gp->_is_circle = true;
                break;
            case 'D':
                _gp->_drawOnly = true;
                break;
            case 'I':
                _gp->_is_ignore_initial = true;
                if (optarg != NULL)
                    _gp->_initial_sol_file_name = optarg;
                break;
            case 'b':
                _gp->_bit_diff = atoi(optarg);
                break;

            case 'c':
                _gp->_is_cache = true;
                break;
            case 'd':
                if (optarg != NULL)
                    _gp->_debug_level = atoi(optarg);
                else _gp->_debug_level = INT_MAX;
                _gp->_is_verbose = true;
                break;
            // case 'e':
            //    _gp->_is_exact = true;
            //    break;
            // case 'g': {
            //    std::stringstream ss(optarg);
            //    ss >> ug >> lg;
            //    _gp->_upperLayer->setGMax(ug > 0 || _gp->_is_exact ? ug : lg);
            //    _gp->_lowerLayer->setGMax(lg > 0 || _gp->_is_exact ? lg : ug);
            //    break;
            // }
            case 'i':
                _gp->_is_input_file = true;
                if (optarg != NULL)
                    _gp->_input_file_name = optarg;
                break;
            // case 'k': {
            //    std::stringstream ss(optarg);
            //    ss >> uk >> lk;
            //    _gp->_upperLayer->setKMax(uk > 0 || _gp->_is_exact ? uk : lk);
            //    _gp->_lowerLayer->setKMax(lk > 0 || _gp->_is_exact ? lk : uk);
            //    break;
            // }
            case 'l':
                _gp->_is_load_seeds = true;
                break;
            case 'n':
                _gp->_is_neighbor_cache = true;
                break;
            case 'o':
                _gp->_is_output_file = true;
                if (optarg != NULL)
                    _gp->_output_file_name = optarg;
                break;
            case 'q':
                _gp->_is_qt = true;
                break;
            case 'r': {
                int r = atoi(optarg);
                if (r < 0) r = 0;
                _gp->_random_times = r;
                break;
            }
            case 't':
                _gp->_is_timed = true;
                break;
            case 's':
                _gp->_is_save_seeds = true;
                break;
            case 'v':
                _gp->_is_verbose = true;
                break;
            case 'h':
            default:
                printUsage(argv[0]);
                exit(EXIT_FAILURE);
                break;
        }
    }
}

void Parser::initLayers()
{
    if (_gp->_upperLayer == nullptr)
        _gp->_upperLayer = new Layer(_gp, "upperLayer", true, true);

    if (_gp->_lowerLayer == nullptr)
        _gp->_lowerLayer = new Layer(_gp, "lowerLayer", false, false);
}

void Parser::randomGenerateOrAdjustGlobalParameters()
{
    if (_gp->_cp[SUPER] <= 0) {
        _gp->_cp[SUPER] = _gp->_rr.Ur(10, 100);
        if (_gp->isVerbose())
            info_msg("randomed cp[SUPER] = " + std::to_string(_gp->_cp[SUPER]));
    }
    if (_gp->_cp[HYPER] <= 0) {
        _gp->_cp[HYPER] = _gp->_rr.Ur(0.8, 1.0)*_gp->_cp[SUPER];
        if (_gp->isVerbose())
            info_msg("randomed cp[HYPER] = " + std::to_string(_gp->_cp[HYPER]));
    }
    if (_gp->_pa[SUPER] != 0) {
        _gp->_pa[SUPER] = 0;
        if (_gp->isVerbose())
            info_msg("pa[SUPER] = 0");
    }
    if (_gp->_pa[HYPER] <= 0) {
        _gp->_pa[HYPER] = _gp->_rr.Ur(1.0);
        if (_gp->isVerbose())
            info_msg("randomed pa[HYPER] = " + std::to_string(_gp->_pa[HYPER]));
    }
    if (_gp->_bp[POOR] <= 0) {
        _gp->_bp[POOR] = _gp->_rr.Ur(1.0);
        if (_gp->isVerbose())
            info_msg("randomed bp[POOR] = " + std::to_string(_gp->_bp[POOR]));
    }
    if (_gp->_bp[MEDIUM] <= 0) {
        _gp->_bp[MEDIUM] = _gp->_rr.Ur(1.1, 3.0);
        if (_gp->isVerbose())
            info_msg("randomed bp[MEDIUM] = " + std::to_string(_gp->_bp[MEDIUM]));
    }
    if (_gp->_bp[RICH] <= 0) {
        _gp->_bp[RICH] = _gp->_rr.Ur(3.1, 5.0);
        if (_gp->isVerbose())
            info_msg("randomed bp[RICH] = " + std::to_string(_gp->_bp[RICH]));
    }
    _gp->_maxbp = std::max(_gp->_maxbp, _gp->_bp[RICH]);
    if (_gp->_theta[1] <= 0) {
        _gp->_theta[1] = 1;
        if (_gp->isVerbose())
            info_msg("theta[1] = 1");
    }
    if (_gp->_theta[3] <= 0) {
        _gp->_theta[3] = _gp->_rr.Ur(0.01, 0.05);
        if (_gp->isVerbose())
            info_msg("randomed theta[3] = " + std::to_string(_gp->_theta[3]));
    }
    if (_gp->_theta[4] <= 0) {
        _gp->_theta[4] = _gp->_rr.Ur(0.01, 0.05);
        if (_gp->isVerbose())
            info_msg("randomed theta[4] = " + std::to_string(_gp->_theta[4]));
    }
    // if (theta[2] <= 0) theta[2] = _gp->_rr.Ur(1);
    if (_gp->_theta[2] <= 0) {
        _gp->_theta[2] = _gp->_rr.Ur(1.0)*(_gp->_theta[4]*100)/log(1+_gp->_pa[HYPER]);
        if (_gp->isVerbose())
            info_msg("randomed theta[2] = " + std::to_string(_gp->_theta[2]));
    }
    for (auto& a1 : _gp->_at)
        for (auto& a2 : a1)
            for (auto& a3 : a2)
                if (a3 <= 0) {
                    a3 = _gp->_rr.Ur(10);
                    if (_gp->isVerbose())
                        info_msg("randomed att = " + std::to_string(a3));
                }
    if (_gp->numFamily() <= 0) {
        _gp->_numFamily = _gp->_ri.Ui(50, 250);
        // _gp->_numFamily = _gp->_ri.Ui(500, 2500);
        if (_gp->isVerbose())
            info_msg("randomed #family = " + std::to_string(_gp->_numFamily));
    }
    if ((int) _gp->_vecFamily.size() != _gp->numFamily()) {
        _gp->_vecFamily.resize(_gp->numFamily());
        for (auto& fam : _gp->_vecFamily) {
            if (fam) continue;
            FAMILY_TYPE f = static_cast<FAMILY_TYPE>(_gp->_ri.Ui(POOR, RICH));
            fam = new Family(_gp->_rr.Ur(-100, 100), _gp->_rr.Ur(-100, 100), f, _gp->_bp[f]);
            _gp->setMapVecFamily(f, fam);
        }
    }
    _gp->_upperLayer->setVecFamily(_gp->_vecFamily);
    _gp->_lowerLayer->setVecFamily(_gp->_vecFamily);
}
