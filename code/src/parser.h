#ifndef PARSER_H
#define PARSER_H

#include "util.h"
#include "layer.h"

class Parser
{

    friend std::ostream& operator << (std::ostream& out, const Parser& parser) {
        parser._gp->show();
        return out;
    }

public:
    Parser(GlobalParameters* gp) : _gp(gp) {}
    ~Parser() {}
    void parseCommand(int argc, char** argv);
    void printUsage(const std::string& str = "");
    void readInputParameterFile(const std::string& str = "");
    void writeOutputParameterFile(const std::string& str = "");
    void randomGenerateOrAdjustGlobalParameters();
    void show(std::ostream& out = std::cout) { out << *this << std::endl; }
    bool isQT() const { return _gp->_is_qt; }
    bool isVerbose() const { return _gp->_is_verbose; }
    int debugLevel() const { return _gp->_debug_level; }
    void initLayers();

private:
    GlobalParameters* _gp;
};

#endif // PARSER_H
