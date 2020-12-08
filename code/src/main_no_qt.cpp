#include <bits/stdc++.h>
#include <getopt.h>
#include <boost/dynamic_bitset.hpp>
#include "gurobi_c++.h"
#include "util.h"
#include "layer.h"
#include "parser.h"

using namespace std;

int main(int argc, char** argv) 
{
    Layer *upperLayer = nullptr, *lowerLayer = nullptr;
    Parser p(argc, argv, upperLayer, lowerLayer);
}
