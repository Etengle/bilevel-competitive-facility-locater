#ifndef LOCATION_H
#define LOCATION_H
#include "util.h"

class Location {

    friend std::ostream& operator << (std::ostream& out, const Location& loc) {
        out << "(" << loc._x << "," << loc._y << ")";
        return out;
    }

public:
    Location(const double& x = 0, const double& y = 0) : _x(x), _y(y) {}
    ~Location() {}

    inline double dist(const Location& aloc) const {
        return std::sqrt((_x-aloc._x)*(_x-aloc._x) + (_y-aloc._y)*(_y-aloc._y));
    }
    inline double dist(const Location* aloc) const {
        error_if_not(aloc, "a location is nullptr");
        return std::sqrt((_x-aloc->_x)*(_x-aloc->_x) + (_y-aloc->_y)*(_y-aloc->_y));
    }
    void show(std::ostream& out = std::cout) const { out << *this << std::endl; }
    inline double x() { return _x; }
    inline double y() { return _y; }
private:
    double _x = 0;
    double _y = 0;
};

#endif // LOCATION_H
