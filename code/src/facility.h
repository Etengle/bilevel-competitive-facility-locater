#ifndef FACILITY_H
#define FACILITY_H
#include "location.h"

class Facility : public Location {

    friend std::ostream& operator << (std::ostream& out, const Facility& fac) {
        const Location loc = static_cast<Location>(fac);
        out << "{id: " << fac._id << ", loc: " << loc
            << ", pa: " << fac._pa << ", cp: " << fac._cp
            << ", cost: " << fac._cost << ", estUtil: " << fac._estUtil << "}";
        return out;
    }

public:
    Facility(const double& x = 0, const double& y = 0,
             const double& cost = 0,
             const double& cp = 0,
             const double& pa = 0,
             const int& id = -1)
        : Location(x, y), _cost(cost), _cp(cp), _pa(pa), _id(id) {}
    ~Facility() {}

    void show(std::ostream& out = std::cout) const { out << *this << std::endl; }
    void setCost(const double& cost) {
         warning_if_not(cost > 0, "cost should be a positive number instead of " + std::to_string(cost));
        _cost = cost;
    }
    void setEstUtil(const double& estUtil) {
        warning_if_not(estUtil > 0, "estimated utility should be a positive number instead of " + std::to_string(estUtil));
       _estUtil = estUtil;
    }
    inline const double& cost() const { return _cost; }
    inline const double& cp() const { return _cp; }
    inline const double& pa() const { return _pa; }
    inline const int& id() const { return _id; }
    inline const double& estUtil() const { return _estUtil; }

private:
    double _cost = 0.0; // cost
    double _cp = 0.0; // commodity price
    double _pa = 0.0; // parking attraction
    int _id = -1;
    double _estUtil = 0.0;
};

#endif // FACILITY_H
