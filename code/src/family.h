#ifndef FAMILY_H
#define FAMILY_H
#include "location.h"

class Family : public Location {

    friend std::ostream& operator << (std::ostream& out, const Family& fam) {
        const Location loc = static_cast<Location>(fam);
        out << "{loc: " << loc << ", type: " << _ftName.at(fam._type) << ", bp: " << fam._bp << "}";
        return out;
    }

public:
    Family(const double& x = 0, const double& y = 0, const FAMILY_TYPE type = FAMILY_UNINIT,
           const double& bp = 0)
        : Location(x, y), _type(type), _bp(bp) {}
    ~Family() {}

    void show(std::ostream& out = std::cout) const { out << *this << std::endl; }
    inline const FAMILY_TYPE& type() const { return _type; }
    inline const double& bp() const { return _bp; }

private:
    FAMILY_TYPE _type = FAMILY_UNINIT;
    double _bp = 0.0; // buying power
};

#endif // FAMILY_H
