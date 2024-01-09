#ifndef RESONANCETYPE_HPP
#define RESONANCETYPE_HPP

#include <iostream>
#include "ParticleType.hpp"

class ResonanceType: virtual public ParticleType{
    public:
    
    double GetWidth() const;
    void Print() const;

    ResonanceType(const char* name, const double mass, const int charge,
    const double width);

    private:
    double const fWidth_;
    

};
#endif