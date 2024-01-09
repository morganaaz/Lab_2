#ifndef PARTICLETYPE_HPP
#define PARTICLETYPE_HPP

#include <iostream>
//#include "ResonanceType.hpp"
class ParticleType{

   public:
    ParticleType(const char* name, const double mass, const int charge);
    const char* GetName()const;
    double GetMass()const;
    int GetCharge()const;
    virtual void Print()const;
    virtual double GetWidth() const;

    private:
    const char* fName_;
    double const fMass_;
    int const fCharge_;
    
};
#endif