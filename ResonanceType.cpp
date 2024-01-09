#include "ResonanceType.hpp"
#include "ParticleType.hpp"
#include <iostream>

ResonanceType::ResonanceType(const char* name, const double mass,
                             const int charge, const double width) :  ParticleType(name, mass, charge), fWidth_(width){}

double ResonanceType::GetWidth()const{
    return fWidth_;
}
void ResonanceType::Print()const {
    ParticleType::Print();
    std::cout<<"Particle Width:\n"<<GetWidth()<<'\n';
}