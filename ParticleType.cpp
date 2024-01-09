#include "ParticleType.hpp"
#include <iostream>

//const tutte le variabili
ParticleType::ParticleType(const char* name, const double mass, const int charge) : fName_((char* const)name),
                                     fMass_(mass), fCharge_(charge){}
const char* ParticleType::GetName()const {return fName_;}
double ParticleType::GetMass()const {return fMass_;}
int ParticleType::GetCharge()const {return fCharge_;}
double ParticleType::GetWidth()const {return 0;}
void ParticleType::Print()const  {
    std::cout<<"Particle Name:\n"<<GetName()<<'\n';
    std::cout<<"Particle Mass:\n"<<GetMass()<<'\n';
    std::cout<<"Particle Charge:\n"<<GetCharge()<<'\n';
}