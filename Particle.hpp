#ifndef PARTICLE_HPP
#define PARTICLE_HPP

//#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include <array>
#include <cmath> //

class Particle
{
public:
    Particle(const char *name, double const FPx, double const FPy, double const FPz);
    Particle();

    //static void ArrayPrint();
    void ParticlePrint();

    double getFPx()const;
    double getFPy()const;
    double getFPz()const;
    int getfIndex() const;

    double GetMass() const;
    double GetCharge() const;
    double Energy() const;
    /**/
    double InvMass(Particle &p) const;
    /**/
    void SetP(double px, double py, double pz);
    void SetfIndex(int newIndex);
    void SetfIndex(const char *name);

    //static void AddParticleType(int);
    static void AddParticle(const char* name, double mass,
                            int charge, double width = 0);

    int Decay2body(Particle &dau1, Particle &dau2) const;
    void Boost(double bx, double by, double bz);
    void PrintArray();
    void PrintParticle() const; 

    //void setFPx(double px);
    //void setFPy(double py);
    //void setFPz(double pz);

private:
    static const int fMaxNumParticleType = 10;
    static std::array<ParticleType *, fMaxNumParticleType> fParticleType;
    // definisco un array statico con fMNPT elementi di tipo PT, lo chiamo fPT
    static int fNParticleType; // ciclo while in cpp
    int fIParticle = -1; 
    int fIndex;
    static int FindParticle(const char *name);

    double fPx = 0;
    double fPy = 0;
    double fPz = 0;
};
#endif