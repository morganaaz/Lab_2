#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "Particle.hpp"

#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdlib> //for RAND_MAX

int Particle::fNParticleType = 0;

// essendo static è necessaria una def out of line
std::array<ParticleType *, Particle::fMaxNumParticleType> Particle::fParticleType;

Particle::Particle(const char *name, double const FPx, double const FPy, double const FPz)
    : fPx(FPx), fPy(FPy), fPz(FPz), fIParticle(FindParticle(name))
{

    fIParticle = FindParticle(name);
    if (fIParticle == -1)
    {
        std::cout << "No corresponding type found" << '\n';
    }
}

Particle::Particle() : fIParticle(-2), fPx(0), fPy(0), fPz(0) {}

int Particle::FindParticle(const char *findName)
{
    for (int i = 0; i < fNParticleType; i++)
    {
        if (findName == (fParticleType[i]->GetName()))
            return i;
    }
    return -1;
}

void Particle::SetfIndex(const char *name)
{
    if ((Particle::FindParticle(name)) < 0)
    {
        std::cout << "Name not found" << '\n';
    }
    else
    {
        fIParticle = FindParticle(name);
    }
}

void Particle::SetfIndex(int newIndex)
{
    if (newIndex < 0)
    {
        std::cout << "Index must be positive" << '\n';
    }
    else
    {
        fIParticle = newIndex;
    }

    /*int Particle::FindParticle(const char *name) //
                            //sost name con findn per chiarezza
        {
        int index = 13;
        static int fIndexParticle = FindParticle(findName);
        static int fNParticleType = 0;

        while (fNParticleType < fMaxNumParticleType)
        {
            // if(name==fParticleType[i]->getName()){
            //== compara gli indirizzi di memoria delle variabili, non il contenuto
            if (strcmp(findName, fParticleType[fNParticleType]->getName()) == findName)
            {
                index = fNParticleType;
                break;
            }

            ++fNParticleType;
        }
        */

    // return(index);
    if (fIParticle == 13)
    {
        std::cout << "No corresponding type found" << '\n';
    }
}

// APT prende l'output di FP e se risulta una nuova particella,
// l'aggiunge all'array
// non prende nulla in input

void Particle::AddParticle(const char *newName, double mass,
                           int charge, double width)
{
    if (Particle::FindParticle(newName) < 0)
    {
        if (fNParticleType < fMaxNumParticleType)
        {
            width == 0 ? fParticleType[fNParticleType] = new ParticleType(newName, mass, charge)
                       : fParticleType[fNParticleType] = new ResonanceType(newName, mass, charge, width);

            fNParticleType++;

            /*const char *newName = fParticleType[fIndex]->getName();
                    if (FindParticle(newName) == 13)
                    {
                        fParticleType[fNParticleType++] = new ParticleType(
                            fParticleType[fIndex]->getName(),
                            fParticleType[fIndex]->getMass(),
                            fParticleType[fIndex]->getCharge()
                            // fParticleType[fIndex]->getWidth()
                        );
            */
            std::cout << "New particle added" << '\n';
        }
        else
        {
            std::cout << "Error: max number of particles reached" << '\n';
        }
    }
    else
    {
        std::cout << "Error: particle type already exists" << '\n';
    }
}

int Particle::getfIndex() const { return fIndex; }

double Particle::getFPx() const { return fPx; };
double Particle::getFPy() const { return fPy; };
double Particle::getFPz() const { return fPz; };

double Particle::GetMass() const
{
    return fParticleType[fIndex]->GetMass();
};

double Particle::GetCharge() const
{
    return fParticleType[fIndex]->GetCharge();
};

double Particle::Energy() const
{
    double energy = sqrt(GetMass() * GetMass() + pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2));
    return energy;
}

/**/double Particle::InvMass(Particle &p)const
    {
        return sqrt(pow(Energy() + p.Energy(), 2) - pow((fPx + p.getFPx()), 2)
                                                  + pow((fPy + p.getFPy()), 2)
                                                  + pow((fPz + p.getFPz()), 2));
    };
/**/

void Particle::SetP(double px, double py, double pz)
{
    fPx = px;
    fPy = py;
    fPz = pz;
}

// aggiunta da file Virtuale controllare variabili
int Particle::Decay2body(Particle &dau1, Particle &dau2) const
{
    if (GetMass() == 0.0)
    {
        printf("Decayment cannot be preformed if mass is zero\n");
        return 1;
    }

    double massMot = GetMass();
    double massDau1 = dau1.GetMass();
    double massDau2 = dau2.GetMass();

    if (fIParticle > -1)
    { // add width effect

        // gaussian random numbers

        float x1, x2, w, y1;

        double invnum = 1. / RAND_MAX;
        do
        {
            x1 = 2.0 * rand() * invnum - 1.0;
            x2 = 2.0 * rand() * invnum - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);

        w = sqrt((-2.0 * log(w)) / w);
        y1 = x1 * w;

        massMot += fParticleType[fIParticle]->GetWidth() * y1;
    }

    if (massMot < massDau1 + massDau2)
    {
        printf("Decayment cannot be preformed because mass is too low in this channel\n");
        return 2;
    }

    double pout = sqrt((massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) * (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) / massMot * 0.5;

    double norm = 2 * M_PI / RAND_MAX;

    double phi = rand() * norm;
    double theta = rand() * norm * 0.5 - M_PI / 2.;
    dau1.SetP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi), pout * cos(theta));
    dau2.SetP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi), -pout * cos(theta));

    double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

    double bx = fPx / energy;
    double by = fPy / energy;
    double bz = fPz / energy;

    dau1.Boost(bx, by, bz);
    dau2.Boost(bx, by, bz);

    return 0;
}
void Particle::Boost(double bx, double by, double bz)
{

    double energy = Energy();

    // Boost this Lorentz vector
    double b2 = bx * bx + by * by + bz * bz;
    double gamma = 1.0 / sqrt(1.0 - b2);
    double bp = bx * fPx + by * fPy + bz * fPz;
    double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

    fPx += gamma2 * bp * bx + gamma * bx * energy;
    fPy += gamma2 * bp * by + gamma * by * energy;
    fPz += gamma2 * bp * bz + gamma * bz * energy;
}

void Particle::PrintArray()
{
    for (int e = 0; e < fNParticleType; ++e)
    {
        std::cout << "Particle: " << e << '\n'
                  << "Index: " << fIParticle << '\n'
                  << "Name: " << fParticleType[e]->GetName() << '\n'
                  << "Mass: " << fParticleType[e]->GetMass() << '\n'
                  << "Charge: " << fParticleType[e]->GetCharge() << '\n';
    }
};

void Particle::PrintParticle() const
{
    for (int e = 0; e < fNParticleType; ++e)
    {
        std::cout << "Particle: " << e << '\n'
                  << "Index: " << fIParticle << '\n'
                  << "Name: " << fParticleType[e]->GetName() << '\n'
                  << "Px: " << fPx << '\n'
                  << "Py: " << fPy << '\n'
                  << "Pz: " << fPz << '\n';
    }
};