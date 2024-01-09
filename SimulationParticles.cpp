#include <iostream>
#include "/home/morganaartemisia/root/include/TH1.h"
#include "/home/morganaartemisia/root/include/TMath.h"
#include "/home/morganaartemisia/root/include/TRandom.h"
#include "/home/morganaartemisia/root/include/TFile.h"
#include "Particle.hpp"

/*
    compilare con: g++ -o SimulationParticles.exe -Wall -Wextra Simulation_particles.cpp ParticleType.cpp ResonanceType.cpp Particle.cpp `root-config --cflags --libs`

    eseguire con: ./SimulationParticles.exe
*/

int main()
{
    std::cout<<'\n';

    Int_t nGen = 100000;
    Int_t nParticlesPerGen = 100;
    Int_t dimArrayParticles = nParticlesPerGen + (Int_t)(0.2 * nParticlesPerGen);
    Double_t theta, phi, p, x;
    gRandom-> SetSeed();

    Particle::AddParticle("pi+",0.13957,1);
    Particle::AddParticle("pi-",0.13957,-1);
    Particle::AddParticle("k+",0.49367,1);
    Particle::AddParticle("k-",0.49367,-1);
    Particle::AddParticle("p+",0.93827,1);
    Particle::AddParticle("p-",0.93827,-1);
    Particle::AddParticle("k*",0.89166,0,0.050);

    //Creo gli istogrammi
    TH1F *hParticleType_distribution = new TH1F ("hParticleType_distribution", "Distribution of particle types",7, 0, 7);
    TH1F *hTheta_distribution = new TH1F ("hTheta_distribution", "Theta istribution (polar coordinate)",1000, 0, TMath::Pi());
    hTheta_distribution->SetMinimum(0);    
    hTheta_distribution->SetMaximum(nGen*.2);
    TH1F *hPhi_distribution = new TH1F ("hPhi_distribution", "Phi distribution (azimutal coordinate)",1000, 0, 2*TMath::Pi());
    hTheta_distribution->SetMinimum(0);    
    hTheta_distribution->SetMaximum(nGen*.2);
    TH1D *hMomentum_distribution = new TH1D("hMomentum_distribution", "Momentum distribution", 1000, 0, 10);
    TH1D *hTrasverseMomentum_distribution = new TH1D ("hTrasverseMomentum_distribution", "Transverse momentum distribution", 1000, 0, 10);
    TH1D *hEnergy_distribution = new TH1D ("hEnergy_distribution", "Energy distribution", 1000, 0, 10);
    TH1D *hTotInvariantMass_distribution = new TH1D ("hTotInvariantMass_distribution", "Invariant mass distribution of all generated particles", 80, 0, 2);
    TH1D *hInvariantMassD_distribution = new TH1D ("hInvariantMassD_distribution", "Invariant mass distribution of all particles with discordant charge", 80, 0, 2);
    TH1D *hInvariantMassDPiK_distribution = new TH1D ("hInvariantMassDPiK_distribution", "Invariant mass distribution of all pi-k particles with discordant charge", 80, 0, 2);
    TH1D *hInvariantMassC_distribution = new TH1D ("hInvariantMassC_distribution", "Invariant mass distribution of all particles with concordant charge", 80, 0, 2);
    TH1D *hInvariantMassCPiK_distribution = new TH1D ("hInvariantMassCPiK_distribution", "Invariant mass distribution of all pi-k particles with concordant charge", 80, 0, 2);
    TH1D *hInvariantMassDecay_distribution = new TH1D ("hInvariantMassDecay_distribution", "Invariant mass distribution of k* particles decay", 80, 0, 2);
    hInvariantMassDecay_distribution -> Sumw2();
    Particle simulation[dimArrayParticles];
    std::cout << "\nLoading simulation\n";
    
    for(Int_t i = 0; i < nGen; i++)
    {
        Int_t resonance = 0;
        for(Int_t j = 0; j < nParticlesPerGen; j++)
        {
            x = gRandom -> Rndm();
            if(x < 0.4)
                simulation[j].SetfIndex(0);
            else if(x < 0.8)
                simulation[j].SetfIndex(1);
            else if(x < 0.85)
                simulation[j].SetfIndex(2);
            else if(x < 0.90)
                simulation[j].SetfIndex(3);
            else if(x < 0.945)
                simulation[j].SetfIndex(4);
            else if(x < 0.99)
                simulation[j].SetfIndex(5);
            else
            {
                simulation[j].SetfIndex(6);
                if (x < 0.995)
                {
                    simulation[nParticlesPerGen + resonance].SetfIndex(2);
                    ++resonance;
                    simulation[nParticlesPerGen + resonance].SetfIndex(1);
                    ++resonance;
                    
                }
                else
                {
                    simulation[nParticlesPerGen + resonance].SetfIndex(3);
                    ++resonance;
                    simulation[nParticlesPerGen + resonance].SetfIndex(0);
                    ++resonance;                   
                }
                simulation[j].Decay2body(simulation[nParticlesPerGen + resonance - 2], simulation[nParticlesPerGen + resonance - 1]);
            }
            hParticleType_distribution -> Fill(simulation[j].getfIndex());
            theta = gRandom -> Uniform(0, TMath::Pi());
            hTheta_distribution -> Fill(theta);
            phi = gRandom -> Uniform(0, 2*TMath::Pi());
            hPhi_distribution -> Fill(phi);
            p = gRandom -> Exp(1);
            hMomentum_distribution -> Fill(p);
            simulation[j].SetP(p*TMath::Sin(theta)*TMath::Cos(phi),
                               p*TMath::Sin(theta)*TMath::Sin(phi),
                               p*TMath::Cos(theta));
            hTrasverseMomentum_distribution -> Fill(sqrt(pow(simulation[j].getFPx(),2)
                                                    + pow(simulation[j].getFPy(),2)));
            hEnergy_distribution -> Fill(simulation[j].Energy());            
        }
        for (Int_t particle_1 = 0; particle_1 < nParticlesPerGen + resonance; particle_1++)
        {
            for (Int_t particle_2 = particle_1 + 1; particle_2 < nParticlesPerGen + resonance; particle_2++)
            {
                hTotInvariantMass_distribution -> Fill(simulation[particle_1].InvMass(simulation[particle_2]));
                if (simulation[particle_1].GetCharge()*simulation[particle_2].GetCharge() < 0.0)
                {
                    hInvariantMassD_distribution -> Fill(simulation[particle_1].InvMass(simulation[particle_2]));
                    if ((simulation[particle_1].getfIndex() <= 1 &&
                         simulation[particle_2].getfIndex() > 1 &&
                         simulation[particle_2].getfIndex() <= 3) ||
                        (simulation[particle_2].getfIndex() <= 1 &&
                         simulation[particle_1].getfIndex() > 1 &&
                         simulation[particle_1].getfIndex() <= 3))

                    hInvariantMassDPiK_distribution -> Fill(simulation[particle_1].InvMass(simulation[particle_2]));
                }
                else if (simulation[particle_1].GetCharge()*simulation[particle_2].GetCharge() > 0.0)
                {
                    hInvariantMassC_distribution -> Fill(simulation[particle_1].InvMass(simulation[particle_2]));
                    if ((simulation[particle_1].getfIndex() <= 1 &&
                         simulation[particle_2].getfIndex() > 1 &&
                         simulation[particle_2].getfIndex() <= 3) ||
                        (simulation[particle_2].getfIndex() <= 1 &&
                         simulation[particle_1].getfIndex() > 1 &&
                         simulation[particle_1].getfIndex() <= 3))
                    hInvariantMassCPiK_distribution -> Fill(simulation[particle_1].InvMass(simulation[particle_2]));
                }
            }
        }
        for (Int_t particle_1 = 0; particle_1 < resonance; particle_1 = particle_1 + 2)
        {
            hInvariantMassDecay_distribution -> Fill(simulation[nParticlesPerGen + particle_1].InvMass(simulation[nParticlesPerGen + 1 + particle_1]));
        }
    }
    TFile * theFile = new TFile("Particles_simulation.root","RECREATE");

    //Graphics
    hParticleType_distribution -> GetXaxis() -> SetTitle("Particle types ");
    hParticleType_distribution -> GetXaxis() -> SetBinLabel(1,"pi+");
    hParticleType_distribution -> GetXaxis() -> SetBinLabel(2,"pi-");
    hParticleType_distribution -> GetXaxis() -> SetBinLabel(3,"k+");
    hParticleType_distribution -> GetXaxis() -> SetBinLabel(4,"k-");
    hParticleType_distribution -> GetXaxis() -> SetBinLabel(5,"p+");
    hParticleType_distribution -> GetXaxis() -> SetBinLabel(6,"p-");
    hParticleType_distribution -> GetXaxis() -> SetBinLabel(7,"k*");
    hParticleType_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hParticleType_distribution -> Write();

    hTheta_distribution -> GetXaxis() -> SetTitle("Theta (rad)");
    hTheta_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hTheta_distribution -> Write();

    hPhi_distribution -> GetXaxis() -> SetTitle("Phi (rad)");
    hPhi_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hPhi_distribution -> Write();

    hEnergy_distribution -> GetXaxis() -> SetTitle("Energy (GeV)");
    hEnergy_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hEnergy_distribution -> Write();

    hMomentum_distribution -> GetXaxis() -> SetTitle("Momentum (GeV)");
    hMomentum_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hMomentum_distribution -> Write();

    hTrasverseMomentum_distribution -> GetXaxis() -> SetTitle("Trasverse momentum (GeV)");
    hTrasverseMomentum_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hTrasverseMomentum_distribution -> Write();

    hTotInvariantMass_distribution -> GetXaxis() -> SetTitle("Invariant mass (GeV/c^2)");
    hTotInvariantMass_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hTotInvariantMass_distribution -> Write();

    hInvariantMassD_distribution -> GetXaxis() -> SetTitle("Invariant mass (GeV/c^2)");
    hInvariantMassD_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hInvariantMassD_distribution -> Write();

    hInvariantMassDPiK_distribution -> GetXaxis() -> SetTitle("Invariant mass (GeV/c^2)");
    hInvariantMassDPiK_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hInvariantMassDPiK_distribution -> Write();

    hInvariantMassC_distribution -> GetXaxis() -> SetTitle("Invariant mass (GeV/c^2)");
    hInvariantMassC_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hInvariantMassC_distribution -> Write();

    hInvariantMassCPiK_distribution -> GetXaxis() -> SetTitle("Invariant mass (GeV/c^2)");
    hInvariantMassCPiK_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hInvariantMassCPiK_distribution -> Write();

    hInvariantMassDecay_distribution -> GetXaxis() -> SetTitle("Invariant mass (GeV/c^2)");
    hInvariantMassDecay_distribution -> GetYaxis() -> SetTitle("Occurrences");
    hInvariantMassDecay_distribution -> Write();

    theFile -> Close();

    std::cout << "\nGeneration completed successfully\n";

}