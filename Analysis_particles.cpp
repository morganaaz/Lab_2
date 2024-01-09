#include <iostream>
#include "/home/morganaartemisia/root/include/TFile.h"
#include "/home/morganaartemisia/root/include/TH1F.h"
#include "/home/morganaartemisia/root/include/TStyle.h"
#include "/home/morganaartemisia/root/include/TCanvas.h"
#include "/home/morganaartemisia/root/include/TPaveStats.h"

/*
    compilare con: g++ -o Analysis_particles.exe -Wall -Wextra Analysis_particles.cpp `root-config --cflags --libs`

    eseguire con: ./Analysis_particles.exe
*/

int main()
{
    Int_t nGen = 100000;
    TFile *theFile = new TFile("Particles_simulation.root");
    TH1F* hParticleTypes_distribution = (TH1F*) theFile -> Get("hParticleType_distribution");
    TH1F* hTheta_distribution = (TH1F*) theFile -> Get("hTheta_distribution");
    TH1F* hPhi_distribution = (TH1F*) theFile -> Get("hPhi_distribution");
    TH1F* hEnergy_distribution = (TH1F*) theFile -> Get("hEnergy_distribution");
    TH1F* hMomentum_distribution = (TH1F*) theFile -> Get("hMomentum_distribution");
    TH1F* hTrasverseMomentum_distribution = (TH1F*) theFile -> Get("hTrasverseMomentum_distribution");
    TH1F* hTotInvaiantMass_distribution = (TH1F*) theFile -> Get("hTotInvaiantMass_distribution");
    TH1F* hInvaiantMassD_distribution = (TH1F*) theFile -> Get("hInvaiantMassD_distribution");
    TH1F* hInvaiantMassDPiK_distribution = (TH1F*) theFile -> Get("hInvaiantMassDPiK_distribution");
    TH1F* hInvaiantMassC_distribution = (TH1F*) theFile -> Get("hInvaiantMassC_distribution");
    TH1F* hInvaiantMassCPiK_distribution = (TH1F*) theFile -> Get("hInvaiantMassCPiK_distribution");
    TH1F* hInvariantMassDecay_distribution = (TH1F*) theFile -> Get("hInvariantMassDecay_distribution");

    std::cout << "Types proportions\n" << "Expected entries\n"
              << "pi+ (1) " << 0.4*nGen*100 << '\n'
              << "pi- (2) " << 0.4*nGen*100 << '\n'
              << "k+ (3) " << 0.05*nGen*100 << '\n'
              << "k- (4) " << 0.05*nGen*100 << '\n'
              << "p+ (5) " << 0.045*nGen*100 << '\n'
              << "p- (6) " << 0.045*nGen*100 << '\n'
              << "k* (7) " << 0.01*nGen*100 << '\n'
              << "Entries of each type \n";
        for (Int_t i = 1; i < 8; i++)
        {
            std::cout<< "Bin " << i << " (" << hParticleTypes_distribution -> GetBinContent(i)
                                        << " ± " << hParticleTypes_distribution -> GetBinError(i) << ")\n";
        }

        std::cout << "Entries of each histogram \n"
        << " Particle Types " << hParticleTypes_distribution -> GetEntries() << '\n'
        << " Theta " << hTheta_distribution -> GetEntries() << '\n'
        << " Phi " << hPhi_distribution -> GetEntries() << '\n'
        << " Momentum " << hMomentum_distribution -> GetEntries() << '\n'
        << " Trasversal Momentum " << hTrasverseMomentum_distribution -> GetEntries() << '\n'
        << " Energy " << hEnergy_distribution -> GetEntries() << '\n'
        << " Total invariant mass " << hTotInvaiantMass_distribution -> GetEntries() << '\n'
        << " Invariant mass of total particles with discordant charge " << hInvaiantMassD_distribution -> GetEntries() << '\n'
        << " Invariant mass of Pi and K particles with discordant charge " << hInvaiantMassDPiK_distribution -> GetEntries() << '\n'
        << " Invariant mass of total particles with concordant charge " << hInvaiantMassC_distribution -> GetEntries() << '\n'
        << " Invariant mass of Pi and K particles with concordant charge " << hInvaiantMassCPiK_distribution -> GetEntries() << '\n'
        << " Invariant mass decay " << hInvariantMassDecay_distribution -> GetEntries() << '\n';

        gStyle -> SetOptStat(1111);
        gStyle -> SetOptFit(1);

        TCanvas *cDistributions = new TCanvas("cDistributions", "Distributions");
        cDistributions -> Divide(2,2);
        cDistributions -> cd(1);
        hParticleTypes_distribution -> Draw("LEGO2z");
        std::cout<< "Theta and Phi fit:\n";
        cDistributions -> cd(2);
        hMomentum_distribution -> Fit("expo");
        hMomentum_distribution -> Draw();
        cDistributions -> cd(3);
        hTheta_distribution -> Fit("pol0");
        hTheta_distribution -> Draw();
        cDistributions -> cd(4);
        hPhi_distribution -> Fit("pol0");
        hPhi_distribution -> Draw();
        cDistributions -> Print("Particle_distributions_analysis.pdf");
        cDistributions -> Print("Particle_distributions_analysis.C");
        cDistributions -> Print("Particle_distributions_analysis.root");

        TCanvas *cMomentum = new TCanvas("cMomentum", "Momentum");
        cMomentum -> Divide(2,1);
        std::cout << "\nMomentum fit\n";
        cMomentum -> cd(1);
        hMomentum_distribution -> Fit("expo");
        hMomentum_distribution -> Draw();
        std::cout<< "\nExpected mean 1 GeV\n" << "Data mean (" << hMomentum_distribution -> GetMean()
                 << " ± " << hMomentum_distribution -> GetMeanError() << ") GeV\n";
        cMomentum -> cd(2);
        hTrasverseMomentum_distribution -> Draw();
        cMomentum -> Print("Momentum_distributions_analysis.pdf");
        cMomentum -> Print("Momentum_distributions_analysis.C");
        cMomentum -> Print("Momentum_distributions_analysis.root");

        TCanvas *cInvariantMass = new TCanvas("cInvariantMass", "Invariant mass");
        cInvariantMass -> Divide(3,2);
        cInvariantMass -> cd(1);
        hTotInvaiantMass_distribution -> Draw();
        cInvariantMass -> cd(2);
        hInvaiantMassD_distribution -> Draw();
        cInvariantMass -> cd(3);
        hInvaiantMassDPiK_distribution -> Draw();
        cInvariantMass -> cd(4);
        hInvaiantMassC_distribution -> Draw();
        cInvariantMass -> cd(5);
        hInvaiantMassCPiK_distribution -> Draw();
        cInvariantMass -> cd(6);
        hInvariantMassDecay_distribution -> Draw();
        cInvariantMass -> Print("Invariant_mass_analysis.pdf");

        TH1F *hDifferencePiK_distribution = new TH1F("hDifferencePiK_distribution",
        "Invariant mass difference between combinations of discordant and concordant pi-k particles", 80, 0, 2);
        hDifferencePiK_distribution -> Sumw2();
        hDifferencePiK_distribution -> Add(hInvaiantMassDPiK_distribution,hInvaiantMassCPiK_distribution,1,-1);
        hDifferencePiK_distribution -> SetMinimum(0);
        hDifferencePiK_distribution -> GetXaxis() -> SetTitle("Invariant mass (GeV/c^2)");
        hDifferencePiK_distribution -> GetYaxis() -> SetTitle("Occurrences");

        TH1F *hDifferenceCD_distribution = new TH1F("hDifferenceCD_distribution",
        "Invariant mass difference between combinations of discordant and concordant particles", 80, 0, 2);
        hDifferenceCD_distribution -> Sumw2();
        hDifferenceCD_distribution -> Add(hInvaiantMassD_distribution,hInvaiantMassC_distribution,1,-1);
        hDifferenceCD_distribution -> SetMinimum(0);
        hDifferenceCD_distribution -> GetXaxis() -> SetTitle("Invariant mass (GeV/c^2)");
        hDifferenceCD_distribution -> GetYaxis() -> SetTitle("Occurrences");
        TCanvas *cDecay = new TCanvas("cDecay", "k* decay");
        cDecay -> Divide(1,3);
        cDecay -> cd(1);
        std::cout<< "\nInvariant mass decay\n";
        hInvariantMassDecay_distribution -> Fit("gaus");
        hInvariantMassDecay_distribution -> SetStats(0);
        hInvariantMassDecay_distribution -> SetStats(1);
        hInvariantMassDecay_distribution -> Draw();
        cDecay -> cd(2);
        std::cout<< "\nInvariant mass difference between combinations of discordant and concordant particles\n";
        hDifferenceCD_distribution -> Fit("gaus");
        hDifferenceCD_distribution -> SetEntries(hDifferenceCD_distribution -> Integral());
        hDifferenceCD_distribution -> Draw();
        cDecay -> cd(3);
        std::cout<< "\nInvariant mass difference between combinations of discordant and concordant pi-k particles\n";
        hDifferencePiK_distribution -> Fit("gaus");
        hDifferencePiK_distribution -> SetEntries(hDifferencePiK_distribution -> Integral());
        hDifferencePiK_distribution -> Draw();
        std::cout << '\n';
        cDecay -> Print("Invariant_mass_decay_analysis.pdf");
        cDecay -> Print("Invariant_mass_decay_analysis.C");
        cDecay -> Print("Invariant_mass_decay_analysis.root");
        
        theFile -> Close();
    
}