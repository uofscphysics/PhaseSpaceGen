#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include "TGenPhaseSpace.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "constants.h"
#include "physics.h"
//#include "radcorr.h"
#ifdef PLOTS
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#endif

int main(int argc, char *argv[]) {
  if (argc < 6) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0]
              << " std::string file_name, int gen_num = 10000, float energy=10.6, float q2_min=0.0, float q2_max=12.0"
              << std::endl;
    exit(1);
  }
  std::string file_name = argv[1];
  long long gen_num = atoll(argv[2]);
  float energy = atof(argv[3]);
  float q2_min = atof(argv[4]);
  float q2_max = atof(argv[5]);

#ifdef PLOTS
  TH1D *W_hist = new TH1D("w", "w", 500, 0.0, energy);
  TH1D *missMass = new TH1D("MM", "MM", 500, 0.0, energy);

  TH2D *WvsQ2 = new TH2D("wvsq2", "wvsq2", 500, 0.0, energy, 500, q2_min, q2_max);
#endif

  std::ofstream myfile(file_name);
  if (!myfile.is_open()) {
    std::cerr << "Did not open file " << file_name << std::endl;
    exit(1);
  }

  TLorentzVector target(0.0, 0.0, 0.0, MASS_P);
  TLorentzVector beam(0.0, 0.0, energy, energy);
  TLorentzVector cms = beam + target;

  //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses[3] = {MASS_E, MASS_PIP, MASS_N};

  // Get random seed to set randomness in TRandom3 for TGenPhaseSpace
  std::mt19937_64 prng;
  auto seed = std::random_device{}();
  prng.seed(seed);
  delete gRandom;
  auto TRandSeed = gRandom = new TRandom3(prng());
  auto event = std::unique_ptr<TGenPhaseSpace>(new TGenPhaseSpace());

  event->SetDecay(cms, 3, masses);
  int n = 0;
  int total = 0;
  while (n < gen_num) {
    Double_t weight = event->Generate();
    auto Eprime = event->GetDecay(0);
    auto Pip = event->GetDecay(1);
    auto Neut = event->GetDecay(2);

    double W = physics::W_calc(beam, *Eprime);
    double Q2 = physics::Q2_calc(beam, *Eprime);

    if (Q2 > q2_min && Q2 < q2_max) {
#ifdef PLOTS
      W_hist->Fill(W);
      WvsQ2->Fill(W, Q2);
      TLorentzVector *mm = new TLorentzVector(cms);
      *mm -= *Eprime;
      *mm -= *Pip;
      missMass->Fill(mm->M());
#endif

      if (n++ % 1000 == 0) std::cout << "\t" << n << "\r" << std::flush;
      myfile << "\t3 0.93827231 1 0 1 11 " << energy << " 2212 0 " << weight << std::endl;
      myfile << "1 0 1 11 0 0 " << Eprime->Px() << " " << Eprime->Py() << " " << Eprime->Pz() << " " << Eprime->E()
             << " " << Eprime->M() << " 0 0 0" << std::endl;
      myfile << "2 0 1 211 0 0 " << Pip->Px() << " " << Pip->Py() << " " << Pip->Pz() << " " << Pip->E() << " "
             << Pip->M() << " 0 0 0" << std::endl;
      myfile << "3 0 1 2112 0 0 " << Neut->Px() << " " << Neut->Py() << " " << Neut->Pz() << " " << Neut->E() << " "
             << Neut->M() << " 0 0 0" << std::endl;
    }

    if (total++ > 5 * gen_num) {
      std::cerr << "[" << __FUNCTION__ << "] Ended with break";
      break;
    }
  }
  myfile << std::endl;
  myfile.close();
  std::cout << gen_num << " " << total << " " << n << std::endl;

#ifdef PLOTS
  auto f = new TFile(Form("%s.root", argv[0]), "RECREATE");
  f->cd();
  W_hist->Write();
  missMass->Write();
  WvsQ2->SetOption("COLZ");
  WvsQ2->Write();
  f->Write();
#endif

  exit(0);
}
