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
  TH1D *W_rad = new TH1D("w_rad", "w_rad", 500, 0.0, energy);

  TH2D *WvsQ2 = new TH2D("wvsq2", "wvsq2", 500, 0.0, energy, 500, q2_min, q2_max);
  TH2D *WvsQ2_rad = new TH2D("wvsq2_rad", "wvsq2_rad", 500, 0.0, energy, 500, q2_min, q2_max);
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
  Double_t masses[2] = {MASS_E, MASS_P};

  // Get random seed to set randomness in TRandom3 for TGenPhaseSpace
  std::mt19937_64 prng;
  auto seed = std::random_device{}();
  prng.seed(seed);
  delete gRandom;
  auto TRandSeed = gRandom = new TRandom3(prng());
  auto event = std::make_unique<TGenPhaseSpace>();

  event->SetDecay(cms, 2, masses);
  bool radiative_corrections = false;
  int n = 0;
  int total = 0;
  while (n < gen_num) {
    Double_t weight = event->Generate();
    auto Eprime = event->GetDecay(0);
    auto Proton = event->GetDecay(1);

    double W = physics::W_calc(beam, *Eprime);
    double Q2 = physics::Q2_calc(beam, *Eprime);

    if (Q2 > q2_min && Q2 < q2_max) {
#ifdef PLOTS
      W_hist->Fill(W);
      WvsQ2->Fill(W, Q2);
#endif

      if (n++ % 1000 == 0) std::cout << "\t" << n << "\r" << std::flush;
      myfile << "\t2 0.93827231 1 0 1 11 " << energy << " 2212 0 " << weight << std::endl;
      myfile << "1 0 1 11 0 0 " << Eprime->Px() << " " << Eprime->Py() << " " << Eprime->Pz() << " " << Eprime->E()
             << " " << Eprime->M() << " 0 0 0" << std::endl;
      myfile << "2 0 1 2212 0 0 " << Proton->Px() << " " << Proton->Py() << " " << Proton->Pz() << " " << Proton->E()
             << " " << Proton->M() << " 0 0 0" << std::endl;
    }
    if (radiative_corrections) {
      std::cerr << "Not implemented yet" << std::endl;
      exit(1);
      // auto rad = std::make_shared<RadCorr>(energy, Eprime, 0.001);
      // Eprime = rad->corrected();
      double W = physics::W_calc(beam, *Eprime);
      double Q2 = physics::Q2_calc(beam, *Eprime);

#ifdef PLOTS
      W_rad->Fill(W);
      WvsQ2_rad->Fill(W, Q2);
#endif

      if (Q2 > q2_min && Q2 < q2_max) {
        myfile << "\t2 0.93827231 1 0 1 11 " << energy << " 2212 0 " << weight << std::endl;
        myfile << "1 0 1 11 0 0 " << Eprime->Px() << " " << Eprime->Py() << " " << Eprime->Pz() << " " << Eprime->E()
               << " " << Eprime->M() << " 0 0 0" << std::endl;
        myfile << "2 0 1 2212 0 0 " << Proton->Px() << " " << Proton->Py() << " " << Proton->Pz() << " " << Proton->E()
               << " " << Proton->M() << " 0 0 0" << std::endl;
      }
    }

    if (total++ > 20 * gen_num) {
      std::cerr << "[" << __FUNCTION__ << "] Ended with break";
      break;
    }
  }
  myfile << std::endl;
  myfile.close();
  std::cout << gen_num << " " << total << " " << n << std::endl;

#ifdef PLOTS
  auto f = new TFile("ElasticGen.root", "RECREATE");
  f->cd();
  W_hist->Write();
  WvsQ2->Write();
  if (radiative_corrections) {
    W_rad->Write();
    WvsQ2_rad->SetOption("COLZ");
    WvsQ2_rad->Write();
  }
  f->Write();
#endif

  exit(0);
}
