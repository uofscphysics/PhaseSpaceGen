#include <chrono>
#include <fstream>
#include <iostream>
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
  TH1D *W_hist = new TH1D("w", "w", 500, 0.0, energy / 2.0);
  TH2D *WvsQ2 = new TH2D("wvsq2", "wvsq2", 500, 0.0, energy / 2.0, 500, q2_min, q2_max);
#endif

  std::ofstream myfile(file_name);

  TLorentzVector target(0.0, 0.0, 0.0, MASS_P);
  TLorentzVector beam(0.0, 0.0, energy, energy);
  TLorentzVector cms = beam + target;

  //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses[4] = {MASS_E, MASS_PIP, MASS_PIM, MASS_PI0};

  auto event = std::make_shared<TGenPhaseSpace>();
  event->SetDecay(cms, 4, masses);
  int n = 0;
  int total = 0;
  while (n < gen_num) {
    Double_t weight = event->Generate();
    auto Eprime = event->GetDecay(0);
    auto Pip = event->GetDecay(1);
    auto Pim = event->GetDecay(2);
    auto Pi0 = event->GetDecay(3);

    double W = physics::W_calc(beam, *Eprime);
    double Q2 = physics::Q2_calc(beam, *Eprime);

    if (Q2 > q2_min && Q2 < q2_max) {
#ifdef PLOTS
      W_hist->Fill(W, weight);
      WvsQ2->Fill(W, Q2, weight);
#endif

      if (n++ % 1000 == 0) std::cout << "\t" << n << "\r" << std::flush;
      myfile << "\t3 0.93827231 1 0 1 11 " << energy << " 2212 0 " << weight << std::endl;
      myfile << "1 0 1 11 0 0 " << Eprime->Px() << " " << Eprime->Py() << " " << Eprime->Pz() << " " << Eprime->E()
             << " " << Eprime->M() << " 0 0 0" << std::endl;
      myfile << "2 0 1 211 0 0 " << Pip->Px() << " " << Pip->Py() << " " << Pip->Pz() << " " << Pip->E() << " "
             << Pip->M() << " 0 0 0" << std::endl;
      myfile << "3 0 1 -211 0 0 " << Pim->Px() << " " << Pim->Py() << " " << Pim->Pz() << " " << Pim->E() << " "
             << Pim->M() << " 0 0 0" << std::endl;
      myfile << "4 0 1 111 0 0 " << Pi0->Px() << " " << Pi0->Py() << " " << Pi0->Pz() << " " << Pi0->E() << " "
             << Pi0->M() << " 0 0 0" << std::endl;
    }

    if (total++ > 5 * gen_num) {
      std::cerr << "[" << __FUNCTION__ << "] Ended with break";
      break;
    }
  }
  myfile << std::endl;
  myfile.close();
  std::cout << n << " " << total << " generated " << 100 * n / total << "\% accepted " << std::endl;

#ifdef PLOTS
  auto f = new TFile(Form("%s.root", argv[0]), "RECREATE");
  f->cd();
  W_hist->Write();
  WvsQ2->SetOption("COLZ");
  WvsQ2->Write();
  f->Write();
#endif

  exit(0);
}
