#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include "TGenPhaseSpace.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TRandomGen.h"
#include "constants.h"
#include "physics.h"
#ifdef PLOTS
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#endif

int main(int argc, char *argv[]) {
  std::string file_name = "";
  long genNum = 10000;
  long seed;
  float energy = 10.6;
  float q2_min = 0.0;
  float q2_max = 45.0;
  float w_min = 0.0;
  float w_max = 7.0;
  bool print_help;
  bool clas12MCgen = false;

  auto cli =
      (clipp::option("-h", "--help").set(print_help) % "print help",
       (clipp::option("--docker").set(clas12MCgen) % "For clas12-mcgen, nothing to do with running in docker.",
        clipp::option("-N", "--trig") &
            clipp::value("genNum", genNum) % "Number of events to generate [Default 10,000]",
        clipp::option("--seed") & clipp::value("seed", seed) % "Random seed number",
        clipp::option("-E", "--energy") & clipp::value("energy", energy) % "Set Beam Energy [Default 10.6]",
        clipp::option("-q2_min", "--q2_min") & clipp::value("q2_min", q2_min) % "Min Q^2 value [Default: 0.0]",
        clipp::option("-q2_max", "--q2_max") & clipp::value("q2_max", q2_max) % "Max Q^2 value [Default: 12.0]",
        clipp::option("-w_min", "--w_min") & clipp::value("w_min", w_min) % "Min W value [Default: 0.0]",
        clipp::option("-w_max", "--w_max") & clipp::value("w_max", w_max) % "Max W value [Default: 12.0]",
        clipp::value("file_name.dat", file_name) % "Filename **Overwritten if using docker option to TwoPi.dat**"));

  clipp::parse(argc, argv, cli);
  if (clas12MCgen) {
    std::cout << "Using clas12-mcgen settings\n";
  } else if (print_help || file_name == "") {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(1);
  }

  // Get random seed to set randomness in TRandom3 for TGenPhaseSpace
  size_t randomSeed;
  if (!clas12MCgen) {
    std::mt19937_64 prng;
    auto deviceSeed = std::random_device{}();
    prng.seed(deviceSeed);
    randomSeed = prng();
  } else {
    file_name = "TwoPi.dat";
    randomSeed = seed;
  }

#ifdef PLOTS
  // If plots is defined let's make the plots too
  TH1D *W_hist = new TH1D("w", "w", 500, 0.0, w_max);
  TH2D *WvsQ2 = new TH2D("wvsq2", "wvsq2", 500, 0.0, w_max, 500, q2_min, q2_max);
#endif

  // Setup target 4 vector
  TLorentzVector target(0.0, 0.0, 0.0, MASS_P);
  // Setup beam 4 vector, E^2 = P^2 + M^2 => P = sqrt(E^2 - M^2)
  TLorentzVector beam(0.0, 0.0, sqrtf(energy * energy - MASS_E * MASS_E), energy);
  TLorentzVector cms = beam + target;

  //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses[4] = {MASS_E, MASS_PIP, MASS_PIM, MASS_P};

  // Open file and if it fails quit with error message
  std::ofstream myfile(file_name);
  if (!myfile.is_open()) {
    std::cerr << "Did not open file " << file_name << std::endl;
    exit(1);
  }

  // Open file and if it fails quit with error message
  std::string csv_name = "twopi.csv";
  std::ofstream csvfile(csv_name);
  if (!csvfile.is_open()) {
    std::cerr << "Did not open file csvfile" << std::endl;
    exit(1);
  }

  csvfile
      << "energy,weight,W,Q2,ep_px,ep_py,ep_pz,ep_e,pip_px,pip_py,pip_pz,pip_e,pim_px,pim_py,pim_pz,pim_e,p_px,p_py,p_"
         "pz,p_e\n";

  // Delete roots default random generator and replace with the
  delete gRandom;
  // Use random based on std::mt19937_64
  auto TRandSeed = gRandom = new TRandomMT64(randomSeed);

  // Make phase space generator
  auto event = std::unique_ptr<TGenPhaseSpace>(new TGenPhaseSpace());
  // Set up phase space conditions

  // Setup number generated and total number accepted to 0's
  size_t acc = 0;
  size_t total = 0;
  size_t per = genNum / 100;
  // Loop until the accepted number of events is greater than number needed aka genNum
  while (acc < genNum) {
    event->SetDecay(cms, 4, masses);
    // Generate new event
    Double_t weight = event->Generate();
    // Get output 4 vectors for the event
    auto Eprime = event->GetDecay(0);
    auto Pip = event->GetDecay(1);
    auto Pim = event->GetDecay(2);
    auto P = event->GetDecay(3);

    // Calcultae W and Q2 for acceptance [and plots]
    double W = physics::W_calc(beam, *Eprime);
    double Q2 = physics::Q2_calc(beam, *Eprime);

    // If it's inside the kinematic region add it to the file
    if (Q2 > q2_min && Q2 < q2_max && W > w_min && W < w_max) {
#ifdef PLOTS
      W_hist->Fill(W, weight);
      WvsQ2->Fill(W, Q2, weight);
#endif

      if (acc++ % per == 0)
        std::cout << "\t" << std::setw(5) << acc / per << "%"
                  << "\r" << std::flush;

      csvfile << std::scientific;
      csvfile << energy << ",";
      csvfile << weight << ",";
      csvfile << W << ",";
      csvfile << Q2 << ",";
      csvfile << Eprime->Px() << ",";
      csvfile << Eprime->Py() << ",";
      csvfile << Eprime->Pz() << ",";
      csvfile << Eprime->E() << ",";
      csvfile << Pip->Px() << ",";
      csvfile << Pip->Py() << ",";
      csvfile << Pip->Pz() << ",";
      csvfile << Pip->E() << ",";
      csvfile << Pim->Px() << ",";
      csvfile << Pim->Py() << ",";
      csvfile << Pim->Pz() << ",";
      csvfile << Pim->E() << ",";
      csvfile << P->Px() << ",";
      csvfile << P->Py() << ",";
      csvfile << P->Pz() << ",";
      csvfile << P->E() << "\n";

      myfile << "\t4 0.93827231 1 0 1 11 " << energy << " 2212 0 " << weight << "\n";
      myfile << "1 0 1 11 0 0 " << Eprime->Px() << " " << Eprime->Py() << " " << Eprime->Pz() << " " << Eprime->E()
             << " " << Eprime->M() << " 0 0 0"
             << "\n";
      myfile << "2 0 1 211 0 0 " << Pip->Px() << " " << Pip->Py() << " " << Pip->Pz() << " " << Pip->E() << " "
             << Pip->M() << " 0 0 0"
             << "\n";
      myfile << "3 0 1 -211 0 0 " << Pim->Px() << " " << Pim->Py() << " " << Pim->Pz() << " " << Pim->E() << " "
             << Pim->M() << " 0 0 0"
             << "\n";
      myfile << "4 0 1 2212 0 0 " << P->Px() << " " << P->Py() << " " << P->Pz() << " " << P->E() << " " << P->M()
             << " 0 0 0"
             << "\n";
    }

    // Check to see if the total number of events is 100x greater than generated
    // If we have 100x gen number quit because something is probably wrong
    if (total++ > 100 * genNum) {
      std::cerr << "[" << __FUNCTION__ << "] Ended with break\n";
      break;
    }
  }
  // Close the file
  myfile << "\n";
  myfile.close();

  // Print report at the end
  std::cout << acc << " accepted, " << total << " generated " << 100 * acc / total << "\% accepted " << std::endl;

#ifdef PLOTS
  auto f = new TFile(Form("%s.root", file_name.c_str()), "RECREATE");
  f->cd();
  W_hist->Write();
  WvsQ2->SetOption("COLZ");
  WvsQ2->Write();
  f->Write();
#endif

  // return 0, the way to make sure unix knows that the program ended SUCCESFULLY, 1 means FAILURE
  return 0;
}
