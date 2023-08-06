#ifndef gasNESTdet_hh
#define gasNESTdet_hh

#include "../Detectors/VDetector.hh"
using namespace std;

//NOTES: best g1 for DD 0.1193, but for tritium 0.1146; S2 noise 1.9, 7.5%; g1_gas 0.1019, 0.1012
//s2fano 3.6, 0.9; eField in gas 6.25, 6.2; e- life 650, 750 us; fid vol 80-130, 38-305 us; gasGap 4.25, 4.5 mm
//DISCLAIMER: Slight differences from official published values due to private LUX algorithms

class gasNESTdet: public VDetector {

public:

  gasNESTdet() {
   cerr << "*** Detector definition message ***" << endl;
   cerr << "You are currently using the gasNESTdet detector." << endl << endl;

    // Call the initialization of all the parameters
    Initialization();
  };
  ~gasNESTdet() override = default;

  // Do here the initialization of all the parameters that are not varying as a function of time
  void Initialization() override {

    // Primary Scintillation (S1) parameters
    g1 = 0.073;  // phd per S1 phot at dtCntr (not phe). Divide out 2-PE effect
    sPEres = 0.58;  // single phe resolution (Gaussian assumed)
    sPEthr = 0.35;  // POD threshold in phe, usually used IN PLACE of sPEeff
    sPEeff = 1.00;  // actual efficiency, can be used in lieu of POD threshold
    noiseBaseline[0] = 0.0;  // baseline noise mean in PE (Gaussian)
    noiseBaseline[1] = 0.0;  // baseline noise width in PE (Gaussian)
    noiseBaseline[2] = 0.0;  // baseline noise mean in e- (for grid wires)
    noiseBaseline[3] = 0.0;  // baseline noise width in e- (for grid wires)
    P_dphe = 0.2;  // chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

    coinWind = 100;  // S1 coincidence window in ns
    coinLevel = 2;   // how many PMTs have to fire for an S1 to count
    numPMTs = 89;    // For coincidence calculation

    OldW13eV = true;  // for matching EXO-200's W measurement
    // the "Linear noise" terms as defined in Dahl thesis and by Dan McK
    noiseLinear[0] = 3e-2;  // S1->S1 Gaussian-smeared with noiseL[0]*S1
    noiseLinear[1] = 3e-2;  // S2->S2 Gaussian-smeared with noiseL[1]*S2

    // Ionization and Secondary Scintillation (S2) parameters
    g1_gas = .0655;  // phd per S2 photon in gas, used to get SE size
    s2Fano = 3.61;   // Fano-like fudge factor for SE width
    s2_thr = 300.;  // the S2 threshold in phe or PE, *not* phd. Affects NR most
    E_gas = 12.;    // field in kV/cm between liquid/gas border and anode //// Though, never ever used here, let's hope. EC, 6-May-2022.
    eLife_us = 2200.;  // the drift electron mean lifetime in micro-seconds

    // Thermodynamic Properties
    inGas = true;
    T_Kelvin = 273.; //1910.04211
    p_bar = 1.0; //1910.04211

    // Data Analysis Parameters and Geometry
    dtCntr = 40.;  // center of detector for S1 corrections, in usec.
    dt_min = 20.;  // minimum. Top of detector fiducial volume
    dt_max = 60.;  // maximum. Bottom of detector fiducial volume

    radius = 50.;  // millimeters (fiducial rad)
    radmax = 50.;  // actual physical geo. limit

    TopDrift = 150.;  // mm not cm or us (but, this *is* where dt=0)
    // a z-axis value of 0 means the bottom of the detector (cathode OR bottom
    // PMTs)
    // In 2-phase, TopDrift=liquid/gas border. In gas detector it's GATE, not
    // anode!
    anode = 152.5;  // the level of the anode grid-wire plane in mm
    // In a gas TPC, this is not TopDrift (top of drift region), but a few mm
    // above it
    gate = 147.5;  // mm. This is where the E-field changes (higher)
    // in gas detectors, the gate is still the gate, but it's where S2 starts
    cathode = 1.00;  // mm. Defines point below which events are gamma-X

    // 2-D (X & Y) Position Reconstruction
    PosResExp = 0.015;     // exp increase in pos recon res at hi r, 1/mm
    PosResBase = 70.8364;  // baseline unc in mm, see NEST.cpp for usage
  }

  // S1 PDE custom fit for function of z
  // s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binom dist) e.g.
  double FitS1(double /*xPos_mm*/, double /*yPos_mm*/, double /*zPos_mm*/,
               LCE /*map*/) override {
    return 1.;  // unitless, 1.000 at detector center
  }

  // Drift electric field as function of Z in mm
  // For example, use a high-order poly spline
  double FitEF ( double /*xPos_mm*/, double /*yPos_mm*/, double /*zPos_mm*/ ) override { // in V/cm
      return 440.0;
  }

  std::vector<double> FitDirEF(double xPos_mm, double yPos_mm, double /*zPos_mm*/) override {
    //    std::vector<double> field_dir = {xPos_mm, yPos_mm, 0};
    std::vector<double> field_dir = {0, 0, 1.0};
    return field_dir;
  }

    // S2 PDE custom fit for function of r
  // s2polA + s2polB*r[mm] + s2polC*r^2+... (QE included, for binom dist) e.g.
  double FitS2(double /*xPos_mm*/, double /*yPos_mm*/, LCE /*map*/) override {
    return 1.;  // unitless, 1.000 at detector center
  }

  vector<double> FitTBA(double /*xPos_mm*/, double /*yPos_mm*/,
                        double /*zPos_mm*/) override {
    vector<double> BotTotRat(2);

    BotTotRat[0] = 0.6;  // S1 bottom-to-total ratio
    BotTotRat[1] = 0.4;  // S2 bottom-to-total ratio, typically only used for
                         // position recon (1-this)

    return BotTotRat;
  }
};
  
#endif
