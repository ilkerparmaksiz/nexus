// ----------------------------------------------------------------------------
// nexus | UniformElectricDriftField.cc
//
// This class defines a homogeneuos electric drift field with constant
// drift lines from cathode to anode and parallel to a cartesian axis
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "UniformElectricDriftField.h"
#include "SegmentPointSampler.h"

#include <Randomize.hh>

#include <math.h>
#include "CLHEP/Units/SystemOfUnits.h"


namespace nexus {

  using namespace CLHEP;


  UniformElectricDriftField::UniformElectricDriftField
  (G4double anode_position, G4double cathode_position, EAxis axis):
    BaseDriftField(),
    axis_(axis), anode_pos_(anode_position), cathode_pos_(cathode_position),
    drift_velocity_(0.), transv_diff_(0.), longit_diff_(0.),  light_yield_(0.)
  {
    // initialize random generator with dummy values
    rnd_ = new SegmentPointSampler(G4LorentzVector(0.,0.,0.,-999.),
                                   G4LorentzVector(0.,0.,0.,-999.));
    steplimitCount_=0;
    tempAnodePos_=0;

  }



  UniformElectricDriftField::~UniformElectricDriftField()
  {
    delete rnd_;
  }



  G4double UniformElectricDriftField::Drift(G4LorentzVector& xyzt)
  {
    // If the origin is not between anode and cathode,
    // the charge carrier, obviously, doesn't move.
    if (!CheckCoordinate(xyzt[axis_]))
      return 0.;

    // Stopping the drift


    // Set the offset according to relative anode-cathode pos
    G4double secmargin = -1. * micrometer;
    if (anode_pos_ > cathode_pos_) secmargin = -secmargin;

    // Calculate drift time and distance to anode
    //G4double drift_length = fabs(xyzt[axis_] - anode_pos_);
    G4double drift_length;
    if (steplimit_==0) {
        drift_length= fabs(xyzt[axis_] - anode_pos_);
        tempAnodePos_=anode_pos_;

    }
    else {
        // This portion is added for making sure that electrons will stop and killed if they are in  a physical volume
        if( xyzt.t()==0) steplimitCount_=0;

        if(steplimitCount_<steplimit_) {
            steplimitCount_ += steps_;
            tempAnodePos_ = xyzt[axis_] - steps_;
        }else{
            tempAnodePos_=anode_pos_;
        }
        drift_length=fabs(xyzt[axis_] - tempAnodePos_);
    }
    G4double drift_time = drift_length / drift_velocity_;
    //G4cout<<"Drift_Length --> "<<drift_length<<G4endl;
    // Calculate longitudinal and transversal deviation due to diffusion
    G4double transv_sigma = transv_diff_ * sqrt(drift_length);
    G4double longit_sigma = longit_diff_ * sqrt(drift_length);
    G4double time_sigma = longit_sigma / drift_velocity_;

    G4ThreeVector position;
    G4double time;

    for (G4int i=0; i<3; i++) {
      if (i != axis_)  {     // Transverse coordinate
        position[i] = G4RandGauss::shoot(xyzt[i], transv_sigma);
      }
      else { // Longitudinal coordinate
        //position[i] = anode_pos_ + secmargin;
        position[i] = tempAnodePos_ + secmargin;
        G4double deltat = G4RandGauss::shoot(0, time_sigma);
        time = xyzt.t() + drift_time + deltat;
        if (time < 0.) time = xyzt.t() + drift_time;
      }
    }

    // Calculate step length as euclidean distance between initial
    // and final positions
    G4ThreeVector displacement = position - xyzt.vect();
    G4double step_length = displacement.mag();

    // Set the new time and position of the drifting charge
    xyzt.set(time, position);
    return step_length;
  }



  G4LorentzVector UniformElectricDriftField::GeneratePointAlongDriftLine(const G4LorentzVector& origin, const G4LorentzVector& end)
  {

    // if (origin != rnd_->GetPrePoint()) {
      // G4LorentzVector end(origin);
      // Drift(end);

    rnd_->SetPoints(origin, end);
    //   }


    return rnd_->Shoot();
  }



  G4bool UniformElectricDriftField::CheckCoordinate(G4double coord)
  {
    G4double max_coord = std::max(anode_pos_, cathode_pos_);
    G4double min_coord = std::min(anode_pos_, cathode_pos_);
    //   G4cout << "max = " << max_coord << ", min = " << min_coord << G4endl;
    return !((coord > max_coord) || (coord < min_coord));
  }


} // end namespace nexus
