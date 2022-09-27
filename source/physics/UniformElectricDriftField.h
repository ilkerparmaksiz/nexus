// ----------------------------------------------------------------------------
// nexus | UniformElectricDriftField.h
//
// This class defines a homogeneuos electric drift field with constant
// drift lines from cathode to anode and parallel to a cartesian axis
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef UNIFORM_ELECTRIC_DRIFT_FIELD_H
#define UNIFORM_ELECTRIC_DRIFT_FIELD_H

#include "BaseDriftField.h"
#include <geomdefs.hh>

class G4Material;


namespace nexus {

  class SegmentPointSampler;

  class UniformElectricDriftField: public BaseDriftField
  {
  public:
    /// Constructor providing position of anode and cathode,
    /// and axis parallel to the drift lines
    UniformElectricDriftField(G4double anode_position=0.,
                              G4double cathode_position=0.,
                              EAxis axis=kZAxis);

    /// Destructor
    ~UniformElectricDriftField();

    /// Calculate final state (position, time, drift time and length, etc.)
    /// of an ionization electron
    G4double Drift(G4LorentzVector& xyzt);

    G4LorentzVector GeneratePointAlongDriftLine(const G4LorentzVector&, const G4LorentzVector&);

    // Setters/getters

    void SetAnodePosition(G4double);
    G4double GetAnodePosition() const;

    void SetCathodePosition(G4double);
    G4double GetCathodePosition() const;

    void SetDriftVelocity(G4double);
    G4double GetDriftVelocity() const;

    void SetLongitudinalDiffusion(G4double);
    G4double GetLongitudinalDiffusion() const;

    void SetTransverseDiffusion(G4double);
    G4double GetTransverseDiffusion() const;

    void SetAttachment(G4double);
    G4double GetAttachment() const;

    void SetLightYield(G4double);
    virtual G4double LightYield() const;

    void SetNumberOfPhotons(G4double);
    G4double GetNumberOfPhotons() const;
    virtual G4double GetStepLimit() const;
    void SetStepLimit(G4double,G4double ) ;


  private:
    /// Returns true if coordinate is between anode and cathode
    G4bool CheckCoordinate(G4double);



  private:

    EAxis axis_; ///< Axis parallel to field lines

    G4double anode_pos_;   ///< Anode position in axis axis_
    G4double cathode_pos_; ///< Cathode position in axis axis_

    G4double drift_velocity_; ///< Drift velocity of the charge carrier
    G4double transv_diff_;    ///< Transverse diffusion
    G4double longit_diff_;    ///< Longitudinal diffusion
    G4double attachment_;
    G4double light_yield_;
    G4double num_ph_;

    SegmentPointSampler* rnd_;
    G4double steplimit_;
    G4double steplimitCount_;
    G4double tempAnodePos_;
    G4double steps_;

  };


  // inline methods ..................................................

  inline void UniformElectricDriftField::SetAnodePosition(G4double p)
  { anode_pos_ = p; }

  inline G4double UniformElectricDriftField::GetAnodePosition() const
  { return anode_pos_; }

  inline void UniformElectricDriftField::SetCathodePosition(G4double p)
  { cathode_pos_ = p; }

  inline G4double UniformElectricDriftField::GetCathodePosition() const
  { return cathode_pos_; }

  inline void UniformElectricDriftField::SetDriftVelocity(G4double dv)
  { drift_velocity_ = dv; }

  inline G4double UniformElectricDriftField::GetDriftVelocity() const
  { return drift_velocity_; }

  inline void UniformElectricDriftField::SetLongitudinalDiffusion(G4double ld)
  { longit_diff_ = ld; }

  inline G4double UniformElectricDriftField::GetLongitudinalDiffusion() const
  { return longit_diff_; }

  inline void UniformElectricDriftField::SetTransverseDiffusion(G4double td)
  { transv_diff_ = td; }

  inline G4double UniformElectricDriftField::GetTransverseDiffusion() const
  { return transv_diff_; }

  inline void UniformElectricDriftField::SetAttachment(G4double a)
  { attachment_ = a; }

  inline G4double UniformElectricDriftField::GetAttachment() const
  { return attachment_; }

  inline void UniformElectricDriftField::SetLightYield(G4double ly)
  { light_yield_ = ly; }

  inline G4double UniformElectricDriftField::LightYield() const
  { return light_yield_; }

  inline void UniformElectricDriftField::SetNumberOfPhotons(G4double nph)
  { num_ph_ = nph; }

  // This funcitons are added to so that we will have fine steps in 5mm gap this way if electrons get git a volume they are discarded
  inline void UniformElectricDriftField::SetStepLimit(G4double stlm,G4double steps)
  {
      steplimit_ = stlm;
      steps_ = steps;
  }

  inline G4double UniformElectricDriftField::GetStepLimit() const { return steplimit_; }




} // end namespace nexus

#endif
