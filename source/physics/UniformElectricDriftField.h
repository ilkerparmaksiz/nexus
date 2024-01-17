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
    UniformElectricDriftField(double anode_position=0.,
                              double cathode_position=0.,
                              EAxis axis=kZAxis);

    /// Destructor
    ~UniformElectricDriftField();

    /// Calculate final state (position, time, drift time and length, etc.)
    /// of an ionization electron
    double Drift(G4LorentzVector& xyzt);

    G4LorentzVector GeneratePointAlongDriftLine(const G4LorentzVector&, const G4LorentzVector&);

    // Setters/getters

    void SetAnodePosition(double);
    double GetAnodePosition() const;

    void SetCathodePosition(double);
    double GetCathodePosition() const;

    void SetDriftVelocity(double);
    double GetDriftVelocity() const;

    void SetLongitudinalDiffusion(double);
    double GetLongitudinalDiffusion() const;

    void SetTransverseDiffusion(double);
    double GetTransverseDiffusion() const;

    void SetAttachment(double);
    double GetAttachment() const;

    void SetLightYield(double);
    virtual double LightYield() const;

    void SetNumberOfPhotons(double);
    double GetNumberOfPhotons() const;
    virtual double GetStepLimit() const;
    void SetStepLimit(double,double ) ;
    void SetELGap(double);
    virtual double GetELGap() const;	
  private:
    /// Returns true if coordinate is between anode and cathode
    G4bool CheckCoordinate(double);



  private:

    EAxis axis_; ///< Axis parallel to field lines

    double anode_pos_;   ///< Anode position in axis axis_
    double cathode_pos_; ///< Cathode position in axis axis_

    double drift_velocity_; ///< Drift velocity of the charge carrier
    double transv_diff_;    ///< Transverse diffusion
    double longit_diff_;    ///< Longitudinal diffusion
    double attachment_;
    double light_yield_;
    double num_ph_;
    double el_gap;

    SegmentPointSampler* rnd_;
    double steplimit_;
    double steplimitCount_;
    double tempAnodePos_;
    double steps_;

  };


  // inline methods ..................................................

  inline void UniformElectricDriftField::SetAnodePosition(double p)
  { anode_pos_ = p; }

  inline double UniformElectricDriftField::GetAnodePosition() const
  { return anode_pos_; }

  inline void UniformElectricDriftField::SetCathodePosition(double p)
  { cathode_pos_ = p; }

  inline double UniformElectricDriftField::GetCathodePosition() const
  { return cathode_pos_; }

  inline void UniformElectricDriftField::SetDriftVelocity(double dv)
  { drift_velocity_ = dv; }

  inline double UniformElectricDriftField::GetDriftVelocity() const
  { return drift_velocity_; }

  inline void UniformElectricDriftField::SetLongitudinalDiffusion(double ld)
  { longit_diff_ = ld; }

  inline double UniformElectricDriftField::GetLongitudinalDiffusion() const
  { return longit_diff_; }

  inline void UniformElectricDriftField::SetTransverseDiffusion(double td)
  { transv_diff_ = td; }

  inline double UniformElectricDriftField::GetTransverseDiffusion() const
  { return transv_diff_; }

  inline void UniformElectricDriftField::SetAttachment(double a)
  { attachment_ = a; }

  inline double UniformElectricDriftField::GetAttachment() const
  { return attachment_; }

  inline void UniformElectricDriftField::SetLightYield(double ly)
  { light_yield_ = ly; }

  inline double UniformElectricDriftField::LightYield() const
  { return light_yield_; }

  inline void UniformElectricDriftField::SetNumberOfPhotons(double nph)
  { num_ph_ = nph; }

  // This funcitons are added to so that we will have fine steps in 5mm gap this way if electrons get git a volume they are discarded
  inline void UniformElectricDriftField::SetStepLimit(double stlm,double steps)
  {
      steplimit_ = stlm;
      steps_ = steps;
  }

  inline double UniformElectricDriftField::GetStepLimit() const { return steplimit_; }
 inline double UniformElectricDriftField::GetELGap() const {return el_gap;}
 inline void  UniformElectricDriftField::SetELGap(double t) {el_gap=t;}


} // end namespace nexus

#endif
