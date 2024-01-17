// ----------------------------------------------------------------------------
// nexus | BaseDriftField.h
//
// Abstract base class for the description of drift fields.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef BASE_DRIFT_FIELD_H
#define BASE_DRIFT_FIELD_H

#include <G4VUserRegionInformation.hh>
#include <G4LorentzVector.hh>

class G4Material;


namespace nexus {

  /// This is an abstract base class for the description of electric 
  /// (or electromagnetic) drift fields. It inherits from 
  /// G4VUserRegionInformation so that it can be passed to the drift
  /// process attached to a G4Region.
  
  class BaseDriftField: public G4VUserRegionInformation
  {
  public:

    /// Constructor
    BaseDriftField();
    /// Destructor
    virtual ~BaseDriftField();

    /// Calculates the final position and time of a charge carrier
    /// drifting under the influence of the field. Returns the step length.
    virtual double Drift(G4LorentzVector&) = 0;

    /// Returns a random 4D point (space and time) along a drift line
    virtual G4LorentzVector 
      GeneratePointAlongDriftLine(const G4LorentzVector&, const G4LorentzVector&) = 0;

    virtual double LightYield() const;
    virtual double GetStepLimit() const;
    virtual double GetELGap() const;


  private:
    void Print() const;

  };

  // INLINE DEFINITIONS //////////////////////////////////////////////
  
  inline BaseDriftField::BaseDriftField() {}
  
  inline BaseDriftField::~BaseDriftField() {}

  inline double BaseDriftField::LightYield() const {return 0.;}

  inline void BaseDriftField::Print() const {}

  // This funcitons are added to so that we will have fine steps in 5mm gap this way if electrons get git a volume they are discarded
  inline double BaseDriftField::GetStepLimit() const {return 0; }
  inline double BaseDriftField::GetELGap() const {return 0;}	

} // end namespace nexus

#endif
