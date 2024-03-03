// ----------------------------------------------------------------------------
// nexus | PhysicsUtils.cc
//
// Commonly used functions to generate random quantities.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------
#include "PhysicsUtils.h"

#include <G4Material.hh>
#include <G4SystemOfUnits.hh>

namespace nexus {

 void BuildThePhysicsTable(G4PhysicsTable* &theFastIntegralTable) {
  if (theFastIntegralTable) return;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  // create new physics table

  if(!theFastIntegralTable)
    theFastIntegralTable = new G4PhysicsTable(numOfMaterials);

  for (G4int i=0 ; i<numOfMaterials; i++) {

  	G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
  	  new G4PhysicsOrderedFreeVector();

  	// Retrieve vector of scintillation wavelength intensity for
  	// the material from the material's optical properties table.

  	G4Material* material = (*theMaterialTable)[i];

    G4MaterialPropertiesTable* mpt = material->GetMaterialPropertiesTable();

    if (mpt) {

  	  G4MaterialPropertyVector* theFastLightVector =
  	    mpt->GetProperty("ELSPECTRUM");

  	  if (theFastLightVector) {
        ComputeCumulativeDistribution(*theFastLightVector, *aPhysicsOrderedFreeVector);
		  }
  	}

  	// The scintillation integral(s) for a given material
  	// will be inserted in the table(s) according to the
  	// position of the material in the material table.

  	theFastIntegralTable->insertAt(i,aPhysicsOrderedFreeVector);
  }
}

void ComputeCumulativeDistribution(const G4PhysicsOrderedFreeVector& pdf, G4PhysicsOrderedFreeVector& cdf) {
  G4double sum = 0.;
  cdf.InsertValues(pdf.Energy(0), sum);

  for (unsigned int i=1; i<pdf.GetVectorLength(); ++i) {
    G4double area =
      0.5 * (pdf.Energy(i) - pdf.Energy(i-1)) * (pdf[i] + pdf[i-1]);
    sum = sum + area;
    cdf.InsertValues(pdf.Energy(i), sum);
  }
}

void GetPhotonPol(G4ThreeVector &momentum, G4ThreeVector &polarization){
  // Generate a random direction for the photon
  // (EL is supposed isotropic)
  G4double cos_theta = 1. - 2.*G4UniformRand();
  G4double sin_theta = sqrt((1.-cos_theta)*(1.+cos_theta));

  G4double phi = twopi * G4UniformRand();
  G4double sin_phi = sin(phi);
  G4double cos_phi = cos(phi);

  G4double px = sin_theta * cos_phi;
  G4double py = sin_theta * sin_phi;
  G4double pz = cos_theta;

  momentum = G4ThreeVector(px, py, pz);

  // Determine photon polarization accordingly
  G4double sx = cos_theta * cos_phi;
  G4double sy = cos_theta * sin_phi;
  G4double sz = -sin_theta;

  polarization = G4ThreeVector(sx, sy, sz);

  G4ThreeVector perp = momentum.cross(polarization);

  phi = twopi * G4UniformRand();
  sin_phi = sin(phi);
  cos_phi = cos(phi);

  polarization = cos_phi * polarization + sin_phi * perp;
  polarization = polarization.unit();

  return;
}

}
