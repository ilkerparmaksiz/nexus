#!/bin/bash
## Usage
## sh GenerateSingleMacro.sh
#General variables
## Functions
ExistOrCreate () {
  if [ ! -d $1 ]
    then
    mkdir -p $1
    cd $1
  else
      cd $1
  fi
}
## Variables
## Geometry
Pressure=${Pressure:-5}  # in terms of bar
S1Yield=${S1Yield:=-25510}
## Used for uniform field
DriftField=${DriftField:-438} # V/cm
ELField=${ELField:-7857.14} # V/cm
## For Dealing with Methane
## We need to reduce yields according to data
GainMean=${GainMean:-1}
GainSdev=${GainSdev:-1}
CH4Concentration=${Concentration:-1}

## Needles
Needle=${Needle:-"Needle_9cm"}
## COMSOL
useCOMSOL=${useCOMSOL:-false}

ComsolPath=${ComsolPath:-"/home/argon/Projects/Ilker/CRAB_COMSOL/With_Needles/13.5k_8k_7k/5bar/"}

## Garfield
GasFile=${GasFile:-"/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe100_5bar.gas"}

MacroFile=${MacroFile:-"Garfield"}

SimPath=${SimPath:-"/home/argon/Projects/Ilker/NewNexus/macros/CRAB_Diffusion/test"}

## Output file
OutPutPath=${OutPutPath:-"/home/argon/Projects/Ilker/CRAB_COMSOL/With_Needles/13.5k_8k_7k/5bar/RootFiles/"}
ExistOrCreate "$OutPutPath"
OutPutFileName=${OutPutFileName:-"G4Test_${MacroFile}_${Pressure}bar_${Needle}"}

NexusPath=${NexusPath:-"/home/argon/Projects/Ilker/NewNexus"}
# Temporary solution
 ExistOrCreate "$SimPath/data"
if [ -f "$SimPath/data/ELProfile.txt" ]; then
  echo "ELProfile.txt exist"
else
  cp $NexusPath/data/ELProfile.txt $SimPath/data
  cp $NexusPath/data/Needle_* $SimPath/data
fi


## Copy Macro Files
ConfigMacro="${SimPath}/CRAB_${MacroFile}_${Needle}_${Pressure}bar.config.mac"
InitMacro="${SimPath}/CRAB_${MacroFile}_${Needle}_${Pressure}bar.init.mac"

echo "Starting to generate files"
echo "ConfigMacro --> ${ConfigMacro}"
echo "InitMacro   --> ${InitMacro}"
ExistOrCreate "$SimPath"

cp "${NexusPath}/macros/CRAB_Diffusion/CRAB_Garfield.init.mac" ${InitMacro}
cp "${NexusPath}/macros/CRAB_Diffusion/CRAB_Garfield.config.mac" ${ConfigMacro}


## Modify Init file
sed -i "s#.*RegisterMacro.*#/nexus/RegisterMacro CRAB_${MacroFile}_${Needle}_${Pressure}bar.config.mac#" ${InitMacro}

## Modify the Config file

sed -i "s#.*GasPressure.*#/Geometry/CRAB0/GasPressure ${Pressure} bar#" ${ConfigMacro}
sed -i "s#.*ScintYield.*#/Geometry/CRAB0/ScintYield ${S1Yield}#" ${ConfigMacro}
sed -i "s#.*fieldDrift.*#/Geometry/CRAB0/fieldDrift ${DriftField}#" ${ConfigMacro}
sed -i "s#.*fieldEL.*#/Geometry/CRAB0/fieldEL ${ELField}#" ${ConfigMacro}
sed -i "s#.*GainReduction.*#/Geometry/CRAB0/GainReduction ${GainMean} ${GainSdev} ${CH4Concentration} #" ${ConfigMacro}
sed -i "s#.*GasFile.*#/Geometry/CRAB0/GasFile ${GasFile}#" ${ConfigMacro}
sed -i "s#.*useCOMSOL.*#/Geometry/CRAB0/useCOMSOL ${useCOMSOL}#" ${ConfigMacro}
sed -i "s#.*ComsolPath.*#/Geometry/CRAB0/ComsolPath ${ComsolPath}#" ${ConfigMacro}
sed -i "s#.*Mode.*#/Generator/SPModified/Mode ${Needle}#" ${ConfigMacro}
sed -i "s#.*output_file.*#/nexus/persistency/output_file ${OutPutPath}${OutPutFileName}#" ${ConfigMacro}
echo "Finished Creating Files .."
