#!/bin/bash
## Generate Macros files for Multiple Needles with specific
## Usage
## sh GenerateMultipleMacro.sh
#General variables
export GainMean=1
export GainSdev=1
export Concentration=1

export Pressure=5

Path="/media/argon/5TB2_rooks/CRAB/Sim/CRAB_Diffusion"
export SimPath="${Path}/Xenon/Xenon_5bar"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe100_5bar.gas"
export ComsolPath="/home/argon/Projects/Ilker/CRAB_COMSOL/With_Needles/13.5k_8k_7k/5bar/"
RunNeedles (){
  export MacroFile=$1
  export useCOMSOL=$2
  echo "Running Multi Macro generator"
  export Needle="Needle_4cm" # Specific Variables
  sh GenerateSingleMacro.sh
  export Needle="Needle_9cm"
  sh GenerateSingleMacro.sh
  export Needle="Needle_14cm"
  sh GenerateSingleMacro.sh
}
## Generating for 5 bar
# Create events with comsol
RunNeedles Comsol true
# Create events with uniformfield
RunNeedles Uniform false

### Generate Folder for Xe_6bar
export Pressure=5.89

export SimPath="${Path}/Xenon/Xenon_6bar"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe100_5_89bar.gas"
export ComsolPath="/home/argon/Projects/Ilker/CRAB_COMSOL/With_Needles/13.5k_8k_7k/5.89bar/"


# Create events with comsol
RunNeedles Comsol true
# Create events with uniformfield
RunNeedles Uniform false



### Generate Folder for Xe_8bar
export Pressure=8.09
export SimPath="${Path}/Xenon/Xenon_8bar"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe100_8.09bar.gas"
export ComsolPath="/home/argon/Projects/Ilker/CRAB_COMSOL/With_Needles/13.5k_8k_7k/8bar/"


# Create events with comsol
RunNeedles Comsol true
# Create events with uniformfield
RunNeedles Uniform false

### Generate Folder for CH4 0.6  For 6bar
export GainMean=0.13
export GainSdev=0.13
export Concentration=0.338

#export S1Yield=25510
export S1Yield=3316.3
#export ELField=
export Pressure=5.89
export SimPath="${Path}/XeMethane6bar"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe99.409_CH4_0.591_5.89bar.gas"
export ComsolPath="/home/argon/Projects/Ilker/CRAB_COMSOL/With_Needles/13.5k_8k_7k/5.89bar/"


# Create events with comsol
RunNeedles Comsol true
# Create events with uniformfield
RunNeedles Uniform false
