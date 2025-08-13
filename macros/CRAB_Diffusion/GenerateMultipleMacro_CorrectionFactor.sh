#!/bin/bash
## Generate Macros files for Multiple Needles with specific
## Usage
## sh GenerateMultipleMacro.sh
#General variables
export GainMean=1
export GainSdev=1
export Concentration=1
export Pressure=5
export S1Yield=25510
export StepDistance=0.003

Path="/media/argon/5TB2_rooks/CRAB/Sim/CRAB_Diffusion/CorrectionFactor"

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

export SimPath="${Path}/XeMethane_101"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe99.899_CH4_0.101_5bar.gas"
#export ComsolPath="/home/argon/Projects/Ilker/CRAB_COMSOL/With_Needles/13.5k_8k_7k/5bar/"
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
# Create events with comsol
RunNeedles Comsol true
# Create events with uniformfield
RunNeedles Uniform false




#export ELField=
export MacroFile="Comsol"
export useCOMSOL=true
export SimPath="${Path}/XeMethane_174"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe99.826_CH4_0.174_5bar.gas"


# Create events with comsol
RunNeedles Comsol true
# Create events with uniformfield
RunNeedles Uniform false


#export ELField=
export MacroFile="Comsol"
export useCOMSOL=true
export SimPath="${Path}/XeMethane_338"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe99.662_CH4_0.338_5bar.gas"
export StepDistance=0.003

# Create events with comsol
RunNeedles Comsol true
# Create events with uniformfield
RunNeedles Uniform false
