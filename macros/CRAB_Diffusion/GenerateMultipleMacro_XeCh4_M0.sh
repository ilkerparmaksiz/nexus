#!/bin/bash
## Generate Macros files for Multiple Needles with specific
## Usage
## sh GenerateMultipleMacro.sh
#General variables
export GainMean=0.55
export GainSdev=0.09
export Concentration=0.101

export Pressure=5
#export S1Yield=25510
export S1Yield=11989.7 # (25510*0.47)
#export ELField=
Path="/media/argon/5TB2_rooks/CRAB/Sim/CRAB_Diffusion/Diffusion_Mode0/XenonMethane"
export SimPath="${Path}/XeMethane_101"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe99.899_CH4_0.101_5bar_M0.gas"
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
# Create events with comsol
RunNeedles Comsol true
# Create events with uniformfield
RunNeedles Uniform false

### Generate Folder for CH4 0.174
export GainMean=0.37
export GainSdev=0.06
export Concentration=0.174

#export S1Yield=25510
export S1Yield=10969.3 # (25510*0.47)
#export ELField=
export MacroFile="Comsol"
export useCOMSOL=true
export SimPath="${Path}/XeMethane_174"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe99.826_CH4_0.174_5bar_M0.gas"


# Create events with comsol
RunNeedles Comsol true
# Create events with uniformfield
RunNeedles Uniform false


### Generate Folder for CH4 0.338
export GainMean=0.37
export GainSdev=0.13
export Concentration=0.338

#export S1Yield=25510
export S1Yield=3316.3
#export ELField=
export MacroFile="Comsol"
export useCOMSOL=true
export SimPath="${Path}/XeMethane_338"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe99.662_CH4_0.338_5bar_M0.gas"
export StepDistance=0.003

# Create events with comsol
RunNeedles Comsol true
# Create events with uniformfield
RunNeedles Uniform false
