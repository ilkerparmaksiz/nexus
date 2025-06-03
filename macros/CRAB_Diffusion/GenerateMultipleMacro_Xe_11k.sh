## Usage
## sh GenerateMultipleMacro.sh
#General variables
export GainMean=1
export GainSdev=1
export Concentration=1

export Pressure=8.09
export ELField=10857.14
Path="/media/argon/5TB2_rooks/CRAB/Sim/CRAB_Diffusion"
export SimPath="${Path}/Xenon/Xenon_8bar_11k"
export OutPutPath="$SimPath/H5files/"
export GasFile="/home/argon/Projects/Ilker/CRAB_Diffusion/CRAB/GasFiles/Xe100_8.09bar.gas"
export ComsolPath="/home/argon/Projects/Ilker/CRAB_COMSOL/With_Needles/15.6k_8k_7k/8bar/"
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

