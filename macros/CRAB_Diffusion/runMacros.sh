#!/bin/bash
### example running
## evt=50 source runMacros.sh
#set -x (For Debugging)
# functions
counter=0
increment(){
  counter=$((counter+1))
}

ExistOrCreate () {
  if [ ! -d $1 ]
    then
    mkdir -p $1
    cd $1
  else
      cd $1
  fi
}

export Photon=${Photon:-2000000}
Events=${evt:-1}

NexusPath=/home/argon/Projects/Ilker/NewNexus
BuildFolder=build_all
#SimPath=/media/argon/5TB2_rooks/CRAB/Sim/CRAB_Diffusion/Xenon/Xenon_6bar
SimPath=/media/argon/5TB2_rooks/CRAB/Sim/CRAB_Diffusion/Xenon/Xenon_8bar
#SimPath=/media/argon/5TB2_rooks/CRAB/Sim/CRAB_Diffusion/XeMethane_101
#SimPath=/home/argon/Projects/Ilker/NewNexus/macros/CRAB_Diffusion/XeMethane_174
#SimPath=/home/argon/Projects/Ilker/NewNexus/macros/CRAB_Diffusion/XeMethane_338
LogPath=${SimPath}/logs
source $NexusPath/OpticksRun
#Scan Directories and simulate using nohup
if [ -d "$SimPath" ]; then
  ExistOrCreate "$LogPath"
  cd $SimPath
  #find init mac files
  find "$SimPath" -type f -name "*.init.mac" | while read File; do
    export GEOM=CRAB_$counter
    export G4CXOpticks__setGeometry_saveGeometry=$HOME/.opticks/GEOM/$GEOM
    echo "GEOM enviroment variable is set to $GEOM"
    FileName=$(basename "$File")
    echo "Processing $FileName"
    cd $SimPath
    nohup "${NexusPath}/${BuildFolder}/nexus" -n "$Events" "$FileName" > "$LogPath/log_$counter.txt" 2>&1 &
    sleep 120
    increment
  done
  cd ../../
else
  echo "There is no init file in $SimPath"
fi
echo "Executation of simulation is completed ... "
#set +x (For Debuggin)
