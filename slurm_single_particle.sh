#!/bin/bash                                                                                                                                                                                                                            

#SBATCH -J CRAB0      # A single job name for the array
#SBATCH -p run         # Partition
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=10GB         # Memory request (Mb)
#SBATCH -t 5-10:08           # Maximum execution time (D-HH:MM)
#SBATCH -o /dev/null        # Standard output
#SBATCH -e /dev/null        # Standard error

##SBATCH -o %A_%a.out        # Standard output
##SBATCH -e %A_%a.err        # Standard error

## Job options 
JOBNUMBER=${SLURM_ARRAY_TASK_ID}
SEED=`echo "scale = 2;  $JOBNUMBER" | bc`

## Path to output files and G4Executable
#outputDir="/media/argon/Ilker/CRAB/Sim/new/S1"
outputDir="/media/argon/Ilker/CRAB/Sim/Isotropic"
PathToG4Executable="/home/argon/Projects/Ilker/latest_nexus/build/nexus"

#outputDir="/media/ilker/Ilker/SimResults/Sim/S1"
#PathToG4Executable="/home/ilker/Projects/latest_CRAB0/build/nexus"


## Events
NumberOfEvents=3613132


### Options
drift=false
EL=false
cluster=false

Pressure=9.7
Run=S1

###Info
ANumber=82
Mass=210
Syield=25510   ## 1/MeV
ELYield=970    ## photons/cm*electron		
ELifeTime=1000 
ELGap=7
## Source Position
pos="0 0 -10.43"
## Naming Macro and Root Files
NAME="${Run}_${SEED}"
OUT="${NAME}_${Pressure}_bar"
InitMACRO="${NAME}_${Pressure}_bar_init.mac"
configMACRO="${NAME}_${Pressure}_bar_config.mac"
PhotonFile="Photon_${NAME}.txt"

## Paths to the Filles

RootFile="${outputDir}/output/${OUT}"

## InitMacro

## making the macro
init_MACRO="${outputDir}/macros/${InitMACRO}"

## making the macro
config_MACRO="${outputDir}/macros/${configMACRO}"

## Counts 
PathToCounts="${outputDir}/counts"
echo "Removing Older Files ..."
rm $init_MACRO
rm $config_MACRO
rm $PathToCounts/$PhotonFile
rm "$RootFile.h5"

echo $init_MACRO

###P
#Physics
echo "/PhysicsList/RegisterPhysics G4EmStandardPhysics_option4"  >>${init_MACRO}
echo "/PhysicsList/RegisterPhysics G4DecayPhysics"  >>${init_MACRO}
echo "/PhysicsList/RegisterPhysics G4RadioactiveDecayPhysics"  >>${init_MACRO}
echo "/PhysicsList/RegisterPhysics G4OpticalPhysics"  >>${init_MACRO}
echo "/PhysicsList/RegisterPhysics NexusPhysics"  >>${init_MACRO}
echo "/PhysicsList/RegisterPhysics G4StepLimiterPhysics"  >>${init_MACRO}

#Geometry
echo "/nexus/RegisterGeometry CRAB"  >>${init_MACRO}
#Generator
echo "/nexus/RegisterGenerator SingleParticleGenerator"  >>${init_MACRO}


#File Saving
echo "/nexus/RegisterPersistencyManager PersistencyManager"  >>${init_MACRO}

#Actions
echo "/nexus/RegisterTrackingAction DefaultTrackingAction"  >>${init_MACRO}
echo "/nexus/RegisterEventAction DefaultEventAction"  >>${init_MACRO}
echo "/nexus/RegisterRunAction DefaultRunAction"  >>${init_MACRO}
echo "/nexus/RegisterSteppingAction CRABAnalysisSteppingAction"  >>${init_MACRO}
echo "/nexus/RegisterMacro ${config_MACRO}"  >>${init_MACRO}


## Remove the Config File if it exists
echo $config_MACRO

## Configurations
echo  "/run/verbose 1"     >>${config_MACRO}
echo  "/event/verbose 1" >>${config_MACRO}
echo  "/tracking/verbose 0" >>${config_MACRO}
echo  "/process/em/verbose 1" >>${config_MACRO}

## Randomizing
echo "/nexus/random_seed ${SEED}" >> ${config_MACRO}


echo "/PhysicsList/Nexus/clustering ${cluster}"  >>${config_MACRO}
echo "/PhysicsList/Nexus/drift ${drift}"  >>${config_MACRO}
echo "/PhysicsList/Nexus/electroluminescence ${EL}"  >>${config_MACRO}


#GEOMETRY

echo "/Geometry/CRAB/gas_pressure ${Pressure} bar"  >>${config_MACRO}
echo "/Geometry/CRAB/chamber_diam  15. cm"  >>${config_MACRO}
echo "/Geometry/CRAB/chamber_length 43.18 cm"  >>${config_MACRO}
echo "/Geometry/CRAB/chamber_thickn 2. mm"  >>${config_MACRO}
echo "#/Geometry/CRAB/SourcePosition -4.5 -4.5 -4.5 cm"  >>${config_MACRO}
echo "/Geometry/CRAB/SourcePosition ${pos} cm"  >>${config_MACRO}
echo "#/Geometry/CRAB/SourcePosition 0 0 0 cm"  >>${config_MACRO}

echo "/Geometry/CRAB/scinYield ${Syield}  1/MeV"  >>${config_MACRO}
echo "/Geometry/CRAB/ELYield ${ELYield} 1/cm"  >>${config_MACRO}


echo "/Geometry/CRAB/ElecLifTime ${ELifeTime} ms"  >>${config_MACRO}
echo "/Geometry/CRAB/ELGap ${ELGap} mm"  >>${config_MACRO}  

###  Active
echo "/Geometry/CRAB/Active_diam 8.5 cm"  >>${config_MACRO}
echo "/Geometry/CRAB/Active_length 42 cm"  >>${config_MACRO}

### SourceEncloser
echo "/Geometry/CRAB/SourceEn_diam 10 mm"  >>${config_MACRO}
echo "/Geometry/CRAB/SourceEn_holedi 5 mm"  >>${config_MACRO}
echo "/Geometry/CRAB/SourceEn_offset 5.7 cm"  >>${config_MACRO}
echo "/Geometry/CRAB/PMT1_Pos 2.00 cm"  >>${config_MACRO}
echo "/Geometry/CRAB/PMT3_Pos 3.52 cm"  >>${config_MACRO}






#ACTIONS
echo "/Actions/CRABAnalysisSteppingAction/FileSave true"  >>${config_MACRO}
echo "/Actions/CRABAnalysisSteppingAction/FileName ${PhotonFile}"  >>${config_MACRO}
echo "/Actions/CRABAnalysisSteppingAction/FilePath ${PathToCounts}"  >>${config_MACRO}



# GENERATOR
### Isotropic Photon  Emmission



echo "/Generator/SingleParticle/region FIELDCAGE "  >>${config_MACRO}
echo "/Generator/SingleParticle/particle opticalphoton "  >>${config_MACRO}
echo "/Generator/SingleParticle/min_energy 7.295 eV "  >>${config_MACRO}
echo "/Generator/SingleParticle/max_energy 7.295 eV "  >>${config_MACRO}
echo "/Generator/SingleParticle/min_costheta -1 "  >>${config_MACRO}
echo "/Generator/SingleParticle/max_costheta 1 "  >>${config_MACRO}
echo "/Generator/SingleParticle/min_phi -3.14 "  >>${config_MACRO}
echo "/Generator/SingleParticle/max_phi 3.14 "  >>${config_MACRO}





#Actions
echo "#/Actions/DefaultEventAction/energy_threshold 100 eV"  >>${config_MACRO}
echo "#/Actions/DefaultEventAction/max_energy 1.00 GeV"  >>${config_MACRO}
echo "#/Actions/DefaultEventAction/min_energy 100 eV"  >>${config_MACRO}

# PERSISTENCY
echo "/nexus/persistency/eventType background"  >>${config_MACRO}
echo "/nexus/persistency/outputFile ${RootFile}"  >>${config_MACRO}
echo "#/run/printProgress 100"  >>${config_MACRO}



${PathToG4Executable} -n ${NumberOfEvents} ${init_MACRO}

#mv "${OUTFILE}" "${OUTPUTDIR}/"
