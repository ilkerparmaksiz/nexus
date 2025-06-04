#!/bin/bash
#SBATCH -J CRAB # A single job name for the array
#SBATCH -c 1 # Number of cores
#SBATCH -p node2
#SBATCH --mem 30000 # Memory request (6Gb)
#SBATCH -t 10-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o /dev/null # Standard output
#SBATCH -e /dev/null # Standard error

start=`date +%s`

# Set the configurable variables
JOBNAME="10bar_G4"

TYPE="Reflections"
N_EVENTS=33
#N_EVENTS=20
Reflections=true
BuildFolder=buildv2
NexusPath=/home/argon/Projects/Ilker/NewNexus
SimPath=/home/argon/Projects/Ilker/NewNexus/out/June_4_2025
alias mc='(cd "${NexusPath}")'
#echo "CRABPATH is $CRABPATH"
## if the folder does nt exist , this will create it
ExistOrCreate () {
  if [ ! -d $1 ]
    then
    mkdir -p $1
    cd $1
  else
      echo "Folder Exist !"
      cd $1
  fi
}

# Create the directory
#source "/home/argon/Projects/Ilker/gxsim/CRAB/macros/run.sh test"

ExistOrCreate "$SimPath"
ExistOrCreate "$SimPath/$JOBNAME/$TYPE/jobid_${SLURM_ARRAY_TASK_ID}"
#ExistOrCreate "$JOBNAME/$TYPE/jobid_${SLURM_ARRAY_TASK_ID}"
# Copy the macro file
cp ${NexusPath}/macros/CRAB.init.mac .
cp ${NexusPath}/macros/CRAB.config.mac .
filePath=$SimPath/$JOBNAME/$TYPE/jobid_"${SLURM_ARRAY_TASK_ID}"/alpha

# Setup nexus and run
echo "Setting Up Code" 2>&1 | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt
#source /home/argon/Projects/Krishan/gxsim/CRAB/setup_cluster.sh

# Calculate the unique seed number	
SEED=$((${N_EVENTS}*(${SLURM_ARRAY_TASK_ID} - 1) + ${N_EVENTS}))
echo "The seed number is: ${SEED}" 2>&1 | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt


# Replace the number of events in the file as well as the event index
sed -i "s#.*random_seed.*#/nexus/random_seed ${SEED}#" CRAB.config.mac
sed -i "s#.*SteelReflect.*#/Geometry/CRAB0/SteelReflect ${Reflections}#" CRAB.config.mac
sed -i "s#.*output_file.*#/nexus/persistency/output_file ${filePath}#" CRAB.config.mac
sed -i "s#.*RegisterMacro.*#/nexus/RegisterMacro CRAB.config.mac #" CRAB.init.mac
# NEXUS
echo "Running Nexus" 2>&1 | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt
${NexusPath}/${BuildFolder}/nexus -n ${N_EVENTS} CRAB.init.mac 2>& 1 | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt

echo; echo; echo;

echo "FINISHED....EXITING" 2>&1 | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt

end=`date +%s`
let deltatime=end-start
let hours=deltatime/3600
let minutes=(deltatime/60)%60
let seconds=deltatime%60
printf "Time spent: %d:%02d:%02d\n" $hours $minutes $seconds | tee -a log_crab"${SLURM_ARRAY_TASK_ID}".txt
