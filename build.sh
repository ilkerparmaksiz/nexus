#!/bin/bash
echo "Welcome to the build script"
echo "First 3 Arguments are $1,$2,$3"
sleep 2
echo "Cleaning tmp File /tmp/argon/opticks"
rm -rf /tmp/argon/opticks
sleep 2
if [ "$1" == "run" ]; then
	source OpticksRun
	buildFolder=build
	runNumber=30
	## Assignt the build folder 
	if [ -n "$2" ]; then
		buildFolder="$2"
	fi
	## Assign the Number of runs
	if [ -n "$3" ]; then
		runNumber=$3
	fi
	## First Compile
	echo "Compiling the code folder is ${buildFolder}"
	sleep 1
	cd "$buildFolder" && cmake .. && make -j4 && cd ..
	## Then run
	echo " Running the code total events are ${runNumber} "
	sleep 1
	./"${buildFolder}"/nexus -n "${runNumber}" macros/CRAB.init.mac

elif [ "$1" == "debug" ]; then
  source OpticksDebug
  buildFolder=build_G4OpticksTest
  runNumber=30
  ## Assignt the build folder
  if [ -n "$2" ]; then
    buildFolder="$2"
  fi
  ## Assign the Number of runs
  if [ -n "$3" ]; then
    runNumber=$3
  fi
  ## First Compile
  echo "Compiling the code"
  sleep 1
	cd "$buildFolder" && cmake .. && make -j4 && cd ..
  ## Then run
  echo " Running the code "
  sleep 1
  ./"${buildFolder}"/nexus -n "${runNumber}" macros/debug/CRAB_Debug.init.mac
elif [ "$1" == "build" ]; then
  source OpticksEnv
  buildFolder=build_G4OpticksTest
  if [ -n "$2" ]; then
   buildFolder=$2
  fi
  echo "Compiling the code"
  sleep 1
  cd "$buildFolder" && cmake .. && make -j4 && cd ..
  echo " Running the code "
  sleep 1
  ./"${buildFolder}"/nexus -i macros/CRAB.init.mac
fi
