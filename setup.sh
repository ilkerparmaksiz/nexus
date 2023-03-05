# Set the path to the Geant4 Installation
export G4INSTALL=/home/argon/Software/Geant4/geant4-v11.0.3/install;
export PATH=$G4INSTALL/bin:$PATH;
export LD_LIBRARY_PATH=$G4INSTALL/lib:$LD_LIBRARY_PATH;

source $G4INSTALL/bin/geant4.sh;


# Path to ROOT
export ROOTSYS=~/Software/ROOT/build_cxx17/;
export PATH=$ROOTSYS/bin/:$PATH;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib;

# Garfield
#export GARFIELD_INSTALL=/home/argon/Projects/Krishan/garfieldpp/install
#export GARFIELD_HOME=/home/argon/Projects/Krishan/garfieldpp/
#export CMAKE_PREFIX_PATH=/home/argon/Projects/Krishan/garfieldpp/install:$CMAKE_PREFIX_PATH
#export HEED_DATABASE=$GARFIELD_INSTALL/share/Heed/database
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GARFIELD_INSTALL/lib

# DEGRAD
#export DEGRAD_HOME=/home/argon/Projects/Krishan/Degrad

# NEST
#export NEST_INCLUDE_DIRS=/home/argon/Projects/Krishan/NEST/install/include/NEST
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/argon/Projects/Krishan/NEST/install/lib

#export CRABPATH=/home/argon/Projects/Krishan/gxsim/CRAB/
