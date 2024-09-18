# Set the path to the Geant4 Installation
export Prog=$HOME/Programs
export G4INSTALL=$Prog/GEANT4/geant4-v11.1.1/install
export PATH=$G4INSTALL/bin:$PATH
export DYLD_LIBRARY_PATH=$G4INSTALL/lib:$DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH=$G4INSTALL/lib:$LD_LIBRARY_PATH

source $G4INSTALL/bin/geant4.sh


# Path to ROOT
export ROOTSYS=$Prog/ROOT/build_cxx17
export PATH=$ROOTSYS/bin/:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ROOTSYS/lib

# Garfield
export GARFIELD_INSTALL=$Prog/garfieldpp/install
export GARFIELD_HOME=$Prog/garfieldpp
source $GARFIELD_HOME/build/setupGarfield.sh

export CMAKE_PREFIX_PATH=$GARFIELD_INSTALL:$CMAKE_PREFIX_PATH
export HEED_DATABASE=$GARFIELD_INSTALL/share/Heed/database
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GARFIELD_INSTALL/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$GARFIELD_INSTALL/lib



# DEGRAD
export DEGRAD_HOME=$Prog/Degrad
#export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/opt/homebrew/Cellar/gcc/11.3.0_2/lib;
#export PATH=/opt/homebrew/bin:$PATH

# NEST
export NEST_INCLUDE_DIRS=$Prog/install/include/NEST
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$Prog/nest/install/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/argon/Programs/nest/install/lib;

export gcem_DIR=$Prog/nest/install/lib/cmake/gcem
#export CRABPATH=$HOME/Projects/Ilker/gxsim/CRAB

# Add the crab exe to the path
export PATH=$CRABPATH/build:$PATH;
