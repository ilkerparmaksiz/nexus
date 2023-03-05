if [-d "build"]
then 
  rm -r build
  echo "Directory Removed"
fi
source setup.sh
mkdir build
cd build
cmake ..
make -j8
cd ..
