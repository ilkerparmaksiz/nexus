source setup.sh
[ ! -d "build" ] && mkdir build
cd build && cmake .. && make -j8 && cd ..
./build/nexus -i macros/CRAB/CRAB_SingleParticle_init.mac
