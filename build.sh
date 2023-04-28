#source setup.sh
[ ! -d "build" ] && mkdir build
cd build && cmake .. && make -j8 && cd ..
./build/nexus -i macros/CRAB/Camera/CRAB_Diffusion_Setup_init.mac 
