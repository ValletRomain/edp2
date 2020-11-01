# A litle script to execute the code for user

# Compilation
mkdir build
cd build
cmake ..
make
cd ..

# output is the directory wich contains the result
mkdir output

./build/godunov input/solve_transport_1d_1 output/  # resolve of equation of transport
./build/godunov input/solve_burgers_1d_1 output/    # resolve of equation of burgers
./build/godunov_error input/error_transport_1d_1 output/    # compute of error and time for equation of transport
./build/godunov_error input/error_burgers_1d_1 output/      # compute of error and time for equation if burgers


