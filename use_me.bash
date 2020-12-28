# A litle script to execute the code for user

# Compilation
echo Creation du dossier build et déplacement dans build
mkdir build
cd build
echo

echo Chargement de cmake
cmake ..
echo

echo Compilation de solve et compute_error
make
cd ..
echo

# output is the directory wich contains the result
echo Creation de output
mkdir output
echo

echo Resolution de l equation de transport
./build/solve -grme input/solve_transport_1d/solve_transport_1d_1 output/  # resolve of equation of transport with Godunov, Rusanov and MUSCL method
echo

echo Resolution de l equation de burgers
./build/solve -grme input/solve_burgers_1d/solve_burgers_1d_1 output/    # resolve of equation of burgers with Godunov, Rusanov and MUSCL method
echo

echo Calcul des erreurs et du temps pour l equation de transport
./build/compute_error -grm input/error_transport_1d/error_transport_1d_1 output/    # compute of error and time for equation of transport with Godunov, Rusanov and MUSCL method
echo

echo Calcul des erreurs et du temps pour l equation de burgers
./build/compute_error -grm input/error_burgers_1d/error_burgers_1d_1 output/      # compute of error and time for equation if burgers with Godunov, Rusanov and MUSCL method
echo

