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
cd output
mkdir solve
cd ..
echo

echo Resolution de l equation de transport
./build/solve -gr input/solve/test output/solve/  # resolve of equation of Saint Venant with Godunov and Rusanov
echo
