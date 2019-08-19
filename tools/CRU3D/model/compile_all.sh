#Compile tetgen library
cd mesh/tetgen1.4.3
make clean
make tetlib

#Compile mesh
cd ..
make clean
make

#Compile linescan
cd ../linescan
make clean
make

#Compile model (cru3d)
cd ..
make clean
make
