# ClothSimulation
To build the project, 
```
mkdir build
cd build
cmake ..
cmake .. -DCMAKE_BUILD_TYPE=Release 
make 
```

To run the simulation, 
```
mode with cloth hanging:
./final_project 

mode with cloth's corners being helded:
./final_project 1 
```

Here are the keys to toggle between different settings of coefficients:
```
   D ... Toggle Damp coefficient 
   B ... Toggle Bend coefficient 
   S ... Toggle Stretch coefficient
   R ... Toggle Shear coefficient
```
