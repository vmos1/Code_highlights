


## How to run this code: 

## Compilation
g++ -lm -std=c++11 mgrid_ntl.cpp -o a1 

## Running code: 
### L=32, iters per layer = 3, block size =2, Generate near-null vectors=1, coupling m=-0.02, levels =3, non-telescoping(ntl) flag=1, ntl copies=2
./a1 32 3 2 1 -0.02 3 1 2
