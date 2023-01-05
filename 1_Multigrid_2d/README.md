
## General idea

### Multigrid solvers 
![Schematic]()

The general idea of multigrid is 

### Non-telescoping approach
Since the matrix size decreases at lower levels of the MG procedure, this can result in inefficient utilization of GPUs. One possible solultion to this is 
the non-telescoping approach, which invovles using multiple blocking schemes at lower levels and combining the information from these copies to get an improved estimate of the solution. 
Here, we explore this at the lowest 2 levels of Multigrid for the Gauge Laplace and Wilson (in progress) operators.


## Instructions to run this code: 
- This code needs the package [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- The code to run and its dependencies are in the directory `code`

## Compilation
`g++ -lm -std=c++11 mgrid_ntl.cpp -o a1 `

## Running code: 
`./a1 32 3 2 1 -0.02 3 1 2`
L=32, iters per layer = 3, block size =2, Generate near-null vectors=1, coupling m=-0.02, levels =3, non-telescoping(ntl) flag=1, ntl copies=2

## Moving results data: 
You need to manually move the data using 
`mv results*.txt ../results_data`

## Analysis: 
Run the notebooks in the folder **analysis_code**
