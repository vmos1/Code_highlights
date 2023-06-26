
## General idea
The main goal is to solve a linear system of the form A x = b
Conventionally, one uses efficient Krylov solvers such as Conjugate Gradient to obtain x given *A and *b.

### Multigrid solvers 
#### Schematic
![Schematic](https://github.com/vmos1/Code_highlights/blob/main/1_Multigrid_2d/images/MGrid_schematic.png)

As shown in the figure above, a typical Multigrid solver uses solutions of coarsened forms of the matrix to construct improved estimates of the operator at finer levels. The power of the Multigrid approach is shown in the figure blow. The iterations taken to arrive at the correct solution decreases dramatically as one goes down MG levels, to coarser matrices.
![MG_performance](https://github.com/vmos1/Code_highlights/blob/main/1_Multigrid_2d/images/MG_iterations_vs_levels.png)

### Non-telescoping approach
#### Schematic
![Non-telescoping schematic](https://github.com/vmos1/Code_highlights/blob/main/1_Multigrid_2d/images/ntl_schematic.png)
Since the matrix size decreases at lower levels of the MG procedure, this decreases the number of available parallel process, resulting in idling of GPU threads and hence inefficient utilization of GPUs. One possible solultion to this is the non-telescoping approach, which invovles using multiple blocking schemes at lower levels and combining the information from these copies using a minimal residual procedure to get an improved estimate of the solution. Here, we explore this at the lowest 2 levels of Multigrid for the Gauge Laplace and Wilson (in progress) operators.

![Results for Laplace operator](https://github.com/vmos1/Code_highlights/blob/main/1_Multigrid_2d/images/ntl_laplace.png)
The plot above shows the iterations for convergence as a function of the number of non-telescoping copies. There is a significant performance gain due to the minimal residual as seen by the speed-up obtained even for 1 copy.
We are currently working on developing such a scheme for the gauge Wilson operator.


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
