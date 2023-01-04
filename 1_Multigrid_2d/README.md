


## Instructions to run this code: 
- This code needs the package [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- Enter the directory `code`

## Compilation
`g++ -lm -std=c++11 mgrid_ntl.cpp -o a1 `

## Running code: 
 L=32, iters per layer = 3, block size =2, Generate near-null vectors=1, coupling m=-0.02, levels =3, non-telescoping(ntl) flag=1, ntl copies=2
`./a1 32 3 2 1 -0.02 3 1 2`

## Moving results data: 
You need to manually move the data using 
`mv results*.txt ../results_data`

## Analysis: 
Run the notebooks in the folder **analysis_code**
