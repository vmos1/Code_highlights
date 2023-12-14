# Introduction
To study strong forces such as the nuclear force, one performs Markov chain Monte-Carlo simulations on large lattices to obtain configurations that are then used to compute physical observables. These simulations produce large amounts of data. Reading from them directly and analyzing them in real-time is inconvenient. Databases offer a convenient way to store processed data that can be used for real-time analysis.

# Data
The model being studied is an $SU(4)$ gauge theory with one fermion flavor, which is of interest as a composite dark matter model. The simulations are done with Mobius-Domain-wall fermions on lattices of size $16 ^ 3 \times 8 \times 16$.
The data is obtained by running Monte-Carlo simulations using the [Grid](https://github.com/paboyle/Grid) library. The resulting dataset is very large (~100GB). The output is then parsed to obtain data for analysis.


### Locating the phase transition
The main goal is to locate the *confinement phase transition* of this theory.
Below is a figure showing the values of one observable : the polyakov loop as a function of the parameter $\beta$. The discontinuity in the region around $\beta \sim 10.8$ indicates the location of this transition.

![Polyakovloop discontinuity](https://github.com/vmos1/Code_highlights/blob/main/4_Monte_Carlo_analysis_sql_dbase/images/polyakov_loop.png)

### Variation in Monte-Carlo time
To further understand the behavior, we plot the variation of observables in Monte-Carlo time.
Time_series variation of polyakov loop | Scatter plot of Polyakov loop  |
:-------------:|:---------------:
![Time_series variation of polyakov loop](https://github.com/vmos1/Code_highlights/blob/main/4_Monte_Carlo_analysis_sql_dbase/images/time_series_polyakov.png) |![Scatter plot of Polyakov loop](https://github.com/vmos1/Code_highlights/blob/main/4_Monte_Carlo_analysis_sql_dbase/images/scatter_polyakov.png)

In the left figure, we see the variation of the Polyakov loop for two values of $\beta$ on either side of the transition. For the larger $\beta$ value, there is larger transient region before the observable attains steady-state values. To the right, we see a scatter plot of the complex Polyakov loop along the real-imaginary axis. As we cross the phase transition, the polyakov loop moves away from the origin (0,0) to a non-zero real value.

### Repository information
The Table below describes the important codes and their locations

| Name | Description |
| --- | ---|
| [1_HMC_write_to_database.ipynb](https://github.com/vmos1/Code_highlights/blob/main/4_Monte_Carlo_analysis_sql_dbase/code/1_HMC_write_to_database.ipynb) | Code to write to sql database files |
| [4_Monte_Carlo_analysis_sql_dbase/code/2_analysis.ipynb](https://github.com/vmos1/Code_highlights/blob/main/4_Monte_Carlo_analysis_sql_dbase/code/2_analysis.ipynb) | Notebook to visualize general behavior and Monte-Carlo time variation |

## Summary: 
The results show that we can identify the approximate location of the confinement transition using the variation of physical observables. In this region, we need to perform measurements to extract other physical quantities.
