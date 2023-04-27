# Fit the antipodal 4 pt function to the expected form of a tower of exponentials
The goal is to fit the antipodal four-point function to the expected form to extract coefficient and exponents to compute physical quantities.
Such exponential fits are quite challenging and unstable. 
We deal with this using the methods of Bayesian model averaging (outlined in [arxiv:2008.01069](https://arxiv.org/abs/2008.01069)). We first perform a variety of fits, allowing variations in the dataset and the model complexity. We then perform a weighted average of the fit parameters obtained from these different fits using the Akaike Information Criterion (AIC) which penalizes fits with less input data and fewer free parameters. 
Using this procedure we obtain more robust estimates of the physical quantities.


For the test data given here, which is for the free theory, the form is a series of integer exponents. 
This folder is self-contained with all the required data. All notebooks can be run directly after cloning the repository.

## Interactive fitting
### Fitting a single correlation function: 

$$ C_1 = a_{11} \left[ e^{-E_1 \ c \ t} + e^{-E_1 \ c \ (L_t-t)} \right] + a_{12} \left[ e^{-E_2 \ c \ t} + e^{-E_2 \ c \ (L_t-t)} \right] + + a_{13} \left[ e^{-E_3 \ c \ t} + e^{-E_3 \ c \ (L_t-t)} \right]  + + a_{14} \left[ e^{-E_4 \ c \ t} + e^{-E_4 \ c \ (L_t-t)} \right]  + \ldots $$

![](https://github.com/vmos1/Code_highlights/blob/main/2_Correlated_Fits_QFE/images/fit_img1.png)

We also implement the model averaging procedure described in [arxiv:2008.01069](https://arxiv.org/abs/2008.01069)

### Combined fit of 4 different correlators

$$ C_0 = a_{11} \left[ e^{-E_1 \ c \ t} + e^{-E_1 \ c \ (L_t-t)} \right] + a_{12} \left[ e^{-E_3 \ c \ t} + e^{-E_3 \ c \ (L_t-t)} \right] + a_{13} \left[ e^{-E_5 \ c \ t} + e^{-E_5 \ c \ (L_t-t)} \right]  + a_{14} \left[ e^{-E_7 \ c \ t} + e^{-E_7 \ c \ (L_t-t)} \right]  + \ldots $$

$$ C_2 = a_{21} \left[ e^{-E_2 \ c \ t} + e^{-E_2 \ c \ (L_t-t)} \right] + a_{22} \left[ e^{-E_3 \ c \ t} + e^{-E_3 \ c \ (L_t-t)} \right] + a_{23} \left[ e^{-E_5 \ c \ t} + e^{-E_5 \ c \ (L_t-t)} \right]  + a_{24} \left[ e^{-E_7 \ c \ t} + e^{-E_7 \ c \ (L_t-t)} \right]  + \ldots $$

$$ C_4 = a_{31} \left[ e^{-E_3 \ c \ t} + e^{-E_3 \ c \ (L_t-t)} \right] + a_{32} \left[ e^{-E_5 \ c \ t} + e^{-E_5 \ c \ (L_t-t)} \right] + a_{33} \left[ e^{-E_7 \ c \ t} + e^{-E_7 \ c \ (L_t-t)} \right]  + a_{34} \left[ e^{-E_9 \ c \ t} + e^{-E_9 \ c \ (L_t-t)} \right]  + \ldots $$

$$ C_6 = a_{41} \left[ e^{-E_1 \ c \ t} + e^{-E_1 \ c \ (L_t-t)} \right] + a_{42} \left[ e^{-E_3 \ c \ t} + e^{-E_3 \ c \ (L_t-t)} \right] + a_{43} \left[ e^{-E_6 \ c \ t} + e^{-E_6 \ c \ (L_t-t)} \right]  + a_{44} \left[ e^{-E_7 \ c \ t} + e^{-E_7 \ c \ (L_t-t)} \right]  + \ldots $$

![](https://github.com/vmos1/Code_highlights/blob/main/2_Correlated_Fits_QFE/images/fit_img2.png)

## Description

The Table below describes the location of codes and data.

| Name | Description |
| --- | ---|
| [fit_4pt_function/code/1_fit_4pt_function.ipynb](https://github.com/vmos1/Code_highlights/tree/main/2_Correlated_Fits_QFE/fit_4pt_function/code/1_fit_4pt_function.ipynb) | Notebook to fit a single 4pt function |
|[fit_4pt_function/code/3_combined_fit_4pt_fcn.ipynb](https://github.com/vmos1/Code_highlights/tree/main/2_Correlated_Fits_QFE/fit_4pt_function/code/3_combined_fit_4pt_fcn.ipynb) | Notebook to perform a combined fit of 4pt functions with different 'l's |
| [fit_4pt_function/data/free_theory/](https://github.com/vmos1/Code_highlights/tree/main/2_Correlated_Fits_QFE/fit_4pt_function/data/free_theory) |Location of the free theory data |
| [fit_4pt_function/code/2_model_avg_single_corr/main.py](https://github.com/vmos1/Code_highlights/tree/main/2_Correlated_Fits_QFE/fit_4pt_function/code/2_model_avg_single_corr/main.py)|Code to perform model averaging |
