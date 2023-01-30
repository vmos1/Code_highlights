

## Input data: 
Correlators for analysis are stored in `data/stored_data`

## Analyzing data

- Fit a single 4pt function using [code/1_fit_4pt_function.ipynb](https://github.com/vmos1/Code_highlights/tree/main/2_Correlated_Fits_QFE/fit_4pt_function/code/1_fit_4pt_function.ipynb) 

- To perform a combined fit of 4pt functions with different 'l's use [code/3_combined_fit_4pt_fcn.ipynb](https://github.com/vmos1/Code_highlights/tree/main/2_Correlated_Fits_QFE/fit_4pt_function/code/3_combined_fit_4pt_fcn.ipynb)

## Perform model averaging 
1. `conda activate \<conda env name\>`
2. `python code/2_model_avg_single_corr/main.py`
3. Results are stored in `data/stored_results/`

## Analyze model averaged results

Run the notebook [Combine_model-avg_and_troubleshoot.ipynb](https://github.com/vmos1/Code_highlights/blob/main/2_Correlated_Fits_QFE/fit_4pt_function/code/2_model_avg_single_corr/Combine_model-avg_and_troubleshoot.ipynb)
