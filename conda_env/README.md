## General info
These files contain information on build a custom conda environment


## Instructions for building the conda environment

1. cd <conda_env_location>
2. conda env create -prefix ./env_py3 --file environment.yml
3. conda activate env_py3
4. Run each command in the file `additions.txt`
5. Create a jupyter kernel using `python -m ipykernel install --user --name env_py3 --display-name env_py3`
6. `conda deactivate`



