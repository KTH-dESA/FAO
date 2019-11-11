# FAO
This repository contains the working code for the FAO project. It inherits the 
version of the code developed for the NWSAS project which has the following functions: 
  1) Estimation of the irrigation water demand based on irrigated area, crop 
  type, crop calendar, effective rainfall ..etc.
  2) Estimation of the energy demand for pumping ground water. 
  3) Estimation of the energy demand for brackish water desalination.
  4) Comparison between different electricity supply options to computer the 
  least cost technology based on (LCOE). 

The model is contained on a python package with all required methods 
to run all four steps previously mentioned. Three `Jupyter notebooks` are 
provided which allow to run the NWSAS case study from begining to end, calculating the 
water demand (`water_demand_runer.ipynb`), energy demand (`energy_demand_runner.ipynb`)
and least-cost energy supply options (`least_cost_runner.ipynb`).

## Installation
To install the required dependencies, the easiest way is to install the 
miniconda or conda package manager and create a conda environment running 
`conda env create -n <name-of-environment> -f environment.yml` in the conda 
shell or git bash (replace \<name-of-environment\> by a custom name for your 
environment). Afterwards, activate the environment with `conda activate 
<name-of-environment>` and install the `pyeto` package with `pip install 
git+https://github.com/woodcrafty/PyETo.git@8ca15bd6f70eda39f42b40581e7d5cab4bf76fbc#egg=PyETo`

## Running the model
To run the model, first activate the previously created conda environment by 
running `conda acivate <name-of-environment>` and then `jupyter notebook`. 
Alternatively, you can start the Anaconda navigator, select the previously 
created environment and start a Jupyter notebook session. Open any of the 
runner files and follow the steps.