# Climate, Land and Energy Analytical study of Potential Nexus Issues in Jordan and Morocco
This repository contains the working code for the “Implementing the 2030 Agenda 
for water efficiency/productivity and water sustainability in NENA countries” project. 
The project is funded by the UN Food and Agricultural Organization, KTH is 
providing analyses and support on decision making across the climate, land, 
energy and water spheres in Morocco and Jordan.

## Model background: NWSAS Nexus analysis
The bases of the model implemented for the Jordan and Morocco cases, 
are inherited from previous work developed for the NWSAS basin. The version of 
the NWSAS code included in this repository has the following functions: 
  1) Estimation of the irrigation water demand based on irrigated area, crop 
  type, crop calendar, effective rainfall ..etc.
  2) Estimation of the energy demand for pumping ground water. 
  3) Estimation of the energy demand for brackish water desalination.
  4) Comparison between different electricity supply options to compute the 
  least cost technology based on LCOE. 

From the original NWSAS model a python package `Nexus_tool` was developed 
containing all required methods to run all four steps previously mentioned. 
Three `Jupyter notebooks` are provided allowing to run the NWSAS case study 
from begining to end, calculating the water demand (`NWSAS model/water_demand_runer.ipynb`), 
energy demand (`NWSAS model/energy_demand_runner.ipynb`) and least-cost energy 
supply options (`NWSAS model/least_cost_runner.ipynb`).

## The Jordan and Moroccan Souss-Massa basin cases
The Nexus analysis developed for the Jordan and the Moroccan Souss-Massa basin 
cases, softlinks a water balance model developed by SEI and the countries 
stakeholders using the Water Evaluation and Planning system (WEAP), and the 
GIS-based energy model developed by KTH. The model uses WEAP for estimating water 
supplies based on climate-driven hydrological routines that calculate rainfall 
runoff and groundwater recharge. It estimates usage patterns for the main water 
sectors, evaluates the productivity of cropping systems under different climate 
futures and assess their impact on energy and water systems. Furthermore, the 
energy component implements GIS-based methodologies to estimate energy 
requirements for groundwater and surface water pumping, new water desalination 
projects and major wastewater treatment plants. Finally, least-cost generation 
technologies to supply electricity for water and agricultural sectors are identified.

## Installation
To install the required dependencies, install the miniconda or conda package 
manager and create two conda environments running 
`conda env create -n <name-of-environment> -f envs/notebook.yml` and 
`conda env create -n <name-of-environment> -f envs/snakemake.yml` in the conda 
shell or git bash (replace `<name-of-environment>` by a custom name for your 
environment). Afterwards, activate the environment with `conda activate 
<name-of-environment>`.

The `notebook` environment is intended to run the standalone jupyter notebooks, 
provided for each step of the model. While the `snakemake` environment is intended 
for running the entire automated workflow of the project.

## Running the model
### Jupyter notebooks
To run the model, first activate the previously created conda environment for 
notebooks by running `conda acivate <name-of-environment>` and then run 
`jupyter notebook`. Alternatively, you can start the Anaconda navigator, 
select the previously created environment and start a Jupyter notebook session. 
Open any of the runner files of either model and follow the steps.

### Snakemake
Activate the previously created conda environment for snakemake by running 
`conda acivate <name-of-environment>`. Then open the snakemake file of the 
model you wish to run and set the user inputs to the ones desired:
```python
######################### User defined parameters #############################
dash_folder = "../Morocco dashboard" #set the dashboard folder path to save the results

SCENARIOS = ['Reference', 'Desalination','Irrigation Intensification'] #define the list of scenarios to run
CLIMATE = ['Trend', 'Climate Change'] #define the list of climates to run
W_RATE = [0, 0.3, 0.5, 0.7] #define the list of wind technology cost reduction rate
PV_RATE = [0, 0.3, 0.5, 0.7] #define the list of solar PV technology cost reduction rate
GRID_RATE = [0, -0.3, -0.5, -0.7] #define the list of grid price increment rate
###############################################################################
```
Run `snakemake -s Moroccan\ model/snakefile -n` for a dry run (of the Moroccan 
case). This will list all of the different jobs the workflow would run, without 
actually running anything.

Finally, run `snakemake -s Moroccan\ model/snakefile -j` to run all of the listed jobs.

## Visualazing results
Results for te case studies of Jordan and the Souss-Massa river basin, can be 
explored through interactive dashboards found in:
* Jordan: [https://jordan-nexus-model.herokuapp.com/](https://jordan-nexus-model.herokuapp.com/)
* Moroco: [https://souss-massa-nexus-model.herokuapp.com/](https://souss-massa-nexus-model.herokuapp.com/)

For exploration of results locally, you can run the visualization from your 
machine in a local host by running `python Morocco\ dashboard/app.py` 
(for the Moroccan case) and going to http://localhost:8080/ in your browser.