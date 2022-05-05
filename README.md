# Climate, Land and Energy Analytical study of Potential Nexus Issues in Jordan and Morocco
This repository contains the working code for the “Implementing the 2030 Agenda 
for water efficiency/productivity and water sustainability in NENA countries” project. 
The project is funded by the UN Food and Agricultural Organization, KTH is 
providing analyses and support on decision making across the climate, land, 
energy and water spheres in Morocco and Jordan.

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
projects and major wastewater treatment plants. Finally, solar PV pumping is 
evaluated as a clean option to supply electricity for the agricultural sector.

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
notebooks by running `conda activate <name-of-environment>` and then type 
`jupyter notebook`. Alternatively, you can start the Anaconda navigator, 
select the previously created environment and start a Jupyter notebook session. 
Open any of the runner files of either model and follow the steps.

### Snakemake
Activate the previously created conda environment for snakemake by running 
`conda activate <name-of-environment>`. Then open the snakemake file of the 
model you wish to run and set the user inputs to the ones desired:
```python
######################### User defined parameters #############################
dash_folder = "dashboard" #set the dashboard folder

CLIMATE = ['Trend', 'Climate Change'] #define the list of climates to run
SCENARIOS = ['Reference', 'Desalination','Irrigation Intensification', 'Reference Wastewater Reuse', 'Desalination Wastewater Reuse']
W_RATE = [0]
BT_RATE = [0]
PV_RATE = [0]
GRID_RATE = [0]
PV_LEVELS = [10, 20, 50]
BUTANE_PHASEOUT_SCENARIOS = [None, 2040, 2030]
###############################################################################
```
Run `snakemake -s Morocco\ model/snakefile -n` for a dry run (of the Moroccan 
case). This will list all of the different jobs the workflow would run, without 
actually running anything.

Finally, run `snakemake -s Morocco\ model/snakefile -j` to run all of the listed jobs.

## Visualazing results
Results for te case studies of Jordan and the Souss-Massa river basin, can be 
explored through interactive dashboards found in:
* Jordan: [https://jordan-nexus-model.herokuapp.com/](https://jordan-nexus-model.herokuapp.com/)
* Moroco: [https://souss-massa-nexus.herokuapp.com/](https://souss-massa-nexus.herokuapp.com/)

For exploration of results locally, you can run the visualization from your 
machine in a local host by running `python Morocco\ model/dashboard/index.py` 
(for the Moroccan case) and going to http://localhost:8080/ in your browser. 
Notice that you would need to create and environment using one of the `environment.yml` 
files found inside each case dashboard folder, and activate it before running the 
dashboard locally.

## Credits

**Conceptualization:** [Youssef Almulla](https://www.kth.se/profile/almulla), [Camilo Ramirez](https://www.kth.se/profile/camilorg), [Brian Joyce](https://www.sei.org/people/brian-joyce/), [Annette Huber-Lee](https://www.sei.org/people/annette-huber-lee/) <br /> and [Francesco Fuso-Nerini](https://www.kth.se/profile/ffn) <br />
**Methodology:** [Youssef Almulla](https://www.kth.se/profile/almulla) and [Camilo Ramirez](https://github.com/camiloramirezgo) on the energy model, [Youssef Almulla](https://www.kth.se/profile/almulla) on the decarbonization strategies of the agricultural sector, [Camilo Ramirez](https://github.com/camiloramirezgo) on the softlinking of models, [Brian Joyce](https://www.sei.org/people/brian-joyce/) on the water-agriculture (WEAP-MABIA) model, all on the participatory approach and overall Nexus approach <br />
**Software:** [Camilo Ramirez](https://github.com/camiloramirezgo) & [Youssef Almulla](https://github.com/JZF07) <br />
**Management and Advisory support:** [Annette Huber-Lee](https://www.sei.org/people/annette-huber-lee/) & [Francesco Fuso-Nerini](https://www.kth.se/profile/ffn)<br />

**Acknowledgements** <br />
We would like to acknowledge the Food and Agricultural Organization (FAO), the Jordan Ministry of Water and Irrigation, the Regional Office for Agricultural Development in Souss-Massa (ORMVA), and a special thanks goes to Domitille Vallee and Jiro Ariyama (FAO), and Prof. Lahcen Kenny (FAO consultant) for their valuable contribution to this work by providing many of the input data and validating modelling assumptions.
