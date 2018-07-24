# README #

Repositorio para datos del proyecto PACA de síntesis planetaria

## General use 

The following sections explain how to use the scripts for processing simulation data.

### Using filtered.zip

- In filtered.zip there will be 3 csv files (no perturbation, low perturbation, and high perturbation) with the following columns:


| Column Name   | Meaning       | Importance (*)  | Initial parameter? |
| ------------- |:-------------:| :-----:| ----------|
| ident     	| unique system identifier | high | Y |
| it   			| iteration index      |   none | Y |
| t 			| simulation time (yr)      | low    | Y |
| a(i)			| planet semi-major axis (AU)| high		| |
| emegas(i)		| planet gas mass (M_earth) | medium | |
| emepla(i)/emet| planet solid mass (M_earth)| high | |
| rplanet(i)/radtie| planet solid radius | low | |
| emestar		| stellar mass (g) | high | Y |
| rc            | disk outer cutoff radius (AU) | high | Y |
| qest | Toomre Q at min radius | none | Y |
| sigmag_0 | maximum dust surface density (g/cm^3) | high | Y |
| emed | disk mass (M_sun) | high | Y |
| gama | surface density power law exponent | none |Y 
| apert | perturbation amplitude | high|Y 
| fpert | perturbation length scale | none |Y 
| constmigI | type I migration rate | none |Y 
| emetal | metallicity wrt solar | high |Y 
| taugas | gas dissipation timescale (yr) | high | Y 

(*) For data in filtered.zip

- The master distributions for systems in filtered are:
![](distributions.png?raw=true)
- All three files in filtered.zip were generated using the same initial parameters except for the perturbation parameter Apert (= 0, 0.1, 0.3).

### If using repo locally

- Clone this repo (planetas_data)
- Choose a .zip in the following list: `properyam.zip pab*.zip yamqc*.zip`. Other zip files in the repo are not formatted for data reduction pipeline, so you can just unzip them and look inside.
- Type in your terminal `sh pipeline.sh properyam` or whatever your chosen .zip filename prefix is. Note: do not include the .zip portion.
- Your processed data and simulation info is now in the `/properyam/` folder (or whatever the name of your .zip filename prefix is).
- Raw data will be in finalresults.csv (containing info on all the planets).
- Reduced data (containing consolidated info on each system) will be in `tr_finalresults.csv` and `gt_finalresults.csv`. See below for more information.
- Simulation data is in some .txt files as described in the output to `pipeline.sh`

## Specific use - Analysis tools

- `quality.py` creates some diagnostic .txt files. It is meant to be used as `python quality.py createdfiles.txt totalsims.txt finished.txt` within a folder that also contains the `results/` folder. 
- `reduce.py` is meant to be used with `finalresults.csv` files created by `pipeline.sh`. This python script creates two data files, `tr_finalresults.csv` and `gt_finalresults.csv`. The reduced data files contain Center of Mass, Number of Planets, Mass budget for planets, Mass efficiency (mass budget/disk mass), Sigmag0 (maximum solid surface density), Disk Mass, Characteristic radius, Stellar Mass, Metallicity and Gas depletion timescale information for each simulation. `tr_finalresults.txt` has this information for planets below 10 earth masses, and `gt_finalresults.txt` for planets above 10 earth masses. It should be used as `python reduce.py finalresults.csv`. 

## References
- Planet formation model (Y. Miguel) https://arxiv.org/pdf/1106.3281.pdf
- Models of population synthesis (W. Benz et al) https://arxiv.org/pdf/1402.7086.pdf
- Perturbation model (P. Pinilla) https://arxiv.org/pdf/1112.2349.pdf
