# RME-model
Matlab Codes for Adjustment of Radar-Gauge Rainfall Discrepancy Due to Raindrop Drift and Evaporation

## Description
### Introduction 
It must be fully considered that raindrops observed by weather radars cannot fall vertically to the ground without changing in size. In the process of rainfall, raindrop location changes due to wind drift and raindrop size changes due to evaporation. Our study proposes a fully formulated scheme to numerically simulate both raindrop drift and evaporation in the air and reduces the uncertainties of radar rainfall estimation, which contributes to the improvement of quantitative precipitation estimation from radar polarimetry and allows a better understanding of precipitation processes.

We provide here some matlab codes for Discrepancy Adjustment Model. All these codes are used in our experiments.
	
### Runfiles
The files in the runfiles/directory are used to set parameters and run the Discrepancy Adjustment Model. You can create your own runfiles by copying an existing one and then run with different parameters and runfile combination.

These runfiles are included:
1.	**`run_files.m`**: Script file that calls for other functions of the model
2.	**`ComputeRainDSD (path)`**: Compute DSD of rainfall
3.	**`RME_SetupModel(path)`**: STEP 1 of RME model. Get the permanent info to run the model, including domain_net and layers.
4.	**`RME_InitModel (path, events)`**: STEP 2 of RME model. Get the event basic info to run the model.
5.	**`RME_RunModel (path, events)`**: STEP 3 of RME model. Simulate the raindrop microphysical process for each event.
6.	**`RME_PostModel (path, events)`**: STEP 4 of RME model. Process the result of the raindrop microphysical simulation.

## Publication
We have published an academic paper which titled “Adjustment of Radar‐Gauge Rainfall Discrepancy Due to Raindrop Drift and Evaporation Using the Weather Research and Forecasting Model and Dual‐Polarization Radar”.

Please visit <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019WR025517> for more details and results of experiments.

Citation: Dai Q, Yang Q, Han D, Rico‐Ramirez M A, Zhang S. Adjustment of radar‐gauge rainfall discrepancy due to raindrop drift and evaporation using the Weather Research and Forecasting model and dual‐polarization radar[J]. Water Resources Research, 2019, 55(11): 9211-9233.

## Contact
If you have some problems or find some bugs in the codes, please email: qd_gis@163.com.
