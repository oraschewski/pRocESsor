# pRocESsor

## For layer-optimized SAR processing of mobile pRES profiles (v1.0)

pRocESsor was developed to handle profiling data acquired with a mobilized Autonomous phase-sensitive Radio Echo Sounder (ApRES), in the following referred to as mobile pRES. Currently, the main application of pRocESsor is to improve the imaging of englacial stratigraphy recorded with a mobile pRES using Layer-optimized SAR processing. pRocESsor is implemented with a modular approach, enabling to include other applications and acquisition methods of the ApRES.

The processing method is detailed in:  
  Oraschewski et al. (2024) Layer-optimized SAR processing with a mobile phase-sensitive radar: a proof of concept for detecting the deep englacial stratigraphy of Colle Gnifetti, Switzerland/Italy. The Cryosphere [preprint: https://doi.org/10.5194/egusphere-2023-2731].

**Author:**  
Falk Oraschewski  
University of Tübingen  
04.07.2024

The ApRES was developed by the British Antarctic Survey (BAS) and the University College London (UCL). 

## Run pRocESsor

Using `run_pRocESsor.m`, the input data is processed as defined in the selected configuration file. Example config files are povided under `./config`.

Raw ApRES `.dat`-files have to be stored in `dirRaw`. Intermediate pre-processed data after applying the standard FMCW signal processing is stored in `dirPreProc/filePreProc` and the final output data in `dirProc/fileProc` as `.mat`-files.  
In addition, GPS data and firn density data can be provided as `.csv`-files in `dirGPS/fileGPS` and `dirSup/fileDensity`.

Processed data can be loaded quickly by running:  
    `addpath(genpath(path_to_pRocESsor))`  
    `cfg = ConfigHandler(config_file);`  
    `data = load(cfg.fileProc);`  

