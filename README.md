# pRocESsor
Processor for pRES surveys (v0.1; for mobile pRES profiling)

The pRocESsor was developed to handle profiling data acquired with the (A)pRES radar, developed by the British Antarctic Survey (BAS). It is implemented with a modular approach to be expandable for further applications and acquisitions methods.

To run the pRocESsor, a configuration file as stored under './config' needs to be created and run using 'run_pRocESsor.m'. The processed data are stored as defined by the settings 'dirProc'/'fileProc'.

The processed data can be obtained by running
  addpath(genpath(path_to_pRocESsor))
  cfg = ConfigHandler(config_file);
  data = load(cfg.fileProc);

The implemented pRES profile processing is detailed in
Oraschewski et al., [preprint], Layer-optimized SAR processing with a mobile phase-sensitive radar for detecting the deep englacial stratigraphy of Colle Gnifetti, Switzerland/Italy, https://doi.org/10.5194/egusphere-2023-2731.
