function cfg = config_T1_LOSAR()
% Configuration file for pRocESsor.
%
% Mobile pRES profile
% Transect 1, 2021 Colle Gnifetti campaign, Uni Tübingen
% Processing: Layer-optimized SAR
%
% Falk Oraschewski, 01.08.2023

%% Mandatory settings
cfg.typeSurvey      = 'profile';

%% File settings
% Root directory of project data
cfg.dirProject      = '..'; % Path to project directory, can be absolute or relative.

% Relative directories to different types of project data
cfg.dirRaw      = 'data/raw_data/transect1'; % Subpath to raw data (.dat-files)
cfg.dirPreProc  = 'data/preproc_data'; % Subpath to pre-processed data
cfg.dirProc     = 'data/proc_data_mat'; % Subpath to final processed data output
cfg.dirGPS      = 'data/raw_data/gps'; % Subpath to GPS data
cfg.dirSup      = 'data/supplement_data'; % Subpath to any supplementary data

% Filenames
cfg.filePreProc = 'T1_p8_md100.mat'; % Pre-processed data file (after fmcw_range is applied).
cfg.fileProc    = 'T1_LOSAR_filtered.mat'; % Processed output data file
cfg.filePos     = 'T1_pos_smooth.mat'; % Pre-processed positioning data file
cfg.fileGPS     = 'gpsAll.csv'; % Input file for GPS data (currently only .csv files are supported with the fields 'filename', 'latitude', 'longitude' and 'elevation').
cfg.fileDensity = 'densityKCC_LinearlyExtrapolated.csv'; % Density data (currently only .csv files are supported with the fields 'range' and 'density').

% File update settings
cfg.newPreProc  = false; % Option to update cfg.filePreProc
cfg.newPos      = true; % Option to update cfg.filePos
cfg.newSlopes   = true;
cfg.newProc     = true; % Option to update cfg.fileProc


%% Processing settings
% General settings
cfg.lowestNoise     = true;  % Option to select burst with lowest noise floor in each trace
cfg.correctDensity  = true;  % Option to correct for density 
cfg.cropCables      = true;  % Option to correct for cable length
cfg.matchShape      = false; % Option to apply a shape matching filter (following Rahman et al., 2014)

% Survey type specific
cfg.methodSAR       = 'layer';          % Options: 'layer', 'interpolation'
cfg.lengthSAR       = 5;                % Synthetic aperture length (m)
cfg.methodPos       = 'smooth';         % Position processing method; Options: 'smooth'
cfg.methodSlope     = 'linefitting';
cfg.filterSlopes    = 1;
cfg.paramSmoothing  = 0.3;
cfg.profileSpacing  = 0.1; % Profile 
cfg.useAntennaPos   = true;
cfg.defineLines     = true;
cfg.numLines        = {1};%, 2};
cfg.indLines        = {[1:1265]};%, [1702:2414]};

cfg.plotProcessing = true;

%% Survey setup settings
cfg.maxDepth    = 100; % MaxDepth for fmcw_range
cfg.padding     = 8; % Padding factor for fmcw_range
cfg.selInd      = 1:1265; % Select subset of input files (Can be used to exclude faulty traces or parts of the data).
cfg.skipFirst   = 0; % Skip files at start (performed after selInd)
cfg.skipLast    = 0; % Skip files at end (performed after selInd)

% Radar setup
cfg.cableTx     = 5; % Cable length for transmitting antenna (m)
cfg.cableRx     = 5; % Cable length for receiving antenna (m)
cfg.distAntenna = 2.7; % Antenna seperation (m)

% Position settings
cfg.useGPS = true; % Use GPS data
cfg.typeGPS = 'csv'; % GPS input type, options: 'csv'
cfg.surveyEPSG = 2056; % EPSG code for Cartesian coordinate system for the project


%cfg.multilooking = false;
%cfg.multilooking_pixels = [3 2];  % In along- and across track direction.



%cfg.signal_correction = true;
%cfg.signal_bins = 250; % Bins at the bottom of the image used for noise computation.


end