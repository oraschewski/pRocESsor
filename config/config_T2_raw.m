function cfg = config_T2_raw()
% Configuration file for pRocESsor.
%
% Mobile pRES profile
% Transect 2, 2021 Colle Gnifetti campaign, Uni Tübingen
% Processing: Interpolation -> Data are only rescaled to equidistant grid.
%
% Falk Oraschewski, 01.08.2023

%% Define survey type
cfg.typeSurvey      = 'profile';

%% File settings
% Root directory of project data
cfg.dirProject      = '..'; % Path to project directory, can be absolute or relative.

% Relative directories to different types of project data
cfg.dirRaw      = 'data/raw_data/transect2';    % Subpath to raw data (.dat-files)
cfg.dirPreProc  = 'data/preproc_data';          % Subpath to pre-processed data
cfg.dirProc     = 'data/proc_data_mat';         % Subpath to final processed data output
cfg.dirGPS      = 'data/raw_data/gps';          % Subpath to GPS data
cfg.dirSup      = 'data/supplement_data';       % Subpath to any supplementary data

% Filenames
cfg.fileProc    = 'T2_interp.mat';              % File for processed output data
cfg.filePreProc = 'T2_p8_md100.mat';            % File for pre-processed data (after FMCW signal processing is applied).
cfg.filePos     = 'T2_pos_smooth.mat';          % File for pre-processed positioning data
cfg.fileGPS     = 'gpsAll.csv';                 % File for GPS data     (currently only .csv files are supported with the fields 'filename', 'latitude', 'longitude' and 'elevation').
cfg.fileDensity = 'densityKCC.csv';             % file for density data (currently only .csv files are supported with the fields 'range' and 'density').

% File update settings
cfg.newPreProc  = false;            % Option to update cfg.filePreProc
cfg.newPos      = false;            % Option to update cfg.filePos
cfg.newProc     = true;             % Option to update cfg.fileProc

%% General processing settings
% Processing options
cfg.lowestNoise     = true;         % Option to select burst with lowest noise floor in each trace
cfg.cropCables      = true;         % Option to correct for cable length
cfg.correctDensity  = true;         % Option to correct for firn density 
cfg.matchShape      = false;        % Option to apply a shape matching filter (following Rahman et al., 2014)

% Position settings
cfg.useGPS          = true;         % Load GPS data associated with .dat-files
cfg.typeGPS         = 'csv';        % GPS input type, options: 'csv'
cfg.surveyEPSG      = 2056;         % EPSG code for Cartesian coordinate system for the project

% FMCW signal processing settings
cfg.maxDepth        = 100;          % MaxDepth for fmcw_range [m]
cfg.padding         = 8;            % Padding factor for fmcw_range

% File selection settings
cfg.selInd          = 1:713;        % Select subset of input files (Can be used to exclude faulty traces or parts of the data).
cfg.skipFirst       = 0;            % Skip files at start (performed after selInd)
cfg.skipLast        = 0;            % Skip files at end (performed after selInd)

% Radar setup
cfg.cableTx         = 5;            % Cable length for transmitting antenna [m]
cfg.cableRx         = 5;            % Cable length for receiving antenna [m]
cfg.distAntenna     = 2.7;          % Antenna seperation [m]

%% Survey type specific settings
% Processing options
cfg.methodProc      = 'interp';     % Options: 'losar', 'interpolation', 'movmean'
cfg.methodPos       = 'smooth';     % Position processing method; Options: 'smooth', 'org'

% Processing settings
cfg.lengthSAR       = 5;            % Synthetic aperture length [m]
cfg.profileSpacing  = 0.1;          % Equidistant spacing for processed profile [m]
cfg.posSmoothing    = 0.3;          % Smoothing parameter for profile positioning

% Linewise processing
% cfg.defineLines     = false;        % Option to process data linewise
% cfg.numLines        = {1};          % Line numbers
% cfg.indLines        = {[1:713]};   % Indices associated with line numbers
end
