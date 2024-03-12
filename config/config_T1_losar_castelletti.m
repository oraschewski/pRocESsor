function cfg = config_T1_losar_castelletti()
% Configuration file for pRocESsor.
%
% Mobile pRES profile
% Transect 1, 2021 Colle Gnifetti campaign, Uni Tübingen
% Processing: Layer-optimized SAR with horizontal summation
%
% Falk Oraschewski, 11.03.2024

%% Mandatory settings
cfg.typeSurvey      = 'profile';

%% File settings
% Root directory of project data
cfg.dirProject      = '..'; % Path to project directory, can be absolute or relative.

% Relative directories to different types of project data
cfg.dirRaw      = 'data/raw_data/transect1';    % Subpath to raw data (.dat-files)
cfg.dirPreProc  = 'data/preproc_data';          % Subpath to pre-processed data
cfg.dirProc     = 'data/proc_data_mat';         % Subpath to final processed data output
cfg.dirGPS      = 'data/raw_data/gps';          % Subpath to GPS data
cfg.dirSup      = 'data/supplement_data';       % Subpath to any supplementary data

% Filenames
cfg.fileProc    = 'T1_losar_filtered_castelletti.mat'; % Processed output data file
cfg.filePreProc = 'T1_p8_md100.mat';            % Pre-processed data file (after fmcw_range is applied).
cfg.filePos     = 'T1_pos_smooth.mat';          % Pre-processed positioning data file
cfg.fileGPS     = 'gpsAll.csv';                 % Input file for GPS data (currently only .csv files are supported with the fields 'filename', 'latitude', 'longitude' and 'elevation').
cfg.fileDensity = 'densityKCC_LinearlyExtrapolated.csv'; % Density data (currently only .csv files are supported with the fields 'range' and 'density').

% File update settings
cfg.newPreProc  = false; % Option to update cfg.filePreProc
cfg.newPos      = false; % Option to update cfg.filePos
cfg.newProc     = true; % Option to update cfg.fileProc
cfg.newSlopes   = false; % Option to avoid reprocessing the slope


%% Processing settings
% General settings
cfg.lowestNoise     = true;  % Option to select burst with lowest noise floor in each trace
cfg.cropCables      = true;  % Option to correct for cable length
cfg.correctDensity  = true;  % Option to correct for density 
cfg.matchShape      = false; % Option to apply a shape matching filter (following Rahman et al., 2014)

% Survey type specific
cfg.methodSAR       = 'losar_castelletti';      % Options: 'losar', 'interpolation', 'movmean'
cfg.lengthSAR       = 5;            % Synthetic aperture length (m)
cfg.methodPos       = 'smooth';     % Position processing method; Options: 'smooth', 'org'
cfg.methodSlope     = 'linefit';    % Slope estimation methods; options: 'linefit', 'compute' 
cfg.filterSlopes    = true;         % Option to filter slopes using a 2d moving median filter
cfg.paramSmoothing  = 0.3;          % Smoothing parameter for positioning
cfg.profileSpacing  = 0.1;          % Equidistant spacing for processed profile 
cfg.useAntennaPos   = false;        % Estimate antenna positions (currently this has no effect)
cfg.defineLines     = false;        % Option to process data linewise
cfg.numLines        = {1};          % Line numbers
cfg.indLines        = {[1:1265]};   % Indices associated with line numbers

%% Survey setup settings
cfg.maxDepth    = 100;      % MaxDepth for fmcw_range
cfg.padding     = 8;        % Padding factor for fmcw_range
cfg.selInd      = 1:1265;   % Select subset of input files (Can be used to exclude faulty traces or parts of the data).
cfg.skipFirst   = 0;        % Skip files at start (performed after selInd)
cfg.skipLast    = 0;        % Skip files at end (performed after selInd)

% Radar setup
cfg.cableTx     = 5;    % Cable length for transmitting antenna (m)
cfg.cableRx     = 5;    % Cable length for receiving antenna (m)
cfg.distAntenna = 2.7;  % Antenna seperation (m)

% Position settings
cfg.useGPS = true;      % Use GPS data
cfg.typeGPS = 'csv';    % GPS input type, options: 'csv'
cfg.surveyEPSG = 2056;  % EPSG code for Cartesian coordinate system for the project
end