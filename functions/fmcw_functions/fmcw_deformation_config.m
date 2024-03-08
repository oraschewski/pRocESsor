function cfg = fmcw_deformation_config

% Default configuration file for fmcw_process SN049 - R14
%
% Edit this and rename for custom processing settings.

% Close previous run
%close all
%clear

% Default config
cfg.notes = 'default config file';

% Use test data
cfg.useTestFiles = 0;
data_dir = '..\..\..\..\..\data\foraschewski\data_large\ApRES\HammarryggenIceRise';
cfg.filename1 = fullfile(data_dir,'p0_Survey_2019-01-07_132649_VV.dat');
cfg.filename2 = fullfile(data_dir,'p0_Survey_2019-12-04_170339_VV.dat');

% Pre-processing
cfg.fRange = [2e8 4e8]; % Tx frequency range to use - full range
%cfg.fRange = [2.5e8 3.5e8];  % half bandwidth
cfg.doManualChirpSelect = 0;
cfg.fchirplist = [1:100];
cfg.gchirplist = [57:100];
cfg.doClean = 0;
cfg.noisePowerLimit = 0.0015; % allowable percent difference in power from quadtratic fit to chirpnum-power trend
cfg.nNoisest = 2;

% phase-processing
cfg.p = 20; % padfactor (interpolation factor gi generating profile)
% note: higher pad factor increases the likelihood of estimating the
% correct bin  lag.
cfg.maxRange = 700; % maximum bed range
cfg.winFun = @blackman;
%cfg.winFun = @blackmanharris; % less spectral leakage so better to pick up near bed internals

% Find bed (or max range for vsr)
cfg.bedMethod = 'xcor'; % 'ampThresh' 'maxAmp' 'xcor' ???
cfg.bedMaxShift = 1; % m
cfg.ampThreshdB = -50; % dB - minimum bed strength
cfg.bedSearchRange = [500 600]; % bed search range (m)

% Bed alignment
cfg.doAlignBed = 1; %

% Range error estimate
cfg.errorMethod = 'emperical'; % 'emperical' 'assumedNoiseFloor'
cfg.noiseFloordB = -100; % Assumed level of noise


% Chunk lag matching (depth-dependent co-registration)
cfg.minDepth = 10; % to avoid breakthrough and cables etc (cables 2m in 2013, 5m in 2013).
cfg.bedBuffer = 10; % m buffer to exclude spectral leakage from bed return
cfg.coarseStep = 4;
cfg.coarseChunkWidth = 8; % long segments more uniquely define lag - except if there is high strain
cfg.maxStrain = 0.005; % maximum strain magnitude to search for
cfg.minAmpCor = 0.95; % Minimum amplitude correlation to use
cfg.minAmpCorProm = 0.05; % Minimum difference between max correlation and next best local maximum

% Chunk phase difference
cfg.fineStep = 4;
cfg.fineChunkWidth = 8; % between 4 to 8 is a good compromise
cfg.doUseCoarseOffset = 1; % uses coarse offset determined above to specify rough lag for fine offset
cfg.CoarseOffsetMethod = 'interpolate'; % Options: 'polynomial', 'exponential', 'interpolate'
cfg.polyOrder = 1;

% Strain fitting
cfg.doStrainPointwise = 1; % Determine individual strain rates between points.
cfg.doStrainFit = 1; % Fit strain rates
cfg.minCohereCoarse = 0.5; % minimum correlation on lag offset (i.e. confidence we have the right lag)
cfg.minCohereFine = 0.5; % minimum coherence to use in strain estimate
cfg.minCoherePhase = 0.5; % phase coherence over win...
cfg.movingMean = 5; % number of reflector pairs used for determining the moving mean.
cfg.firnDepth = 100; % used to select depth range for vertical strain estimate
cfg.fitMethod = 'robust'; % 'menke' 'robust' 'regress'
cfg.fit.type = 'linear'; % 'linear' 'quadratic'
cfg.fit.errorMetric = 'standardError'; % 'standardError' 'alpha'
cfg.minPointsToFit = 3; % minimum number of depth offset to use in estimate

% Bulk lag matching (co-registration)
cfg.doBulkAllignment = 0; %
cfg.bulkAlignRange = [40 80];
cfg.maxOffsetM = 10; % 10m recoverable offset near surface

% Output
cfg.verbose = 1; % results to screen
cfg.doSaveOutput = 0; % save mat file

% Plots
cfg.doPlotAll = 0; % plot lots of other stuff

% Individual control of all plots
cfg.doPlotAlignBulk = 0;
cfg.doPlotAlignCoarse = 0;
cfg.doPlotAlignFine = 1;
cfg.doPlotVSR = 0;
cfg.doPlotFit = 0;
cfg.doPlotBedShift = 0;
cfg.doPlotMelt = 0;
cfg.doPlotResid = 0;

