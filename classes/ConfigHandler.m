classdef ConfigHandler
    % CONFIGHANDLER create a config object that includes all settings for
    % processing an ApRES survey.
    
    properties
        % Survey type
        typeSurvey = 'profile'; % Survey type (Currently only option: 'profile')
        
        % File settings
        dirProject              % Root directory for project
        
        dirRaw                  % Directory for raw/input files
        dirPreProc              % Directory for pre-processed pRES data
        dirProc                 % Directory for processed output
        dirGPS                  % Directory for GPS files
        dirSup                  % Directory for supplementary inputs
        
        filePreProc             % Pre-processed data file
        fileProc                % Processed data file
        filePos                 % Processed position data file
        fileGPS                 % GPS data file
        fileDensity             % Density data file
        
        % Update settings
        newPreProc = false;     % Flag to re-perform pre-processing
        newPos = false;         % Flag to re-process positioning
        newSlopes = true;       % Flag to re-process slopes
        newProc = false;        % Flag to re-process final output
        
        % ApRES settings
        T = 1;                  % Pulse duration [s]
        B = 2e8;                % Bandwidth [Hz]
        f0 = 2e8;               % Start frequency [Hz]
        fc = 3e8;               % Center frequency [Hz]
        K                       % Sweep rate (K = 2pi*B/T)
        icePermittivity = 3.17; % Ice permittivity
        
        % Survey setup settings
        maxDepth = 1000;        % Maximum depth
        padding = 1;            % Padding
        skipFirst = 0;          % Skip first data points
        skipLast = 0;           % Skip last data points
        selInd = [];            % Selected indices
        
        % Radar settings
        useAttenuation = 1      % Selected attenuation settings
        useTx = 1               % Selected transmitter (currently has no effect)
        useRx = 1               % Selected receiver (currently has no effect)
        cableTx                 % Transmitter cable length
        cableRx                 % Receiver cable length
        distAntenna             % Distance between antennas
        
        % Positioning settings
        useGPS = true;          % Flag indicating usage of GPS
        typeGPS                 % Type of GPS
        surveyEPSG double       % Survey EPSG code
        
        % Optional settings
        lowestNoise = false     % Flag to select burst with lower noise floor
        burstMean = false       % Option to compute mean of all chirps in burst
        correctDensity = false  % Flag to use density correction
        cropCables = false      % Flag to use cable cropping
        matchShape = false      % Flag to use signal matching
        
        % Survey type specific settings
        methodProc                  % Processing method (options are survey specific)
        lengthSAR                   % SAR length
        methodPos                   % Positioning method
        methodSlope = 'linefit';    % Slope method
        filterSlopes = true         % Flag to use of slope filtering
        slopeMax = 10;              % Slope range boundaries for 'linefit' [degrees]
        slopeInterval = 0.1;        % Slope range interval for 'linefit' [degrees]
        slopeFilterWin = [10 10];   % Slope filtering window [bins]

        posSmoothing = 0.3;     % Smoothing parameter
        profileSpacing = 0.1;   % Profile spacing [m]
        useAntennaPos           % Flag to use antenna position (currently this has no effect)
        
        defineLines             % Flag to define seperate lines in data
        numLines                % Number of lines
        indLines                % Indices of lines
    end
    
    methods
        function obj = ConfigHandler(configFile)
            % CONFIGHANDLER Constructor for ConfigHandler class.
            %   Initializes the ConfigHandler object using the provided
            %   configuration file.
            %
            %   configFile: A structure containing configuration settings.
            %               Fields should correspond to class properties.
            %
            %   Note: At the moment the config file needs to be written as 
            %   a Matlab function, but switching to .json file or the like
            %   is intended in the future.
            
            if nargin > 0
                % Update properties with values from the config structure
                fields = fieldnames(configFile);
                allowedFields = properties(obj);
                for i = 1:length(fields)
                    if ismember(fields{i}, allowedFields)
                        obj.(fields{i}) = configFile.(fields{i});
                    else
                        error('Invalid field name in config: %s. \nAllowed fields are: %s', fields{i}, strjoin(allowedFields, ', '));
                    end
                end
            end

            % Derived properties
            obj.K = 2 * pi * obj.B / obj.T;

            % Obtain absolute file paths
            obj = obj.obtainFullfiles();
        end
    end

    methods (Access = private)
        function obj = obtainFullfiles(obj)
            % OBTAINFULLFILES Obtain absolute file paths for configHandler
            %
            % Description:
            %   For directiory and file properties in the configHandler 
            %   object, absolute paths are obtained based on the root
            %   directory of the project.
            
            % Check if the project directory exists and ensure it is an absolute path
            if exist(obj.dirProject, 'dir') == 7
                obj.dirProject = fullfile(obj.dirProject);
            else
                % If the directory does not exist, assume obj.dirProject is
                % a relative path from a general data folder set by a
                % system variable 'dirData'.
                % Note from Falk: I use this to access data in my local
                % root data directory.
                obj.dirProject = fullfile(dirData, obj.dirProject);
            end

            % Combine sub-directories with project directory
            dirProperties = ["dirRaw", "dirPreProc", "dirProc", "dirGPS", "dirSup"];
            for n = 1:length(dirProperties)
                if ~isempty(obj.(dirProperties(n)))
                    obj.(dirProperties(n)) = fullfile(obj.dirProject, obj.(dirProperties(n)));
                end
            end

            % Associate files with directories
            fileProperties  = ["filePreProc", "filePos",    "fileProc", "fileGPS", "fileDensity"];
            dirRelated      = ["dirPreProc",  "dirPreProc", "dirProc",  "dirGPS",  "dirSup"];

            for n = 1:length(fileProperties)
                if ~isempty(obj.(fileProperties(n)))
                    if ~isempty(obj.(dirRelated(n)))
                        obj.(fileProperties(n)) = fullfile(obj.(dirRelated(n)), obj.(fileProperties(n)));
                    else
                        error("The file variable '%s' was set, but the corresponding directory variable '%s' is missing. Please include it in the config file.", fileProperties(n), dirRelated(n))
                    end
                end
            end
        end
    end
end