classdef ConfigHandler
    % CONFIGHANDLER: A class that handles configuration settings for the pRES survey processing.
    %   This class encapsulates various properties representing different aspects
    %   of the survey setup and processing settings.
    
    properties
        % Survey type
        typeSurvey  % Type of survey (option: profile)
        
        % File settings
        dirProject  % Root project directory
        
        dirRaw      % Directory for raw/input files
        dirPreProc  % Directory for pre-processed pRES data
        dirProc     % Directory for processed output
        dirGPS      % Directory for GPS files
        dirSup      % Directory for supplementary inputs
        
        filePreProc % Pre-processed data file
        filePos     % Processed position data file
        fileProc    % Processed data file
        fileGPS     % GPS data file
        fileDensity % Density data file
        
        % Update settings
        newPreProc = false; % Flag indicating to (re-processed data
        newPos = false;     % Flag indicating new position data
        newSlopes = true;  % Flag indicating new slope data
        newProc = false;    % Flag indicating new output data
        
        % Default pRES settings
        T = 1;               % Pulse duration
        B = 2e8;             % Bandwidth
        f0 = 2e8;            % Start frequency
        fc = 3e8;            % Center frequency
        K                    % Sweep rate (K = 2pi*B/T)
        icePermittivity = 3.18; % Ice permittivity
        useAttenuators      % Flag indicating usage of attenuators
        
        % Survey setup settings
        maxDepth = 1000;     % Maximum depth
        padding = 1;         % Padding
        skipFirst = 0;       % Skip first data points
        skipLast = 0;        % Skip last data points
        selInd = [];         % Selected indices
        
        % Radar settings
        useAttenuator = 1   % Flag indicating usage of attenuator
        useTx = 1           % Flag indicating usage of transmitter
        useRx = 1           % Flag indicating usage of receiver
        cableTx             % Transmitter cable length
        cableRx             % Receiver cable length
        distAntenna         % Distance between antennas
        
        % Positioning settings
        useGPS = true;      % Flag indicating usage of GPS
        typeGPS             % Type of GPS
        surveyEPSG double  % Survey EPSG code
        
        % Optional settings
        lowestNoise = false      % Option to select burst with lower noise floor
        burstMean = false        % Option to compute mean of all chirps in burst
        correctDensity = false  % Flag indicating usage of density correction
        cropCables = false      % Flag indicating usage of cable cropping
        matchShape = false      % Flag indicating usage of signal matching
        
        % Survey type specific settings
        lengthSAR           % SAR length
        methodPos           % Positioning method
        methodSlope = 'linefit'; % Slope method
        filterSlopes = true  % Flag indicating usage of slope filtering
        slopeMax = 10;       % Slope range boundaries for 'linefit'
        slopeInterval = 0.1; % Slope range interval for 'linefit'
        methodProc                  % Processing method (options are survey specific)
        
        paramSmoothing      % Smoothing parameter
        profileSpacing      % Profile spacing
        useAntennaPos       % Flag indicating usage of antenna position
        
        defineLines         % Definition of lines
        numLines            % Number of lines
        indLines            % Indices of lines
    end
    
    methods
        function obj = ConfigHandler(configFile)
            % CONFIGHANDLER: Constructor for ConfigHandler class.
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
            % OBTAINFULLFILES: Obtains absolute file paths for directories and files.
            %   Updates directory and file properties to contain absolute
            %   paths based on the project directory.
            
            % Check if the project directory exists and ensure it is an absolute path
            if exist(obj.dirProject, 'dir') == 7
                obj.dirProject = fullfile(obj.dirProject);
            else
                % If the directory does not exist, assume obj.dirProject is
                % a relative path from a general data folder set by a
                % system variable 'dirData'.
                % Note, Falk: This is an approach that I use to access data
                % in my local root data directory.
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