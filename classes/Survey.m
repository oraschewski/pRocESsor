classdef Survey
    % SURVEY is the base class for all possible ApRES surveys.
    %   The ApRES can be employed in various ways for different
    %   applications. While many of the required processing steps for the
    %   different survey types differ, many are also the same.
    %
    %   The survey object is the base class for this project which entails
    %   all the processing steps that the different survey types have in
    %   common.

    properties
        % Surve type
        type string

        % Data Properties 
        range double
        rangeCor double
        twtt double
        spec double
        specCor double
        lat double = NaN;
        lon double = NaN;
        x double = NaN;
        y double = NaN;
        elev double = NaN;
        pos double = NaN;
        distGPS double = NaN;
        nSubBursts uint16
        nBurstSel double
        nAttenuators uint16
        attenuator1 int16
        afGain int16
        txAnt int16
        rxAnt int16
        temp1 double
        temp2 double
        batteryVoltage double
        time double
        timestamp datetime
    	dtBin double

        % Radar settings
        T double
        f0 double
        f1 double
        fc double
        B double
        K double
        dt double
        permittivity double = 3.17;
        ci double
        lambdaC double
        padding double = 2;
        
        % Radar setup
        cableTx double = NaN;
        cableRx double = NaN;
        distAntenna double = NaN;

        % File properties
        numFiles uint16
        fileName string
        
        % Metadata
        flags = [];
        config
    end
    
    properties(Constant)
       version = 1.0;
    end

    methods (Static)
        function obj = initiateSurvey(config)
            % INITIATESURVEY creates or loads survey object.
            
            % Check if a file for the pre-processed survey data already
            % exsists and shall not be overwritten.
            if ~config.newPreProc && exist(config.filePreProc, 'file') == 2
                % Load existing Survey object.
                pRESdata = load(config.filePreProc,'-mat');

                % Unfold survey object.
                obj = pRESdata.pRESdata;

                % Update configurations in Survey object.
                obj.config = config;
            else
                % Otherwise, create new survey object depending on the
                % survey type.
                switch lower(config.typeSurvey) % Make survey typeSurvey is case insensitive.
                    case 'profile'
                        obj = Profile(config);
                    %case 'timeseries'
                    %    obj = Timeseries(config);
                    %case 'repeat'
                    %    obj = Repeat(config);
                    otherwise
                        error('Invalid survey type.');
                end

                %%% Saving the new survey object.
                % Create survey object
                pRESdata = obj;

                % Exclude local folder structure from config field for file sharing.
                cs = obj.config;
                [cs.dirProject, cs.dirRaw, cs.dirPreProc, cs.dirProc, cs.dirGPS, cs.dirSup, cs.filePreProc, cs.fileProc, cs.filePos, cs.fileGPS, cs.fileDensity] = deal([]);
                % Update config field
                pRESdata.config = cs;
                
                % Save survey object
                save(config.filePreProc,'pRESdata');
                clear pRESdata
            end
        end
    end


    methods
        function obj = Survey(config)
            % SURVEY Constructor for Survey class
            %   This method creates the base survey object, providing the
            %   pre-processed pRES data. It is applicable to all survey
            %   types.
            
            % Set survey type
            obj.type = config.typeSurvey;

            % Find all .dat-files in raw data directory.
            rawFiles = dir(fullfile(config.dirRaw, '**\*.dat'));
            if isempty(rawFiles)
                rawFiles = dir(fullfile(config.dirRaw, '**\*.DAT'));
            end
            
            % Use all files if no subset is selected.
            if isempty(config.selInd)
                config.selInd = 1:length(rawFiles);
            end

            % Potentially skip few first and last files.
            config.selInd = config.selInd(1 + config.skipFirst:end - config.skipLast);

            % Create wait bar for reading files.
            numFiles = length(config.selInd);
            disp(['Reading ',num2str(numFiles),' files from ',config.dirRaw, newline]);
            wb = waitbar(0,'Reading Files');

            % Iterate over all selected files.
            for ii = 1:numFiles
                % Obtain file name
                iiSelInd = config.selInd(ii);
                fileName = fullfile(rawFiles(iiSelInd).folder, rawFiles(iiSelInd).name);
                %disp(fileName)

                % Load the file.
                vdat = fmcw_load(fileName);
                
                % ToDo: Include functionality to handle MIMO files and
                % multiple attenuations.
                    % Number of combinations that are used
                    %numUseVdats = config.useAttenuator * config.useTx * config.useRx;

                % At the moment only the first antenna combination can be
                % used.
                % Select attenuation setting
                if config.useAttenuation > 1
                    vdats = fmcw_split_all(vdat);
                    vdat = vdats(config.useAttenuation);
                end

                % Handle multiple bursts:
                if config.lowestNoise
                    % Select burst with lowest noise floor.
                    [~,~,~,spec] = fmcw_range(vdat,config.padding,4*config.maxDepth,@blackman);
                    [vdat,selBurst] = fmcw_burst_select(vdat,spec);
                elseif config.burstMean
                    % Compute the mean of the burst.
                    vdat = fmcw_burst_mean(vdat);
                else
                    % Alternatively, use first burst
                    selBurst = 1;
                    vdat.vif =  vdat.vif(selBurst,:);
                    vdat.Startind =  vdat.Startind(selBurst);
                    vdat.Endind =  vdat.Endind(selBurst);
                    vdat.chirpNum =  vdat.chirpNum(selBurst);
                    vdat.chirpAtt =  vdat.chirpAtt(selBurst);
                    vdat.chirpTime = vdat.chirpTime(selBurst);
                end

                % Range compresion of FMCW data (following Brennan et al., 2014)
                [range,~,specCor,spec] = fmcw_range(vdat,config.padding,config.maxDepth,@blackman);

                % Initalize data structure for survey.
                if ii == 1
                    obj.range = range';
                    obj.twtt = obj.range / (physconst('lightspeed') / sqrt(vdat.er) / 2);
                    obj.dtBin = obj.twtt(2) - obj.twtt(1);

                    % Preset survey variables
                    obj.spec = NaN(length(range),numFiles);
                    obj.specCor = NaN(length(range),numFiles);
                    obj.nSubBursts = NaN(1,numFiles);
                    obj.nBurstSel = NaN(1,numFiles);
                    obj.nAttenuators = NaN(size(obj.nSubBursts));
                    obj.attenuator1 = NaN(1,numFiles); %NaN(4,numFiles);
                    obj.afGain = NaN(1,numFiles); %NaN(4,numFiles);
                    obj.txAnt = NaN(1,numFiles); %NaN(8,numFiles);
                    obj.rxAnt = NaN(1,numFiles); %NaN(8,numFiles);
                    obj.temp1 = NaN(size(obj.nSubBursts));
                    obj.temp2 = NaN(size(obj.nSubBursts));
                    obj.batteryVoltage = NaN(size(obj.nSubBursts));
                    obj.time = NaN(size(obj.nSubBursts));
                    obj.timestamp = NaT(size(obj.nSubBursts));
                    obj.T = vdat.T;
                    obj.f0 = vdat.f0;
                    obj.f1 = vdat.f1;
                    obj.fc = vdat.fc;
                    obj.B = vdat.B;
                    obj.K = vdat.K;
                    obj.dt = vdat.dt;
                    obj.permittivity = vdat.er;
                    obj.ci = vdat.ci;
                    obj.lambdaC = vdat.lambdac;
                    obj.numFiles = numFiles;
                    obj.fileName = strings(1, numFiles);
                end

                % Fill in trace data
                obj.spec(:,ii) = spec;
                obj.specCor(:,ii) = specCor;
                obj.nSubBursts(ii) = vdat.SubBurstsInBurst;
                obj.nBurstSel(ii) = selBurst;
                obj.nAttenuators(ii) = vdat.NAttenuators;
                obj.attenuator1(ii) = vdat.Attenuator_1;
                obj.afGain(ii) = vdat.Attenuator_2;
                obj.txAnt(ii) = vdat.TxAnt;
                obj.rxAnt(ii) = vdat.RxAnt;
                obj.temp1(ii) = vdat.Temperature_1;
                obj.temp2(ii) = vdat.Temperature_2;
                obj.batteryVoltage(ii) = vdat.BatteryVoltage;
                obj.time(ii) = vdat.chirpTime;
                obj.timestamp(ii) = datetime(obj.time(ii), 'ConvertFrom', 'datenum');
                obj.fileName(ii) = rawFiles(iiSelInd).name;

                % update waitbar
                waitbar(ii/(numFiles-1.0),wb);
            end
            % Close waitbar
            close(wb);

            % Add optional fields from config
            fields = {'cableTx', 'cableRx', 'distAntenna', 'padding'};
            for i = 1:length(fields)
                if isprop(config, fields{i})
                    obj.(fields{i}) = config.(fields{i});
                end
            end

            obj.config = config;

            if config.useGPS
                obj = loadGPS(obj);
            end
        end


        function obj = loadGPS(obj)
            % LOADGPS loads GPS data.
            %
            %   Options:
            %   csv: comma-seperated csv file with the fields
            %       filename, latitude, longitude and elevation
            
            switch lower(obj.config.typeGPS) % Make typeGPS case insensitive.
                case 'csv'
                    % Open GPS table
                    gpsData = readtable(obj.config.fileGPS,Delimiter=',');
                    
                    % Extract GPS data
                    filenames = gpsData.('filename');
                    latitude = gpsData.('latitude');
                    longitude = gpsData.('longitude');
                    elevation = gpsData.('elevation');

                    % Add .dat to file names if not existing
                    for indFile = 1:numel(filenames)
                        if ~contains(filenames{indFile},'.dat')
                            filenames{indFile} = [filenames{indFile} '.dat'];
                        end
                    end

                    % Assign GPS data by matching with filename
                    for indObj = 1:obj.numFiles
                        indGPS = find(obj.fileName(indObj)==filenames);
                        if isempty(indGPS)
                            error(['Survey/loadGPS: filename "' char(obj.fileName(indObj)) '" is missing in GPS file.'])
                        end
                        obj.lat(indObj) = latitude(indGPS);
                        obj.lon(indObj) = longitude(indGPS);
                        obj.elev(indObj) = elevation(indGPS);
                    end
                otherwise
                    error('Invalid GPS input type.');
            end
            
            % Turn lat/lon into x/y coordinates.
            proj = projcrs(obj.config.surveyEPSG);
            [obj.x, obj.y] = projfwd(proj, obj.lat, obj.lon);
            obj.pos = [obj.x; obj.y; obj.elev];

            % Compute distance array from distances between GPS points.
            if numel(obj.x)>2
                distgps = cumsum(sqrt(diff(obj.x).^2 + diff(obj.y).^2));
                obj.distGPS = [0 distgps];
            end
        end

        
        function obj = processSurvey(obj,config)
            %PROCESSSURVEY runs all general pRES survey processing steps
            %   Process survey runs the processing steps that are
            %   potentially applicable to all survey types.
            
            % Option to crop data for cable length
            if config.cropCables
                obj = obj.cropCables();
            end

            % Option to correct for firn density
            if config.correctDensity
                obj = obj.correctDensity();
            else
                obj.rangeCor = obj.range;
            end

            % Option to apply a shape matching filter (following Rahman et al., 2014)
            if config.matchShape
                obj = obj.matchShape();
            end
        end

        
        function obj = cropCables(obj)
            %CROP_CABLES crops the bins that correspond to the cables.
            %   The first range bins of pRES data only contain noise that
            %   is recorded before the signals can have travelled through
            %   the cables. This function cuts these range bins.
            
            % Find airwave at first amplitude maximum
            [~,airwaveInd] = findpeaks(mean(abs(obj.spec),2),'Npeaks',20,'SortStr','descend');
            airwaveInd = min(airwaveInd);
            
            % Compute the velocity of the wave in the cable
            vfCable = (obj.cableTx + obj.cableRx)/(airwaveInd * obj.dtBin * physconst('lightspeed') - obj.distAntenna);
            
            % Compute the traveltime in the cable
            cableTime = (obj.cableTx + obj.cableRx) / (physconst('lightspeed') * vfCable);
            
            % Compute the bins that are recorded while the signal is
            % traveling through the cable.
            cableBin = round(cableTime ./ obj.dtBin);
            
            flag = {['cropCables: cropped first ' num2str(cableBin) ' bins.']};
            obj.flags = [obj.flags flag];
            disp([flag{1} newline])
            
            % Crop data (i.e. to remove cable length)
            obj.range = obj.range(1:end-cableBin+1);
            obj.twtt = obj.twtt(1:end-cableBin+1);
            obj.spec = obj.spec(cableBin:end,:);
            obj.specCor = obj.specCor(cableBin:end,:);
            
            % Correct for phase of surface wave to be 0.
            obj.spec = obj.spec .*  repmat(exp(-1i*(angle(obj.spec(1,:)))), size(obj.spec,1), 1);
            obj.specCor = obj.specCor .*  repmat(exp(-1i*(angle(obj.specCor(1,:)))), size(obj.specCor,1), 1);
        end


        function obj = correctDensity(obj)
            %DENSITY_CORRECTION computes the density corrected range of
            %each bin.
            %   This method applies a density correction based on the
            %   density data given in config.fileDensity. This data needs
            %   to be provided as a table with the entries 'range' and
            %   'density'.
            %
            %   The script matches the two-way traveltimes with the
            %   expected two-way traveltimes for a given depth as given by
            %   the density data.

            % Load density dada
            densityData = readtable(obj.config.fileDensity);
            rangeData = densityData.('range');
            density = densityData.('density');
            
            % Compute permittivity
            permittivityIce = 3.17;
            densityIce = 917;
            permittivityCore = (density/densityIce * (permittivityIce^(1/3) - 1) + 1).^3;
            
            % Compute twtt.
            cfirn = physconst('lightspeed')./sqrt(permittivityCore); % Compute wave velocity in firn.
            t2firn = 2*diff(rangeData)./cfirn(1:end-1); % Compute twtt per bin.
            t2firn = [0; cumsum(t2firn)]; % Compute full twtt.
            
            % Derive corresponding depth from twtt of ApRES data.
            range_apres = zeros(size(obj.twtt));
            for n = 1:numel(obj.twtt)
                % Find coordinate of next lowest bin in t2firn.
                nt2 = find(t2firn <= obj.twtt(n), 1, 'last');
                % Compute additional time after nt2.
                t2diff = obj.twtt(n) - t2firn(nt2); 
                % Compute additional range after nt2 by wave velocity in the layer.
                rdiff = t2diff * cfirn(nt2)/2; 
                % Add up the range of nt2 with the additional component.
                range_apres(n) = rangeData(nt2) + rdiff;
            end
            obj.rangeCor = range_apres;

            flag = {['correctDensity: twtt was corrected for the firn density using the file: ', obj.config.fileDensity]};
            obj.flags = [obj.flags flag];
            disp([flag{1} newline])
        end


        function obj = matchShape(obj)
            %MATCHSHAPE Matches the recorded with expected shape repsonses
            %of reflections.
            %   This script convolutes the recorded data with a
            %   filter_shape given by the expected response of specular
            %   reflections to improve their detectability.
            %
            %   This method is implemented following:
            %   Rahman et al. (2014), Shape matching algorithm for
            %   detecting radar hidden targets to determine internal layers
            %   in Antarctic ice shelves, 2014 IEEE Radar Conference
            %   (RadarCon).

            % Create the filter
            % (Defined based on typical signal responses).
            win_len = 4 * obj.padding; % Padding dependened window length of convolution
            filter_range = (-win_len:win_len) * obj.dtBin*1e9;
            filter_shape = normpdf(filter_range,0,4.75).*obj.dtBin*1e9;

            size_spec = size(obj.spec);
            spec_abs = abs(obj.spec);
            
            thr = 0.05;
            
            sm = zeros(size_spec);
            
            for n = 1:size_spec(2)
                for m = 1:size_spec(1)
                    fs_start = max(1, win_len + 2 - m);
                    fs_end = min(numel(filter_shape), win_len + 1 + size_spec(1) - m);
                    fs = filter_shape(fs_start:fs_end);
                    
                    sn = spec_abs(:,n)';
                    gm_start = max(1, m - win_len);
                    gm_end = min(numel(sn), m + win_len);
                    gm = sn(gm_start:gm_end);
                    
                    sc = max(fs)./max(gm);
                    ssc = gm*sc;
                    
                    rank = 1 - abs(ssc - fs)./thr;
                    rank = (rank + abs(rank))./2;
                    rank_cum = mean(rank);

                    sm(m,n) = rank_cum;
                end
            end
            
            % Create shape matched signal by combining with phase from original data
            obj.spec = sm .* obj.spec./spec_abs;
            obj.specCor = sm .* obj.specCor./spec_abs;
            
            flag = {'matchShape: filtered signal by shape matching with typical signal response.'};
            obj.flags = [obj.flags flag];
            disp([flag{1} newline])
        end
        

        function obj = loadSurvey(obj_org,filename)
            % LOADSURVEY loads a survey object .mat-file
            pRESdata = load(filename,'-mat');
            obj = pRESdata.pRESdata;

            obj.config = obj_org.config;
        end


        function saveSurvey(obj, filename)
            % LOADSURVEY saves a survey object .mat-file after excluding
            % local folder structure settings.

            % Create structure to save
            pRESdata = obj;

            % Exclude local folder structure
            cs = obj.config;
            [cs.dirProject, cs.dirRaw, cs.dirPreProc, cs.dirProc, cs.dirGPS, cs.dirSup, cs.filePreProc, cs.fileProc, cs.filePos, cs.fileGPS, cs.fileDensity] = deal([]);
            pRESdata.config = cs;

            % Save file.
            save(filename,'pRESdata');
        end
    end
end

