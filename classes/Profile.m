classdef Profile < Survey
    %PROFILE is the survey subclass for profile surveys.
    %   This survey subclass includes processing of profile surveys that
    %   were generated using a mobile pRES setup.
    %
    %   For details, see:
    %   Oraschewski et al. (2024) Layer-optimized SAR processing with a
    %   mobile mobile phase-sensitive radar: a proof of concept for
    %   detecting the deep englacial stratigraphy of Colle Gnifetti,
    %   Switzerland/Italy. The Cryosphere.

    properties
        % Profile domain and positioning
        profilePos double
        profileDist
        profileSize double
        profileLine double
        profileSpacing double
        distFit double
        dist
        lineNumber
        txPos double
        rxPos double

        % Data variables
        imgProc double
        imgProcIncoherent double
        imgSlope double
        imgSlopeFiltered double
    end


    methods
        function obj = Profile(config)
            % PROFILE Constructor for Profile class which is a subclass
            % of the Survey class
            %   This method creates the Profile object, by taking the base
            %   class properties from the Survey class. This gives the
            %   pre-processed ApRES data for the profile.
            
            % Get base class properties
            obj = obj@Survey(config);% Get base class properties
            % Include lineNumber property
            obj.lineNumber = ones([1, obj.numFiles]);
        end
        
        
        function obj = processProfile(obj,config)
            %PROCESSPROFILE runs the profile processing as set in the
            %configuration
            %   This method executes the processing of the Profile object
            %   as set in the configuration. First, the general survey
            %   processing steps are executed which are potentially
            %   applicable for all survey types. Afterwards, the processing
            %   specifically for profile objects is executed.
            
            % General survey processing
            obj = obj.processSurvey(config);

            % Estimate antenna position (Currently, this has no effect)
            if config.useAntennaPos
                obj = obj.estimateAntennaPos();
            end

            % Obtain lines as defined in config.
            if config.defineLines
                obj = obj.obtainLines();
            end

            % (Re-)Compute profile position (if not deactivated)
            if config.newPos || ~isfile(config.filePos)
                obj = obj.positionProfile();
                PosData.distFit = obj.distFit;
                PosData.profilePos = obj.profilePos;
                PosData.profileLine = obj.profileLine;
                PosData.profileSize = obj.profileSize;
                PosData.profileDist = obj.profileDist;

                save(config.filePos, '-struct', 'PosData');
            else
                % Or load existing position data
                PosData = load(config.filePos);
                obj.distFit = PosData.distFit;
                obj.profilePos = PosData.profilePos;
                obj.profileLine = PosData.profileLine;
                obj.profileSize = PosData.profileSize;
                obj.profileDist = PosData.profileDist;
            end

            % (Re-)run processing
            if config.newProc || ~isfile(config.fileProc)
                obj = obj.processSAR();
                obj.saveSurvey(config.fileProc)
            else
                % Or load existing processed data
                obj = obj.loadSurvey(config.fileProc);
            end
        end


        function obj = obtainLines(obj)
            %OBTAINLINES obtains linenumbers as defined in config file.
            %
            % This function allows to handle several lines in the radar
            % data seperately

            % Number of lines
            nLines = numel(obj.config.numLines);

            % Indices associated with lines
            inds = [obj.config.indLines{:}];
            obj.numFiles = numel(inds);
            
            % Cut data outside lines
            vars = {'spec', 'specCor', 'lat', 'lon', 'x', 'y', 'elev', ...
                'pos', 'distGPS', 'nSubBursts', 'nAttenuators', ...
                'attenuator1', 'afGain', 'txAnt', 'rxAnt', 'temp1', ...
                'temp2', 'batteryVoltage', 'time', 'timestamp', ...
                'partNumber', 'partName', 'lineNumber', 'fileName', ...
                'txPos', 'rxPos', 'profileSpacing', 'dist'};
            for m = 1:numel(vars)
                var = vars{m};
                try %#ok<TRYNC> 
                    obj.(var) = obj.(var)(:,inds);
                end
            end

            % Assign line numbers
            trace_start = 1;
            for n = 1:nLines
                trace_end = trace_start + numel(obj.config.indLines{n}) - 1;
                obj.lineNumber(trace_start:trace_end) = obj.config.numLines{n};
                trace_start = trace_end + 1;
            end
        end


        function obj = positionProfile(obj)
            %POSITIONPROFILE Compute position for processed profile
                        
            % Update Coordinate system (in case the EPSG was changed in
            % config without reprocessing the positioning)
            proj = projcrs(obj.config.surveyEPSG);
            [obj.x, obj.y] = projfwd(proj, obj.lat, obj.lon);
            obj.pos = [obj.x; obj.y; obj.elev];

            % Compute position
            switch lower(obj.config.methodPos) % Make survey types case insensitive.
                case 'org'
                    obj.profilePos = obj.pos;
                    obj.profileLine = obj.lineNumber;
                    obj.profileSize = [length(obj.range), length(obj.lineNumber)]; % ToDo: Split this into different line components.
                case 'smooth'
                    datKeys = unique(obj.lineNumber);
                    nLines = numel(datKeys);
                    
                    for n = 1:nLines        
                        indOrg = (obj.lineNumber == datKeys(n));
                        lnPosOrg = obj.pos(:,indOrg);

                        % Obtain differences in x and y direction
                        xDiff = lnPosOrg(1,end) - lnPosOrg(1,1);
                        yDiff = lnPosOrg(2,end) - lnPosOrg(2,1);
                        
                        % Access orientation of section
                        [dir1, dir2] = deal(1); % Define main and secondary direction of line, where its potential values 1 and 2 refer to x and y direction.
                        if abs(xDiff) > abs(yDiff)
                            dir2 = 2;
                            %dirSign = sign(xDiff);
                        else
                            dir1 = 2;
                            %dirSign = sign(yDiff);
                        end
                        
                        % Run spline interpolation
                        f = fit(lnPosOrg(dir1,:)',lnPosOrg(dir2,:)','smoothingspline','SmoothingParam',obj.config.posSmoothing);
                        
                        % Create fine spaced array along the interpolated
                        % smoothed profile
                        dir1SF = linspace(lnPosOrg(dir1,1), lnPosOrg(dir1,end), 100*length(lnPosOrg(1,:)));
                        dir2SF = f(dir1SF)';

                        % Compute full length of the interpolated profile
                        dir1SFdiff = diff(dir1SF);
                        dir2SFdiff = diff(dir2SF);
                        diffSF = [0 sqrt(dir1SFdiff.^2 + dir2SFdiff.^2)];
                        distSF = cumsum(diffSF);

                        % Create equally spaced array with profile length
                        obj.profileDist = 0:obj.config.profileSpacing:distSF(end);

                        % Find x and y coordinates of equally spaced points
                        dir1Pos = interp1(distSF, dir1SF, obj.profileDist);
                        dir2Pos = interp1(distSF, dir2SF, obj.profileDist);
                        
                        % Assign dir1/dir2 to corresponding coordinate.
                        if dir1 == 1
                            x = dir1Pos;
                            y = dir2Pos;
                        else
                            x = dir2Pos;
                            y = dir1Pos; 
                        end

                        % Find distance of trace points on the smooth line.
                        [~,indSF] = min(sqrt((dir1SF - lnPosOrg(dir1,:)').^2 + (dir2SF - lnPosOrg(dir2,:)').^2),[],2);
                        [indSF, indindSF] = unique(indSF);
                        obj.distFit = distSF(indSF);

                        % Obtain z positions.
                        z = interp1(obj.distFit, lnPosOrg(3,indindSF), obj.profileDist,'linear');
                        
                        profilePosNew = [x; y; z];
                        if n==1
                            obj.profilePos = profilePosNew;
                            obj.profileLine = ones(1, length(profilePosNew)) * datKeys(n);
                            obj.profileSize = [length(obj.range), length(profilePosNew)];
                        else
                            obj.profilePos = [obj.profilePos, profilePosNew];
                            obj.profileLine = [obj.profileLine, ones(1, length(profilePosNew)) * datKeys(n)];
                            obj.profileSize = [obj.profileSize; [length(obj.range), length(profilePosNew)]];
                        end
                    end
                otherwise
                    error(['positionProfile: positioning method "' obj.config.methodPos '" is unknown.']);
            end

            flag = {['positionProfile: positioning method "' obj.config.methodPos '".']};
            obj.flags = [obj.flags flag];
            disp([flag{1} newline])
        end


        function obj = processSAR(obj)
            %PROCESSSAR Execute the selected SAR processing method
            switch lower(obj.config.methodProc) % Make survey types case insensitive.
                case {'losar', 'losar_castelletti'}
                    obj = obj.losar();
                % case 'backprojection' % Excluded for v1.0, testing
                % required.
                %    obj = obj.backprojection();
                case {'interpolation', 'interp'}
                    obj = obj.interpolation();
                case 'movmean'
                    obj = obj.movmean();
                otherwise
                    error(['processSAR: SAR processing method "' obj.config.methodProc '" is unknown.']);
            end
        end


        function obj = losar(obj)
            %LOSAR Executes Layer-optimized SAR processing
            flag = {'Profile processing method: layer-optimized SAR'};
            obj.flags = [obj.flags flag];
            disp([flag{1} newline])
            
            % Create empty matrizies for processed images
            obj.imgProc = zeros([obj.profileSize(1,1) sum(obj.profileSize(:,2))]);
            obj.imgProcIncoherent = zeros(size(obj.imgProc));
            obj.imgSlope = zeros(size(obj.imgProc));
            obj.imgSlopeFiltered = zeros(size(obj.imgProc));

            % Find number of individual lines. 
            datKeys = unique(obj.lineNumber);
            sarKeys = unique(obj.profileLine);
            nLines = numel(datKeys);

            % Iterate over individual lines.
            for k = 1:nLines
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%% Find layer slopes from phase %%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

                if obj.config.newSlopes || ~isfile(obj.config.fileProc)
                    % Obtain slopes using one of the implemented methods
                    switch lower(obj.config.methodSlope)
                        case 'compute'
                            obj = obj.slopeCompute(k);
                        case 'linefit'
                            obj = obj.slopeLineFitting(k);
                        %case 'radon' Excluded for v1.0, testing required.
                        otherwise
                            error(['losar: Slope finding method "' obj.config.methodSlope '" is unknown.']);
                    end
                    % Save survey with updated slopes
                    obj.saveSurvey(obj.config.fileProc)
                else
                    % Load survey with previously computed slopes
                    obj = obj.loadSurvey(obj.config.fileProc);
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% LoSAR processing: Sum coherently along slope layers %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                % Obtain indizes of lines in original and processed image
                indOrg = (obj.lineNumber == datKeys(k)); % Indizes of line in orignal radargram
                lnPosOrg = obj.pos(:,indOrg);

                indSAR = (obj.profileLine == sarKeys(k)); % Indizes of line in processed radargram
                lnPosSAR = obj.profilePos(:,indSAR);

                % Load original data
                imgOrg    = obj.specCor(:,indOrg); % Corrected data

                % Load derived layer slopes for current line
                if obj.config.filterSlopes
                    imgSlopeSAR = obj.imgSlopeFiltered(:,indSAR);
                else
                    imgSlopeSAR = obj.imgSlope(:,indSAR);
                end

                % Create data arrays for SAR processed image
                imgSAR = zeros(obj.profileSize(k,:));
                imgSARIncoherent = zeros(obj.profileSize(k,:));

                % Process profile point by point
                nRange = 1:obj.profileSize(k,1);
                mRange = 1:obj.profileSize(k,2);

                for m = mRange
                    % Find closest points in SAR range
                    mDistAll = sqrt((lnPosOrg(1,:) - lnPosSAR(1,m)).^2 + (lnPosOrg(2,:) - lnPosSAR(2,m)).^2); % Distance of SAR point to all traces
                    indOrgInRange = find(mDistAll < obj.config.lengthSAR/2); % Indizes of datapoints in SAR range
                    distInRange = mDistAll(indOrgInRange); % Reduce dist to points in range
                    posInRange = lnPosOrg(:,indOrgInRange); % Reduce position vector to points in range
                        
                    % Check if SAR point is before or after closest trace
                    [~, indMinDist] = min(distInRange);
                        
                    if m>1 && m< obj.profileSize(k,2)
                        prev_dist = sqrt((posInRange(1,indMinDist) - lnPosSAR(1,m-1)).^2 + (posInRange(2,indMinDist) - lnPosSAR(2,m-1)).^2);
                        next_dist = sqrt((posInRange(1,indMinDist) - lnPosSAR(1,m+1)).^2 + (posInRange(2,indMinDist) - lnPosSAR(2,m+1)).^2);
                        
                        % Check if data point is closer to previous or next point on SAR array.
                        if prev_dist < next_dist
                            distInRange(1:indMinDist) = -distInRange(1:indMinDist);
                        else
                            distInRange(1:indMinDist-1) = -distInRange(1:indMinDist-1);
                        end
                    elseif m == obj.profileSize(k,2)
                        % At end of line, make distInRange negative.
                        distInRange(1:indMinDist) = -distInRange(1:indMinDist);
                        % At beginning do nothing.
                    end
    
                    % Sort reference points according to their distance
                    [~, idxSort] = sort(distInRange);
                    distInRange = distInRange(idxSort);
                    indOrgInRange = indOrgInRange(idxSort);
    
                    if strcmpi(obj.config.methodProc, 'losar')
                        for n = nRange
                            % Perform operations on smaller window to increase performance.
                            isRange = [max([1,n-40]), min([obj.profileSize(k,1),n+40])]; % Vertical range for small image
                            imgSmall = imgOrg(isRange(1):isRange(2), indOrgInRange);
                            imgSmallTraces = 1:size(imgSmall,2);
        
                            slopeFilt =  imgSlopeSAR(n,m);
                            % drn = dr(n);
                            drn = median(diff(obj.range));
                            y = n - distInRange * tand(slopeFilt)/drn - isRange(1) + 1;        
                            
                            yMax = isRange(2) - isRange(1) + 1;
        
                            % Obtain floor and ceil bin values and ratios
                            % for its dominance
                            [yFl,yCe,ratioFl,ratioCe,yFlInRange,yCeInRange] = getFloorCeil(y, yMax);
                            
                            % Obtain phase values at these bins
                            yFlInd = sub2ind(size(imgSmall),yFl,imgSmallTraces(yFlInRange));
                            yFlPower = imgSmall(yFlInd);
                            yCeInd = sub2ind(size(imgSmall),yCe,imgSmallTraces(yCeInRange));
                            yCePower = imgSmall(yCeInd);
        
                            % Combine phase points of yFl and yCe in one
                            % array and scale for the dominance ratios.
                            PowerAll = [ratioFl .* yFlPower, ratioCe .* yCePower];
        
                            % Sum in the phase space to find maximal
                            % coherence in the phase.
                            imgPower = sum(PowerAll);
                            imgPowerIncoherent = sum(abs(PowerAll));
                            
                            % Save number of considered points to scale
                            % results.
                            imgNumPoints = sum([ratioFl, ratioCe]);
        
                            imgSAR(n,m) = imgPower./imgNumPoints;
                            imgSARIncoherent(n,m) = imgPowerIncoherent./imgNumPoints;
                        end
                    elseif strcmpi(obj.config.methodProc, 'losar_castelletti')
                        % This method approximates the LO-SAR
                        % implementation as presented by
                        % Castelletti et al. (2019) Layer optimized SAR
                        % processing and slope estimation in radar sounder
                        % data. J. Glaciol. 65, 983–988.
                        % https://doi.org/10.1017/jog.2019.72
                        %
                        % This method was only implemented for comparision
                        % and does not recompute the slopes by maximizing
                        % the SNR after the phase correction, but takes the
                        % slopes estimated with our LO-SAR approach and
                        % performs the coherent summation along range bins.

                        % Extract data in range                        
                        BinInRange = imgOrg(:, indOrgInRange);
    
                        % Load slope
                        slopeFilt =  imgSlopeSAR(:,m);
                        
                        % Compute phase correction
                        dPhi = tand(slopeFilt) * distInRange * 4 * pi / obj.lambdaC;
                        
                        % Apply phase correction
                        BinInRange = BinInRange .* exp(1i*(dPhi));
    
                        % Sum power over each range bin
                        sumPower = sum(BinInRange, 2);
                        sumPowerIncoherent = sum(abs(BinInRange),2);
                        
                        % Find number of considered traces to scale results
                        imgNumPoints = numel(indOrgInRange);
    
                        % Scale power and assign to SAR trace
                        imgSAR(:,m) = sumPower./imgNumPoints;
                        imgSARIncoherent(:,m) = sumPowerIncoherent./imgNumPoints;
                    end
                    disp(['Layer-optimized SAR processing, trace: ' num2str(m) '/' num2str(numel(mRange))]);
                end
                obj.imgProc(:,indSAR) = imgSAR;
                obj.imgProcIncoherent(:,indSAR) = imgSARIncoherent;
            end
        end


        function obj = interpolation(obj)
            %INTERPOLATION assigns to each profile position the data value of
            %the nearest recorded data point.
            %
            % Falk Oraschewski, 19.10.2022
            
            disp('Profile processing method: Interpolation')
    
            % Find number of individual lines. 
            datKeys = unique(obj.lineNumber);
            sarKeys = unique(obj.profileLine);
            nLines = numel(datKeys);
    
            imgOrg = abs(obj.specCor);
            obj.imgProc = zeros([obj.profileSize(1,1) sum(obj.profileSize(:,2))]);
                
            traceLen = obj.profileSize(1,1);
                
            for k = 1:nLines
                indOrg = (obj.lineNumber == datKeys(k)); % Indizes of line in orignal radargram
                indSAR = (obj.profileLine == sarKeys(k)); % Indizes of line in sar processed radargram
                    
                linePosOrg = obj.pos(:,indOrg);
                lineImgOrg = imgOrg(:,indOrg);
                
                linePosSAR = obj.profilePos(:,indSAR);
                lineLenSAR = length(linePosSAR);
                lineImgSAR = zeros(traceLen, lineLenSAR);
    
                % Find main direction of SAR profile
                xDiff = linePosSAR(1,end) - linePosSAR(1,1);
                yDiff = linePosSAR(2,end) - linePosSAR(2,1);
                
                if abs(xDiff) > abs(yDiff)
                    [lineXorg, xind, ~] = unique(linePosOrg(1,:),'stable');
                    lineImgOrg = lineImgOrg(:,xind);
                    
                    for kk = 1:traceLen
                        lineImgSAR(kk,:) = interp1(lineXorg, lineImgOrg(kk,:), linePosSAR(1,:),'nearest');
                    end
                else
                    [lineYorg, yind, ~] = unique(linePosOrg(2,:),'stable');
                    lineImgOrg = lineImgOrg(:,yind);
                    
                    for kk = 1:traceLen
                        lineImgSAR(kk,:) = interp1(lineYorg, lineImgOrg(kk,:), linePosSAR(2,:),'nearest');
                    end
                end
                
                obj.imgProc(:,indSAR) = lineImgSAR;
            
            end
        end


        function obj = movmean(obj)
            %MOVMEAN computes the moving mean along the pRES profile.
            %   This function computes the moving mean along a mobile pRES
            %   profile over the SAR length. Essentially, this corresponds
            %   to an unfocused SAR computation.
            %
            % Falk Oraschewski, 11.03.2022
            
            disp('Profile processing method: Moving mean')

            % Create empty matrizies for processed images
            obj.imgProc = zeros([obj.profileSize(1,1) sum(obj.profileSize(:,2))]);
            obj.imgProcIncoherent = zeros(size(obj.imgProc));

            % Find number of individual lines. 
            datKeys = unique(obj.lineNumber);
            sarKeys = unique(obj.profileLine);
            nLines = numel(datKeys);

            % Iterate over individual lines.
            for k = 1:nLines
                % Obtain indizes of lines in original and processed image
                indOrg = (obj.lineNumber == datKeys(k)); % Indizes of line in orignal radargram
                lnPosOrg = obj.pos(:,indOrg);

                indSAR = (obj.profileLine == sarKeys(k)); % Indizes of line in processed radargram
                lnPosSAR = obj.profilePos(:,indSAR);

                lineLenSAR = obj.profileSize(k,2);

                lineImg = zeros(obj.profileSize(1,1), lineLenSAR);
                lineImgIncoherent = zeros(obj.profileSize(1,1), lineLenSAR);

                mRange = 1:lineLenSAR;
                for m = mRange
                    % Find points in SAR range
                    mDistAll = sqrt((lnPosOrg(1,:) - lnPosSAR(1,m)).^2 + (lnPosOrg(2,:) - lnPosSAR(2,m)).^2); % Distance of SAR point to all traces
                    indOrgInRange = find(mDistAll < obj.config.lengthSAR/2); % Indizes of datapoints in SAR range

                    %  Compute moving mean of the points in range
                    lineImg(:,m) = abs(mean(obj.specCor(:,indOrgInRange), 2)); % Coherent summation
                    lineImgIncoherent(:,m) = mean(abs(obj.specCor(:,indOrgInRange)), 2); % Incoherent summation
                end
                % Assign line images to full output matrix
                obj.imgProc(:,indSAR) = lineImg;
                obj.imgProcIncoherent(:,indSAR) = lineImgIncoherent;
            end
        end


        function obj = slopeCompute(obj, numLine)
            %SLOPECOMPUTE compute slopes from along-track phase differences
            %   Detailed explanation goes here
            %
            %   26.05.2023, Falk Oraschewski

            % Find number of individual lines. 
            datKeys = unique(obj.lineNumber);
            sarKeys = unique(obj.profileLine);

            % Obtain indizes of lines in original and processed image
            indOrg = (obj.lineNumber == datKeys(numLine)); % Indizes of line in orignal radargram
            lnPosOrg = obj.pos(:,indOrg);

            indSAR = (obj.profileLine == sarKeys(numLine)); % Indizes of line in processed radargram
            lnPosSAR = obj.profilePos(:,indSAR);

            % Load original data
            imgPhiRaw = angle(obj.spec(:,indOrg)); % Raw phase

            % Compute phase differences
            imgPhiRawDiff = diff(imgPhiRaw,1,2);
            % Correct for jumps between +pi and -pi
            imgPhiRawDiff(imgPhiRawDiff > pi) = imgPhiRawDiff(imgPhiRawDiff > pi) - 2 * pi;
            imgPhiRawDiff(imgPhiRawDiff < -pi) = imgPhiRawDiff(imgPhiRawDiff < -pi) + 2 * pi;

            distMid = obj.distFit(1:end-1) + diff(obj.distFit)/2;
            imgPhiRawDiffFilt = imgPhiRawDiff;
            for n = 1:obj.profileSize(numLine,1)
                f = fit(distMid',imgPhiRawDiff(n,:)','smoothingspline','SmoothingParam',0.9);
                imgPhiRawDiffFilt(n,:) = f(distMid);

%                 imagesc(imgPhiRawDiffFilt)
%                 colormap(redblue)
%                 clim([-pi pi])
%                 drawnow;
            end

            % Compute change of reflector depths
            imgDeltaR = imgPhiRawDiff * obj.lambdaC/(4 * pi);
            imgDeltaRFilt = imgPhiRawDiffFilt * obj.lambdaC/(4 * pi);
            
            % Compute point distance parallel to path:
            deltaStep = diff(obj.distGPS);
            deltaStepMode = mode(round(deltaStep,2)); % Most frequent step size
            deltaStep(deltaStep<0.5*deltaStepMode) = deltaStepMode;

            % Find orientation of smoothed radar line
            profileOrientation = atand(gradient(lnPosSAR(2,:))./gradient(lnPosSAR(1,:)));
            % Find orientation of steps between traces
            stepOrientation = atand(diff(lnPosOrg(2,:))./diff(lnPosOrg(1,:)));
            % Find closest point on smoothed line for every
            % trace step
            posStep = movmean(lnPosOrg,2,2);
            posStep = posStep(:,2:end);
            allDist = sqrt((posStep(1,:) - lnPosSAR(1,:)').^2 + (posStep(2,:) - lnPosSAR(2,:)').^2);
            [~, indClosest] = min(allDist);
            % Find orientation mismatch
            dOrientation = abs(stepOrientation - profileOrientation(indClosest));
            dOrientation(dOrientation>45) = 45;
            % Correct step length parallel to radar track
            deltaStep = deltaStep.*abs(cosd(dOrientation));

            % Compute layer slopes
            imgSlopeOrg = atand(-imgDeltaR./deltaStep);
            % Apply median filter to suppress extrema
            imgSlopeFilt = atand(-imgDeltaRFilt./deltaStep);
            imgSlopeFilt = medfilt2(imgSlopeFilt,[20,20]);


            % Interpolate slope onto SAR profile points
            distMid = obj.distFit + [0, diff(obj.distFit)/2];
            imgSlopeSAR = zeros([obj.profileSize(numLine,1), obj.profileSize(numLine,2)]);
            imgSlopeFiltSAR = zeros([obj.profileSize(numLine,1), obj.profileSize(numLine,2)]);
            imgSlopeOrg = [imgSlopeOrg(:,1), imgSlopeOrg];
            imgSlopeFilt = [imgSlopeFilt(:,1), imgSlopeFilt];
            for n = 1:obj.profileSize(numLine,1)
                imgSlopeFiltSAR(n,:) = interp1(distMid, imgSlopeFilt(n,:), obj.profileDist);
                imgSlopeSAR(n,:) = interp1(distMid, imgSlopeOrg(n,:), obj.profileDist);
            end
            obj.imgSlope(:,indSAR) = imgSlopeSAR;
            obj.imgSlopeFiltered(:,indSAR) = imgSlopeFiltSAR;
        end


        function obj = slopeLineFitting(obj, numLine)
            %SLOPELINEFITTING slopes from along-track phase differences
            %   Detailed explanation goes here
            %
            %   30.05.2023, Falk Oraschewski

            
            % Find number of individual lines. 
            datKeys = unique(obj.lineNumber);
            sarKeys = unique(obj.profileLine);

            % Obtain indizes of lines in original and processed image
            indOrg = (obj.lineNumber == datKeys(numLine)); % Indizes of line in orignal radargram
            lnPosOrg = obj.pos(:,indOrg);

            indSAR = (obj.profileLine == sarKeys(numLine)); % Indizes of line in processed radargram
            lnPosSAR = obj.profilePos(:,indSAR);
            
            % Define empty slope matrix.
            imgSlopeSAR = zeros([obj.profileSize(numLine,1), obj.profileSize(numLine,2)]);

            % Load phase corrected data
            imgOrg = obj.specCor(:,indOrg);
            imgPhiCor = imgOrg./abs(imgOrg); % Raw phase
                        
            % Obtain range differences
            dr = diff(obj.range);
            dr = [dr; dr(end)];

            % Find orientation of smoothed radar line
            profileOrientation = atand(gradient(lnPosSAR(2,:))./gradient(lnPosSAR(1,:)));

            % Find length of line section 
            lnLenSAR = length(lnPosSAR);
            traceLen = obj.profileSize(numLine,1);

            % Initiate counter for intermediate saving.
            nSave = 0;
            
            % Iterate over traces
            mRange = 1:lnLenSAR;
            for m = mRange
                tic
                % Find data points within SAR length
                allDist = sqrt((lnPosOrg(1,:) - lnPosSAR(1,m)).^2 + (lnPosOrg(2,:) - lnPosSAR(2,m)).^2);
                
                % Find points in range of the SAR length
                irIndOrg = find(allDist < obj.config.lengthSAR/2); % Indizes of datapoints in range
                irDistDir = allDist(irIndOrg); % Reduce dist to points in range
                irPos = lnPosOrg(:,irIndOrg); % Reduce position vector to points in range
                               
                % Find orientation mismatch
                irPointOrientation = atand((irPos(2,:) - lnPosSAR(2,m))./(irPos(1,:) - lnPosSAR(1,m)));     
                dOrientation = abs(irPointOrientation - profileOrientation(m));
                dOrientation(dOrientation>45) = 45;

                % Correct distance to point parallel to radar track
                irDist = irDistDir.*abs(cosd(dOrientation));

                % Check if sar point is before or after closest datapoint
                [~, indMinDist] = min(irDist);
                
                if m>1 && m<lnLenSAR
                    prev_dist = sqrt((irPos(1,indMinDist) - lnPosSAR(1,m-1)).^2 + (irPos(2,indMinDist) - lnPosSAR(2,m-1)).^2);
                    next_dist = sqrt((irPos(1,indMinDist) - lnPosSAR(1,m+1)).^2 + (irPos(2,indMinDist) - lnPosSAR(2,m+1)).^2);
                    
                    % Check if data point is closer to previous or next point on
                    % SAR array.
                    if prev_dist < next_dist
                        irDist(1:indMinDist) = -irDist(1:indMinDist);
                    else
                        irDist(1:indMinDist-1) = -irDist(1:indMinDist-1);
                    end
                elseif m == lnLenSAR
                    % At end of line, make ir_dist negative.
                    irDist(1:indMinDist) = -irDist(1:indMinDist);
                    % At beginning do nothing.
                end
                
                % Avoid points with same distance
                [~, idxUnique] = unique(irDist);
                irDist = irDist(idxUnique);
                irIndOrg = irIndOrg(idxUnique);

                % Sort reference points according to their distance
                [~, idxSort] = sort(irDist);
                irDist = irDist(idxSort);
                irIndOrg = irIndOrg(idxSort);
                
                testSlopes = -obj.config.slopeMax:obj.config.slopeInterval:obj.config.slopeMax;

                for n = 1:traceLen
                % parfor n = 1:traceLen % Use parallel computing to speed up code
                    imgPhase = zeros(size(testSlopes));
                    imgNumPoints = zeros(size(testSlopes));
                    
                    drn = dr(n);
            
                    % ToDo: Make is_range dependent on SAR width
                    % Perform operations on smaller window to increase performance.
                    isRange = [max([1,n-40]), min([traceLen,n+40])]; % Vertical range for small image
                    imgSmallPhase = imgPhiCor(isRange(1):isRange(2), irIndOrg);
                    
                    imgSmallTraces = 1:size(imgSmallPhase,2);
                    
                    yMax = isRange(2) - isRange(1) + 1;
                    for o = 1:length(testSlopes)
                        % Turn slope difference into bin value
                        y = n - irDist * tand(testSlopes(o))/drn - isRange(1) + 1;
                        
                        % Obtain floor and ceil bin values and ratios
                        % for its dominance
                        [yFl,yCe,ratioFl,ratioCe,yFlInRange,yCeInRange] = getFloorCeil(y, yMax);
                        
                        % Obtain phase values at these bins
                        yFlInd = sub2ind(size(imgSmallPhase),yFl,imgSmallTraces(yFlInRange));
                        yFlPhase = imgSmallPhase(yFlInd);
                        yCeInd = sub2ind(size(imgSmallPhase),yCe,imgSmallTraces(yCeInRange));
                        yCePhase = imgSmallPhase(yCeInd);

                        % Combine phase points of yFl and yCe in one
                        % array and scale for the dominance ratios.
                        PhaseCorAll = [ratioFl .* yFlPhase, ratioCe .* yCePhase];

                        % Sum in the phase space to find maximal
                        % coherence in the phase.
                        imgPhase(o) = sum(PhaseCorAll);
                        
                        % Save number of considered points to scale
                        % results.
                        imgNumPoints(o) = sum([ratioFl, ratioCe]);
                    end
                    % Find slope with maximal coherence
                    [~, ipha] = max(abs(imgPhase)./imgNumPoints);
                    slope = testSlopes(ipha);

                    imgSlopeSAR(n,m) = slope;
                end

                % Intermediate saving (manual interaction with code needed
                % to recover only partially estimated slopes).
                nSave = nSave + 1;                
                if nSave == 100
                    obj.imgSlope(:,indSAR) = imgSlopeSAR;
                    obj.saveSurvey(obj.config.fileProc);
                    nSave = 0;
                end
                disp(['Layer-optimized SAR processing, trace: ' num2str(m) '/' num2str(numel(mRange))]);
                toc
            end

            % Filter slopes using median filter
            imgSlopeFiltSAR = medfilt2(imgSlopeSAR, obj.config.slopeFilterWin);

            % Assign slope variables
            obj.imgSlope(:,indSAR) = imgSlopeSAR;
            obj.imgSlopeFiltered(:,indSAR) = imgSlopeFiltSAR;

            % Intermediate plots for testing
            % figure()
            % imagesc(obj.imgSlope)
            % colorbar()
            % colormap(redblue)
            % clim([-30 30])
            % 
            % figure()
            % imagesc(obj.imgSlopeFiltered)
            % colorbar()
            % colormap(redblue)
            % clim([-30 30])
        end


        function obj = estimateAntennaPos(obj)
            %ESTIMATEANTENNAPOS estimates the position of the antennas at
            %each trace on the profile.
            %   This method estimates the antenna positions for each trace
            %   on the profile, based on the positions of the previous and
            %   following trace locations.
            %   
            %   Note: Currently, the estimated antenna positions are not
            %   used further and the method is only included for potential
            %   future use.

            disp('Antenna positions are estimated. Currently this has no effect.')
            
            dx = diff(obj.x);
            dy = diff(obj.y);
            dz = diff(obj.elev);
            dthrsh = 10;
            antenna_averaging = 3;
            dBig = any([abs(dx) > dthrsh*mean(abs(dx)); abs(dy) > dthrsh*mean(abs(dy))]);
            dx(dBig) = NaN;
            dy(dBig) = NaN;
            dz(dBig) = NaN;    
            dx = movmean(dx, antenna_averaging,'omitnan');
            dy = movmean(dy, antenna_averaging,'omitnan');
            dz = movmean(dz, antenna_averaging,'omitnan');
            dx = [dx dx(end)];
            dy = [dy dy(end)];
            dz = [dz dz(end)];
            
            % Calculate normalised delta
            deltaAntenna = [dx; dy; dz];
            deltaAntenna = deltaAntenna ./ vecnorm(deltaAntenna);
        
            % Estimated Antenna positions
            obj.txPos = obj.pos + deltaAntenna * obj.config.distAntenna / 2;
            obj.rxPos = obj.pos - deltaAntenna * obj.config.distAntenna / 2;
        end
    end
end

