function site = fmcw_deformation(filename1,filename2,cfgi,fBurst,gBurst)

% site = fmcw_deformation(filename1,filename2,cfg)
% FMCW Radar Processor for vertical strain rate from a pair of bursts
%
% args:
% filename1: file 1 to process
% filename2: file 2 to process
% cfg = configuration settings structure (or only the fields you wish to
% change from defaults)
%
% Craig Stewart
% 2014

global cfg

%% Load processing settings
cfg = fmcw_process_config_deformation; % Load default configuration
% if nargin >2
%     % overwrite defaults
%     fieldnames = fields(cfgi);
%     fieldlist = [];
%     for ii = 1:length(fieldnames);
%         thisfield = fieldnames{ii};
%         cfg = setfield(cfg,thisfield,getfield(cfgi,thisfield));
%         fieldlist = [fieldlist ', ' thisfield];
%     end
%     fieldlist = fieldlist(3:end);
%     cfg.notes = [cfg.notes '. User modified fields: ' fieldlist '.'];
% end

%% Define Constants
daysPerYear = 365.25;

%% Select and load data
if nargin == 0
    if cfg.useTestFiles
        % test data
        filename1 = cfg.filename1;
        filename2 = cfg.filename2;
    else
        [file,path] = uigetfile({'*.dat;*.DAT;*.000;*.mat','Radar files: .dat, .DAT, .000 and .mat'},'Choose first (or both) radar files to calculate strain rate','multiselect','on');
        if isa(file,'double') % no files chosen
            return
        end
        if ischar(file) % then we've only selected one file as second is in a different dir. repeat
            filename1 = [path file];
            [file,path] = uigetfile({'*.dat;*.DAT;*.000;*.mat','Radar files: .dat, .DAT, .000 and .mat'},'Choose second radar files to calculate strain rate');
            filename2 = [path file];
        else
            filename1 = [path file{1}];
            filename2 = [path file{2}];
        end
    end
end
[path1,name1,ext1] = fileparts(filename1);
[path2,name2,ext2] = fileparts(filename2);
fBurst = 1;
gBurst = 1;
% Load
f = fmcw_load(filename1,fBurst); 
g = fmcw_load(filename2,gBurst);

%% Print Output
Disp(' ')
Disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
Disp(['Radar processor ' mfilename])
Disp(' ')
Disp('Using files:')
Disp('Filename                                  Date/Time')
Disp([name1 '   ' datestr(f.TimeStamp)])
Disp([name2 '   ' datestr(g.TimeStamp)])
dt = g.TimeStamp-f.TimeStamp; % time difference days
Disp(['Time interval: ' num2str(dt) ' days'])
Disp(' ')
%Disp('Processing settings:')
%Disp(cfg)
%Disp(' ')
% Display comand line to repeat this test
% Disp('To repeat this processing paste the following command from the clipboard')
% cmd = ['fmcw_melt(''' name1 ext1 ''',''' name2 ext2 ''')'];
% clipboard('copy',cmd)
% Disp(cmd)


%% Clean shots
% Chirps
if cfg.doManualChirpSelect
    Disp('Using chirp subset defined in config')
    f = fmcw_burst_subset(f,cfg.fchirplist);
    g = fmcw_burst_subset(g,cfg.gchirplist);
end

% Frequency range
if cfg.fRange(1) > 2e8 || cfg.fRange(2) < 4e8
    % Cull f.specCor range
    Disp(['Cropping frequency range to: ' num2str(cfg.fRange(1)/1e6) '-' num2str(cfg.fRange(2)/1e6) ' MHz'])
    Disp(' ')
    f = fmcw_cull_freq(f,cfg.fRange);
    g = fmcw_cull_freq(g,cfg.fRange);
end

% Keep only one attenuator setting (otherwise phase standard deviation and
% error estimate is wrong. (as attenuators have different phase delays)
attSetList = unique(f.chirpAtt,'stable');
if numel(attSetList)>1
    disp(['Warning: multiple attenuator settings found file ' name1 '. Keeping first set only.'])
    fs = fmcw_burst_split_by_att(f);
    f = fs(1); % keeing results from first attenuator setting only.
end
attSetList = unique(g.chirpAtt,'stable');
if numel(attSetList)>1
    disp(['Warning: multiple attenuator settings found file ' name2 '. Keeping first set only.'])
    gs = fmcw_burst_split_by_att(g);
    g = gs(1);
end

if cfg.doClean
    % Cull gross contaminated chirps
    [f,nBad] = fmcw_cull_bad(f,cfg.noisePowerLimit);
    Disp(['Removed ' int2str(nBad) ' contaminated chirps from ' name1])
    [g,nBad] = fmcw_cull_bad(g,cfg.noisePowerLimit);
    Disp(['Removed ' int2str(nBad) ' contaminated chirps from ' name2])
    
    % Cull noisey remaining
    Disp(['Culling noisest ' int2str(cfg.nNoisest) ' shots'])
    f = fmcw_cull_noisey(f,cfg.nNoisest);
    g = fmcw_cull_noisey(g,cfg.nNoisest);
end

%% Phase sensitive processing
% Average all shots in burst and phase process this average shot
fm = fmcw_burst_mean(f);
%[f.rangeCoarse,f.rangeFine,f.specCor,f.specRaw] = fmcw_range(fm,cfg.p,cfg.maxRange,cfg.winFun); % note - should this be changed to a weighted mean for cases where the attenuation changes within the burst
[f.rangeCoarse,f.rangeFine,f.specCor,f.specRaw] = fmcw_range(fm,cfg.p,cfg.maxRange,cfg.winFun); % note - should this be changed to a weighted mean for cases where the attenuation changes within the burst
gm = fmcw_burst_mean(g);
%[g.rangeCoarse,g.rangeFine,g.specCor,g.specRaw] = fmcw_range(gm,cfg.p,cfg.maxRange,cfg.winFun);
[g.rangeCoarse,g.rangeFine,g.specCor,g.specRaw] = fmcw_range(gm,cfg.p,cfg.maxRange,cfg.winFun);
if f.rangeCoarse~=g.rangeCoarse
    error('ranges are different')
end
dr = mean(diff(f.rangeCoarse));

%% Set max depth or relate it to bed
switch cfg.maxDepthMethod
    case 'auto'
        % Use auto detected bed
        Disp(['Searching for bed using method: ' cfg.bedMethod])
        f.bn = fmcw_findbed(f.rangeCoarse,abs(f.specCor),cfg.bedSearchRange,cfg.bedMethod,cfg.ampThreshdB);
        f.bedDepth = f.rangeCoarse(f.bn) + f.rangeFine(f.bn);
        g.bn = fmcw_findbed(g.rangeCoarse,abs(g.specCor),cfg.bedSearchRange,cfg.bedMethod,cfg.ampThreshdB);
        g.bedDepth = g.rangeCoarse(g.bn) + g.rangeFine(g.bn);
        sr_fit.maxDepth = min([f.bedDepth g.bedDepth]) - cfg.bedBuffer;
        Disp(['Shot 1: bed found at ' num2str(f.bedDepth) ' m'])
        Disp(['Shot 2: bed found at ' num2str(g.bedDepth) ' m'])
        %bed.dh = g.bedDepth - f.bedDepth;
        %Disp(['Range difference: ' num2str(bed.dh) ' m'])
        
    case 'config'
        % Use maxDepth from config
        sr_fit.maxDepth = cfg.maxDepthConfig;
        Disp(['Using user defined max fit depth ' num2str(cfg.maxDepthConfig)])
        
    case 'manual'
        % Manually define max depth (graphical)

        % Plot Amplitudes vs range
        figure
        s1ahl = plot(f.rangeCoarse,20*log10(abs(f.specCor)),'r');
        hold on
        s2ahl = plot(g.rangeCoarse,20*log10(abs(g.specCor)),'b');
        ylabel('Vrms (dB)')
        xlabel('range (m)')
        legend([s1ahl(1) s2ahl(1)],{name1,name2},'Location','SouthWest','interpreter','none')
        
        % Graphically select range
        [sr_fit.maxDepth,~] = ginput(1);
        y = get(gca,'ylim');
        plot([sr_fit.maxDepth sr_fit.maxDepth],y,'k');
        set(gca,'xlim',[0 1.5*sr_fit.maxDepth]);
        h = fillXRange([cfg.firnDepth sr_fit.maxDepth],'facecol',[0.6 0.6 0.6],'facealpha',0.3,'edgecol','none');
        pause(1)
        close
        Disp(['Using user defined max fit depth ' num2str(sr_fit.maxDepth)])
end

%% ALIGN BED: Align profiles at bed

if cfg.doAlignBed
    % Locate the bed (in case not done above)
    if ~isfield(f,'bn')
        f.bn = fmcw_findbed(f.rangeCoarse,abs(f.specCor),cfg.bedSearchRange,cfg.bedMethod,cfg.ampThreshdB);
        g.bn = fmcw_findbed(g.rangeCoarse,abs(g.specCor),cfg.bedSearchRange,cfg.bedMethod,cfg.ampThreshdB);
    end
    f.range = f.rangeCoarse + f.rangeFine; % total range shot 1
    g.range = g.rangeCoarse + g.rangeFine; % total range shot 2
    f.bedDepth = f.range(f.bn);
    g.bedDepth = g.range(g.bn);
    
    Disp(['Alinging profiles at bed. Initially the bed can be found at:'])
    Disp([num2str(f.bedDepth) ' m, for shot 1.'])
    Disp([num2str(g.bedDepth) ' m, for shot 2.'])
    
    bnDiff = g.bn-f.bn;
    
    if abs(bnDiff*dr)<1
        g.specRawUnshifted = g.specCor; % keep a copy of g.specCor for plotting
        g.specCorUnshifted = g.specCor; % keep a copy of g.specCor for plotting
        g.specRaw = circshift(g.specRaw,[0 -bnDiff]); % lagg offset
        g.specCor = circshift(g.specCor,[0 -bnDiff]); % lagg offset
        g.rangeFine = circshift(g.rangeFine,[0 -bnDiff]); % lagg offset
        g.bn = f.bn;

        Disp(['Shifting profile 2, ' int2str(bnDiff) ' steps left to align bed. (' num2str(bnDiff*dr) 'm)'])

        f.range = f.rangeCoarse + f.rangeFine; % total range shot 1
        g.range = g.rangeCoarse + g.rangeFine; % total range shot 2
        f.bedDepth = f.range(f.bn);
        g.bedDepth = g.range(g.bn);

        Disp(['Shot 1: new bed found at ' num2str(f.bedDepth) ' m'])
        Disp(['Shot 2: new bed found at ' num2str(g.bedDepth) ' m'])
    else
        Disp(['A bed shift of ' int2str(bnDiff) ' steps, respectively ' num2str(bnDiff*dr) 'm, was found.'])
        Disp(['As this exceeds the expactation limits, no bed shift is undertaken.'])
    end
    % ToDo: Maybe include a plot.
    
    if 0
        switch cfg.bedShiftMethod
            case 'rangeDiff'
                % Just use the (coarse + fine) range difference at peak
                % now use differences across a range of points just before the
                % bed (minimal influence from other reflectors)
                bed.winInds = round(cfg.rangeBedWin(1)/dr):round(cfg.rangeBedWin(2)/dr);
                f.bedInds = f.bn+bed.winInds;
                f.bedDepths = f.range(f.bedInds);
                %bed.winInds = bed.winInds + 3 % test offset windows incorrectly
                g.bedInds = g.bn+bed.winInds;
                g.bedDepths = g.range(g.bedInds);
                bed.dhs = g.bedDepths-f.bedDepths; % range difference at a range of depths near bed peak
                % Select the true offset from the group
                if rem(length(bed.dhs),2)==0
                    bed.dhs = bed.dhs(2:end); % want odd length for median to choose an existing value
                end
                bed.dh = median(bed.dhs);
                bed.dhsDiff = bed.dhs-bed.dh;
                % Check that there's a decent consensus on which way to go - 
                % This assumes that we've correctly aligned the shot segments
                % to less than lambdac/4... or about 0.14m, if we're offset by
                % some integer half-wavelength then we won't know.
                bed.dhsIsWrapped = round(abs(bed.dhsDiff)/(f.lambdac/2));
                bed.fractionWrapped = sum(bed.dhsIsWrapped)/length(bed.dhs);
                if bed.fractionWrapped < cfg.bed.maxWrapFraction
                    % Errors
                    bed.dhsPhaseGrad = std(bed.dhs(~bed.dhsIsWrapped)); % from phase difference gradient across the reflector
                    bed.dheEmperical = sqrt(f.rangeError(f.bn).^2 + g.rangeError(g.bn).^2); % error in bedshift
                    bed.dhe = max([bed.dhsPhaseGrad bed.dheEmperical]);
                else
                    disp('Phase wrapping going on...')
                    figure
                    hist(bed.dhs)
                    figure, plot(bed.dhs)
                    hold on
                    plot(1,median(bed.dhs),'r.')
                    % Define bedshift and shift error based on stats of group
                    bed.dh = mean(bed.dhs);
                    bed.dhe = std(bed.dhs); % we've got to assume a half wavelength error possible
                end

                % This is the old method which just takes the range difference
                % at a single bed reflector.
                %bed.dh = g.bedDepth - f.bedDepth;
                %bed.dhe = sqrt(f.rangeError(f.bn).^2 + g.rangeError(g.bn).^2); % error in bedshift
                % Note this currently ignores the phase gradient error... and bin lag error... to do

            case 'xcorr'
                % Use xcorr as for internals
                f.bedInds = find((f.rangeCoarse>=min(f.rangeCoarse(f.bn) + cfg.xcorBedWin) & f.rangeCoarse<max(f.rangeCoarse(f.bn) + cfg.xcorBedWin))); % depth bins to use (f)
                %bed.maxlag = ceil(cfg.maxBedOffset/dr); % max bin lags
                bed.coarseShift = g.rangeCoarse(g.bn) - f.rangeCoarse(f.bn); % coarseRange between bed peaks
                bed.xcorSearchRange = bed.coarseShift + f.lambdac*[-1 1]; % coarse bedshift +/- 1 wavelength margin
                bed.xcorSearchBinLags = ceil(bed.xcorSearchRange/dr); % lags to use
                [bed.RANGEIND,bed.AMPCOR,bed.COR,bed.LAGS,bed.PE,bed.PSE] = fmcw_xcorr(f.specCor,g.specCor,f.bedInds,bed.xcorSearchBinLags,f.phaseStdError,g.phaseStdError,cfg.p);
                %[~,bed.mci] = max(bed.AMPCOR);
                [~,bed.mci] = max(bed.COR); % max (complex) coherence
                bed.cor = bed.COR(bed.mci);
                bed.RANGE = interp1(1:numel(f.rangeCoarse),f.rangeCoarse,bed.RANGEIND);
                bed.range = bed.RANGE(bed.mci); % bin centre range (weighted by amplitude of f.specCor.*g.specCor)
                bed.lags = bed.LAGS(bed.mci);
                bed.ampCor = bed.AMPCOR(bed.mci);
                bed.coherence = bed.COR(bed.mci);
                %bed.phaseCor = bed.PHASECOR(bed.mci);
                bed.pe = bed.PE(bed.mci);
                bed.pse = bed.PSE(bed.mci);

                % Calculate the total depth shift from the integer bin lags and the phase shift
                bed.lagdh = bed.lags*dr; %
                bed.phasedh = -angle(bed.cor)*(f.lambdac/(4*pi)); % note: phase changes in opposite sense to range
                % Warn if we're near wrapping....
                if abs(angle(bed.cor)) > (2/3)*pi
                    disp(['Warning: phase diff over bed window = ' num2str(angle(bed.cor)) ' radians.'])
                    %keyboard
                end
                %bed.phasedh = -angle(bed.cor)./((4*pi/f.lambdac) - (4*bed.RANGE*f.K/f.ci^2)); % this is the full equation including the term generated by the last term in (13)
                bed.dh = bed.lagdh + bed.phasedh; % range change between shots % Brennan et al. eq 14 and 15;
                bed.dhe = bed.pse*(f.lambdac/(4*pi)); % using the standard error of the phase difference estimated across the range bin
                g.bedInds = f.bedInds + bed.lags;
        end
        Disp(['Bed range change estimated using method: ' cfg.bedShiftMethod])
        %bed.refLevel = fit.zeroShiftDepth;
        %Disp(['Bedchange relative to reference depth: ' num2str(bed.refLevel) ' m'])
        Disp(['Bed range change: ' num2str(bed.dh) '+/- ' num2str(bed.dhe) ' m'])

        if cfg.doPlotBedShift || cfg.doPlotAll
            figure
            subplot(1,2,1)
            [thetaf,rhof] = cart2pol(real(f.specRaw),imag(f.specRaw));
            [thetag,rhog] = cart2pol(real(g.specRaw),imag(g.specRaw));
            polar(thetaf(f.bedInds),rhof(f.bedInds),'r')
            hold on
            polar(thetag(g.bedInds),rhog(g.bedInds),'b')
            legend('f','g')
            title('Raw')

            subplot(1,2,2)
            [thetaf,rhof] = cart2pol(real(f.specCor),imag(f.specCor));
            [thetag,rhog] = cart2pol(real(g.specCor),imag(g.specCor));
            polar(thetaf(f.bedInds),rhof(f.bedInds),'r')
            hold on
            polar(thetag(g.bedInds),rhog(g.bedInds),'b')
            legend('f','g')
            title('Phase corrected to bin centre')
        end
    end
end

%% Estimate error
switch cfg.errorMethod
    case 'emperical'
        % Process each shot separately to get phase standard deviation at each depth
        % shot 1 (f.specCor)
        [~,~,F] = fmcw_range(f,cfg.p,cfg.maxRange,cfg.winFun);
        f.phaseStdDev = std(F)./abs(f.specCor); % phase standard deviation of the burst
        f.phaseStdError = f.phaseStdDev/sqrt(size(f.vif,1)); % phase standard error of the mean shot - using sqrt(n)
        % shot 2 (g.specCor)
        [~,~,G] = fmcw_range(g,cfg.p,cfg.maxRange,cfg.winFun);
        g.phaseStdDev = std(G)./abs(g.specCor); % phase standard deviation of the burst
        g.phaseStdError = g.phaseStdDev/sqrt(size(g.vif,1)); % phase standard error of the mean shot - using sqrt(n)
    case 'assumedNoiseFloor'
        % Error estimate by assuming a noise level and calculating phase
        % noise from this
        noiseFloor = 10.^(cfg.noiseFloordB/20);
        f.phaseStdError = noiseFloor./abs(f.specCor); % phase error shot 1
        g.phaseStdError = noiseFloor./abs(g.specCor); % phase error shot 2
end
%f.rangeError = angle(f.specCor)./((4*pi/f.lambdac) - (4*f.rangeCoarse*f.K/f.ci^2));
%g.rangeError = angle(g.specCor)./((4*pi/g.lambdac) - (4*g.rangeCoarse*g.K/g.ci^2));
f.rangeError = fmcw_phase2range(f.phaseStdError,f.lambdac,f.rangeCoarse,f.K,f.ci);
g.rangeError = fmcw_phase2range(g.phaseStdError,g.lambdac,g.rangeCoarse,g.K,g.ci);

%% ALIGN COARSE: Co-register profile segments (amplitude xcorr)
% Cross correlate internals gi segments to get vertical shift as a function of range
% xcor in big chunks to get minimum chance of wrapping errors.
AC.maxOffset = cfg.maxStrain*(sr_fit.maxDepth-cfg.minDepth); %
maxlag = ceil(AC.maxOffset/dr); % max bin lags
AC.stepSizeM = 5; % cfg.coarseChunkWidth/2; %cfg.coarseChunkWidth/2;
binStart = [cfg.minDepth:AC.stepSizeM:sr_fit.maxDepth-cfg.coarseChunkWidth]; % measure offset over a wider range to plot
[AC.range,AC.dh,AC.lags,AC.ampCor,AC.ampCorProm] = deal(zeros(size(binStart)));
for ii = 1:numel(binStart)
    depthRange = [binStart(ii) binStart(ii)+cfg.coarseChunkWidth];
    fi = find((f.rangeCoarse>=min(depthRange) & f.rangeCoarse<max(depthRange))); % depth bins to use (f)
    [AC.RANGEIND(ii,:),AC.AMPCOR(ii,:),~,AC.LAGS(ii,:)] = fmcw_xcorr(f.specCor,g.specCor,fi,maxlag);
    AC.RANGE(ii,:) = interp1(1:numel(f.rangeCoarse),f.rangeCoarse,AC.RANGEIND(ii,:));
    [~,mci] = max(AC.AMPCOR(ii,:));
    AC.range(ii) = AC.RANGE(ii,mci); % bin centre range (weighted by amplitude of f.specCor.*g.specCor)
    AC.lags(ii) = AC.LAGS(ii,mci);
    AC.ampCor(ii) = AC.AMPCOR(ii,mci);
    AC.dh(ii) = dr*AC.lags(ii); % Range offset (m) between segments
    
    % Quality checks on best correlation
    % Check whether correlation is limited by maxlag
    if mci == 1 || mci == size(AC.LAGS,2) 
        AC.ampCor(ii) = 0;
    end
    % Check prominence of peak (how much better than the next match)
    [cpk,~] = findpeaks(AC.AMPCOR(ii,:),'sortstr','descend','npeaks',2);
    if isempty(cpk) % no peaks!
        AC.ampCorProm(ii) = 0;
    elseif numel(cpk)==1
        AC.ampCorProm(ii) = 1; % this is the only maximum
    else
        AC.ampCorProm(ii) = cpk(1) - cpk(2); % Absolute prominence
    end
end
AC.isGood = AC.ampCor>=cfg.minAmpCor & AC.ampCorProm>=cfg.minAmpCorProm;

% Now fit a polnomial through the lags
AC.P = polyfit(AC.range(AC.isGood),AC.lags(AC.isGood),cfg.polyOrder);

if cfg.doPlotAlignCoarse || cfg.doPlotAll
    figure
    clear ax
    ax(1) = subplottight(3,1,1);
    ii = f.rangeCoarse<sr_fit.maxDepth*1.1;
    plotAmp(f,g,ii);
    title('Profile co-registration - amplitude cross-correlation')
    
    ax(2) = subplottight(3,1,2);
    %sh = surf(AC.RANGE',AC.LAGS',AC.AMPCOR','edgecol','none');
    sh = pcolor(AC.RANGE',AC.LAGS',AC.AMPCOR'); % ,'edgecol','none'
    shading interp
    view(0,90)
    caxis([0.85 1])
    ylabel('bin lag')
    ch = colorbar('East');
    ylabel(ch,'amplitude correlation')
    set(gca,'ydir','normal')
    hold on
    plot3(AC.range,AC.lags,ones(size(AC.range)),'w.') % all lags
    plot3(AC.range(AC.isGood),AC.lags(AC.isGood),ones(1,sum(AC.isGood)),'g.') % only good lags
    AC.lagsPoly = polyval(AC.P,f.rangeCoarse(ii)); % generate smoothed lags
    plot3(f.rangeCoarse(ii),AC.lagsPoly,ones(1,sum(ii)),'g') % poly fit to lags
    
    ax(3) = subplottight(3,1,3);
    plot(AC.range,AC.ampCor)
    ylabel('correlation')
    xlabel('depth')
    
    %set(gcf,'pos',[232 554 560 420])
    linkaxes(ax,'x')
    
    %keyboard
end

%% ALIGN FINE: Estimate phase shift between profile segments (complex xcor)
switch cfg.phaseDiffMethod
    case 'peakDiff'
        % Directly difference phase using value at peak and one bin either
        % side for error
        [peaks,peaki] = findpeaks(abs(f.specCor).*double(f.rangeCoarse<sr_fit.maxDepth));
        for ii = 1:length(peaks)
            % Get coarse offsets at bin centres
            if cfg.doPolySmoothCoarseOffset
                peaklag = round(polyval(AC.P,f.rangeCoarse(peaki(ii)))); % generate smoothed lags
            else
                peaklag = round(interp1(AC.range,AC.lags,f.rangeCoarse(peaki(ii)),'linear','extrap')); % offset between shots at this peak (m)
            end
%            peaklag = round(AC.dhInterp/dr); % bin lag
            jj = peaki(ii)-1:peaki(ii)+1;
            fg = f.specCor(jj).*conj(g.specCor((jj)+peaklag)); % one point either side of peak
            AF.phasedh = -angle(fg)*(f.lambdac/(4*pi)); % note: phase changes in opposite sense to range
            AF.lagdh = peaklag*dr; %
            AF.DH(ii,:) = AF.lagdh + AF.phasedh; % range change between shots % Brennan et al. eq 14 and 15;
            
            if cfg.doPlotAlignFine  || cfg.doPlotAll
                if mod(ii,23) == 0 % only plot some
                
                    figure
                    clear ax
                    ax(1) = subplottight(2,1,1);
                    wi = peaki(ii)-100:peaki(ii)+100;
                    %jj = peaki(ii)-1:peaki(ii)+1;
                    plot(f.rangeCoarse(wi),20*log10(abs(f.specCor(wi))),'r');
                    hold on
                    plot(f.rangeCoarse(jj),20*log10(abs(f.specCor(jj))),'ro');
                    plot(f.rangeCoarse(peaki(ii)),20*log10(abs(f.specCor(peaki(ii)))),'r+');
                    plot(g.rangeCoarse(wi),20*log10(abs(g.specCor(wi+peaklag))),'b');
                    plot(f.rangeCoarse(jj),20*log10(abs(g.specCor(jj+peaklag))),'bo');
                    plot(f.rangeCoarse(peaki(ii)),20*log10(abs(g.specCor(peaki(ii)+peaklag))),'b+');
                    ylabel('Vrms (dB)')
                    %xlabel('range (m)')
                    %legend({name1,name2},'Location','SouthWest','interpreter','none')
                    title('Profile phase difference - peaks only (peak 1)')
                    
                    ax(2) = subplottight(2,1,2);
                    plot(f.rangeCoarse(wi),angle(f.specRaw(wi)),'r');
                    hold on
                    plot(g.rangeCoarse(wi),angle(g.specRaw(wi+peaklag)),'b');
                    ylabel('phase (rad)')
                    linkaxes(ax,'x')
                    
                    %keyboard
                end
            end
        end
        dh  = AF.DH(:,2);
        dhe = (AF.DH(:,3)-AF.DH(:,1))/2; % average of upper and lower error...
        range = transpose(f.rangeCoarse(peaki));
        AF.coherence = ones(size(range));
        
    case 'xcor'
        % complex correlation: xcor in small chunks to get good depth resolution
        % this also gives us the AF.coherence of the segments
        stepSizeM = cfg.chunkWidth; % cfg.chunkWidth/2; cfg.chunkWidth;
        % binStart = [cfg.minDepth:stepSizeM:fit.maxDepth-cfg.chunkWidth]; % measure offset over a wider range to plot
        binStart = [cfg.minDepth:cfg.chunkStep:sr_fit.maxDepth-cfg.chunkWidth]; % Noted: added by Falk, take intermediate steps to consider good reflectors at chunk-boundaries
        %OffsetRange = AC.maxOffset;
        for ii = 1:numel(binStart)
            depthRange = [binStart(ii) binStart(ii) + cfg.chunkWidth];
            binDepth = mean(depthRange);
            maxlag = ceil(AC.maxOffset/dr); % max bin lags
            fi = find((f.rangeCoarse>=min(depthRange) & f.rangeCoarse<max(depthRange))); % depth bins to use (f)
            [AF.RANGEIND(ii,:),AF.AMPCOR(ii,:),AF.COR(ii,:),AF.LAGS(ii,:),AF.PE(ii,:),AF.PSE(ii,:)] = fmcw_xcorr(f.specCor,g.specCor,fi,maxlag,f.phaseStdError,g.phaseStdError,cfg.p);
            AF.RANGE(ii,:) = interp1(1:numel(f.rangeCoarse),f.rangeCoarse,AF.RANGEIND(ii,:));
            if cfg.doUseCoarseOffset % Define the bin lag from the coarse correlation
                % Get coarse offsets at bin centres
                if cfg.doPolySmoothCoarseOffset
                    AC.dhInterp = dr*polyval(AC.P,binDepth); % generate smoothed lags
                else
                    AC.dhInterp = interp1(AC.range,AC.dh,binDepth,'linear','extrap'); %
                end
                [~,AF.mci(ii)] = min(abs(AC.dhInterp/dr-AF.LAGS(ii,:))); % bin lags index
            else
                [~,AF.mci(ii)] = max(AF.AMPCOR(ii,:)); % use best lag from fine cor
            end
            AF.ampCor(ii) = AF.AMPCOR(ii,AF.mci(ii));
            AF.cor(ii) = AF.COR(ii,AF.mci(ii)); % complex correlation at best amp correlation point
            % will be overwritten later...
            range(ii) = AF.RANGE(ii,AF.mci(ii)); % bin centre range (weighted by amplitude of f.specCor.*g.specCor)
            % will be overwritten later...
        end
        AF.PHASECOR = abs(AF.COR)./AF.AMPCOR;
        AF.lagvec = AF.LAGS(1,:);
        
        % Unwrap
        if cfg.doSmartUnwrap
            % redefine best lags as those closest to zero phase diff
            dphi = 2*pi*f.fc/(f.B*cfg.p); % phase change over 1 bin (brennan eq 16)
            mci_pm = AF.mci - fix(angle(AF.cor)/dphi); % lag to give smallest phase diff near good correlation point
            % Note Falk: Why fix()? numWavelengthError smaller with round()
            mci_pm(mci_pm<1) = mci_pm(mci_pm<1) + round(2*pi/dphi); % deal with edge effects
            mci_pm(mci_pm>size(AF.COR,2)) = mci_pm(mci_pm>size(AF.COR,2)) - round(2*pi/dphi); % deal with edge effects
            %mci_pmu = round(unwrap(mci_pm*dphi)/dphi); % lag with no large steps
            
            % Try following phase difference minimum
            % start from maximum correlation point(s)
%             if ~cfg.doBulkAllignment
                [~,si] = max(AF.ampCor); % index of best correlation point
                mci_pmu(si) = mci_pm(si);
%             else
%                 % the above method is suseptable to starting in the wrong phase
%                 % catchment (integer wavelength error) if binWidth is small. Better to
%                 % use the 0 lag at the centre of the bulk alignment region as
%                 % this has used a much larger xcor window
%                 [~,si] = min(abs(range-mean(cfg.bulkAlignRange))); % range bin closest to centre of bulk alignment range
%                 mci_pmu(si) = find(AF.lagvec==0);
%             end
            % Forwards to end
            for ii = si+1:length(AF.mci)
                [~,pki] = findpeaks(-abs(angle(AF.COR(ii,:)))); % index of phase difference minimum
                [~,bpki] = min(abs(mci_pmu(ii-1)-pki)); % closest minimum to last
                mci_pmu(ii) = pki(bpki);
                %mci_pmu(ii) = mci_pmu(ii-1) + round(angle(c(mci_pm(ii-1),ii))/dphi); % lag to give smallest phase diff
            end
            % Backwards to start
            for ii = si-1:-1:1
                [~,pki] = findpeaks(-abs(angle(AF.COR(ii,:)))); % index of phase difference minimum
                [~,bpki] = min(abs(mci_pmu(ii+1)-pki)); % closest minimum to last
                mci_pmu(ii) = pki(bpki);
                %mci_pmu(ii) = mci_pmu(ii-1) + round(angle(c(mci_pm(ii-1),ii))/dphi); % lag to give smallest phase diff
            end
            AF.mciu = mci_pmu;
            
            % Now check that starting at the bulk alignment centre got the
            % same interger pathlength/wavelength as most of the local amp
            % correlations.
            binLagError = mci_pm - mci_pmu;
            numWavelengthError = sum(abs(binLagError)>0);
        end
        
        % Extract values at chosen offsets
        if cfg.doSmartUnwrap
            igood = AF.mciu;
        else
            igood = AF.mci;
        end
        for ii = 1:length(igood)
            AF.cor(ii) = AF.COR(ii,igood(ii));
            range(ii) = AF.RANGE(ii,igood(ii)); % bin centre range (weighted by amplitude of f.specCor.*g.specCor)
            AF.lags(ii) = AF.LAGS(ii,igood(ii));
            AF.ampCor(ii) = AF.AMPCOR(ii,igood(ii));
            AF.coherence(ii) = AF.COR(ii,igood(ii));
            AF.phaseCor(ii) = AF.PHASECOR(ii,igood(ii));
            AF.pe(ii) = AF.PE(ii,igood(ii));
            AF.pse(ii) = AF.PSE(ii,igood(ii));
        end
        
        % Calculate the total depth shift from the integer bin lags and the phase shift
        AF.lagdh = AF.lags*dr; % 
        AF.phasedh = -angle(AF.cor)*(f.lambdac/(4*pi)); % note: phase changes in opposite sense to range
        %AF.phasedh = -angle(AF.cor)./((4*pi/f.lambdac) - (4*AF.RANGE*f.K/f.ci^2)); % this is the full equation including the term generated by the last term in (13)
        dh = AF.lagdh + AF.phasedh; % range change between shots % Brennan et al. eq 14 and 15;
        dhe = AF.pse*(f.lambdac/(4*pi)); % using the standard error of the phase difference estimated across the range bin
end

if cfg.doPlotAlignFine || cfg.doPlotAll
    figure
    
    %sh = surf(transpose(AF.RANGE),transpose(AF.LAGS),transpose(angle(AF.COR)),'edgecol','none'); % ,transpose(AF.AMPCOR)
    sh = surf(transpose(AF.RANGE),transpose(AF.LAGS),transpose(abs(AF.COR)),transpose(angle(AF.COR)),'edgecol','none'); % ,transpose(AF.AMPCOR)
    set(sh,'alphadata',transpose(AF.AMPCOR))
    
    view(0,90)
    ylabel('bin lag')
    xlabel('depth (m)')
    colormap jet
    set(gca,'clim',[-pi pi])
    ch = colorbar('East');
    ylabel(ch,'phase diff')
    set(gca,'ydir','normal')
    hold on
    plot3(range,AF.lagvec(AF.mci),pi*ones(size(range)),'k.','markersize',30); % best amp
    h(2) = plot3(range,AF.lagvec(AF.mci),pi*ones(size(range)),'w.','markersize',20); % best amp
    if cfg.doSmartUnwrap
        plot3(range,AF.lagvec(AF.mciu),pi*ones(size(range)),'k.','markersize',18); % unwrapped
        h(3) = plot3(range,AF.lagvec(AF.mciu),pi*ones(size(range)),'g.','markersize',10); % unwrapped
        plot3(range(si),AF.lagvec(AF.mci(si)),pi,'k.','markersize',40); % best amp
        h(1) = plot3(range(si),AF.lagvec(AF.mci(si)),pi,'m.','markersize',30); % best amp
    end
    
    if cfg.doSmartUnwrap
        legend(h,{'start (best amp-cor)','best amp-cor lag','unwrapped'})
    else
        %legend(h,{'start (best amp-cor)','best amp-cor lag'})
    end
    %set(gca,'xticklabel',[])
    title('Phase difference from x-cor')
%    set(gcf,'pos',[792 527 565 447])
    
    %keyboard
end

if cfg.doSmartUnwrap
    if numWavelengthError>(0.5*numel(range))
        disp('Warning: Integer half-wavelength ambiguity detected - check manually')
        %keyboard
    elseif numWavelengthError>(0.1*numel(range))
        disp(['Warning ' int2str(numWavelengthError) ' wavelength errors detected in ' int2str(numel(range)) ' points,'])
    end
end

%% FIT: Estimate internal strain rate
% linear fit to range diff of internals to estimate vertical strain

% Apply criteria to select good reflectors
AC.ampCorInterp = interp1(AC.range,AC.ampCor,range,'nearest','extrap'); % resample coarse amplitude correlation at fine range points
sr_pt.gi = find(range>=10 & ...                %dhe<1 & ...
                range<=sr_fit.maxDepth & ...
                abs(AF.coherence)>=cfg.minCohereFine & ...
                AC.ampCorInterp>=cfg.minCohereCoarse & ...
                AF.phaseCor>= cfg.minCoherePhase); %  all points in range with acceptable quality ( can add other criteria)

sr_fit.gi = sr_pt.gi(range(sr_pt.gi)>=cfg.firnDepth); % cut range at firn depth.
            
% Extract vertical velocity      
if cfg.doStrainPointwise 
    gn = sr_pt.gi(diff(sr_pt.gi)==1); % First of neighboring pair of good reflectors.
    sr_pt.gn = gn;
    sr_pt.range_gn = (range(gn+1)+range(gn))./2;
    sr_pt.range_gn = (range(gn+1)+range(gn))./2;
    dL = dh(gn+1)-dh(gn);
    L = range(gn+1)-range(gn);
    sr_pt.vs = dL./L; % vertical strain
    sr_pt.vsr = sr_pt.vs*daysPerYear/dt; % vertical strain rate (per year)
    sr_pt.vse = sqrt((1./(L.^2)).*(dhe(gn+1).^2+dhe(gn).^2)); %add +(-dL./(L.^2)).^2 to 1. term for errors on range
    sr_pt.vsre = sr_pt.vse*daysPerYear/dt; 


    Nmean = cfg.movingMean;
    sr_pt.range_gn_mm = movmean(sr_pt.range_gn,Nmean);
    sr_pt.vsr_mm = movmean(sr_pt.vsr,Nmean);
    sr_pt.vsr_mm_std = movstd(sr_pt.vsr,Nmean);
    sr_pt.vsre_mm = sqrt(movsum(sr_pt.vsre*2,Nmean))/Nmean;
    disp(mean(sr_pt.vsr));

    strain_P = polyfit(sr_pt.range_gn,sr_pt.vsr,cfg.polyOrder);
    all_p = polyval(strain_P,range);


    %N_gn = 4;
    %gn_alt = fit.gi(((fit.gi(1+N_gn:end)-fit.gi(1:end-N_gn))==4));
    %L_alt = range(gn_alt+N_gn)-range(gn_alt);
    %dL_alt = dh(gn_alt+N_gn)-dh(gn_alt);
    %vs_alt = dL_alt./L_alt; % vertical strain
    %vsr_alt = vs_alt*daysPerYear/dt; % vertical strain rate (per year)
    %vse_alt = sqrt((1./(L_alt.^2)+(-dL_alt./(L_alt.^2)).^2).*(dhe(gn_alt+N_gn).^2+dhe(gn_alt).^2)); %add +(-dL./(L.^2)).^2 to 1. term for errors on range
    %vsre_alt = vse_alt*daysPerYear/dt; 
    %range_gn_alt = (range(gn_alt+N_gn)+range(gn_alt))./2;
    %strain_P_alt = polyfit(range_gn_alt,vsr_alt,cfg.polyOrder);
    %all_p_alt = polyval(strain_P_alt,range);

    %disp(mean(vsr_alt));

    %[fit.M,fit.Me] = menke_fit(predictors,dh(fit.gi),dhe(fit.gi));



    if cfg.doPlotVSR || cfg.doPlotAll
       figure
       clear ax
       ax(1) = subplot(2,1,1);
       plotAmp(f,g);
       h = fillXRange([cfg.firnDepth sr_fit.maxDepth],'facecol',[0.6 0.6 0.6],'facealpha',0.3,'edgecol','none');
       set(gca,'xlim',[0 1.2*sr_fit.maxDepth])
       box on
       title('Strain fitting')

       ax(2) = subplot(2,1,2);
       %plot(range_gn_alt,vsr_alt,'b.'); % best fit line
       hold on
       plot(sr_pt.range_gn,sr_pt.vsr,'k.'); % best fit line
       plot(sr_pt.range_gn_mm,sr_pt.vsr_mm,'b.'); % best fit line
       %errorbar(range_gn_mean,vsr_mean,vsre_mean,-vsre_mean,'r.'); % errors
       %errorbar(sr_pt.range_gn_mm,sr_pt.vsr_mm,sr_pt.vsr_mm_std,-sr_pt.vsr_mm_std,'g.'); % errors
       %errorbar(sr_pt.range_gn,sr_pt.vsr,sr_pt.vsre,-sr_pt.vsre,'k.'); % errors
       hold on
       %plot(range, all_p_alt,'r') % poly fit to lags
       %plot(range, all_p,'g') % poly fit to lags
       %plot(range(fit.gi),dh(fit.gi),'g.','markerSize',20) % used points
       %h = fillXRange([cfg.firnDepth fit.maxDepth],'facecol',[0.6 0.6 0.6],'facealpha',0.3,'edgecol','none');
       ylabel('vertical strain rate (yr^{-1})')
       box on

       linkaxes(ax,'x')
    %   set(gcf,'pos',[1360 554 560 420])
    end
end


if cfg.doStrainFit
    % Loop through fitting to different depth ranges
    if length(sr_fit.gi) < cfg.minPointsToFit
        sr_fit.M = [nan nan];
        sr_fit.Me = [nan nan];
        disp(['warning: only ' int2str(length(sr_fit.gi)) ' good matches in depth range'])
        error('Not enough good correlations in depth range')
    else
        %setup things to fit here - i.e. offset, linear and optionally bending
        sr_fit.midRange = mean(range(sr_fit.gi)); % mean range used in fit
        switch cfg.fit.type
            case 'linear'
                predictors = [ones(numel(range(sr_fit.gi)),1) transpose(range(sr_fit.gi)-sr_fit.midRange)];
                % Note: its important to de-mean range, otherwise errors are wrong.
            case 'quadratic' % for bending
                predictors = [ones(numel(range(sr_fit.gi)),1) (range(sr_fit.gi)-sr_fit.midRange).' ((range(sr_fit.gi)-sr_fit.midRange).^2).'];
    %             % See Jenkins 2006
    %            not yet implemented - we have to change lots below for this...
        end
        switch cfg.fitMethod
            case 'menke'
                [sr_fit.M,sr_fit.Me] = menke_fit(predictors,dh(sr_fit.gi),dhe(sr_fit.gi));% linear fit with error
            case 'robust'
                [sr_fit.M,sr_fit.STATS] = robustfit(predictors(:,2:end),dh(sr_fit.gi)); % linear fit - downweighting outliers
                sr_fit.Me = sr_fit.STATS.se; % ???

                switch cfg.fit.errorMetric
                    case 'standardError'

                    case 'alpha'
                        % Custom level confidence interval
                        sr_fit.alpha = 0.05;
                        sr_fit.Mci = nlparci(sr_fit.M,sr_fit.STATS.resid,'cov',sr_fit.STATS.covb,'alpha',sr_fit.alpha);
                        sr_fit.Me_alpha = diff(sr_fit.Mci')/2;
                        sr_fit.Me(2) = sr_fit.Me_alpha(2);
                end
            case 'regress'
                [sr_fit.M,bint,sr_fit.resid,rint,sr_fit.stats] = regress(dh(sr_fit.gi),predictors);
                sr_fit.Me = diff(bint,2)/4; % trying to get 1 std from 5-95% ci

            case 'fitnlm' % non-linear
                modelFun = @(b,x) b(1)+b(2).*x; % just linear fit
                start = [0 0.001];
                wnlm = fitnlm(range(sr_fit.gi),dh(sr_fit.gi),modelFun,start,'Weight',w);
                line(range(sr_fit.gi),predict(wnlm,xx),'color','b')
                % not finished... the point here was to use weights based on
                % phase error.
        end
    end
    sr_fit.M(1) = sr_fit.M(1)-sr_fit.midRange*sr_fit.M(2); % correct for the fact we removed mean range before fit
    sr_fit.vsr = sr_fit.M(2)*daysPerYear/dt; % vertical strain rate (per year)
    sr_fit.vsre = sr_fit.Me(2)*daysPerYear/dt; % vertical strain rate (per year)
    sr_fit.dh = sr_fit.M(1)+range.*sr_fit.M(2); % linear fit to internal strain
    sr_fit.zeroShiftDepth = -sr_fit.M(1)/sr_fit.M(2); % depth at which fit offset is zero;
    sr_fit.dhe = sqrt(sr_fit.Me(1).^2 + ((range-sr_fit.midRange).*sr_fit.Me(2)).^2); % linear fit to internal strain
    [sr_fit.R,sr_fit.P] = corrcoef(range(sr_fit.gi),dh(sr_fit.gi)); % correlation coefficient
    sr_fit.resid = dh-sr_fit.dh;

    % calculate significance at 95% level

    % Check stats on residule errors from linear fit???
    Disp(' ')
    Disp('Estimating vertical strain using linear fit to phase offset')
    Disp(['fitting using method: ' cfg.fitMethod])
    Disp(['Vertical strain rate: ' num2str(1000*sr_fit.vsr,'%5.2f') ' +/- ' num2str(1000*abs(sr_fit.vsre),'%5.2f') ' millistrain/year'])
    Disp(['R^2 = ' num2str(sr_fit.R(2).^2)])
    Disp(['P = ' num2str(sr_fit.P(2))])
    Disp(' ')

    if cfg.doPlotFit || cfg.doPlotAll
       figure
       clear ax
       ax(1) = subplot(2,1,1);
       plotAmp(f,g);
       h = fillXRange([cfg.firnDepth sr_fit.maxDepth],'facecol',[0.6 0.6 0.6],'facealpha',0.3,'edgecol','none');
       set(gca,'xlim',[0 1.2*sr_fit.maxDepth])
       box on
       title('Strain fitting')

       ax(2) = subplot(2,1,2);
       plotrange = [0 sr_fit.midRange 1.2*sr_fit.maxDepth];
       plotdhe = sqrt(sr_fit.Me(1).^2 + ((plotrange-sr_fit.midRange).*sr_fit.Me(2)).^2);
       plotdh = sr_fit.M(1)+plotrange.*sr_fit.M(2);
       plotdhe_high = plotdh + plotdhe;
       plotdhe_low = plotdh - plotdhe;
       patch([plotrange fliplr(plotrange)],[plotdhe_high fliplr(plotdhe_low)],[0.6 0.6 0.6],'edgeColor','none') % best fit error band
       hold on
       plot(plotrange,plotdh,'g'); % best fit line
       errorbar(range,dh,dhe,-dhe,'k'); % errors
       plot(range(sr_fit.gi),dh(sr_fit.gi),'g.','markerSize',20) % used points
       %h = fillXRange([cfg.firnDepth fit.maxDepth],'facecol',[0.6 0.6 0.6],'facealpha',0.3,'edgecol','none');
       ylabel('range difference (m)')
       box on

       linkaxes(ax,'x')
    %   set(gcf,'pos',[1360 554 560 420])
    end
end

if cfg.doPlotResid || cfg.doPlotAll
    figure
    % Plot residuls from linear
    stem(range,1000*sr_fit.resid,'filled','col',[0.6 0.6 0.6],'linewidth',2)
    hold on
    errorbar(range(sr_fit.gi),1000*sr_fit.resid(sr_fit.gi),1000*dhe(sr_fit.gi),-1000*dhe(sr_fit.gi),'r');
    grid on
    ylabel('Range residual from linear fit (mm)')
    xlabel('range (m)')
    
%    set(gcf,'pos',[1358 -2 560 420])
end

%% Measure bedshift
if 0 %cfg.doMeltEstimate
    % Locate the bed (in case not done above)
    if ~isfield(f,'bn')
        f.bn = fmcw_findbed(f.rangeCoarse,abs(f.specCor),cfg.bedSearchRange,cfg.bedMethod,cfg.ampThreshdB);
        g.bn = fmcw_findbed(g.rangeCoarse,abs(g.specCor),cfg.bedSearchRange,cfg.bedMethod,cfg.ampThreshdB);
    end
    f.range = f.rangeCoarse + f.rangeFine; % total range shot 1
    g.range = g.rangeCoarse + g.rangeFine; % total range shot 2
    f.bedDepth = f.range(f.bn);
    g.bedDepth = g.range(g.bn);
    
    Disp(['Shot 1: bed found at ' num2str(f.bedDepth) ' m'])
    Disp(['Shot 2: bed found at ' num2str(g.bedDepth) ' m'])
    
    switch cfg.bedShiftMethod
        case 'rangeDiff'
            % Just use the (coarse + fine) range difference at peak
            % now use differences across a range of points just before the
            % bed (minimal influence from other reflectors)
            bed.winInds = round(cfg.rangeBedWin(1)/dr):round(cfg.rangeBedWin(2)/dr);
            f.bedInds = f.bn+bed.winInds;
            f.bedDepths = f.range(f.bedInds);
            %bed.winInds = bed.winInds + 3 % test offset windows incorrectly
            g.bedInds = g.bn+bed.winInds;
            g.bedDepths = g.range(g.bedInds);
            bed.dhs = g.bedDepths-f.bedDepths; % range difference at a range of depths near bed peak
            % Select the true offset from the group
            if rem(length(bed.dhs),2)==0
                bed.dhs = bed.dhs(2:end); % want odd length for median to choose an existing value
            end
            bed.dh = median(bed.dhs);
            bed.dhsDiff = bed.dhs-bed.dh;
            % Check that there's a decent consensus on which way to go - 
            % This assumes that we've correctly aligned the shot segments
            % to less than lambdac/4... or about 0.14m, if we're offset by
            % some integer half-wavelength then we won't know.
            bed.dhsIsWrapped = round(abs(bed.dhsDiff)/(f.lambdac/2));
            bed.fractionWrapped = sum(bed.dhsIsWrapped)/length(bed.dhs);
            if bed.fractionWrapped < cfg.bed.maxWrapFraction
                % Errors
                bed.dhsPhaseGrad = std(bed.dhs(~bed.dhsIsWrapped)); % from phase difference gradient across the reflector
                bed.dheEmperical = sqrt(f.rangeError(f.bn).^2 + g.rangeError(g.bn).^2); % error in bedshift
                bed.dhe = max([bed.dhsPhaseGrad bed.dheEmperical]);
            else
                disp('Phase wrapping going on...')
                figure
                hist(bed.dhs)
                figure, plot(bed.dhs)
                hold on
                plot(1,median(bed.dhs),'r.')
                % Define bedshift and shift error based on stats of group
                bed.dh = mean(bed.dhs);
                bed.dhe = std(bed.dhs); % we've got to assume a half wavelength error possible
            end
            
            % This is the old method which just takes the range difference
            % at a single bed reflector.
            %bed.dh = g.bedDepth - f.bedDepth;
            %bed.dhe = sqrt(f.rangeError(f.bn).^2 + g.rangeError(g.bn).^2); % error in bedshift
            % Note this currently ignores the phase gradient error... and bin lag error... to do
            
        case 'xcorr'
            % Use xcorr as for internals
            f.bedInds = find((f.rangeCoarse>=min(f.rangeCoarse(f.bn) + cfg.xcorBedWin) & f.rangeCoarse<max(f.rangeCoarse(f.bn) + cfg.xcorBedWin))); % depth bins to use (f)
            %bed.maxlag = ceil(cfg.maxBedOffset/dr); % max bin lags
            bed.coarseShift = g.rangeCoarse(g.bn) - f.rangeCoarse(f.bn); % coarseRange between bed peaks
            bed.xcorSearchRange = bed.coarseShift + f.lambdac*[-1 1]; % coarse bedshift +/- 1 wavelength margin
            bed.xcorSearchBinLags = ceil(bed.xcorSearchRange/dr); % lags to use
            [bed.RANGEIND,bed.AMPCOR,bed.COR,bed.LAGS,bed.PE,bed.PSE] = fmcw_xcorr(f.specCor,g.specCor,f.bedInds,bed.xcorSearchBinLags,f.phaseStdError,g.phaseStdError,cfg.p);
            %[~,bed.mci] = max(bed.AMPCOR);
            [~,bed.mci] = max(bed.COR); % max (complex) coherence
            bed.cor = bed.COR(bed.mci);
            bed.RANGE = interp1(1:numel(f.rangeCoarse),f.rangeCoarse,bed.RANGEIND);
            bed.range = bed.RANGE(bed.mci); % bin centre range (weighted by amplitude of f.specCor.*g.specCor)
            bed.lags = bed.LAGS(bed.mci);
            bed.ampCor = bed.AMPCOR(bed.mci);
            bed.coherence = bed.COR(bed.mci);
            %bed.phaseCor = bed.PHASECOR(bed.mci);
            bed.pe = bed.PE(bed.mci);
            bed.pse = bed.PSE(bed.mci);
            
            % Calculate the total depth shift from the integer bin lags and the phase shift
            bed.lagdh = bed.lags*dr; %
            bed.phasedh = -angle(bed.cor)*(f.lambdac/(4*pi)); % note: phase changes in opposite sense to range
            % Warn if we're near wrapping....
            if abs(angle(bed.cor)) > (2/3)*pi
                disp(['Warning: phase diff over bed window = ' num2str(angle(bed.cor)) ' radians.'])
                %keyboard
            end
            %bed.phasedh = -angle(bed.cor)./((4*pi/f.lambdac) - (4*bed.RANGE*f.K/f.ci^2)); % this is the full equation including the term generated by the last term in (13)
            bed.dh = bed.lagdh + bed.phasedh; % range change between shots % Brennan et al. eq 14 and 15;
            bed.dhe = bed.pse*(f.lambdac/(4*pi)); % using the standard error of the phase difference estimated across the range bin
            g.bedInds = f.bedInds + bed.lags;
    end
    Disp(['Bed range change estimated using method: ' cfg.bedShiftMethod])
    %bed.refLevel = fit.zeroShiftDepth;
    %Disp(['Bedchange relative to reference depth: ' num2str(bed.refLevel) ' m'])
    Disp(['Bed range change: ' num2str(bed.dh) '+/- ' num2str(bed.dhe) ' m'])
    
    if cfg.doPlotBedShift || cfg.doPlotAll
        figure
        subplot(1,2,1)
        [thetaf,rhof] = cart2pol(real(f.specRaw),imag(f.specRaw));
        [thetag,rhog] = cart2pol(real(g.specRaw),imag(g.specRaw));
        polarplot(thetaf(f.bedInds),rhof(f.bedInds),'r')
        hold on
        polarplot(thetag(g.bedInds),rhog(g.bedInds),'b')
        legend('f','g')
        title('Raw')
        
        subplot(1,2,2)
        [thetaf,rhof] = cart2pol(real(f.specCor),imag(f.specCor));
        [thetag,rhog] = cart2pol(real(g.specCor),imag(g.specCor));
        polarplot(thetaf(f.bedInds),rhof(f.bedInds),'r')
        hold on
        polarplot(thetag(g.bedInds),rhog(g.bedInds),'b')
        legend('f','g')
        title('Phase corrected to bin centre')
    end
end

%% Create output structure
site.t1 = f.TimeStamp; % time
site.t2 = g.TimeStamp; % time
site.dt = dt;
site.file1 = filename1;
site.file2 = filename2;
site.range = range;
site.sr_fit = sr_fit;
site.sr_pt = sr_pt; % /year
% if ~cfg.exportVSR % this allows use to calculate a vsr within the firn for estimating melt only - but not export the used vsr
%     fit.vsr = nan;
%     fit.vsre = nan;
%     fit.resid = nan;
% end
%site.vsr = sr_fit.vsr; % /year
%site.vsre = sr_fit.vsre; % /year
%site.resid = sr_fit.resid;
site.cfg = cfg;
site.AC = AC;
AF.dh = dh; % get a few things we didn't put in AF before
AF.dhe = dhe;
site.AF = AF;
%site.bed = bed;

%% Save output
if cfg.doSaveOutput
    [~,name1,~] = fileparts(filename1);
    [~,name2,~] = fileparts(filename2);
    outfile = [name1 name2 '_meltsite.mat'];
    save(outfile,'site');
end

% Display comand line to repeat this test
% Disp('To repeat this processing paste the following command from the clipboard')
% cmd = ['fmcw_melt(''' filename1 ''',''' filename2 ''')'];
% clipboard('copy',cmd)
% Disp(cmd)
%keyboard

function Disp(text) % Only diplay output if cfg.verbose on
global cfg
if cfg.verbose
    disp(text)
end

function h = fillXRange(x,varargin)
y = get(gca,'ylim');
h = patch([x(1) x(2) x(2) x(1)],[y(1) y(1) y(2) y(2)],[0.6 0.6 0.6]);
set(h,varargin{:})

function plotAmp(f,g,ii)
if nargin<3
    ii = 1:numel(f.rangeCoarse);
end
%plot the standard amplitude profile
fah = plot(f.rangeCoarse(ii),20*log10(abs(f.specCor(ii))),'r');
hold on
gah = plot(g.rangeCoarse(ii),20*log10(abs(g.specCor(ii))),'b');
ylabel('Vrms (dB)')
xlabel('range (m)')
legend([fah(1) gah(1)],{'f','g'},'Location','SouthWest')
ylim([-140 -20])





