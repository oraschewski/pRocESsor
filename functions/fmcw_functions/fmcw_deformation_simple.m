function site = fmcw_deformation_simple(filename1,filename2)
%FMCW_DEFORMATION_SIMPLE determines the deformation between two ApRES
%measurements. It is a simplified version of fmcw_deformation
%   Detailed explanation goes here


%% Load settings
% Processing settings
global cfg
cfg = fmcw_deformation_config;

% Define Constants
daysPerYear = 365.25;

% Load data
[~,name1,~] = fileparts(filename1);
[~,name2,~] = fileparts(filename2);
fBurst = 1;
gBurst = 1;
% Load
f = fmcw_load(filename1,fBurst); 
g = fmcw_load(filename2,gBurst);

%% Print Output
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Radar processor ' mfilename])
disp(' ')
disp('Using files:')
disp('Filename                                  Date/Time')
disp([name1 '   ' datestr(f.TimeStamp)])
disp([name2 '   ' datestr(g.TimeStamp)])
dt = g.TimeStamp-f.TimeStamp; % time difference days
disp(['Time interval: ' num2str(dt) ' days'])
disp(' ')

%% Clean shots
% Chirps
if cfg.doManualChirpSelect
    disp('Using chirp subset defined in config')
    f = fmcw_burst_subset(f,cfg.fchirplist);
    g = fmcw_burst_subset(g,cfg.gchirplist);
end

% Frequency range
if cfg.fRange(1) > 2e8 || cfg.fRange(2) < 4e8
    % Cull f.specCor range
    disp(['Cropping frequency range to: ' num2str(cfg.fRange(1)/1e6) '-' num2str(cfg.fRange(2)/1e6) ' MHz'])
    disp(' ')
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
    disp(['Removed ' int2str(nBad) ' contaminated chirps from ' name1])
    [g,nBad] = fmcw_cull_bad(g,cfg.noisePowerLimit);
    disp(['Removed ' int2str(nBad) ' contaminated chirps from ' name2])
    
    % Cull noisey remaining
    disp(['Culling noisest ' int2str(cfg.nNoisest) ' shots'])
    f = fmcw_cull_noisey(f,cfg.nNoisest);
    g = fmcw_cull_noisey(g,cfg.nNoisest);
end

%% Phase sensitive processing
% Average all shots in burst and phase process this average shot
fm = fmcw_burst_mean(f);
[f.rangeCoarse,f.rangeFine,f.specCor,f.specRaw] = fmcw_range(fm,cfg.p,cfg.maxRange,cfg.winFun); % note - should this be changed to a weighted mean for cases where the attenuation changes within the burst
gm = fmcw_burst_mean(g);
[g.rangeCoarse,g.rangeFine,g.specCor,g.specRaw] = fmcw_range(gm,cfg.p,cfg.maxRange,cfg.winFun);
if f.rangeCoarse~=g.rangeCoarse
    error('ranges are different')
end
dr = mean(diff(f.rangeCoarse));

%% Detect and align bed
bed_align_done = false;
bed_iteration = 0;
bedMethods = [cfg.bedMethod "maxAmp" "xcor" "ampThresh"];
while ~bed_align_done
    bed_iteration = bed_iteration + 1;
    bedMethod = bedMethods(bed_iteration);
    
    disp(['Searching for bed using method: ' bedMethod])
    f.bn = fmcw_findbed(f.rangeCoarse,abs(f.specCor),cfg.bedSearchRange, bedMethod, cfg.ampThreshdB);
    f.bedDepth = f.rangeCoarse(f.bn) + f.rangeFine(f.bn);
    g.bn = fmcw_findbed(g.rangeCoarse,abs(g.specCor),cfg.bedSearchRange, bedMethod, cfg.ampThreshdB);
    g.bedDepth = g.rangeCoarse(g.bn) + g.rangeFine(g.bn);
    maxDepth = min([f.bedDepth g.bedDepth]) - cfg.bedBuffer;
    disp(['Alinging profiles at bed. Initially the bed can be found at:'])
    disp([num2str(f.bedDepth) ' m, for shot 1.'])
    disp([num2str(g.bedDepth) ' m, for shot 2.'])

    bnDiff = g.bn-f.bn;


    if abs(bnDiff*dr)<cfg.bedMaxShift
        bed_align_done = true;
        g.specRawUnshifted = g.specCor; % keep a copy of g.specCor for plotting
        g.specCorUnshifted = g.specCor; % keep a copy of g.specCor for plotting
        g.specRaw = circshift(g.specRaw,[0 -bnDiff]); % lagg offset
        g.specCor = circshift(g.specCor,[0 -bnDiff]); % lagg offset
        g.rangeFine = circshift(g.rangeFine,[0 -bnDiff]); % lagg offset
        g.bn = f.bn;

        disp(['Shifting profile 2, ' int2str(bnDiff) ' steps left to align bed. (' num2str(bnDiff*dr) 'm)'])

        f.range = f.rangeCoarse + f.rangeFine; % total range shot 1
        g.range = g.rangeCoarse + g.rangeFine; % total range shot 2
        f.bedDepth = f.range(f.bn);
        g.bedDepth = g.range(g.bn);

        disp(['Shot 1: new bed found at ' num2str(f.bedDepth) ' m'])
        disp(['Shot 2: new bed found at ' num2str(g.bedDepth) ' m'])
    elseif bed_iteration==4
        bed_align_done = true;
        disp(['A bed shift of ' int2str(bnDiff) ' steps, respectively ' num2str(bnDiff*dr) 'm, was found.'])
        disp(['As this exceeds the expectation limits and all methods to find the bed were tested, no bed shift is undertaken.'])
    elseif ~bed_align_done
        disp(['A bed shift of ' int2str(bnDiff) ' steps, respectively ' num2str(bnDiff*dr) 'm, was found.'])
        disp(['As this exceeds the expactation limits, another method to find the bed is tested.'])
    end
end

%% Estimate range error
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
f.rangeError = fmcw_phase2range(f.phaseStdError,f.lambdac,f.rangeCoarse,f.K,f.ci);
g.rangeError = fmcw_phase2range(g.phaseStdError,g.lambdac,g.rangeCoarse,g.K,g.ci);    
    
%% ALIGN COARSE: Co-register profile segments (amplitude xcorr)
% Cross correlate internals gi segments to get vertical shift as a function of range
% xcor in big chunks to get minimum chance of wrapping errors.
AC.maxOffset = cfg.maxStrain*(maxDepth-cfg.minDepth); %
maxlag = ceil(AC.maxOffset/dr); % max bin lags
binStart = [cfg.minDepth:cfg.coarseStep:maxDepth-cfg.coarseChunkWidth]; % measure offset over a wider range to plot
[AC.range,AC.dh,AC.lags,AC.ampCor,AC.ampCorProm] = deal(zeros(size(binStart)));
for ii = 1:numel(binStart)
    depthRange = [binStart(ii) binStart(ii)+cfg.coarseChunkWidth];
    fi = find((f.rangeCoarse>=min(depthRange) & f.rangeCoarse<max(depthRange))); % depth bins to use (f)
    [AC.RANGEIND(ii,:),AC.AMPCOR(ii,:),~,AC.LAGS(ii,:)] = fmcw_xcorr(f.specRaw,g.specRaw,fi,maxlag);
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
AC.Exp = fit(AC.range(AC.isGood)',AC.lags(AC.isGood)','exp1');
            
if cfg.doPlotAlignCoarse || cfg.doPlotAll
    figure
    clear ax
    ax(1) = subplottight(3,1,1);
    ii = f.rangeCoarse<maxDepth*1.1;
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

% complex correlation: xcor in small chunks to get good depth resolution
% this also gives us the AF.coherence of the segments
binStart = [cfg.minDepth:cfg.fineStep:maxDepth - cfg.fineChunkWidth]; % Noted: added by Falk, take intermediate steps to consider good reflectors at chunk-boundaries
%OffsetRange = AC.maxOffset;
for ii = 1:numel(binStart)
    depthRange = [binStart(ii) binStart(ii) + cfg.fineChunkWidth];
    binDepth = mean(depthRange);
    maxlag = ceil(AC.maxOffset/dr); % max bin lags
    fi = find((f.rangeCoarse>=min(depthRange) & f.rangeCoarse<max(depthRange))); % depth bins to use (f)
    [AF.RANGEIND(ii,:),AF.AMPCOR(ii,:),AF.COR(ii,:),AF.LAGS(ii,:),AF.PE(ii,:),AF.PSE(ii,:)] = fmcw_xcorr(f.specRaw,g.specRaw,fi,maxlag,f.phaseStdError,g.phaseStdError,cfg.p);
    AF.RANGE(ii,:) = interp1(1:numel(f.rangeCoarse),f.rangeCoarse,AF.RANGEIND(ii,:));
    if cfg.doUseCoarseOffset % Define the bin lag from the coarse correlation
        % Get coarse offsets at bin centres
        if strcmpi(cfg.CoarseOffsetMethod,'polynomial')
            AC.dhInterp = dr*polyval(AC.P,binDepth); % generate smoothed lags
        elseif strcmpi(cfg.CoarseOffsetMethod,'exponential')
            AC.dhInterp = dr*AC.Exp(binDepth); % generate smoothed lags
        elseif strcmpi(cfg.CoarseOffsetMethod,'interpolate')
            AC.dhInterp = interp1(AC.range,AC.dh,binDepth,'linear','extrap');
        else
            disp('No coarse offset method defined.')
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
% redefine best lags as those closest to zero phase diff
dphi = 2*pi*f.f0/(f.B*cfg.p); % phase change over 1 bin (brennan eq 16)
dphi = dphi/2; % Falk: Added to get actual fit.
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
    if isempty(pki)
        [~,pki] = max(-abs(angle(AF.COR(ii,:))));
    end
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

igood = AF.mciu;

% Extract values at chosen offsets
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

    plot3(range,AF.lagvec(AF.mciu),pi*ones(size(range)),'k.','markersize',18); % unwrapped
    h(3) = plot3(range,AF.lagvec(AF.mciu),pi*ones(size(range)),'g.','markersize',10); % unwrapped
    plot3(range(si),AF.lagvec(AF.mci(si)),pi,'k.','markersize',40); % best amp
    h(1) = plot3(range(si),AF.lagvec(AF.mci(si)),pi,'m.','markersize',30); % best amp

    legend(h,{'start (best amp-cor)','best amp-cor lag','unwrapped'})
    title('Phase difference from x-cor')
%    set(gcf,'pos',[792 527 565 447])
    
    %keyboard
end

if numWavelengthError>(0.5*numel(range))
    disp('Warning: Integer half-wavelength ambiguity detected - check manually')
    %keyboard
elseif numWavelengthError>(0.1*numel(range))
    disp(['Warning ' int2str(numWavelengthError) ' wavelength errors detected in ' int2str(numel(range)) ' points,'])
end

%% FIT: Estimate internal strain rate
% linear fit to range diff of internals to estimate vertical strain

% Apply criteria to select good reflectors
AC.ampCorInterp = interp1(AC.range,AC.ampCor,range,'nearest','extrap'); % resample coarse amplitude correlation at fine range points
sr_pt.gi = find(range>=10 & ...                %dhe<1 & ...
                range<=maxDepth & ...
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
       h = fillXRange([cfg.firnDepth maxDepth],'facecol',[0.6 0.6 0.6],'facealpha',0.3,'edgecol','none');
       set(gca,'xlim',[0 1.2*maxDepth])
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
       h = fillXRange([cfg.firnDepth maxDepth],'facecol',[0.6 0.6 0.6],'facealpha',0.3,'edgecol','none');
       set(gca,'xlim',[0 1.2*maxDepth])
       box on
       title('Strain fitting')

       ax(2) = subplot(2,1,2);
       plotrange = [0 sr_fit.midRange 1.2*maxDepth];
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

f.phaseCor = angle(f.specCor);
f.phaseRaw = angle(f.specRaw);
f.phaseCor = f.phaseCor + pi;
f.phaseRaw = f.phaseRaw + pi;
f.cumphaseCor = cumsum(f.phaseCor);
f.cumphaseRaw = cumsum(f.phaseRaw);

g.phaseCor = angle(g.specCor);
g.phaseRaw = angle(g.specRaw);
g.phaseCor = g.phaseCor + pi;
g.phaseRaw = g.phaseRaw + pi;
g.cumphaseCor = cumsum(g.phaseCor);
g.cumphaseRaw = cumsum(g.phaseRaw);



 figure(6)
plot(f.range, f.phaseCor)
hold on
plot(f.range, f.phaseRaw)
plot(g.range, g.phaseCor)
plot(g.range, g.phaseRaw)
legend('1','2','3','4')

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







