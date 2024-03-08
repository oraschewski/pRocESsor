function [vdat, n_burst] = fmcw_burst_select(vdat, data)

% [vdat, n_burst] = fmcw_burst_select(vdat, data)
%
% Selects the burst with the lower noise level.
%
% Falk Oraschewski
% 17.3.2022

% Check if there are more than one attenuator settings
if length(unique(vdat.chirpAtt)) > 1
    error('Trying to compare data with different attenuator settings')
end

N = floor(0.5*size(data,2));
noise = mean(abs(data(:,N:end)),2);
[~,n_burst] = min(noise);

% Selecr burst with lower noise level.
if size(vdat.vif,1)>1
    vdat.vif =  vdat.vif(n_burst,:);
    vdat.Startind =  vdat.Startind(n_burst);
    vdat.Endind =  vdat.Endind(n_burst);
    vdat.chirpNum =  vdat.chirpNum(n_burst);
    vdat.chirpAtt =  vdat.chirpAtt(n_burst);
    vdat.chirpTime = vdat.chirpTime(n_burst);
    vdat.ChirpsInBurst = 1;
    
    vdat.processing = [vdat.processing {[mfilename ': selected burst number: ' int2str(n_burst)]}];
end
    