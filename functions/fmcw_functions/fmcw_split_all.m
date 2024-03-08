function vdats = fmcw_split_all(vdat)
%FMCW_SPLIT_ALL seperates multiple antenna and attenuator combinations in a vdat.
%   In files that are recorded in MIMO mode (multiple input, multiple
%   output), with multiple attenuator settings and with more than one Chirp
%   per Burst, the information from the various combinations is intricated
%   in the file. This function disentagles them and put them out as one
%   structure with a seperate field for each combinations that in the
%   following can be used as if they were all recorded seperately.
%
%   input: original vdat
%   output: vdat-structure containing multiple elements for each 
%   combination.
%
%   Falk Oraschewski,
%   2022 April 04


% Splits a fmcw radar burst into subgroups to seperate chirps recorded
% between different Antenna combinations and attenuator settings.
%
% Falk Oraschewski
% 03.04.2023

% Get number of Tx-Rx-Att combinations
nTx = numel(vdat.TxAnt); % Number of transmitting antennas
nRx = numel(vdat.RxAnt); % Number of receiving antennas
nAtt = vdat.NAttenuators; % Number of attenuator settings

%nX = nTx * nRx; % Number of antenna combinations
nComb =  nTx * nRx * nAtt; % Total number of combinations

nChirp = vdat.ChirpsInBurst; % Number of Chirps in burst

% Define rules about the order of the applied settings
%whichAtt = repmat((1:nAtt),1,nX);
whichTx  = repelem(repelem((1:nTx),1,nRx),1,nAtt);
whichRx  = repelem(repmat((1:nRx),1,nTx),1,nAtt);

%attList = unique(vdat.chirpAtt,'stable');
vdats(1:nComb) = deal(vdat); % Initiate vdats
for ii = 1:nComb
    chirplist = ii:nComb:nChirp;
    %attCurrent = attList(whichAtt(ii));
    %vdats(ii).processing = [vdat.processing {[mfilename ': keeping chirps with att: ' int2str(real(attCurrent)) '+' int2str(imag(attCurrent)) ]}]; % ', chirps: ' mat2str(vdat.chirpNum(chirplist))
    vdats(ii) = fmcw_burst_subset(vdats(ii),chirplist);
    vdats(ii).TxAnt = 0 * vdats(ii).TxAnt;
    vdats(ii).TxAnt(whichTx(ii)) = 1;
    vdats(ii).RxAnt = 0 * vdats(ii).RxAnt;
    vdats(ii).RxAnt(whichRx(ii)) = 1;
end
