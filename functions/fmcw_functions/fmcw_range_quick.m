function [Rcoarse,spec] = fmcw_range_quick(vdat,p,maxrange,winfun)

% [Rcoarse,Rfine,spec_cor,spec] = fmcw_range(vdat,p,maxrange,winfun)
%
% Phase sensitive processing of FMCW radar data based on Brennan et al. 2013
%
% Based on Paul's scripts but following the nomenclature of:
% "Phase-sensitive FMCW radar imaging system for high precision Antarctic
% ice shelf profile monitoring"
% Brennan, Lok, Nicholls and Corr, 2013
%
% Summary: converts raw FMCW radar voltages into a range for 
%
% input args: 
% vdat = structure containing shot metadata
% vdat.vif = data matrix of vltages as digitised with one chirp per row size(nchirps,nsamples)
% p = pad factor (i.e. level of interpolation to use during fft)
% maxrange = maximum range to crop output to
% winfun = window function handle (defaults to blackman)
%
% outputs:
% Rcoarse = range to bin centres (m)
% Rfine = range to reflector from bin centre (m)
% spec_cor = spectrum corrected. positive frequency half of spectrum with
% ref phase subtracted. This is the complex signal which can be used for
% cross-correlating two shot segements.

% Craig Stewart
% 2013 April 24
% Modified frequencies 10 April 2014
% Removed spec_cor & range_fine, Falk Oraschewski 15 Dec 2022

if nargin < 3
    maxrange = 2000; %m (range to crop output to)
end
if nargin < 4
    winfun = @blackman; % default to blackman window
end

% Ensure odd number of samples per chirp (to get phase centering right)
%vdat = fmcw_burst_make_length_odd(vdat);

% Extract variables from structure to make it readable
%fs = vdat.fs;
%T = vdat.T;
B = vdat.B;
%K = vdat.K;
ci = vdat.ci;
%fc = vdat.fc;
%f0 = vdat.f0;
%lambdac = vdat.lambdac;

% Processing settings
N = size(vdat.vif,2);
xn = round(0.5*(N));
[nchirps,N] = size(vdat.vif);

% Measure the sampled IF signal: FFT to measure frequency and phase of IF
%deltaf = 1/(T*p); % frequency step of FFT
%f = [0:deltaf:fs/2-deltaf]; % frequencies measured by the fft - changed 16 April 2014, was %f = [0:deltaf:fs/2]; 
%Rcoarse = f*ci*T/(2*B); % Range at the centre of each range bin: eq 14 (rearranged) (p is accounted for inf)
%Rcoarse = [0:1/p:T*fs/2-1/p]*ci/(2*B); % Range at the centre of each range bin: eq 14 (rearranged) (p is accounted for inf)

nf = round((p*N)/2 - 0.5); % number of frequencies to recover
%nf = length(f); % changed from above 2014/5/22
%nf = length(Rcoarse); 
win = window(winfun,N); %chebwin(N);  %rectwin(N); %

%% Loop through for each shot in burst
% Calculate phase of each range bin centre for correction
n = (0:nf - 1)';
Rcoarse = transpose(n*ci/(2*B*p));
n = find(Rcoarse<=maxrange,1,'last');
Rcoarse = Rcoarse(1:n); % preallocate
spec = zeros(nchirps,n); % preallocate
for ii = 1:nchirps
    vif = vdat.vif(ii,:);
    vif = vif-mean(vif); % de-mean
    vif = win.*vif.'; % windowed
    %vif = [vif; zeros((p-1)*N,1)]; % zero padded to length p*N
    vifpad = zeros(p*N,1);
    vifpad(1:length(vif)) = vif;
    vifpad = circshift(vifpad,-xn); % signal time shifted so phase centre at start
    %plot(vifpad), keyboard
    fftvif = fft(vifpad); % fft
    fftvif = fftvif(1:n);
    spec(ii,:) = fftvif.*(sqrt(2*p)/length(vifpad))./rms(win); % scale for window, crop for useful range and scale for padding 
end

