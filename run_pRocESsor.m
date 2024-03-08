%%% pRocESsor, for mobile pRES surveys
% This script runs the processing of mobile pRES surveys.
%
% The user only needs to define a configFile which needs to be stored
% under './config'. In this file, the survey type (currently only
% 'profile' is possible), file paths and names and settings for the applied
% processing are defined
% 
% This script loads the config file, initiates the survey object and
% executes the processing.
% 
% Falk Oraschewski, 04.04.2023

% Clear workspace and command window, close open figures
clear
clc
close all

% Add package path
addpath(genpath(pwd))

%% Define configuration file
configFile = config_T1_LOSAR();

%% Execute processing for pRES survey processing
% Handle configuration file
cfg = ConfigHandler(configFile);

% Initiate survey
ProfileSurvey = Survey.createSurvey(cfg);

% Process the survey as set in the configuration.
ProfileSurvey = ProfileSurvey.processProfile(cfg);

%% Plot processe profile
figure()
imagesc(ProfileSurvey.profileDist, ProfileSurvey.rangeCor,20*log10(abs(ProfileSurvey.imgProc)))
xlabel('Distance (m)')
ylabel('Depth (m)')
axis equal
xlim([0, 165.1])
ylim([0, 100])
cb = colorbar;
cb.Label.String = 'Power (dB)';
colormap(gray)
clim([-120, -40])