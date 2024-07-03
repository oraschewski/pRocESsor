%%% pRocESsor
% Processor for mobile pRES surveys
%
% Usage:
% 1. Create a config file for your survey.
%       In this configFile, all information required to process the data
%       can be provided, including file paths and processing settings.
% 2. Define the config file below under configFile.
% 3. Run this script.
% ------------------------------------------------------------------------%
% 4. The processed data is stored as set in 'dirProc/fileProc' and a quick
%    view plot is provided.
% 
% Falk Oraschewski
% University of Tübingen
% 04.04.2023
clear
close all

% Add package path
addpath(genpath(pwd))

%% Set configuration file
configFile = config_T1_losar();

%% Execute pRocESsor
% Create config handle
cfg = ConfigHandler(configFile);

% Initiate survey
ProfileSurvey = Survey.createSurvey(cfg);

% Process the survey as set in the config.
ProfileSurvey = ProfileSurvey.processProfile(cfg);

%% Plot processed profile
figure()
imagesc(ProfileSurvey.profileDist, ProfileSurvey.rangeCor, 20*log10(abs(ProfileSurvey.imgProc)))
xlabel('Distance (m)')
ylabel('Depth (m)')
axis equal
%xlim([0, 165.1])
%ylim([0, 100])
cb = colorbar;
cb.Label.String = 'Power (dB)';
colormap(gray)
clim([-105, -40])
