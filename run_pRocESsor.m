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
ColleProfile = Survey.createSurvey(cfg);

% Process the survey as set in the configuration.
ColleProfile = ColleProfile.processProfile(cfg);
