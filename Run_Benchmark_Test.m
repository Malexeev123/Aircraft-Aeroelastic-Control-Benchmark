%% Wing only and Coupled Full Aircraft Run Script
%
%   --------------------------------------------
%   Script Set Written by Maxim Alexeev
%   Revised:
%   --------------------------------------------
% Body Case 

close all; clc;


addpath(genpath('MatlabFlex\')); addpath(genpath('TestBenchPazy\'));
addpath(genpath('plots\'));

% Run Settings
case_name = 'pazy_krylov_ROM';

body_case   = 'wingOnly'    ; % wingOnly (default) | coupledfull
sim_case    = 'nmhe_nmpc'   ; % openloop (default) | nmhe_nmpc


% Generates the Run Script
sim_init('\home\maxal\Research\TestBenchPazy',  ...
    'body_case',body_case, 'sim_case',sim_case,'overwrite', true)




%%
% Executes the RUn Script
success = sim_run(case_name, body_case,'overwrite', true);




