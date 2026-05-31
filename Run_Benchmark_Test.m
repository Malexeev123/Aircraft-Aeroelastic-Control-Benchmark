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
% note that nmhe_nmpc is actually lqr+nmhe/nmpc for coupledfull

% THIS NEEDS TO BE EDITED WITH ACTUAL FOLDER LOC
% Generates the Run Script
sim_init('\home\maxal\Aircraft-Aeroelastic-Control-Benchmark\TestBenchPazy',  ...
    'body_case',body_case, 'sim_case',sim_case,'overwrite', true)




%%
% Executes the Run Script
[success, sim_hist] = sim_run(case_name, body_case,'overwrite', true);

if success
    fprintf('Simulation Succeeded. See log for plot and sim history directory.\n')
else
    fprintf('Simulation Failed. See log for error description.\n')
end


%% Change Toggle for separate post processing using sim_hist
% (To not have to rerun the entire sim)
separate_postProc = true; 

if separate_postProc
    % Extract relevant data from simulation history for post-processing
    t = sim_hist.t;
    x = sim_hist.x;
    cfg = sim_hist.cfg;
    beam = sim_hist.beam;
    aero = sim_hist.aero;
    base = sim_hist.base;
    log = sim_hist.log;
    close all % Just to not have to retype it 
    out = AeroFlex.sim.postProcess(t, x, cfg, beam, aero, base, log);
end