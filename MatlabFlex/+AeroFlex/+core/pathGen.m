proj  = matlab.project.currentProject;   % works inside any file
projectRoot  = proj.RootFolder;                 % char array with full path
addpath(genpath(projectRoot));                  % recursive, once at start‑up

% addpath(genpath(fullfile(projectRoot,'AeroFlex')));   % packages
addpath(genpath(fullfile(projectRoot,'ctrl/@EstimatorBase')));    % global classes
% savepath
