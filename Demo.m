clear; clc;

% Add code
addpath('./code');

% Set options of MTV registration
opt          = struct;
opt.write    = true; 
opt.speak    = 2;
opt.save_mtv = true;

% Simulate some misaligned BrainWeb images (T1, T2, PD)
DirData       = './data/brainweb';
[Nii3d,Nii2d] = SimulateData('DirRef',DirData, ...
                             'Random',true);
                          
% Do registration                         
spm_mtvcoreg(Nii2d,opt); 