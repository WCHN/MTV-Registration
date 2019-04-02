function BrainWeb

dir_distributed_toolbox = '/home/mbrud/dev/mbrud/code/matlab/distributed-computing';
dir_code                = '../code';
addpath(dir_code);
addpath(dir_distributed_toolbox);

job_mode = 'for';
job_mem  = '6G';
    
%--------------------------------------------------------------------------
% Main settings
%--------------------------------------------------------------------------

opt_gen         = struct;
if strcmp(job_mode,'qsub')
    opt_gen.DirRef  = '/data/mbrud/populations/original/BrainWebNoiseFree/';
    opt_gen.DirTemp = '/data/mbrud/Holly/MTV-reg/';
else
    opt_gen.DirRef  = '../data/brainweb';
    opt_gen.DirTemp = '../Temp/ValidateMTVonBrainWeb';
end
opt_gen.NumRuns  = 2000;
opt_gen.NumChan  = 3;
opt_gen.Speak    = true;
opt_gen.DirPrint = '';

opt_gen.DoMTV   = true;
opt_gen.DoIT    = true;
opt_gen.Run3D   = false;

% MTV flags
opt_mtv         = struct;
opt_mtv.write   = false; 
opt_mtv.bbpad   = 0;
opt_mtv.bbmni   = false;
opt_mtv.mx_tr   = 100;
opt_mtv.mx_rot  = 30;
opt_mtv.tol_scl = 0.5;
opt_mtv.speak   = 2;

% IT flags
opt_it          = struct;
opt_it.sep      = [4 1];
opt_it.tol      = 0.5*[0.02 0.02 0.02 0.001 0.001 0.001];
opt_it.cost_fun = {'mi ','nmi','ecc','ncc'};

if strcmp(job_mode,'qsub')
    opt_mtv.speak = 0;
    opt_gen.Speak = false;
end

if exist(opt_gen.DirTemp,'dir') == 7, rmdir(opt_gen.DirTemp,'s'); end; mkdir(opt_gen.DirTemp); 

%--------------------------------------------------------------------------
% Simulation settings
%--------------------------------------------------------------------------

opt_sim = struct;

opt_sim.do = true; % If simulation should be performed?

% Bias field
opt_sim.bf_mn = 0;
opt_sim.bf_mx = 2; % 2

% Rotation (degrees)
opt_sim.r_mn = 0;
opt_sim.r_mx = 15; % 15
opt_sim.r_mn = opt_sim.r_mn*pi/180;
opt_sim.r_mx = opt_sim.r_mx*pi/180;

% Translation (mm)
opt_sim.t_mn = 0;
opt_sim.t_mx = 50; % 50

% Downsampling
opt_sim.d_mx  = 6; % 6
opt_sim.d_deg = 0;

% Noise
opt_sim.n_mn = 0;
opt_sim.n_mx = 0.5; % 0.3

% Margin
opt_sim.mrg = 20; % 20

% Holly
if strcmp(job_mode,'qsub')
    dir_code_holly = fullfile(opt_gen.DirTemp,'code');
    copyfile(dir_code,dir_code_holly)
else
    dir_code_holly = dir_code;
end
holly = default_holly(job_mode,opt_gen.DirTemp,dir_code_holly,job_mem);
holly = distribute_default(holly);

if strcmp(job_mode,'qsub')
    opt_gen.DirTemp = fullfile(opt_gen.DirTemp,'Runs');
end

opt     = struct;
opt.gen = opt_gen;
opt.sim = opt_sim;
opt.mtv = opt_mtv;
opt.it  = opt_it;

%--------------------------------------------------------------------------
% Start experiment
%--------------------------------------------------------------------------

disp('---------------------------------------------------------------')
fprintf('Experiment: (%s)\n',opt_gen.DirRef)      
 
Runs                             = cell(1,opt_gen.NumRuns);
for r=1:opt_gen.NumRuns, Runs{r} = r; end

[~,Res] = distribute(holly,'RunExperiment','iter',Runs,opt);       

if opt_gen.Run3D 
    save(['Res-N' num2str(opt_gen.NumRuns) '-3d.mat'],'Res');
else
    save(['Res-N' num2str(opt_gen.NumRuns) '-2d.mat'],'Res');
end

%--------------------------------------------------------------------------
% Analyse results
%--------------------------------------------------------------------------

AnalyseResults(Res);