function holly = default_holly(job_mode,dir_out,code2holly,mem)
if nargin < 4, mem = '5G'; end

holly               = struct;
holly.mode          = job_mode;
holly.server.ip     = 'holly';
holly.server.login  = 'mbrud';
holly.client.folder = fullfile(dir_out,'cluster');
holly.server.folder = holly.client.folder;
holly.matlab.bin    = '/share/apps/matlabR2018a';
holly.matlab.addsub = code2holly;
holly.clean         = false;
holly.clean_init    = true;
holly.verbose       = true;
holly.job.mem       = mem;
holly.job.est_mem   = false;
holly.job.sd        = 0.2;
%==========================================================================