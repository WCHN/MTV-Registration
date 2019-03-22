function [dir_output,dir_core,dir_aux] = dir2holly(job_mode,job_name,dir_SegModel_code,dir_aux_toolbox)
if strcmpi(job_mode,'qsub')      
    dir_output = fullfile('/data/mbrud/Validations/segmentation-model',job_name);
    
    % Copy most recent code to holly
    dir_core = fullfile(dir_output,'core');
    if exist(dir_core,'dir'), rmdir(dir_core,'s'); end; mkdir(dir_core);   

    copyfile(dir_SegModel_code,dir_core)
    
    dir_aux = fullfile(dir_output,'aux');
    if exist(dir_aux,'dir'), rmdir(dir_aux,'s'); end; mkdir(dir_aux);   

    copyfile(dir_aux_toolbox,dir_aux)    
else
    dir_output = fullfile('./output',job_name);
    dir_core   = '';
    dir_aux    = '';
end
%==========================================================================