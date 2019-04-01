function Res = RunExperiment(run,opt)

%----------------------------------------------------------------------
% Get parameters
%----------------------------------------------------------------------

DirRef   = opt.gen.DirRef;
DirTemp  = opt.gen.DirTemp;
NumChan  = opt.gen.NumChan;
Run3D    = opt.gen.Run3D;
DoIT     = opt.gen.DoIT;
DoMTV    = opt.gen.DoMTV;
Speak    = opt.gen.Speak;
DirPrint = opt.gen.DirPrint;

do_sim = opt.sim.do;
bf_mn  = opt.sim.bf_mn;
bf_mx  = opt.sim.bf_mx;
r_mn   = opt.sim.r_mn;
r_mx   = opt.sim.r_mx;
t_mn   = opt.sim.t_mn;
t_mx   = opt.sim.t_mx;
n_mn   = opt.sim.n_mn;
n_mx   = opt.sim.n_mx;
d_mx   = opt.sim.d_mx;
d_deg  = opt.sim.d_deg;
mrg    = opt.sim.mrg;

opt_mtv = opt.mtv;
opt_it  = opt.it;

if Speak, disp('---------------------------------------------------------------'); end

% For reproducible results
rng('default') 
rng(run);    

DirTemp = fullfile(DirTemp,num2str(run));
if exist(DirTemp,'dir') == 7, rmdir(DirTemp,'s'); end; mkdir(DirTemp); 

DirSim = fullfile(DirTemp,'sim');
mkdir(DirSim); 

%----------------------------------------------------------------------
% Simulate misaligned images (2D)
%----------------------------------------------------------------------

Idx = randperm(NumChan);
Ref = Idx(1);
Src = Idx(Idx ~= Ref);

Res         = struct;
Res.sim.ix  = Idx;
Res.sim.src = Src;
Res.sim.off = zeros(3);
Res.sim.rot = zeros(3);

% Noise
NoisePrct   = n_mn + rand*(n_mx - n_mn);
Res.sim.noi = NoisePrct;

bf_scl      = bf_mn + rand*(bf_mx - bf_mn);
Res.sim.bf  = bf_scl;

str_sim = sprintf('n=%4.2f, bf=%4.2f | ',NoisePrct,bf_scl);

Offset       = {[0; 0; 0],[0; 0; 0],[0; 0; 0]};
Rotation     = {[0; 0; 0],[0; 0; 0],[0; 0; 0]};
for s = Src

    % Translation
    t                = t_mn + (sign(randn(3,1)).*rand(3,1))*(t_mx - t_mn);  
    Offset{s}        = t;
    Res.sim.off(:,s) = t;
    str_sim          = [str_sim sprintf('t(%i)=[%6.2f,%6.2f,%6.2f], ',s,t(1),t(2),t(3))];

    % Rotation
    r                = r_mn + (sign(randn(3,1)).*rand(3,1))*(r_mx - r_mn);
    Rotation{s}      = r;
    Res.sim.rot(:,s) = r;
    str_sim          = [str_sim sprintf('r(%i)=[%5.2f,%5.2f,%5.2f] | ',s,r(1),r(2),r(3))];   
end

% Down-sampling
ds                      = randi(max(d_mx,1),1);
Res.sim.ds              = ds;
DownSampling            = {[1 1 1],[1 1 1],[1 1 1]};
DownSampling{1}(Idx(1)) = 1/ds;
DownSampling{2}(Idx(2)) = 1/ds;
DownSampling{3}(Idx(3)) = 1/ds;
if ~Run3D
    DownSampling{1}(end) = 1;
    DownSampling{2}(end) = 1;
    DownSampling{3}(end) = 1;
end


Margin            = {[0 0 0],[0 0 0],[0 0 0]};
if mrg > 0
    mrg1              = 2*(randperm(NumChan) - 1);
    Margin{1}(Idx(1)) = mrg;
    Margin{1}(1:3 ~= Idx(1)) = mrg1(1);
    Margin{2}(Idx(2)) = mrg;
    Margin{2}(1:3 ~= Idx(2)) = mrg1(2);
    Margin{3}(Idx(3)) = mrg;
    Margin{3}(1:3 ~= Idx(3)) = mrg1(3);
end

if do_sim
    % Simulate (rotate, translate, etc.)
    [Nii_sim3d,Nii_sim2d,~,~,bb0,Rtrue3d,Rtrue2d] = SimulateData('DirRef',DirRef, ...
                                                                 'NoisePrct',NoisePrct, ...
                                                                 'Offset',Offset, ...
                                                                 'Rotation',Rotation, ...
                                                                 'BiasFieldScl',bf_scl,...
                                                                 'DirSim',DirSim, ...
                                                                 'DownSampling',DownSampling, ...
                                                                 'DownSamplingDeg',d_deg, ...
                                                                 'Margin',Margin, ...
                                                                 'DirPrint',DirPrint);
else
    % Get non-modified input image
    [Nii_sim3d,Nii_sim2d,~,~,bb0,Rtrue3d,Rtrue2d] = SimulateData('DirRef',DirRef,'DirSim',DirSim); fprintf('OBS: No simulation!\n')
end

if Run3D
    Nii_sim = Nii_sim3d;
    Rtrue   = Rtrue3d;
else
    Nii_sim = Nii_sim2d;
    for c=1:NumChan
        bb0{c}(:,end) = 1;
    end
    Rtrue   = Rtrue2d;
end

if Speak, fprintf('r=%i | ds=%i, %s\n',run,ds,str_sim); end

if DoMTV
    %----------------------------------------------------------------------
    % Register with MTV
    %----------------------------------------------------------------------    

    % Do registration
    %----------------------------------------------------------------------
    FileNames = char({Nii_sim(1).dat.fname, ...
                      Nii_sim(2).dat.fname, ...
                      Nii_sim(3).dat.fname}); 

    tstart = tic;

    [~,R_mtv] = spm_mtvcoreg(FileNames,opt_mtv);                

    telap     = toc(tstart);    
    Res.mtv.t = telap;            

    if 0
        spm_check_registration(FileNames);
    end

    % Check error
    %----------------------------------------------------------------------    

    Res.mtv.err = zeros(6,3);            
    Rr          = R_mtv(:,:,Ref);
    str_err     = '';
    
    for s = Src    
        Rsrc = R_mtv(:,:,s);
        Rsrc = Rsrc/Rr;
        Rt   = Rtrue(:,:,s);
        qs   = spm_imatrix(Rsrc);
        qr   = spm_imatrix(Rt);
            
        er               = qr(1:6)' - qs(1:6)';     
        Res.mtv.err(:,s) = er;
        str_err          = [str_err sprintf('e(%i)=[%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f] | ',s,er(1),er(2),er(3),er(4),er(5),er(6))];                
    end

    if Speak, fprintf('r=%i | mtv => t=%7.1f | %s\n',run,telap,str_err); end
    
end

if DoIT
    %----------------------------------------------------------------------
    % Register with IT
    %----------------------------------------------------------------------

    for cf=1:numel(opt_it.cost_fun) % Loop over cost functions

        flags             = opt_it;
        flags.cost_fun    = strtrim(opt_it.cost_fun{cf});
        Res.it.cf(cf).nam = flags.cost_fun;              

        % Do registration
        %------------------------------------------------------------------
        Pr     = Nii_sim(Ref).dat.fname;
        R_it   = zeros(4,4,NumChan);
        tstart = tic;
        for s = Src            
            Ps = Nii_sim(s).dat.fname;

            if Run3D
                q = my_spm_coreg(Pr,Ps,flags);
            else
                q = spm_coreg_2d(Pr,Ps,flags);           
                q = [q(1) q(2) 0 0 0 q(3)];
            end
            
            R_it(:,:,s) = spm_matrix(q(:)');            
        end
        telap           = toc(tstart);    
        Res.it.cf(cf).t = telap;

        if 0
            spm_check_registration(FileNames);
        end

        % Check error
        %------------------------------------------------------------------     

        Res.it.cf(cf).err = zeros(6,3);         
        str_err           = '';
        
        for s = Src        
            Rsrc = R_it(:,:,s);
            Rt   = Rtrue(:,:,s);
            qs   = spm_imatrix(Rsrc);
            qr   = spm_imatrix(Rt);

            er                     = qr(1:6)' - qs(1:6)';           
            Res.it.cf(cf).err(:,s) = er;
            str_err                = [str_err sprintf('e(%i)=[%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f] | ',s,er(1),er(2),er(3),er(4),er(5),er(6))];                
        end

        if Speak, fprintf('r=%i | %s => t=%7.1f | %s\n',run,opt_it.cost_fun{cf},telap,str_err); end

    end % End loop over IT cost functions

end

rmdir(DirTemp,'s');
%==========================================================================