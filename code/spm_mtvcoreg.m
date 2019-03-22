function [q,R,Mmu,Nii_reg] = spm_mtvcoreg(varargin)
% varargin{1}   = Images as Niis or filenames
% varargin{2}   = flags
% varargin{>=3} = Arguments to Powell
%__________________________________________________________________________

if nargin < 3
    if nargin < 2
        flags = struct;
    else
        flags = varargin{2};
    end

    %----------------------------------------------------------------------
    % Parse settings
    %----------------------------------------------------------------------

    % Convergence criteria [tx,ty,tz,rx,ry,rz]
    if ~isfield(flags,'tol'),    flags.tol     = [0.02 0.02 0.02 0.001 0.001 0.001]; end
    % For iteratively applying less smoothing, in unit FWHM, algorithm
    % calls spm_powell numel(flags.fwhm) times
    if ~isfield(flags,'fwhm'),   flags.fwhm    = [12 8 4 0]; end 
    % For iteratively decreasing subsampling of images (f). Needs to be at
    % most as many elements as flags.fwhm (1 == no subsampling)
    if ~isfield(flags,'samp'),   flags.samp    = [3 2 1.5 1]; end
    % Degree to subsample each image
    if ~isfield(flags,'degsamp'),flags.degsamp = 4; end 
    % For iteratively decreasing voxel size of template (mu). Needs to be at
    % most as many elements as flags.fwhm
    if ~isfield(flags,'vxmu'),   flags.vxmu    = [2 1.75 1.5 1]; end
    % Ad-hoc constraints to not move or rotate too much...
    if ~isfield(flags,'mx_tr'),  flags.mx_tr   = 60; end % degrees
    if ~isfield(flags,'mx_rot'), flags.mx_rot  = 15; end % mm
                                 flags.mx_rot  = flags.mx_rot*pi/180;
    % Verbosity levels:
    %                   0 - No verbose
    %                   1 - Displays input and registered for different
    %                   scheds
    %                   2 - 1 + hist fits and subsampled images
    %                   3 - 1 + 2 + realtime view of mtv
    if ~isfield(flags,'speak'),  flags.speak   = 0; end
    % Instead of using maximum bounding-box (bb), uses bb derived from SPM
    % MNI template
    if ~isfield(flags,'bbmni'),  flags.bbmni   = true; end
    % Extra padding of MNI bb
    if ~isfield(flags,'bbpad'),  flags.bbpad   = 0; end 
    % Modality of each input image (e.g., {'MRI','CT'})
    if ~isfield(flags,'mod'),    flags.mod     = {}; end
    % Modify orientation matrices of input images with registration
    % parameters
    if ~isfield(flags,'write'),  flags.write   = false; end 
    % For changing tolerance of algorithm
    if ~isfield(flags,'tol_scl'), flags.tol_scl = 1; end 
end

if nargin > 2
    
    %----------------------------------------------------------------------
    % Step into Powell optimiser
    %----------------------------------------------------------------------
    
    q = optfun(varargin{:});
    return;
end

%--------------------------------------------------------------------------
% Get spm_vol
%--------------------------------------------------------------------------

Nii     = varargin{1};
if ~isa(Nii,'nifti')
    Nii = nifti(Nii);    
end
C       = numel(Nii);
is2d    = numel(Nii(1).dat.dim) == 2;

if isempty(flags.mod)
    flags.mod = repelem({'mri'},1,C);
end

if 0    
    check_obj(Nii,[1 2],12,60,120);
end

%--------------------------------------------------------------------------
% Get lambda from mean brain intensity
%--------------------------------------------------------------------------

lam = zeros(1,C);
nr  = floor(sqrt(C));
nc  = ceil(C/nr);  
for c=1:C
    if flags.speak > 1, figure(664); if c == 1, clf(figure(664)); end; subplot(nr,nc,c); end        

    has_negval = min(Nii(c).dat(:)) < 0;
    
    % Get mean brain intensity                
    if has_negval
        % Fit GMM
        [~,mu_brain] = fit_gmm(Nii(c),flags.speak > 1); 
    else
        % Fit RMM
        [~,mu_brain] = spm_noise_estimate_mod(Nii(c),flags.speak > 1);
    end
    
    lam(c) = 1/double(mu_brain);            
end

%--------------------------------------------------------------------------
% Store final bb parameters
%--------------------------------------------------------------------------

mat0 = zeros(4,4,C);
dm0  = zeros(C,3);
for c=1:C        
    mat0(:,:,c) = Nii(c).mat;
    dmc         = Nii(c).dat.dim;
    dmc         = [dmc 1];
    dm0(c,:)    = dmc(1:3);
end
[Mmu0,dmmu0] = max_bb_orient(mat0,dm0,flags.vxmu(end),flags.bbmni,flags.bbpad);

%--------------------------------------------------------------------------
% Initialise
%--------------------------------------------------------------------------

if is2d
    npar      = 3;
    flags.tol = flags.tol([1 2 6]); % Select x and y translation, and rotation component    
else
    npar      = 6;
end

sc0  = flags.tol(:)';    % Required accuracy
q    = zeros(1,C*npar);

% To make sure that Powell zero centres q
opt_pow          = struct;
opt_pow.mc.do    = true;
opt_pow.mc.npar  = npar;
opt_pow.mc.C     = C;
opt_pow.mc.speak = false;
    
if flags.speak
    % Verbose
    display_results(q,Nii,mat0,Mmu0,dmmu0,lam,npar,1,'Input');
end

%--------------------------------------------------------------------------
% Iterate Powell
%--------------------------------------------------------------------------

for k=1:numel(flags.fwhm)
    
    fwhm = flags.fwhm(k);
    samp = flags.samp(min(k,numel(flags.samp)));
    vxmu = flags.vxmu(min(k,numel(flags.vxmu)));
    
    %----------------------------------------------------------------------
    % Pull out image data, and possibly subsample (if samp > 1) and smooth (if fwhm > 0)
    %----------------------------------------------------------------------
    
    f    = cell(1,C);
    matk = zeros(4,4,C);
    dmk  = zeros(C,3);
    for c=1:C        
        if flags.speak > 1, figure(666); if c == 1, clf(figure(666)); end; subplot(nr,nc,c); end
        
        [fc,matk(:,:,c),dmk(c,:)] = get_img(Nii(c),max(samp,1),fwhm,flags.degsamp);
        f{c}                      = spm_diffeo('bsplinc',fc,[2 2 2 0 0 0]);
        
        if flags.speak > 1, imagesc(fc(:,:,ceil(dmk(c,3)/2))'); axis off xy image; end
    end
    clear fc
    if flags.speak > 1, drawnow; end
        
    %----------------------------------------------------------------------
    % Get mean orientation matrix and dimensions        
    %----------------------------------------------------------------------
    
    [Mmu,dmmu] = max_bb_orient(matk,dmk,vxmu,flags.bbmni,flags.bbpad);
    fovmu      = abs(Mmu*[dmmu'+0.5; 1]-Mmu*[0.5*ones(3,1); 1]);
    fovmu      = fovmu(1:3);

    %----------------------------------------------------------------------
    % Get interpolation grid
    %----------------------------------------------------------------------
    
    x = get_x(dmmu);

    %----------------------------------------------------------------------
    % Initial search values and stopping criterias (changes with decreasing FWHM)
    %----------------------------------------------------------------------
    
    sc            = [];
    for c=1:C, sc = [sc sc0]; end
    iq            = diag(sc*20);
    sc            = max(fwhm/2,1)*sc*flags.tol_scl;
    
    %----------------------------------------------------------------------
    % Start Powell
    %----------------------------------------------------------------------
    
    q = my_spm_powell(q(:),iq,sc,opt_pow,mfilename,f,x,matk,Mmu,dmmu,lam,npar,fovmu,flags,flags.speak > 2);        
    
    if flags.speak
        % Verbose
        display_results(q,Nii,mat0,Mmu0,dmmu0,lam,npar,2,['k=' num2str(k) ', samp=' num2str(samp) ', fwhm=' num2str(fwhm) ', vx=' num2str(vxmu)]);
    end
end

%--------------------------------------------------------------------------
% Prepare output
%--------------------------------------------------------------------------

% Rigid matrices
B = get_rigid_basis;
R = zeros(4,4,C);
for c=1:C
    ix     = (c - 1)*npar + 1:(c - 1)*npar + npar;        
    qc     = q(ix);
    if npar == 3
        qc = [qc(1) qc(2) 0 0 0 qc(3)];
    end
    R(:,:,c)  = spm_dexpm(qc,B); 
end

Nii_reg = nifti;
if flags.write
    % Adjust orientation matrices of input images    
    for c=1:C
        
        FileName = Nii(c).dat.fname;    
        Mf       = Nii(c).mat;
        spm_get_space(FileName, Mmu0\(R(:,:,c)\Mf)); % TODO: double check that this is REALLY the way to set orientation matrix

        Nii_reg(c) = nifti(FileName);
    end    
end
%==========================================================================

%==========================================================================
function o = optfun(q,f,x,mat,Mmu,dmmu,lambda,npar,fovmu,flags,show)
% Computes registration objective function as:
% MTV(a,b,c) - TV(a) - TV(b) - TV(c)
%__________________________________________________________________________

if nargin < 10 || isempty(flags)
    flags.mx_tr  = 50;
    flags.mx_rot = 10*pi/180;
end
if nargin < 11, show = false; end

if outside_fov(fovmu,q,npar) || any(~isfinite(q))
    o = realmax;
    return    
end

% if ~q_ok(q,numel(f),npar,flags)
%     o = 1e10; % TODO: very ad-hoc...
%     return
% end

get_slice = false;
[mtv,tv]  = compute_mtv(q,f,x,mat,Mmu,dmmu,lambda,npar,get_slice,show);
o         = sum(sum(sum(double(mtv)))) + sum( - sum(sum(sum(double(tv),1),2),3));
% o        = mean(mtv(:)) - sum(mean(reshape(tv,[],size(tv,4)),1));
%==========================================================================

%==========================================================================   
function [mtv,tv] = compute_mtv(q,f,x,mat,Mmu,dmmu,lambda,npar,get_slice,show)
if nargin < 9,  get_slice = false; end
if nargin < 10, show      = false; end

B      = get_rigid_basis;
C      = numel(f);
mtv    = single(0);
if nargout > 1 
    tv = zeros([dmmu(1:3) C],'single'); 
end
for c=1:C
    ix = (c - 1)*npar + 1:(c - 1)*npar + npar;
    qc = q(ix);    
    if npar == 3
        qc = [qc(1) qc(2) 0 0 0 qc(3)];
    end
    R  = spm_dexpm(qc,B); 
    
    Mf = mat(:,:,c);
    vx = sqrt(sum(Mf(1:3,1:3).^2));
    M  = Mf\R*Mmu;
    y  = affine_transf(M,x);
        
    if get_slice
        z            = round(dmmu(3)/2);
        [~,Gx,Gy,Gz] = spm_diffeo('bsplins',f{c},y(:,:,z,:),[2 2 2  0 0 0]);  
    else
        [~,Gx,Gy,Gz] = spm_diffeo('bsplins',f{c},y,[2 2 2  0 0 0]);  
    end
    clear y
    
    mtvc = sum((lambda(c)*cat(4,Gx/vx(1),Gy/vx(2),Gz/vx(3))).^2,4);                                                   
    
    if nargout > 1 
        tv(:,:,:,c) = sqrt(mtvc);
    end
    mtv = mtv + mtvc;
end

mtv = sqrt(mtv);

if show
    
    mtv1 = mtv;
    for c=1:C
        mtv1 = mtv1 - tv(:,:,:,c);
    end
    
    figure(667);
    if dmmu(3) == 1
        subplot(111)
        imagesc(mtv1'); axis off xy image;
    else        
        subplot(131)
        mtvc = mtv1(:,:,round(dmmu(3)/2));        
        imagesc(mtvc'); axis off xy image;
        subplot(132)
        mtvc = squeeze(mtv1(:,round(dmmu(2)/2),:));        
        imagesc(mtvc'); axis off xy image;
        subplot(133)
        mtvc = squeeze(mtv1(round(dmmu(1)/2),:,:));        
        imagesc(mtvc'); axis off xy image;
    end
    drawnow;
end
%==========================================================================   

%==========================================================================   
function isoutside = outside_fov(fovmu,q,npar)
isoutside = false;
C         = numel(q)/npar;
margin1   = 0.1;
margin1   = margin1*fovmu;
if npar == 3
    margin1(3) = 0;
end
for c=1:C
    ix     = (c - 1)*npar + 1:(c - 1)*npar + npar;        
    qc     = q(ix);
    if npar == 3
        qc = [qc(1) qc(2) 0 0 0 qc(3)];
    end
    
    if fovmu(1)/2 - qc(1) <= 0 + margin1(1) || fovmu(1)/2 + qc(1) >= fovmu(1) - margin1(1), isoutside = true; break; end
    if fovmu(2)/2 - qc(2) <= 0 + margin1(2) || fovmu(2)/2 + qc(2) >= fovmu(2) - margin1(2), isoutside = true; break; end
    if fovmu(3)/2 - qc(3) <= 0 + margin1(3) || fovmu(3)/2 + qc(3) >= fovmu(3) - margin1(3), isoutside = true; break; end    
end
%==========================================================================   

%==========================================================================
function [mat,dm] = max_bb_orient(amat,adm,vx,bbmni,padding)
% Calculate orientation matrix and dimensions from maximum bounding-box
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 4, bbmni   = true; end
if nargin < 5, padding = 0; end
if numel(vx) == 1, vx = vx*ones([1 3]); end

if numel(padding) == 1, padding = padding*ones(2,3); end

is2d = adm(1,3) == 1;

mn = [ Inf  Inf  Inf]';
mx = [-Inf -Inf -Inf]';
for n=1:size(amat,3)
    dmn = adm(n,:);

    if dmn(3) == 1
        dmn(3) = 0; 
    end

    t = uint8(0:7);
    c = diag(dmn+1)*double([bitshift(bitand(t,bitshift(uint8(1),1-1)),1-1)
                          bitshift(bitand(t,bitshift(uint8(1),2-1)),1-2)
                          bitshift(bitand(t,bitshift(uint8(1),3-1)),1-3)]);
    c  = bsxfun(@plus,amat(1:3,1:3,n)*c,amat(1:3,4,n));
    mx = max(mx,max(c,[],2));
    mn = min(mn,min(c,[],2));
end

mat = spm_matrix(mn-1)*diag([vx 1])*spm_matrix(-[1 1 1]);
dm  = ceil((mat\[mx'+1 1]')');
dm  = dm(1:3);

if bbmni
    V     = struct;
    V.mat = mat;
    V.dim = dm;
    bb0   = world_bb(V);
    mp    = max(bb0) - (max(bb0) - min(bb0))./2;    
    lm    = bb0 - mp;
    
    PthTemplateSPM = fullfile(spm('dir'),'tpm','TPM.nii,');
    V              = spm_vol(PthTemplateSPM);
    nbb            = world_bb(V(1));

    nbb(1,:) = nbb(1,:) - padding(1,:);
    nbb(2,:) = nbb(2,:) + padding(2,:);    
    
    nbb(1,:) = max(nbb(1,:),lm(1,:));
    nbb(2,:) = min(nbb(2,:),lm(2,:));
    
    nbb = mp + nbb;    
    
    if is2d
        nbb(:,3) = 1;
    end
    
    mat = spm_matrix([nbb(1,:) 0 0 0 vx])*spm_matrix([-1 -1 -1]);
    dm  = ceil(mat \ [nbb(2,:) 1]' - 0.1)';
    dm  = dm(1:3);
end

if is2d
    dm(3)    = 1;
    mat(3,4) = 0;
end
%==========================================================================

%==========================================================================
function bb = world_bb(V)
%  world-bb -- get bounding box in world (mm) coordinates

d = V.dim(1:3);
% corners in voxel-space
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
% corners in world-space
tc = V.mat(1:3,1:4)*c;

% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];
%==========================================================================

%==========================================================================
function id = get_id(dm)
id        = cell(3,1);
[id{1:3}] = ndgrid(single(1:dm(1)),single(1:dm(2)),single(1:dm(3)));
id        = cat(4,id{:});
%==========================================================================

%==========================================================================
function x = affine_transf(Affine,x)
dm = size(x);
x  = reshape(x,[prod(dm(1:3)),3]);
x  = x*Affine(1:3,1:3)' + Affine(1:3,4)';
x  = reshape(x,dm);
if dm(3) == 1
    x(:,:,:,end) = 1;
end    
%==========================================================================   

%==========================================================================
function img = smooth_img(img,fwhm,vx)
if nargin < 2, fwhm = 1; end
if nargin < 3, vx   = 1; end

if numel(fwhm) == 1, fwhm = fwhm*ones(1,3); end
if numel(vx) == 1,   vx = vx*ones(1,3); end

if fwhm > 0        
    fwhm = fwhm./vx;            % voxel anisotropy
    s1   = fwhm/sqrt(8*log(2)); % FWHM -> Gaussian parameter

    x  = round(6*s1(1)); x = -x:x; x = spm_smoothkern(fwhm(1),x,1); x  = x/sum(x);
    y  = round(6*s1(2)); y = -y:y; y = spm_smoothkern(fwhm(2),y,1); y  = y/sum(y);
    z  = round(6*s1(3)); z = -z:z; z = spm_smoothkern(fwhm(3),z,1); z  = z/sum(z);

    i  = (length(x) - 1)/2;
    j  = (length(y) - 1)/2;
    k  = (length(z) - 1)/2;

    spm_conv_vol(img,img,x,y,z,-[i,j,k]);   
end
%==========================================================================

%==========================================================================
function display_results(q,Nii,mat,Mmu,dmmu,lambda,npar,plt,title_text)
C = numel(Nii);
f = cell(1,C);
for c=1:C
    f{c} = spm_diffeo('bsplinc',single(Nii(c).dat(:,:,:)),[2 2 2 0 0 0]);
end

x        = get_x(dmmu);
[mtv,tv] = compute_mtv(q,f,x,mat,Mmu,dmmu,lambda,npar);

for c=1:C
    mtv = mtv - tv(:,:,:,c);
end
    
figure(665);
if plt == 1, clf(figure(665)); end

if dmmu(3) == 1
    subplot(1,2,plt);imagesc(mtv'); axis off xy image
    title(title_text)
else
    subplot(2,3,3*(plt - 1) + 1)
    mtvc = mtv(:,:,round(dmmu(3)/2));        
    imagesc(mtvc'); axis off xy image;
    subplot(2,3,3*(plt - 1) + 2)
    mtvc = squeeze(mtv(:,round(dmmu(2)/2),:));        
    imagesc(mtvc'); axis off xy image;
    title(title_text)
    subplot(2,3,3*(plt - 1) + 3)
    mtvc = squeeze(mtv(round(dmmu(1)/2),:,:));        
    imagesc(mtvc'); axis off xy image;
end
drawnow;
%==========================================================================

%==========================================================================
function [noise,mu_brain] = fit_gmm(Nii,speak,K)
if nargin < 2, speak = false; end
if nargin < 3, K     = 2; end
   
% Get image date (vectorised)
f = Nii.dat(:);
f(~isfinite(f)) = [];
f(f == min(f))  = [];
f(f == max(f))  = [];

% Histogram bin voxels
x = min(f):max(f);
c = hist(f(:),x);

% Transpose
x = x';
c = c';

% Intensity distribution hyper-parameters
MU0 = mean(f)*ones([1 K]) + K*randn([1 K]);
b0  = ones([1 K]);
n0  = ones([1 K]);
V0  = ones([1 1 K]);

% Fit GMM
[~,MU,A,PI] = spm_gmm(x,K,c,'BinWidth',1,'GaussPrior',{MU0,b0,V0,n0}, ...
                      'Tolerance',1e-6,'Start','prior','Verbose',0,'Prune',true);

% Compute variance for each class
V = zeros(size(A));
for k=1:size(V,3)
    V(:,:,k) = inv(A(:,:,k));
end

% Get standard deviation of class closest to air (-1000)
sd       = sqrt(squeeze(V));
noise    = min(sd);
[~,ix]   = max(MU);
mu_brain = MU(ix);

if speak
    % Plot histogram + GMM fit
    p = zeros([numel(x) K],'single');
    for k=1:K
        p(:,k) = PI(k)*mvnpdf(x(:),MU(:,k),V(:,:,k));
    end
    sp = sum(p,2) + eps;    
    md = mean(diff(x));
    
    plot(x(:),p,'--',x(:),c/sum(c)/md,'b.',x(:),sp,'r'); 
    drawnow
end
%==========================================================================

%==========================================================================
function [img,mat,dm] = get_img(Nii,samp,fwhm,deg,bc)
if nargin < 2, samp = 1; end
if nargin < 3, fwhm = 0; end
if nargin < 4, deg  = 0; end
if nargin < 5, bc   = 0; end

% Input image properties
img  = Nii.dat(:,:,:);    
mat0 = Nii.mat;
dm0  = size(img);
dm0  = [dm0 1];

if samp == 1 || samp == 0    
    mat = mat0;
    dm  = dm0(1:3);
    
    % Smooth
    vx  = sqrt(sum(mat(1:3,1:3).^2));
    img = smooth_img(img,fwhm,vx);

    img = single(img);
    
    return
end

% Output image properties
vx   = sqrt(sum(mat0(1:3,1:3).^2));
samp = max([1 1 1],ceil(samp*[1 1 1]./vx));
    
D    = diag([1./samp 1]);
mat  = mat0/D;
dm   = floor(D(1:3,1:3)*dm0(1:3)')';

if dm0(3) == 1
    dm(3) = 1;
end

% Make interpolation grid
[x0,y0,z0] = ndgrid(1:dm(1),1:dm(2),1:dm(3));

T = mat0\mat;    

x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

if dm0(3) == 1
    z1 = ones((size(z1)));
end

% Resample
if numel(deg)  == 1, deg  = deg*ones([1 3]);  end
if numel(bc)   == 1, bc   = bc*ones([1 3]);   end

par                 = [deg bc];
C                   = spm_bsplinc(img,par);
img                 = spm_bsplins(C,x1,y1,z1,par);    
img(~isfinite(img)) = 0;

% Smooth
vx  = sqrt(sum(mat(1:3,1:3).^2));
img = smooth_img(img,fwhm,vx);

img = single(img);
%==========================================================================

%==========================================================================
function val = check_obj(Nii,ix,fwhm,mx,stps,samp,deg)
if nargin < 2, ix   = [1 2]; end
if nargin < 3, fwhm = 0;     end
if nargin < 4, mx   = 60;    end
if nargin < 5, stps = 100;   end
if nargin < 6, samp = 1;     end
if nargin < 7, deg  = 0;     end

C    = numel(Nii);
is2d = numel(Nii(1).dat.dim) == 2;

% Get scaling parameters (lambda)
lam = zeros(1,C);
for c=1:C
    [~,mu_brain] = fit_gmm(Nii(c));
    mu           = mean(mu_brain);
    lam(c)       = 1/double(mu);
end

% Get image data (possibly subsample and smooth)
f    = cell(1,C);
matk = zeros(4,4,C);
dmk  = zeros(C,3);
for c=1:C        
    [fc,matk(:,:,c),dmk(c,:)] = get_img(Nii(c),1/max(samp,1),fwhm,deg);
    f{c}                      = spm_diffeo('bsplinc',fc,[2 2 2 0 0 0]);
end
clear fc
    
% Number of registration parameters
if is2d, npar = 3;
else     npar = 6;
end

% Pick two images based on ix
lam = lam(ix);
f   = f(ix);

% Parameter range
par = linspace(-mx,0,stps);
par = [par(1:end - 1) linspace(0,mx,stps)];

% Compute objective for a range of values
val = zeros(1,numel(par));
q   = zeros(1,numel(ix)*npar); 
qr  = zeros(1,6); 
for i=1:numel(par)
    fprintf('%i ',i)
       
    mat = [];
    dm  = [];
    for i1=ix
        mat = cat(3,mat,matk(:,:,i1));
        dm  = [dm; dmk(i1,:)];
    end
    
    qr(1)        = par(i);
    R            = spm_matrix(qr(:)');
    mat(:,:,end) = R*mat(:,:,end);
    
    [Mmu,dmmu] = max_bb_orient(mat,dm,samp);    
    fovmu      = abs(Mmu*[dmmu'+0.5; 1]-Mmu*[0.5*ones(3,1); 1]);
    fovmu      = fovmu(1:3);
    x          = get_x(dmmu);
    
    val(i) = optfun(q,f,x,mat,Mmu,dmmu,lam,npar,fovmu,[],true);
    
    figure(668); plot(par(1:i),val(1:i),'b-x'); drawnow
end
fprintf('Done!\n')
%==========================================================================

%==========================================================================
function ok = q_ok(q,C,npar,flags)
mx_tr  = flags.mx_tr;
mx_rot = flags.mx_rot;
ok     = true;
for c=1:C
    
    ix = (c - 1)*npar + 1:(c - 1)*npar + npar;
    qc = q(ix);    
    if npar == 3
        qc = [qc(1) qc(2) 0 0 0 qc(3)];
    end
    
    if qc(1) > mx_tr  || qc(1) < -mx_tr,  ok = false; return; end
    if qc(2) > mx_tr  || qc(2) < -mx_tr,  ok = false; return; end
    if qc(3) > mx_tr  || qc(3) < -mx_tr,  ok = false; return; end
    if qc(4) > mx_rot || qc(4) < -mx_rot, ok = false; return; end
    if qc(5) > mx_rot || qc(5) < -mx_rot, ok = false; return; end    
    if qc(6) > mx_rot || qc(6) < -mx_rot, ok = false; return; end    
end
%==========================================================================

%==========================================================================
function x = get_x(dmmu)
rng('default') 
rng(1);       
x = get_id(dmmu);  
x = x + rand(size(x),'single');
%==========================================================================    