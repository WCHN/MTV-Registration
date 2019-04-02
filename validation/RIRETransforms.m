function RIRETransforms

%--------------------------------------------------------------------------
% Directories
%--------------------------------------------------------------------------

dir_code = '../code';
addpath(dir_code);

DirData            = '../data/rire/images';
DirTransformations = '../data/rire/transforms';
% DirTemp            = '../Temp/ValidateMTVonRIRE';
DirTemp            = '../Temp/ValidateMTVonRIRE1';

if exist(DirTemp,'dir') == 7, rmdir(DirTemp,'s'); end; mkdir(DirTemp); 

%--------------------------------------------------------------------------
% Settings
%--------------------------------------------------------------------------

DoIT  = true;
DoMTV = true;

% IT flags
OptionsIT.sep      = [4 2]; fprintf('sep = [4 2]\n');
OptionsIT.tol      = 0.5*[0.02 0.02 0.02 0.001 0.001 0.001];
OptionsIT.cost_fun = {'mi','nmi','ecc','ncc'};

% MTV flags
OptionsMTV.speak    = 1;
OptionsMTV.save_mtv = true;

Fix     = 1;
NumChan = 3;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------

Idx = 1:NumChan;
Ref = Fix;
Src = Idx(Idx ~= Ref);

Nii = nifti(spm_select('FPListRec',DirData,'^.*\.nii$'));
Res = struct;

Rtrue = GetRtrue(Nii,Ref,Src,DirTransformations);

if DoIT
    %----------------------------------------------------------------------
    % Register with IT
    %----------------------------------------------------------------------

    for cf=1:numel(OptionsIT.cost_fun) % Loop over cost functions

        flags             = OptionsIT;
        flags.cost_fun    = strtrim(OptionsIT.cost_fun{cf});
        Res.it.cf(cf).nam = flags.cost_fun;              

        % Do registration
        %------------------------------------------------------------------
        Pr     = Nii(Ref).dat.fname;
        R_it   = zeros(4,4,NumChan);
        for s = Src            
            Ps = Nii(s).dat.fname;

            q = my_spm_coreg(Pr,Ps,flags);
            
            R_it(:,:,s) = spm_matrix(q(:)');            
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
            str_err                = [str_err sprintf('e(%i -> %i)=[%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f] | ',s,Ref,er(1),er(2),er(3),er(4),er(5),er(6))];                
        end

        fprintf('%s = %s\n',OptionsIT.cost_fun{cf},str_err);

    end % End loop over IT cost functions
end

if DoMTV
    %----------------------------------------------------------------------
    % Register with MTV
    %----------------------------------------------------------------------    

    % Do registration
    %----------------------------------------------------------------------
    FileNames = char({Nii(1).dat.fname, ...
                      Nii(2).dat.fname, ...
                      Nii(3).dat.fname}); 

%     [~,R_mtv] = spm_mtvcoreg(FileNames,opt_mtv);                         
    var   = load('RegResRIRE-sampdeg0.mat');
    R_mtv = var.RIRE_res.R;        
    clear var
        
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
        str_err          = [str_err sprintf('e(%i -> %i)=[%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f] | ',s,Ref,er(1),er(2),er(3),er(4),er(5),er(6))];                
    end

    fprintf('mtv = %s\n',str_err);    
end

%--------------------------------------------------------------------------
% Look at results
%--------------------------------------------------------------------------

Error = zeros(6,6);

for s=Src
    Rt = Rtrue(:,:,s);    
    qr = spm_imatrix(Rt);
            
    Error(1,:) = Error(1,:) + abs(qr(1:6));
end

for cf=2:5
    for s=Src
        Error(cf,:) = Error(cf,:) + abs(Res.it.cf(cf - 1).err(:,s)');
    end
end

for s=Src
    Error(end,:) = Error(end,:) + abs(Res.mtv.err(:,s)');
end

cf = {'pre','mi ','nmi','ecc','ncc','mtv'};
for i=1:6
    er      = Error(i,:);
    er(4:6) = 180/pi*er(4:6);
    
    str_err = sprintf('[%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f], mn(tr) = %5.2f, mn(rot) = %5.2f', ... 
                      er(1)/2,er(2)/2,er(3)/2,er(4)/2,er(5)/2,er(6)/2,sum(er(1:3))/6,sum(er(4:6))/6);                
    fprintf('%s = %s\n',cf{i},str_err);  
end

return
%==========================================================================

%==========================================================================
function R = GetRtrue(Nii,Ref,Src,DirTransformations)

C = numel(Nii);
R = zeros(4,4,C);

f          = Nii(Ref).dat.fname;
[~,NamFix] = fileparts(f);
ix         = strfind(NamFix,'_');
NamFix     = NamFix(ix(end) + 1:end);

R(:,:,Ref) = eye(4);

for s=Src
    f          = Nii(s).dat.fname;
    [~,NamMov] = fileparts(f);    
    ix         = strfind(NamMov,'_');
    NamMov     = NamMov(ix(end) + 1:end);
    
    FileName = fullfile(DirTransformations,[NamFix '_' NamMov '.standard']);
    R(:,:,s) = TransformationFromFile(FileName);
end
%==========================================================================

%==========================================================================
function T = TransformationFromFile(FileName)
Text            = GetText(FileName);
CornersFromTrue = zeros(3,8);
CornersToTrue   = zeros(3,8);
for i=1:8
    vals = strsplit(Text{15 + i});
    
    CornersFromTrue(1,i) = str2double(vals{2});
    CornersFromTrue(2,i) = str2double(vals{3});
    CornersFromTrue(3,i) = str2double(vals{4});
    
    CornersToTrue(1,i) = str2double(vals{5});
    CornersToTrue(2,i) = str2double(vals{6});
    CornersToTrue(3,i) = str2double(vals{7});
end

X = CornersFromTrue;
Y = CornersToTrue;
X = [X; ones(1,8)];
Y = [Y; ones(1,8)];
T = Y/X;
%==========================================================================

%==========================================================================
function Text = GetText(FileName)
fid  = fopen(FileName,'rt');
Text = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
Text = Text{1};
%==========================================================================