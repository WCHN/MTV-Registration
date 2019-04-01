function RIRE

dir_code = '/home/mbrud/dev/mbrud/code/matlab/MTV-reg/code';
addpath(dir_code);

%--------------------------------------------------------------------------
% Main settings
%--------------------------------------------------------------------------

DirData            = '../data/rire/images';
DirTransformations = '../data/rire/transforms';
DirTemp            = '../Temp/ValidateMTVonRIRE';

if exist(DirTemp,'dir') == 7, rmdir(DirTemp,'s'); end; mkdir(DirTemp); 

DoIT  = false;
DoMTV = true;

% MTV flags
OptionsMTV.speak    = 1;
OptionsMTV.bbpad    = 0;%[0 0 0;0 60 0];
OptionsMTV.degsamp  = 4; 
% OptionsMTV.fwhm     = 0;
OptionsMTV.bbmni    = false;
% OptionsMTV.samp     = 1;
OptionsMTV.mx_tr    = 100;
OptionsMTV.mx_rot   = 30;
OptionsMTV.tol_scl  = 0.5;
OptionsMTV.save_mtv = true;

% NMI flags
OptionsIT.sep      = [4 1];
OptionsIT.tol      = 0.5*[0.02 0.02 0.02 0.001 0.001 0.001];
OptionsIT.cost_fun = {'mi','nmi','ecc','ncc'};

% Set what images are From and To in the registration
% FromImages   = {'ct','pet'};
% ToImages     = {'mr_PD','mr_T1','mr_T2','mr_PD_rectified','mr_T1_rectified','mr_T2_rectified'};
FromImages   = {'mr_PD','mr_T1','mr_T2'};
ToImages     = {'mr_PD','mr_T1','mr_T2'};
MRISequences = {'mr_PD','mr_T1','mr_T2','mr_PD_rectified','mr_T1_rectified','mr_T2_rectified'};

%----------------------------------------------------------------------
% Start RIRE registration
%----------------------------------------------------------------------

[~,DirPatients] = spm_select('FPList',DirData);
NumDirPatients  = size(DirPatients,1);
ResultsIT       = cell(1,NumDirPatients);
ResultsMTV      = cell(1,NumDirPatients);

for d1=1:NumDirPatients
    
    fprintf('----------------------------------------------------------\n');
            
    DirPatient  = strtrim(DirPatients(d1,:));
    PatientName = strsplit(DirPatient,filesep);
    PatientName = PatientName{end};                    
    
    [PthImages,DirImages] = spm_select('FPListRec',DirPatient,'^.*\.nii$');
    NumPthImages          = size(PthImages,1);
    NumDirImages          = size(DirImages,1);            
    
    if DoIT
        %----------------------------------------------------------------------
        % Information theory registration -> Cycle over images
        %----------------------------------------------------------------------    
        
        ResultsIT{d1}  = cell(1,3);

        for d2=1:NumDirImages

            DirFrom     = strtrim(DirImages(d2,:));
            DirFromName = strsplit(DirFrom,filesep);
            DirFromName = DirFromName{end};

            if any(strcmp(FromImages,DirFromName))

                for d3=1:NumDirImages

                    DirTo    = strtrim(DirImages(d3,:));
                    DirToNam = strsplit(DirTo,filesep);
                    DirToNam = DirToNam{end};

                    if any(strcmp(ToImages,DirToNam)) && ~strcmp(DirFromName,DirToNam)

                        PathFrom           = strtrim(spm_select('FPList',DirFrom,'^.*\.nii$'));                    
                        PathTo             = strtrim(spm_select('FPList',DirTo,'^.*\.nii$'));                    
                        PathTransformation = fullfile(DirTransformations,[DirFromName '_' DirToNam(4:end) '.standard']);

                        ResultsCost = cell(1,numel(OptionsIT.cost_fun));

                        for f=1:numel(OptionsIT.cost_fun)

                            ResultsCost{f} = cell(1,3);

                            flags          = OptionsIT;
                            flags.cost_fun = flags.cost_fun{f};

                            fprintf('Subject %s: %s registering from %s to %s... ',PatientName,flags.cost_fun,DirFromName,DirToNam);
                            TimeStart = tic;
                            q         = spm_coreg(PathTo,PathFrom,flags);                   
                            TimeElap  = toc(TimeStart);   
                            fprintf('done!\n')
                                                                           
                            ResultsCost{f}{1} = flags.cost_fun;
                            ResultsCost{f}{2} = error_it(PathFrom,PathTo,q,PathTransformation);
                            ResultsCost{f}{3} = TimeElap;
                        end

                        ResultsIT{d1}{1}{end + 1} = PathFrom;
                        ResultsIT{d1}{2}{end + 1} = PathTo;
                        ResultsIT{d1}{3}{end + 1} = ResultsCost;

                    end
                end        
            end
        end
    end
    
    if DoMTV
        %----------------------------------------------------------------------
        % MTV registration -> Give all images as input
        %----------------------------------------------------------------------

        ResultsMTV{d1} = cell(1,4);
        
        % Collect modalities
        Niis       = nifti(PthImages);
        Modalities = cell(1,NumPthImages);
        for i=1:NumPthImages

            PthImage = strtrim(PthImages(i,:));
            Modality = strsplit(PthImage,filesep);
            Modality = Modality{end - 1};

            if any(strcmp(Modality,MRISequences)) 
                Modalities{i} = 'mri';
            else
                Modalities{i} = Modality;
            end
        end

        % Register
        flags       = OptionsMTV;
        flags.mod   = Modalities;    

        fprintf('Subject %s: registering all with mtv... ',PatientName);
        TimeStart     = tic;
        
        [~,R_mtv,Mmu] = spm_mtvcoreg(Niis,flags);

%         var   = load('R_mtv.mat');
%         R_mtv = var.R_mtv;
%         var   = load('Mmu.mat');
%         Mmu   = var.Mmu;
        
        TimeElap      = toc(TimeStart);   
        fprintf('done!\n')
        
        % Write all registered images to a temp directory
        Dir = fullfile(DirTemp,'All');
        if exist(Dir,'dir') == 7, rmdir(Dir,'s'); end; mkdir(Dir); 

        for c=1:NumDirImages
            PthImage = Niis(c).dat.fname;
            copyfile(PthImage,Dir)
            [~,nam,ext] = fileparts(PthImage);
            nPathFrom   = fullfile(Dir,[nam ext]);
            Mf          = spm_get_space(nPathFrom);
            spm_get_space(nPathFrom, Mmu\(R_mtv(:,:,c)\Mf));
        end
        spm_check_registration(spm_select('FPList',Dir,'^.*\.nii$'));
        
        % Compare with groundtruth
        for d2=1:NumDirImages

            DirFrom     = strtrim(DirImages(d2,:));
            DirFromName = strsplit(DirFrom,filesep);
            DirFromName = DirFromName{end};

            if any(strcmp(FromImages,DirFromName))

                for d3=1:NumDirImages

                    DirTo    = strtrim(DirImages(d3,:));
                    DirToNam = strsplit(DirTo,filesep);
                    DirToNam = DirToNam{end};

                    if any(strcmp(ToImages,DirToNam)) && ~strcmp(DirFromName,DirToNam)

                        PathFrom           = strtrim(spm_select('FPList',DirFrom,'^.*\.nii$'));                    
                        PathTo             = strtrim(spm_select('FPList',DirTo,'^.*\.nii$'));                    
                        PathTransformation = fullfile(DirTransformations,[DirFromName '_' DirToNam(4:end) '.standard']);

                        RFrom = R_mtv(:,:,d2);
                        RTo   = R_mtv(:,:,d3);                                                                                                                     
                        
                        ResultsMTV{d1}{1}{end + 1} = PathFrom;
                        ResultsMTV{d1}{2}{end + 1} = PathTo;                
                        ResultsMTV{d1}{3}{end + 1} = error_mtv(PathFrom,PathTo,RFrom,RTo,Mmu,PathTransformation);
                        ResultsMTV{d1}{4}{end + 1} = TimeElap;
                    end
                end        
            end
        end
    end
end
%==========================================================================

%==========================================================================
function er = error_mtv(PathFrom,PathTo,RFrom,RTo,Mmu,PathTransformation)

% Get ground-truth
[CornersFromTrue,CornersToTrue] = get_true(PathTransformation);

% Parameters
NiiFrom = nifti(PathFrom);
MatFrom = NiiFrom.mat;
DimFrom = NiiFrom.dat.dim;

NiiTo = nifti(PathTo);
MatTo = NiiTo.mat;
DimTo = NiiTo.dat.dim;

% CornersFromEst
IxCorners      = [1 5 3 7 2 6 4 8];
CornersFromEst = [	1    1    1    1
                    1    1    DimFrom(3) 1
                    1    DimFrom(2) 1    1
                    1    DimFrom(2) DimFrom(3) 1
                    DimFrom(1) 1    1    1
                    DimFrom(1) 1    DimFrom(3) 1
                    DimFrom(1) DimFrom(2) 1    1
                    DimFrom(1) DimFrom(2) DimFrom(3) 1]';  
CornersFromEst = CornersFromEst(:,IxCorners);
CornersFromEst = MatFrom*CornersFromEst;

% round(CornersFromEst(1:3,:) - CornersFromTrue,4)

% CornersToEst
r1      = 0;
r2      = 0;
r3      = 180;
MatFlip = spm_matrix([0 0 0 r1*pi/180 r2*pi/180 r3*pi/180]);

CornersToEst = MatFlip*(RTo*(inv(RFrom)*CornersFromEst));

% Error
er = CornersToEst(1:3,:) - CornersToTrue;
%==========================================================================

%==========================================================================
function er = error_it(PathFrom,PathTo,q,PathTransformation)

% Get ground-truth
[CornersFromTrue,CornersToTrue] = get_true(PathTransformation);

% Parameters
NiiFrom = nifti(PathFrom);
MatFrom = NiiFrom.mat;
DimFrom = NiiFrom.dat.dim;
VoxFrom = sqrt(sum(MatFrom(1:3,1:3).^2));

NiiTo = nifti(PathTo);
MatTo = NiiTo.mat;
DimTo = NiiTo.dat.dim;

MatRig = spm_matrix(q);

% CornersFromEst
IxCorners      = [1 5 3 7 2 6 4 8];
CornersFromEst = [	1    1    1    1
                    1    1    DimFrom(3) 1
                    1    DimFrom(2) 1    1
                    1    DimFrom(2) DimFrom(3) 1
                    DimFrom(1) 1    1    1
                    DimFrom(1) 1    DimFrom(3) 1
                    DimFrom(1) DimFrom(2) 1    1
                    DimFrom(1) DimFrom(2) DimFrom(3) 1]';  
CornersFromEst = CornersFromEst(:,IxCorners);
CornersFromEst = MatFrom*CornersFromEst;

% round(CornersFromEst(1:3,:) - CornersFromTrue,4)

% CornersToEst
r1      = 0;
r2      = 0;
r3      = 180;
MatFlip = spm_matrix([0 0 0 r1*pi/180 r2*pi/180 r3*pi/180]);

CornersToEst = MatFlip*(inv(MatRig)*CornersFromEst);

% Error
er = CornersToEst(1:3,:) - CornersToTrue;
%==========================================================================

%==========================================================================
function [CornersFromTrue,CornersToTrue] = get_true(PathTransformation)
fid = fopen(PathTransformation,'rt');
Text = textscan(fid,'%s','Delimiter','\n');
fclose(fid);

Text            = Text{1};
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
%==========================================================================

%==========================================================================
function copy_write_matrix_mtv(PathFrom,PathTo,DirTemp,RFrom,RTo,Mmu,DirFromName,DirToNam)

Dir = fullfile(DirTemp,[DirFromName '_' DirToNam]);

if exist(Dir,'dir') == 7, rmdir(Dir,'s'); end; mkdir(Dir); 

copyfile(PathFrom,Dir)
[~,nam,ext] = fileparts(PathFrom);
nPathFrom   = fullfile(Dir,[nam ext]);
Mf          = spm_get_space(nPathFrom);
spm_get_space(nPathFrom, Mmu\(RFrom\Mf));

copyfile(PathTo,Dir)
[~,nam,ext] = fileparts(PathTo);
nPathTo     = fullfile(Dir,[nam ext]);
Mf          = spm_get_space(nPathTo);
spm_get_space(nPathTo, Mmu\(RTo\Mf));
%==========================================================================