function RIRECorners

dir_code = '../code';
addpath(dir_code);

%--------------------------------------------------------------------------
% Main settings
%--------------------------------------------------------------------------

DirData            = '../data/rire/images';
DirTransformations = '../data/rire/transforms';
% DirTemp            = '../Temp/ValidateMTVonRIRE';
DirTemp            = '../Temp/ValidateMTVonRIRE1';

if exist(DirTemp,'dir') == 7, rmdir(DirTemp,'s'); end; mkdir(DirTemp); 

DoIT  = true;
DoMTV = true;

% IT flags
OptionsIT.sep      = [4 2]; fprintf('sep = [4 2]\n');
OptionsIT.tol      = 0.5*[0.02 0.02 0.02 0.001 0.001 0.001];
OptionsIT.cost_fun = {'mi','nmi','ecc','ncc'};

% MTV flags
OptionsMTV.speak    = 1;
OptionsMTV.save_mtv = true;

% Set what images are From and To in the registration
% FromImages   = {'ct','pet'};
% ToImages     = {'mr_PD','mr_T1','mr_T2','mr_PD_rectified','mr_T1_rectified','mr_T2_rectified'};

FromImages   = {'mr_PD','mr_T1','mr_T2'};
% FromImages   = {'mr_T1'};
ToImages     = {'mr_PD','mr_T1','mr_T2'};
MRISequences = {'mr_PD','mr_T1','mr_T2','mr_PD_rectified','mr_T1_rectified','mr_T2_rectified'};

%----------------------------------------------------------------------
% Start RIRE registration
%----------------------------------------------------------------------

[~,DirPatients] = spm_select('FPList',DirData);
NumDirPatients  = size(DirPatients,1);

% For storing results
Results     = struct;
Results.IT  = cell(1,NumDirPatients);
Results.MTV = cell(1,NumDirPatients);

for pat=1:NumDirPatients % Loop over patients
    
    fprintf('----------------------------------------------------------\n');
            
    DirPatient  = strtrim(DirPatients(pat,:));
    PatientName = strsplit(DirPatient,filesep);
    PatientName = PatientName{end};                    
    
    [PthImages,DirImages] = spm_select('FPListRec',DirPatient,'^.*\.nii$');
    NumPthImages          = size(PthImages,1);
    NumDirImages          = size(DirImages,1);            
    
    if DoIT
        %----------------------------------------------------------------------
        % Information theory registration -> Cycle over images
        %----------------------------------------------------------------------    
        
        Results.IT{pat} = {};
        CntReg          = 1;
        
        for from=1:NumDirImages % Loop over From images

            DirFrom     = strtrim(DirImages(from,:));
            DirFromName = strsplit(DirFrom,filesep);
            DirFromName = DirFromName{end};

            if any(strcmp(FromImages,DirFromName))
                % From is in set to be registered
                
                for to=1:NumDirImages % Loop over To images
                                        
                    DirTo    = strtrim(DirImages(to,:));
                    DirToNam = strsplit(DirTo,filesep);
                    DirToNam = DirToNam{end};

                    if any(strcmp(ToImages,DirToNam)) && ~strcmp(DirFromName,DirToNam)
                        % From-to-To pair exists

                        Results.IT{pat}{CntReg} = cell(1,4);
                        
                        PthFrom           = strtrim(spm_select('FPList',DirFrom,'^.*\.nii$'));                    
                        PthTo             = strtrim(spm_select('FPList',DirTo,'^.*\.nii$'));                    
                        PthTransformation = GetPathTransformation(DirTransformations,DirFromName,DirToNam);

                        ResultsCost = cell(1,numel(OptionsIT.cost_fun));

                        for cf=1:numel(OptionsIT.cost_fun) % Loop over IT cost functions
                            
                            % Set IT flags
                            flags          = OptionsIT;
                            flags.cost_fun = flags.cost_fun{cf};

                            % Do registration
                            fprintf('Subject %s: %s registering from %s to %s... ',PatientName,flags.cost_fun,DirFromName,DirToNam);
                            TimeStart = tic;
                            q         = spm_coreg(PthTo,PthFrom,flags);                   
                            TimeElap  = toc(TimeStart);   
                            fprintf('done!\n')
                                             
                            % Store results
                            ResultsCost{cf}    = cell(1,3);
                            ResultsCost{cf}{1} = flags.cost_fun;
                            ResultsCost{cf}{2} = error_it(PthFrom,PthTo,PthTransformation,q);
                            ResultsCost{cf}{3} = TimeElap;
                        end

                        Results.IT{pat}{CntReg}{1} = PthFrom;
                        Results.IT{pat}{CntReg}{2} = PthTo;
                        Results.IT{pat}{CntReg}{3} = ResultsCost;
                        Results.IT{pat}{CntReg}{4} = error_it(PthFrom,PthTo,PthTransformation);

                        CntReg = CntReg + 1;
                        
                    end % End loop over IT cost functions                    
                end % End loop over To images                   
            end             
        end % End loop over From images
        
    end
    
    if DoMTV
        %----------------------------------------------------------------------
        % MTV registration -> Give all images as input
        %----------------------------------------------------------------------
        
        % Collect modality names
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

        % Set MTV flags        
        flags     = OptionsMTV;
        flags.mod = Modalities;    

        % Register
        fprintf('Subject %s: registering all with mtv... ',PatientName);
        TimeStart = tic;
        
%         [~,R_mtv,Mmu] = spm_mtvcoreg(Niis,flags);

        var   = load('RegResRIRE-sampdeg0.mat');
        R_mtv = var.RIRE_res.R;
        Mmu   = var.RIRE_res.Mmu;
        clear var
        
        TimeElap = toc(TimeStart);   
        fprintf('done!\n')

        if 0
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
        end
        
        % Compare with groundtruth
        CntReg = 1;
        
        for from=1:NumDirImages

            DirFrom     = strtrim(DirImages(from,:));
            DirFromName = strsplit(DirFrom,filesep);
            DirFromName = DirFromName{end};

            if any(strcmp(FromImages,DirFromName))
                               
                for to=1:NumDirImages                   
                    
                    DirTo    = strtrim(DirImages(to,:));
                    DirToNam = strsplit(DirTo,filesep);
                    DirToNam = DirToNam{end};

                    if any(strcmp(ToImages,DirToNam)) && ~strcmp(DirFromName,DirToNam)

                        Results.MTV{pat}{CntReg} = cell(1,5);
                        
                        PthFrom           = strtrim(spm_select('FPList',DirFrom,'^.*\.nii$'));                    
                        PthTo             = strtrim(spm_select('FPList',DirTo,'^.*\.nii$'));       
                        PthTransformation = GetPathTransformation(DirTransformations,DirFromName,DirToNam);                        

                        RFrom = R_mtv(:,:,from);
                        RTo   = R_mtv(:,:,to);                                                                                                                     
                        
                        Results.MTV{pat}{CntReg}{1} = PthFrom;
                        Results.MTV{pat}{CntReg}{2} = PthTo;                
                        Results.MTV{pat}{CntReg}{3} = error_mtv(PthFrom,PthTo,PthTransformation,RFrom,RTo,Mmu);
                        Results.MTV{pat}{CntReg}{4} = TimeElap;
                        Results.MTV{pat}{CntReg}{5} = error_mtv(PthFrom,PthTo,PthTransformation);
                        
                        CntReg = CntReg + 1;
                        
                    end
                end        
            end
        end
    end
    
end % End loop over patients

% Save results
save('ResultsRIRE.mat','Results')

% Analyse results
AnalyseResultsRIRE(Results);
%==========================================================================

%==========================================================================
function er = error_it(PathFrom,PathTo,PathTransformation,q)
if nargin < 4, q = [0 0 0 0 0 0]; end

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
CornersFromEst = [1    1    1    1
                  1    1    DimFrom(3) 1
                  1    DimFrom(2) 1    1
                  1    DimFrom(2) DimFrom(3) 1
                  DimFrom(1) 1    1    1
                  DimFrom(1) 1    DimFrom(3) 1
                  DimFrom(1) DimFrom(2) 1    1
                  DimFrom(1) DimFrom(2) DimFrom(3) 1]';  
CornersFromEst = CornersFromEst(:,IxCorners);
CornersFromEst = MatFrom*CornersFromEst;

if 0
    round(CornersFromEst(1:3,:) - CornersFromTrue,4)
end

CornersToEst = MatRig\CornersFromEst;
% CornersToEst = GetMatFlip*CornersToEst;

% Corner errors
er = CornersToEst(1:3,:) - CornersToTrue;

% Mean absolute error
er = mean(abs(er(:)));
%==========================================================================

%==========================================================================
function er = error_mtv(PathFrom,PathTo,PathTransformation,RFrom,RTo,Mmu)

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

if 0
    round(CornersFromEst(1:3,:) - CornersFromTrue,4)
end

if nargin > 3
    CornersToEst = RTo*(RFrom\CornersFromEst);
else
    CornersToEst = CornersFromEst;
end
% CornersToEst = GetMatFlip*CornersToEst;

% Corners error
er = CornersToEst(1:3,:) - CornersToTrue;

% Mean absolute error
er = mean(abs(er(:)));
%==========================================================================

%==========================================================================
function MatFlip = GetMatFlip
r1      = 0;
r2      = 0;
r3      = 180;
MatFlip = spm_matrix([0 0 0 r1*pi/180 r2*pi/180 r3*pi/180]);
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

%==========================================================================
function PthTransformation = GetPathTransformation(DirTransformations,DirFromName,DirToNam)
if contains(DirFromName,'mr')
    DirFromName = DirFromName(4:end);
end
if contains(DirToNam,'mr')
    DirToNam = DirToNam(4:end);
end
PthTransformation = fullfile(DirTransformations,[DirFromName '_' DirToNam '.standard']);

if ~isfile(PthTransformation)
    error('~isfile(PthTransformation)')
end
%==========================================================================