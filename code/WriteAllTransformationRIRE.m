function WriteAllTransformationRIRE

DirImages          = '/home/mbrud/Data/challenges/RIRE/proc-train/training_001/';
DirTransformations = '/home/mbrud/Data/challenges/RIRE/transformations-mod';
Files              = spm_select('FPList',DirTransformations,'^.*\.standard$');
N                  = size(Files,1);

% Get names of all modalities
AllModalities = {};
for n=1:N    
    FileName = strtrim(Files(n,:));
    Name     = GetModalityName(FileName,'from');
    if ~any(strcmp(AllModalities,Name))
        AllModalities{end + 1} = Name;
    end
    Name     = GetModalityName(FileName,'to');
    if ~any(strcmp(AllModalities,Name))
        AllModalities{end + 1} = Name;
    end
end
    
% Get all transformations
Transforms = containers.Map;
for n=1:N
    FileName         = strtrim(Files(n,:));
    T                = TransformationFromFile(FileName);
    [~,Name]         = fileparts(FileName);    
    Transforms(Name) = T;
    
    Name             = FlipName(Name);
    Transforms(Name) = inv(T);
end

% Defined as from_to
Get = {{'T1_T2','ct_T2','T1_ct'}, ... % T1 to T2
       {'T1_PD','ct_PD','T1_ct'}, ... % T1 to PD
       {'T2_T1','ct_T1','T2_ct'}, ... % T2 to T1
       {'T2_PD','ct_PD','T2_ct'}, ... % T2 to PD
       {'PD_T1','ct_T1','PD_ct'}, ... % PD to T1
       {'PD_T2','ct_T2','PD_ct'}};    % PD to T2
   
R = numel(Get);
for r=1:R       
    GetCorners(DirImages,DirTransformations,Transforms,Get{r});
end
%==========================================================================

%==========================================================================
function [From,To] = GetCorners(DirImages,DirTransformations,Transforms,Get)
Name     = Get{1};
ix       = strfind(Name,'_');
Name     = Name(1:ix - 1);
PathFrom = spm_select('FPListRec',DirImages,['^.*\' Name '.nii$']);

R = Transforms(Get{2})*Transforms(Get{3});  

[From,To] = GetCornersToFrom(PathFrom,R);

FileName = fullfile(DirTransformations,[Get{1} '.standard']);

SetText(FileName,Get{1},From,To);

return
%==========================================================================

%==========================================================================
function SetText(FileName,Names,From,To)

TemplateFileName = '/home/mbrud/Data/challenges/RIRE/transformations-mod/ct_PD.standard';
Text             = GetText(TemplateFileName);

ix       = strfind(Names,'_');
Text{11} = ['From: mr_' Names(1:ix - 1)];
Text{12} = ['To: mr_' Names(ix + 1:end)];

for i=1:8
    Text{15 + i} = sprintf('  %i%13.4f%11.4f%11.4f%11.4f%11.4f%11.4f',i,From(1,i),From(2,i),From(3,i),To(1,i),To(2,i),To(3,i));
end

if isfile(FileName ),delete(FileName); end
fid = fopen(FileName,'w');
for i=1:numel(Text)
    fprintf(fid,Text{i});
    fprintf(fid,'\n');
end
fclose(fid);

return
%==========================================================================

%==========================================================================
function [From,To] = GetCornersToFrom(PathFrom,R)
NiiFrom = nifti(PathFrom);
MatFrom = NiiFrom.mat;
DimFrom = NiiFrom.dat.dim;

r1      = 0;
r2      = 0;
r3      = 180;
MatFlip = spm_matrix([0 0 0 r1*pi/180 r2*pi/180 r3*pi/180]);
% MatFlip = eye(4);

IxCorners = [1 5 3 7 2 6 4 8];
From      = [1    1    1    1
             1    1    DimFrom(3) 1
             1    DimFrom(2) 1    1
             1    DimFrom(2) DimFrom(3) 1
             DimFrom(1) 1    1    1
             DimFrom(1) 1    DimFrom(3) 1
             DimFrom(1) DimFrom(2) 1    1
             DimFrom(1) DimFrom(2) DimFrom(3) 1]';  
From = From(:,IxCorners);
From = MatFlip*(MatFrom*From);

To = R\From;
%==========================================================================

%==========================================================================
function Name = FlipName(Name)
ix   = strfind(Name,'_');
Name = [Name(ix + 1:end) '_' Name(1:ix - 1)];
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

%==========================================================================
function Name = GetModalityName(FileName,Direction)
[~,Name] = fileparts(FileName);
ix       = strfind(Name,'_');
if strcmp(Direction,'from')
    Name = Name(1:ix - 1);
elseif strcmp(Direction,'to')
    Name = Name(ix + 1:end);
end
%==========================================================================