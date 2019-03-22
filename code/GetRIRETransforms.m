PathTransformation = '/home/mbrud/Data/challenges/RIRE/transformations/ct_PD.standard';

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

X=CornersFromTrue;
Y=CornersToTrue;
X = [X; ones(1,8)];
Y = [Y; ones(1,8)];
Y/X