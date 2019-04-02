function AnalyseResultsRIRE(Res)
%
% INPUT:
%
% Res.IT{patient}{reg}{1.PthFrom}
%                     {2.PthTo}
%                     {3.Error}{cf}{1.Name}
%                                  {2.MeanAbsErr}
%                                  {3.Time}
%                     {4.PreMeanAbsErr}
%
% Res.MTV{patient}{reg}{1.PthFrom}
%                      {2.PthTo}
%                      {3.MeanAbsErr}
%                      {4.Time}
%                      {5.PreMeanAbsErr}
%__________________________________________________________________________

ResIT  = Res.IT;
ResMTV = Res.MTV;

% Get errors
Np  = numel(ResIT);
Nc  = numel(ResIT{1}{1}{3});
ser = zeros(1,Nc + 2);
cnt = 0;
for p=1:Np
    Nr = numel(ResIT{p});
    for r=1:Nr
        
        er  = zeros(1,Nc + 2);
        
        er(1) = ResIT{p}{r}{4};
        
        Nc = numel(ResIT{p}{r}{3});
        for c=2:Nc + 1
            er(c) = ResIT{p}{r}{3}{c - 1}{2};
        end
        
        er(c + 1) = ResMTV{p}{r}{3};

        fprintf('p=%i, r=%i        | pre=%7.4f, mi=%7.4f, nmi=%7.4f, ecc=%7.4f, ncc=%7.4f, mtv=%7.4f\n',p,r,er(1),er(2),er(3),er(4),er(5),er(6));
        
        ser = ser + er;
        cnt = cnt + 1;
    end
end

ser = ser/cnt;

fprintf('Results for N=%i | pre=%7.4f, mi=%7.4f, nmi=%7.4f, ecc=%7.4f, ncc=%7.4f, mtv=%7.4f\n',cnt,ser(1),ser(2),ser(3),ser(4),ser(5),ser(6));

return
%==========================================================================