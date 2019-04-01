function AnalyseResults(Res)

var = load('Res-N2000-3d.mat');
Res = var.Res;
clear var

% Implement possibility to limit range of ds,noise,off,rot and bf in
% boxplot

N = numel(Res);
% N = 2000;

val_mtv = get_val(Res,'mtv');
val_mi  = get_val(Res,'mi');
val_nmi = get_val(Res,'nmi');
val_ecc = get_val(Res,'ecc');
val_ncc = get_val(Res,'ncc');

figure(1);
Data = [val_mtv.err.t val_mi.err.t val_nmi.err.t val_ecc.err.t val_ncc.err.t];
subplot(211)
boxplot(Data,'labels',{'MTV','MI','NMI','ECC','NCC'});
set(gca,'yscale','log','ydir','reverse')
title(['Absolute translation error (N=' num2str(N) ')'])
xlabel('method')
ylabel('mm')
grid on

Data = 180/pi*[val_mtv.err.r val_mi.err.r val_nmi.err.r val_ecc.err.r val_ncc.err.r];
subplot(212)
boxplot(Data,'labels',{'MTV','MI','NMI','ECC','NCC'});
set(gca,'yscale','log','ydir','reverse')
title(['Absolute rotation error (N=' num2str(N) ')'])
xlabel('method')
ylabel('degrees')
grid on

% Res{.}.sim.src
% Res{.}.sim.off
% Res{.}.sim.rot
% Res{.}.sim.noi
% Res{.}.sim.bf

% Res{.}.mtv.err
% Res{.}.mtv.t

% Res{.}.it.cf(.).err
% Res{.}.it.cf(.).t
% Res{.}.it.cf(.).nam

% Look at mean off and rot for each subject for different noise levels
% bf on and off as well
% simulate thick-sliced? (6mm)

% %--------------------------------------------------------------------------
% %% One-Sided F-test
% % Test the null hypothesis that the data in x and y comes from 
% % distributions with the same variance, against the alternative that 
% % the population variance of x is greater than that of y.
% %
% % The returned result h = 1 indicates that vartest2 rejects the null 
% % hypothesis at the default 5% significance level, in favor of the 
% % alternative hypothesis that the population variance of x is greater 
% % than that of y.
% 
% for i=1:size(er_mi,1)
%     x = er_mi(i,:)';
%     y = er_mtv(i,:)';
%     vartest2(x,y,'Tail','right')
% end
%==========================================================================

%==========================================================================
function val = get_val(Res,cf)
N         = numel(Res);
val.err.t = [];
val.err.r = [];
val.t     = [];

for n=1:N
    
    src = Res{n}.sim.src;
    
    if strcmp(cf,'mtv')
        val.err.t = [val.err.t; Res{n}.mtv.err([1 2 3],src)];
        val.err.r = [val.err.r; Res{n}.mtv.err([4 5 6],src)];
        val.t     = [val.t      Res{n}.mtv.t];
    elseif strcmp(cf,'mi')
        val.err.t = [val.err.t; Res{n}.it.cf(1).err([1 2 3],src)];
        val.err.r = [val.err.r; Res{n}.it.cf(1).err([4 5 6],src)];
        val.t     = [val.t      Res{n}.it.cf(1).t];
    elseif strcmp(cf,'nmi')
        val.err.t = [val.err.t; Res{n}.it.cf(2).err([1 2 3],src)];
        val.err.r = [val.err.r; Res{n}.it.cf(2).err([4 5 6],src)];
        val.t     = [val.t      Res{n}.it.cf(2).t];        
    elseif strcmp(cf,'ecc')
        val.err.t = [val.err.t; Res{n}.it.cf(3).err([1 2 3],src)];
        val.err.r = [val.err.r; Res{n}.it.cf(3).err([4 5 6],src)];
        val.t     = [val.t      Res{n}.it.cf(3).t];        
    elseif strcmp(cf,'ncc')
        val.err.t = [val.err.t; Res{n}.it.cf(4).err([1 2 3],src)];
        val.err.r = [val.err.r; Res{n}.it.cf(4).err([4 5 6],src)];
        val.t     = [val.t      Res{n}.it.cf(4).t];        
    end
end

val.err.t = abs(val.err.t(:));
val.err.r = abs(val.err.r(:));
val.t     = val.t(:);

val.err.t(val.err.t == 0) = [];
val.err.r(val.err.r == 0) = [];
val.t(val.t == 0)         = [];

% val.err.t = log(val.err.t(:));
% val.err.r = log(val.err.r(:));
%==========================================================================