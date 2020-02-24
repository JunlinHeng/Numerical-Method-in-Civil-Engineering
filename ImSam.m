clear;clc;
cd E:\MatlabCodes\Reliability  % modify the path
addpath(genpath(pwd));
% Importance Sampling using the Design point from FORM

NSIM = 2000;
Funs = FunGen;
DPt = FindDPt([0,0]);
disp('*************** FORM Analysis Results ***************')
fprintf('Iter = %i, ¦Â= %4.2f, Pf = %8.7d\n', DPt{4},DPt{2},DPt{3});
fprintf('Design Point:');disp(DPt{1}');
% disp('*****************************************************')

fprintf('\n')
disp('****************** Crude MCS Results ****************')
x = randn(NSIM,2)';   
g = Funs.gfun(x);
pf = sum(g<0)/NSIM;
beta = -norminv(pf);
fprintf('Sample = %i, ¦Â= %4.2f, Pf = %8.7d\n', NSIM,beta,pf);

fprintf('\n')
disp('*********** Importance MCS at Design Points ***********')
Cov = [1,0;0,1];              % Covariance Matrix
xplus = x + DPt{1};
ig = Funs.gfun(xplus);
weight = ones(1,NSIM);
for i=1:1:NSIM
    JPdf = prod(Funs.fpdf(xplus(:,i)));
    weight(i) = JPdf/Funs.wfun(xplus(:,i),DPt{1},Cov);
end
ifail = (ig<0).*weight;
ipf = sum(ifail)/NSIM;
ibeta = -norminv(ipf);
fprintf('Sample = %i, ¦Â= %4.2f, Pf = %8.7d\n', NSIM,ibeta,ipf);

fprintf('\n')
disp('************** Adaptive Importance MCS  **************')
Abeta = DPt{2};

muY = 1/sqrt(2*pi)/normcdf(-Abeta)*exp(-1*Abeta^2/2);
sigY = sqrt(1+Abeta*muY - muY^2);

Cov = [1,0;0,1].*sigY^2;           % Covariance Matrix
xAdp = x.* sigY + muY;
Ag = Funs.gfun(xAdp);
weight = ones(1,NSIM);
for i=1:1:NSIM
    JPdf = prod(Funs.fpdf(xAdp(:,i)));
    weight(i) = JPdf/Funs.wfun(xAdp(:,i),muY,Cov);
end
Afail = (Ag<0).*weight;
Apf = sum(Afail)/NSIM;
Abeta = -norminv(Apf);
fprintf('Sample = %i, ¦Â= %4.2f, Pf = %8.7d\n', NSIM,Abeta,Apf);



