clear;clc;
cd E:\MatlabCodes\Reliability  % modify the path
addpath(genpath(pwd));
% Reliability Calculation through FORM

Funs = FunGen;
x = [1,0.2];
DPt = FindDPt(x);
xd = DPt{1}
ipf = DPt{3}
ibeta = DPt{2}


DPt = FindDPtNw(x);

Nsim = 1e7;
x1 = randn(Nsim,1).*0.4+1;
x2 = exprnd(0.2,Nsim,1);
g = x1 - x2;
pf = sum(g<=0)/Nsim
beta = -norminv(pf)

ipf / pf - 1
ibeta / beta - 1