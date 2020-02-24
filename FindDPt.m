function DPt = FindDPt(x0)
% Find the Design Point using the Rackwitz-Fuessiler iteration
% Assuming: (x1,x2,..,xn) are Independent Variable
% x0 - Intitial Search Point
MaxI = 1e5;     % Maximum Number of Iterations
Conv = 1e-5;    % Convergence Criteria

Funs = FunGen;
u0 = norminv(Funs.fcdf(x0));
u = u0;

for i=1:1:MaxI
    gu = Funs.gfun(Funs.finv(normcdf(u)));
    dgu = Funs.gdrv(Funs.finv(normcdf(u)));
    lamda = (gu - dgu'*u)/(dgu'*dgu);
    
    u = -lamda*dgu;
    beta = sqrt(u'*u);
    
    if all(abs((u-u0)./u0)<Conv)  % Converged
        break;
    end
    u0 = u;
end
% Iteration Solving

x = Funs.finv(normcdf(u));
pf = normcdf(-beta);
DPt = {x,beta,pf,i};

if i == MaxI
    error('Not Covenverged');
end


end







