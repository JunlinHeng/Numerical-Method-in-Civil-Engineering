function Funs = FunGen
    Funs.gfun = @gfun;      % Limit State Function
    Funs.gdrv = @gdrv;      % Derivate of LSF
    Funs.fpdf = @fpdf;      % PDF of Original Variables
    Funs.fcdf = @fcdf;      % CDF of Original Variables
    Funs.finv = @finv;      % Inverse CDF of Original Variables
    Funs.wfun = @wfun;      % Weight Function for Importance Sampling
end


function g = gfun(x)
% the Limit State Function
    g = 3 - x(1,:) - x(2,:);
end

function drv = gdrv(x)
% Calculate the Partial Derivates from LSF
% x: calculating point
    alph = 1e-5;    % A small value
    drv = 0*x;
    for i =1:1:length(x)
        xinc = x;
        xinc(i) = xinc(i)+alph;
        drv(i) = (gfun(xinc)-gfun(x))./alph;
    end
end

function f = fpdf(x)
% the PDF of the original variables  
%     fx1 = normpdf(x(1),1,0.4);    % x1 ~ N(10,0.4)
%     fx2 = exppdf(x(2),0.2);       % x2 ~ e(0.2)
    fx1 = normpdf(x(1,:),0,1);    % x1 ~ N(0,1)
    fx2 = normpdf(x(2,:),0,1);    % x1 ~ N(0,1)
    f = [fx1;fx2];
end

function f = fcdf(x)
% the CDF of the original variables
%     fx1 = normcdf(x(1),1,0.4);    % x1 ~ N(10,0.4)
%     fx2 = expcdf(x(2),0.2);       % x2 ~ e(0.2)
    fx1 = normcdf(x(1),0,1);    % x1 ~ N(0,1)
    fx2 = normcdf(x(2),0,1);    % x1 ~ N(0,1)
    f = [fx1;fx2];
end

function f = finv(x)
% the inverse CDF of the original variables
%     fx1 = norminv(x(1),1,0.4);    % x1 ~ N(10,0.4)
%     fx2 = expinv(x(2),0.2);       % x2 ~ e(0.2)
    fx1 = norminv(x(1,:),0,1);    % x1 ~ N(0,1)
    fx2 = norminv(x(2,:),0,1);    % x1 ~ N(0,1)
    f = [fx1;fx2];
end

function weight = wfun(x,Mu,Cov)
% General Weight Function for importance sampling (Gassian-based)
% x: variable vector, Mu: Mean
% Cov: Covariance Matrix

weight = 1/(2*pi)^(size(x,1)/2)/sqrt(norm(Cov));
weight = weight*exp(-0.5*(x-Mu)'*(inv(Cov))*(x-Mu));

end











