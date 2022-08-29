function Y = GEDcdf( X, alpha, beta )
%%   GED cumulative density function
%   INPUTS
% X      = value of X, it could be a matrix
% alpha  = alpha parameter
% beta   = beta parameter
%   OUTPUTS
% Y      = GEDcdf(X)

F=@(z) 1/2+sign(z).*gammainc((abs(z)./alpha).^beta,1/beta,'lower')./2; 
%Matlab already use 1/gamma(1/beta) in gammainc(a,b), so we do not need to
%add it 
Y=F(X);
end

