function Y = GEDpdf( X, alpha, beta )
%%   GED probabily density function
%   INPUTS
% X      = value of X, it could be a matrix
% alpha  = alpha parameter
% beta   = beta parameter
%   OUTPUTS
% Y      = GEDpdf(X)

f=@(z) beta*exp(-(abs(z)/alpha).^beta)./(2*alpha*gamma(1/beta));
Y=f(X);
end

