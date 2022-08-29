function  X  = GEDinv(p,alpha,beta)
%   GED inverse cumulative function
%   INPUTS
% p      = probability \in [0,1], it can be a matrix
% alpha  = alpha parameter
% beta   = beta parameter
%   OUTPUTS
% X      = F^{-I}(p)

Inv_1=@(u) -(alpha.^beta*gammaincinv((1-2*u),1/beta)).^(1/beta); % F inverse(u) if u<=0 
Inv_2=@(u)  (alpha.^beta*gammaincinv((2*u-1),1/beta)).^(1/beta); % F inverse(u) if u>0
X=zeros(size(p,1),(size(p,2)));
for w=1:size(p,2)
    for z=1:size(p,1)
        if p(z,w)<=0.5
            X(z,w)=Inv_1(p(z,w));
        else X(z,w)=Inv_2(p(z,w));
        end;
    end;
end;
end

