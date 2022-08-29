function p_erf = prod_exp_erf_all_errs_geo(b,X)

% p_erf = erf(b(1)*prod(X.^b(2:end),2));
p_erf = erf(prod((X(:,1).*(X(:,end).^b(end))).^b(1).* ...
    X(:,2:end-1).^b(2:end-1),2));

if any(isnan(p_erf) | isinf(p_erf))
    p_erf(isnan(p_erf)) = 0;
    p_erf(isinf(p_erf)) = 1;
end
% try
%     
% if ~any(isinf(p_arg)) && ~any(isnan(p_arg))
%     p_erf = erf(p_arg);
% else
%     [idx_inf, jdx_inf] = find(isinf(X.^b(2:end)));
%     u_iis = unique(idx_inf);
%     for ii = 1:numel(u_iis)
%         u_ii = u_iis(ii);
%         j_iis = jdx_inf(idx_inf==u_ii);
%         p_erf(u_ii,1) = all(b(j_iis) > 0 & X(u_ii, j_iis) >=1);
%     end
%     [idx_nan, jdx_nan] = find(isnan(X.^b(2:end)));
%     u_iis = unique(idx_nan);
%     for ii = 1:numel(u_iis)
%         u_ii = u_iis(ii);
%         j_iis = jdx_nan(idx_inf==u_ii);
%         p_erf(u_ii,1) = all(b(j_iis) > 0 & X(u_ii, j_iis) >=1);
%     end
% %     p_erf(idx_nan,1) = b(jdx_nan) > 0 & X(idx_nan, jdx_nan) >=1;   
%     p_erf(~isinf(p_arg) & ~isnan(p_arg),1) = ...
%         erf(p_arg(~isinf(p_arg) & ~isnan(p_arg)));
% end
% 
% if any(isinf(p_erf)) || any(isnan(p_erf))
%     keyboard
% end
% 
% catch
%     
%     keyboard
%     
% end    
%     