i_mdl = 0;

n_cols = size(rte{k_rte}.X,2);

% note the 'extra' zdros in these cols beyond N_pars aren't used in the
% fit function (fitnlm excludes rows with NaNs)
% X(1:2,:) = zeros(2,n_cols);
X = rte{k_rte}.X;

if ~TC_rte
    X(1:N_pars,end) = Exp_e_vals';
else
    X(1:N_pars,end) = Exp_e_vals_TCSC';   
end

% X(1:N_pars,end) = rte{k_rte}.X;
% 
% % dummy response for flags
% Y(1:2,1) = zeros(2,1);
% Y(3:2+n_reps,1) = rte{k_rte}.Y;

all_pars = 1:N_pars;

for i_nps = 1:n_test_pars
    
    b_is = nchoosek(test_pars,i_nps);
    
%     if i_nps == 0 % fit with preset constants and no parameters
%         
%         mdl_0 = fittype(@(n1,n2,sig1,sig2,g1,g2,X) ...
%             prod_erf_updated(n1,n2,sig1,sig2,g1,g2,X), ...
%             'problem',all_par_nm_s); % 
%         
%         mdl{k_rtei_mdl} =  fitnlm(X,rte{k_rte}.Y, ...
%             mdl_0, 'problem',const_vals,'start',[]); %0.5); %
%         rte{k_rte,i_mdl}.AICs(i_opt) = aicbic(rte{k_rte}. ...
%             mdl{i_opt}.LogLikelihood,0,n_spec); 
%                     
%     else
    
    for j_b = 1:size(b_is,1)

        % increment model
        i_mdl = i_mdl +1;
        n_ps(i_mdl) = i_nps; % -1;

        par_ij{i_mdl} = b_is(j_b,:);
        
        idx_bs = b_is(j_b,:);

        const_ijs = all_pars(~ismember(all_pars,idx_bs)); 
        % don't fit constant value vars
        X(const_ijs,end-1) = -1;    
        % do fit other vras
        X(idx_bs,end-1) = 1;    

        if ~TC_rte
            
            mdl{k_rte,i_mdl} =  fitnlm(X,rte{k_rte}.Y,mdl_perf,Exp_e_vals(idx_bs)); 
       
        else
            
            mdl{k_rte,i_mdl} =  fitnlm(X,rte{k_rte}.Y,mdl_perf,Exp_e_vals_TCSC(idx_bs)); 

        end
        
       [~,~,ic] = aicbic(mdl{k_rte,i_mdl}.LogLikelihood,i_nps);
       
        AICs(k_rte,i_mdl) = ic.aicc;

        p_fit{k_rte,i_mdl} = mdl{k_rte,i_mdl}.Fitted; % 
        resids{k_rte,i_mdl} = mdl{k_rte,i_mdl}.Residuals.Raw';
        coeffs{k_rte,i_mdl} = mdl{k_rte,i_mdl}.Coefficients;
        adj_r2(k_rte,i_mdl) = mdl{k_rte,i_mdl}.Rsquared.Adjusted;
                
    end

end

% also include null model with all default values
X(1:N_pars,end-1) = -1;
p_Arr_dft{k_rte} = mdl_perf([],X);
R2_dft(k_rte) = corr(p_Arr_dft{k_rte},rte{k_rte}.Y);

[sort_AICs(k_rte,:), sort_mods(k_rte,:)] = sort(AICs(k_rte,:));
best_mdl_idx = sort_mods(k_rte,1);
% p_fit_best{k_rte} = p_fit{k_rte,best_mdl_idx};
close_mdl_idx{k_rte} = sort_mods(k_rte,(1:find(sort_AICs(k_rte,:) - ...
    sort_AICs(k_rte,1)<2,1,'last')));
best_pars{k_rte} = par_ij{best_mdl_idx};
best_mdl{k_rte} =  mdl{k_rte,best_mdl_idx};
best_R2(k_rte) = adj_r2(k_rte,best_mdl_idx);
best_coeffs{k_rte} = coeffs{k_rte,best_mdl_idx};
% end