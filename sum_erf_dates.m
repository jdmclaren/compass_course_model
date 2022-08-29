function p_Arr = sum_erf_dates(b,X)

% X(1,:) are breadth, X(2,:) are sig_alf, X(3,:) Geo facts
% N_s are possible # steps, p_ArrLats prob arriving at Arr Lat for those
% steps

% X(:,4:5) are N_min_fls and N_max_fls
% X(:,6) is p_ArrLats

% try 
%     
   
N_s =  X(:,4:5);
p_ArrLats =  X(:,6:end);

% p_Arr = zeros(size(X,2));

for i_rep = 1:size(X,1)
    
    N_0 = N_s(i_rep,1);
    
    switch numel(b)
        
        case 1 % no geo factor
            
            p_Arr(i_rep,1) = p_ArrLats(i_rep,1)*erf(sqrt(0.5)*X(i_rep,1).* ...
                N_0.^b(1)./X(i_rep,2));

            for iN = 1:N_s(i_rep,2)-N_0

                p_Arr(i_rep,1) = p_Arr(i_rep,1) + (1-p_Arr(i_rep,1)).* ...
                    p_ArrLats(i_rep,iN+1).*erf(sqrt(0.5)*X(i_rep,1).* ...
                (N_0+iN).^b(1)./X(i_rep,2));

            end
            
        case 2
            
            p_Arr(i_rep,1) = p_ArrLats(i_rep,1)*erf(sqrt(0.5)*X(i_rep,1).* ...
                N_0.^b(1)./X(i_rep,2).*X(i_rep,3).^b(2));

            for iN = 1:N_s(i_rep,2)-N_0

                p_Arr(i_rep,1) = p_Arr(i_rep,1) + (1-p_Arr(i_rep,1)).* ...
                    p_ArrLats(i_rep,iN+1).*erf(sqrt(0.5)*X(i_rep,1).* ...
                (N_0+iN).^b(1)./X(i_rep,2).*X(i_rep,3).^b(2));

            end
            
        case 3 % incl exponent for sigma
            
            p_Arr(i_rep,1) = p_ArrLats(i_rep,1)*erf(sqrt(0.5)*X(i_rep,1).* ...
                N_0.^b(1)./X(i_rep,2).^b(3).*X(i_rep,3).^b(2));

            for iN = 1:N_s(i_rep,2)-N_0

                p_Arr(i_rep,1) = p_Arr(i_rep,1) + (1-p_Arr(i_rep,1)).* ...
                    p_ArrLats(i_rep,iN+1).*erf(sqrt(0.5)*X(i_rep,1).* ...
                (N_0+iN).^b(1)./(X(i_rep,2).*X(i_rep,3).^b(2)).^b(3));

            end
            
        case 4 % incl exponent for sigma & beta
            
%             try
                
            p_Arr(i_rep,1) = p_ArrLats(i_rep,1)*erf(sqrt(0.5)*X(i_rep,1).^b(4).* ...
                N_0.^b(1)./X(i_rep,2).^b(3).*X(i_rep,3).^b(2));

            for iN = 1:N_s(i_rep,2)-N_0

                p_Arr(i_rep,1) = p_Arr(i_rep,1) + (1-p_Arr(i_rep,1)).* ...
                    p_ArrLats(i_rep,iN+1).*erf(sqrt(0.5)*X(i_rep,1).^b(4).* ...
                (N_0+iN).^b(1)./(X(i_rep,2).*X(i_rep,3).^b(2)).^b(3));

            end
            
%             catch
%                 
%                 keyboard
%                 
%             end
    
    end
    
end

% catch
% 
% keyboard
% 


if any(isnan(p_Arr) | isinf(p_Arr))
    p_Arr(isnan(p_Arr)) = 0;
    p_Arr(isinf(p_Arr)) = 1;
end
