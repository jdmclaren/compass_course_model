function [Cf_dst, SE_dst, Cf_lat, SE_lat, Cf_ori, SE_ori] = ...
    det_dst_geo_ori_expnts(mdl,istr,jstr)

    if istr == 2 % includes explicit distance-factor (n_steps)
        Cf_dst = mdl{istr,jstr}.Coefficients.Estimate(3);
        SE_dst = mdl{istr,jstr}.Coefficients.SE(3);
    else % no explicit distance-factor
        Cf_dst = 0;
        SE_dst = NaN;       
    end
        
    if jstr ~= 1 % j == 1 is breadth only 
        
        if jstr ~= 3 % j == 3 is orientn only
            Cf_lat = mdl{istr,jstr}.Coefficients.Estimate(2+istr);
            SE_lat = mdl{istr,jstr}.Coefficients.SE(2+istr);    
            if jstr ~= 2 % j == 2 is geo only 
                Cf_ori = mdl{istr,jstr}.Coefficients.Estimate(3+istr);
                SE_ori = mdl{istr,jstr}.Coefficients.SE(3+istr);      
            else
                Cf_ori = 0;
                SE_ori = NaN;                   
            end
        else
            Cf_ori = mdl{istr,jstr}.Coefficients.Estimate(2+istr);
            SE_ori = mdl{istr,jstr}.Coefficients.SE(2+istr);    
            Cf_lat = 0;
            SE_lat = NaN;               
        end

    else % jstr == 1, no lat or ori effects

        Cf_lat = 0;
        SE_lat = NaN;     
        Cf_ori = 0;
        SE_ori = NaN;          
        
    end
    
end