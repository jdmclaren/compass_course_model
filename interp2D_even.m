function finalPathXY = interp2D_even(pathXY,mult_fact)
% interprets 2D area with multiple number of points fac
    n_points = mult_fact*size(pathXY,1);
    stepLengths = sqrt(sum(diff(pathXY,[],1).^2,2));
    stepLengths = [0; stepLengths]; % add the starting point
    cumulativeLen = cumsum(stepLengths);
    % try doubling total # points (equally spaced in Lat / Lon)
    finalStepLocs = linspace(0,cumulativeLen(end), n_points);
    finalPathXY = interp1(cumulativeLen, pathXY, finalStepLocs);