function [interval] = randinterval(low, high,N_inds)
% returns a random number within specified interval

% previously used for random start locations from polygon

interval = low + ( (high - low) * rand(N_inds,1) );