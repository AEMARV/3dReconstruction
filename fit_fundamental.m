function [F,e1,e2] = fit_fundamental( matches)
%RECONS Summary of this function goes here
%   Detailed explanation goes here
    homMatches = ones(size(matches,1),size(matches,2)+2)';
    homMatches(1:2,:) = matches(:,1:2)';
    homMatches(4:5,:) = matches(:,3:4)';

    % must test the config of matches
    [F,e1,e2] = fundmatrix(homMatches);

end

