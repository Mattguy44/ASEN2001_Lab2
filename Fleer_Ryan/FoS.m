function [F_design] = FoS( Max_load, Uncertainty, P_dsr )
% Determines the allowable design force based on the maximum load joints 
% can take, the uncertainty in their strength, and the desired probabiblity
% of failure.
%
% Input:    Max_load = joint failure load
%           Uncertainty = uncertainty in joint strength
%           P_dsr = desired probability of failure
%
% Output:   F_design = maximum allowable joint load
%
% Author: Ethan Fleer, Oct. 21, 2017
% Modified: Matthew Ryan, Nov. 6, 2017

% Calculate factor of safety
F_dsr = icdf('normal', P_dsr, Max_load, Uncertainty);
n = (Max_load - F_dsr)/Uncertainty;
FOS = 1/(1-n*(Uncertainty/Max_load));

% Calculate max allowable load
F_design = Max_load/FOS;

end

