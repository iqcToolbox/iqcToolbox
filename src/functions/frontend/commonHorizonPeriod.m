function horizon_period = commonHorizonPeriod(horizon_periods)
%% COMMONHORIZONPERIOD function for finding a horizon_period that can 
%  capture all the given horizon_periods
%
%  horizon_period = commonHorizonPeriod(horizon_periods)
%
%  Variables:
%  ---------
%     Input:
%       horizon_periods : n x 2 array of naturals :: [horizons, periods] (in timesteps) of sequences
%     Output:
%       horizon_period : 1 x 2 array of naturals :: common [horizon, period] of sequences
%
%  See also makeNewIndices, Ulft.matchHorizonPeriod

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check validity of inputs
if isempty(horizon_periods)
    horizon_period = horizon_periods;
    return;
else
    validateattributes(horizon_periods,...
                       'numeric',...
                       {'integer',...
                        'nonnegative',...
                        'ncols', 2})
    validateattributes(horizon_periods(:,2), 'numeric', {'positive'})
end
%% Calculate horizon_period
num_hp = size(horizon_periods, 1);
period = 1;
for i = 1:num_hp
    period = lcm(period, horizon_periods(i, 2));
end
horizon = max(horizon_periods(:, 1));

horizon_period = [horizon, period];
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)