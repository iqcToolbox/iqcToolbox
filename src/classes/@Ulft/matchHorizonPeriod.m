function this_lft = matchHorizonPeriod(this_lft, new_horizon_period)
%% MATCHHORIZONPERIOD function to change Ulft object to match a new horizon_period.
%     This will change the horizon_period and associated properties of the lft.
%
%     this_lft = matchHorizonPeriod(this_lft, new_horizon_period)
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft to be modified with a new horizon_period
%         new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
%       Output:
%         this_lft : Ulft object :: the lft with a new horizon_period
%
%     See also Ulft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check correctness and consistency of inputs
validateattributes(new_horizon_period, {'numeric'}, {'size', [1, 2],...
                                                   'integer',...
                                                   'nonnegative'})
[indices, new_horizon_period] = makeNewIndices(this_lft.horizon_period,...
                                               new_horizon_period);

%% Gather properties with new horizon_period
a = this_lft.a(indices);
b = this_lft.b(indices);
c = this_lft.c(indices);
d = this_lft.d(indices);
delta = matchHorizonPeriod(this_lft.delta, new_horizon_period);
performance = matchHorizonPeriod(this_lft.performance,...
                                 new_horizon_period);
disturbance = matchHorizonPeriod(this_lft.disturbance,...
                                 new_horizon_period);

%% Construct new Ulft                             
this_lft = Ulft(a, b, c, d, delta,...
           'horizon_period', new_horizon_period,...
           'performance', performance,...
           'disturbance', disturbance);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)