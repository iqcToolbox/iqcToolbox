function [indices, new_horizon_period] = makeNewIndices(old_horizon_period, new_horizon_period)
%% MAKENEWINDICES function for creating the indices that map sequence A
%  (with horizon_period A) to sequence B (with horizon_period B), such 
%  that B = A(indices)
%
%  indices = makeNewIndices(new_horizon_period, old_horizon_period)
%
%  Variables:
%  ---------
%     Input:
%       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
%       old_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of old sequence
%     Output:
%       indices : array of naturals :: map of indices from sequence in old_horizon_period to
%                                      sequence in new_horizon_period
%       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] of common horizon_period that
%                                                        captures input horizon_periods
%
%  See also commonHorizonPeriod

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check old_horizon_period is valid 
validateattributes(old_horizon_period, 'numeric', {'size', [1,2],...
                                                   'integer',...
                                                   'nonnegative'});
validateattributes(old_horizon_period(2), 'numeric', {'positive'});

%% Check new_horizon_period is valid
validateattributes(new_horizon_period, 'numeric', {'size', [1,2],...
                                                   'integer',...
                                                   'nonnegative'});
assert(old_horizon_period(1) <= new_horizon_period(1),...
       'DeltaBounded:matchHorizonPeriod',...
       ['new_horizon_period(1) must be greater than or equal to ',...
        'current horizon_period(1)']);
assert(old_horizon_period(2) <= new_horizon_period(2),...
       'DeltaBounded:matchHorizonPeriod',...
       ['new_horizon_period(2) must be greater than or equal to ',...
        'old_horizon_period(2)']);

%% Make indices for mapping the old sequence to a new sequence
new_horizon_period(2) = lcm(old_horizon_period(2), new_horizon_period(2));
new_total_time = sum(new_horizon_period);
repeated_indices = old_horizon_period(1) + (1:old_horizon_period(2));
repititions = ceil(new_total_time / old_horizon_period(2));
indices = [1:old_horizon_period(1),...
           repmat(repeated_indices, 1, repititions)];
indices = indices(1, 1:new_total_time);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)