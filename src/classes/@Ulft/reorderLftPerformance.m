function this_lft = reorderLftPerformance(this_lft, new_order)
%% REORDERLFTPERFORMANCE method for reordering the individual performances within a Ulft.
%
%     reordered_lft = reorderLftPerformance(this_lft, new_order)
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft whose performances are to be reordered
%         new_order : array of naturals :: the new ordering of the performance block
%       Output:
%         reordered_lft : Ulft object :: the lft with reordered Performances
%
%     For example, lft = reorderLftPerformance(lft, [2, 1]) will switch the 2nd
%     performance with the 1st performance in lft.performance
%
%     See also Ulft, gatherLft, reorderLftDisturbance, reorderLftDelta

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check correctness of inputs
assert(isa(this_lft, 'Ulft'), 'Ulft:reorderLftPerformance',...
       'First argument must be a Ulft object');
   
number_of_perf = length(this_lft.performance.names);
validateattributes(new_order, {'numeric'},...
                   {'integer', 'positive', 'numel', number_of_perf});
assert(length(unique(new_order)) == number_of_perf,...
       'Ulft:reorderLftPerformance',...
       'Second argument must be an array of nonrepeated integers');

%% Reorder lft
total_time = sum(this_lft.horizon_period);

ending_output_indices = cumsum(this_lft.performance.dim_outs, 1);
starting_output_indices = [ones(1, total_time);
                          ending_output_indices(1:end - 1, :) + 1];

% Loop reordering for each timestep of lft
for i = 1:total_time
    % Generate reordering indices at current timestep
    new_output_indices = [];
    for j = 1:number_of_perf
        new_output_indices = [new_output_indices,...
                             starting_output_indices(new_order(j), i) :...
                             ending_output_indices(new_order(j), i)];
    end
    
    % Reorder lft matrices at current timestep
    this_lft.c{i} = this_lft.c{i}(new_output_indices, :);
    this_lft.d{i} = this_lft.d{i}(new_output_indices, :);
end
this_lft.performance.performances = ...
    this_lft.performance.performances(1, new_order);

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)