function this_lft = reorderLftDisturbance(this_lft, new_order)
%% REORDERLFTDISTURBANCE method for reordering the individual disturbances within a Ulft.
%
%     reordered_lft = reorderLftDisturbance(this_lft, new_order)
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft whose disturbances are to be reordered
%         new_order : array of naturals :: the new ordering of the disturbance block
%       Output:
%         reordered_lft : Ulft object :: the lft with reordered Disturbances
%
%     For example, lft = reorderLftDisturbance(lft, [2, 1]) will switch the 2nd
%     disturbance with the 1st disturbance in lft.disturbance
%
%     See also Ulft, gatherLft, reorderLftDelta, reorderLftPerformance

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check correctness of inputs
assert(isa(this_lft, 'Ulft'), 'Ulft:reorderLftDisturbance',...
       'First argument must be a Ulft object');
   
number_of_dis = length(this_lft.disturbance.names);
validateattributes(new_order, {'numeric'},...
                   {'integer', 'positive', 'numel', number_of_dis});
assert(length(unique(new_order)) == number_of_dis,...
       'Ulft:reorderLftDisturbance',...
       'Second argument must be an array of nonrepeated integers');

%% Reorder lft
total_time = sum(this_lft.horizon_period);

ending_input_indices = cumsum(this_lft.disturbance.dim_ins, 1);
starting_input_indices = [ones(1, total_time);
                          ending_input_indices(1:end - 1, :) + 1];

% Loop reordering for each timestep of lft
for i = 1:total_time
    % Generate reordering indices at current timestep
    new_input_indices = [];
    for j = 1:number_of_dis
        new_input_indices = [new_input_indices,...
                             starting_input_indices(new_order(j), i) :...
                             ending_input_indices(new_order(j), i)];
    end
    
    % Reorder lft matrices at current timestep
    this_lft.b{i} = this_lft.b{i}(:, new_input_indices);
    this_lft.d{i} = this_lft.d{i}(:, new_input_indices);
end
this_lft.disturbance.disturbances = ...
    this_lft.disturbance.disturbances(1, new_order);

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)