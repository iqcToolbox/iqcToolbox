function this_lft = reorderLftDelta(this_lft, new_order)
%% REORDERLFTDELTA method for reordering the individual Deltas within a Ulft
%
%     reordered_lft = reorderLftDelta(this_lft, new_order)
% 
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft whose Deltas are to be reordered
%         new_order : array of naturals :: the new ordering of the Delta block
%       Output:
%         reordered_lft : Ulft object :: the lft with reordered Deltas
%
%     For example, lft = reorderLftDelta(lft, [3, 2]) will switch the 3rd delta
%     with the 2nd delta in lft.delta and ensure that the a, b, and c
%     matrices of lft are appropriately manipulated.
%
%     See also Ulft, gatherLft, reorderLftDisturbance, reorderLftPerformance

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check correctness of inputs
assert(isa(this_lft, 'Ulft'), 'Ulft:reorderLftDelta',...
       'First argument must be a Ulft object');
   
number_of_deltas = length(this_lft.delta.names);
validateattributes(new_order, {'numeric'},...
                   {'integer', 'positive', 'numel', number_of_deltas});
assert(length(unique(new_order)) == number_of_deltas,...
       'Ulft:reorderLftDelta',...
       'Second argument must be an array of nonrepeated integers');

%% Reorder lft
total_time = sum(this_lft.horizon_period);

ending_input_indices = cumsum(this_lft.delta.dim_ins, 1);
starting_input_indices = [ones(1, total_time);
                          ending_input_indices(1:end - 1, :) + 1];
ending_output_indices = cumsum(this_lft.delta.dim_outs, 1);
starting_output_indices = [ones(1, total_time);
                          ending_output_indices(1:end - 1, :) + 1];

% Loop reordering for each timestep of lft
for i = 1:total_time
    % Generate reordering indices at current timestep
    new_input_indices = [];
    new_output_indices = [];
    for j = 1:number_of_deltas
        new_input_indices = [new_input_indices,...
                             starting_input_indices(new_order(j), i) :...
                             ending_input_indices(new_order(j), i)];
        new_output_indices = [new_output_indices,...
                              starting_output_indices(new_order(j), i) :...
                              ending_output_indices(new_order(j), i)];
    end
    
    % Reorder lft matrices at current timestep
    this_lft.a{i} = this_lft.a{i}(new_input_indices, new_output_indices);
    this_lft.b{i} = this_lft.b{i}(new_input_indices, :);
    this_lft.c{i} = this_lft.c{i}(:, new_output_indices);
end
this_lft.delta.deltas = this_lft.delta.deltas(1, new_order);

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)